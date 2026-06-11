#' Bootstrap helper utilities
#'
#' Shared functions for computing bootstrap confidence intervals and
#' plotting bootstrap distributions. Used internally by all bootstrap
#' CI analyses (SingleMeanCI, SinglePropCI, twomeanCI, pairedmeanCI,
#' slopeCI, TwoPropCI).
#'
#' @name bootstrap_utils
#' @import ggplot2
#' @import ggforce
#' @import dplyr
#' @import tibble
NULL

#' Compute a bootstrap confidence interval
#'
#' Calculates either a percentile CI or a bootstrap-SE CI from a
#' bootstrap distribution produced by the \pkg{infer} package.
#'
#' @param boot A data frame with a \code{stat} column containing
#'   bootstrap replicate statistics (as produced by \code{infer::calculate()}).
#' @param obs_stat The observed sample statistic.
#' @param conf_level Confidence level. Values > 1 are treated as percentages
#'   (e.g. 95 is converted to 0.95).
#' @param ci_type Either \code{"bootperc"} for the percentile method or
#'   \code{"bootse"} for the bootstrap-SE method.
#' @param clamp Optional length-2 numeric vector giving lower and upper bounds
#'   for the CI (e.g. \code{c(0, 1)} for proportions).
#'
#' @return A list with components \code{cil}, \code{ciu}, \code{se}, and
#'   \code{zcrit}. For the percentile method, \code{se} and \code{zcrit}
#'   are \code{NULL}.
#'
#' @keywords internal
compute_boot_ci <- function(boot, obs_stat, conf_level, ci_type,
                            clamp = NULL) {
    if (conf_level > 1)
        conf_level <- conf_level / 100

    if (ci_type == "bootperc") {
        ci <- boot |>
            dplyr::summarise(
                cil = quantile(stat, (1 - conf_level) / 2, na.rm = TRUE),
                ciu = quantile(stat, 1 - (1 - conf_level) / 2, na.rm = TRUE)
            )
        se <- NULL
        zcrit <- NULL
    } else {
        zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
        se <- boot |>
            dplyr::summarise(se = sd(stat, na.rm = TRUE)) |>
            dplyr::pull()
        ci <- tibble::tibble(
            cil = obs_stat - zcrit * se,
            ciu = obs_stat + zcrit * se
        )
    }

    cil <- ci |> dplyr::pull(cil)
    ciu <- ci |> dplyr::pull(ciu)

    if (!is.null(clamp)) {
        cil <- max(cil, clamp[1])
        ciu <- min(ciu, clamp[2])
    }

    list(cil = cil, ciu = ciu, se = se, zcrit = zcrit)
}

#' Plot a bootstrap distribution with confidence interval
#'
#' Creates either a dotplot or histogram of bootstrap replicate statistics
#' with dashed vertical lines showing the CI bounds and an explanatory
#' caption describing how the CI was constructed.
#'
#' Bins are centered on the observed statistic with a width adapted to the
#' discreteness of the replicates (see \code{\link{choose_binning}}), so
#' bars stay contiguous and evenly spaced even for discrete statistics
#' like proportions.  When the replicates lie on a lattice, the drawn CI
#' lines are nudged outward to the nearest bin edge so they never slice
#' through a bar (the reported CI bounds are unaffected; this is display
#' only).
#'
#' @inheritParams compute_boot_ci
#' @param dotHist Either \code{"dotplot"} or \code{"histogram"}.
#' @param xlab Label for the x-axis (e.g. \code{"mean"}, \code{"slope"}).
#' @param stat_label Descriptive label for the caption (e.g.
#'   \code{"bootstrap means"}).
#' @param show_lines Logical; if \code{FALSE}, omit the dashed vertical lines
#'   marking the CI bounds. Default \code{TRUE}.
#' @param show_caption Logical; if \code{FALSE}, omit the caption explaining
#'   how the CI was constructed. Default \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
plot_boot_dist <- function(boot, obs_stat, conf_level, ci_type,
                           dotHist = c("dotplot", "histogram"),
                           xlab = "statistic",
                           stat_label = "bootstrap statistics",
                           clamp = NULL,
                           show_lines = TRUE,
                           show_caption = TRUE) {
    dotHist <- match.arg(dotHist)
    ci <- compute_boot_ci(boot, obs_stat, conf_level, ci_type, clamp)
    cil <- ci$cil
    ciu <- ci$ciu

    if (conf_level > 1) conf_level <- conf_level / 100

    if (ci_type == "bootperc") {
        caption <- paste0(
            "CI limits (dashed) are ",
            round((1 - conf_level) / 2 * 100, 1),
            "% and ",
            round((1 - (1 - conf_level) / 2) * 100, 1),
            "%",
            "\n percentiles of ", stat_label
        )
    } else {
        caption <- paste0(
            "CI (dashed) is calculated using SE of ", stat_label,
            "\n (SE = ", round(ci$se, 3), ", z* = ", round(ci$zcrit, 3), ")"
        )
    }

    b <- choose_binning(boot$stat, obs_stat, align = "center")
    idx <- bin_index(boot$stat, obs_stat, b$bw, b$off, 1, "center")

    bars <- data.frame(idx = idx) |>
        dplyr::count(idx)
    xr <- bin_xrange(bars$idx, obs_stat, b$bw, 1, "center")
    bars$xmin <- xr$xmin
    bars$xmax <- xr$xmax

    # For lattice data, nudge the drawn CI lines outward to the nearest bin
    # edge so they never cut through a bar.  Display only; the reported CI
    # is unchanged.
    if (b$lattice) {
        eps <- 1e-7
        cila <- obs_stat + (floor((cil - obs_stat) / b$bw + 0.5 + eps) - 0.5) * b$bw
        ciua <- obs_stat + (ceiling((ciu - obs_stat) / b$bw - 0.5 - eps) + 0.5) * b$bw
    } else {
        cila <- cil
        ciua <- ciu
    }

    pad <- 0.03 * (max(bars$xmax, ciua) - min(bars$xmin, cila))
    xlims <- c(min(bars$xmin, cila) - pad, max(bars$xmax, ciua) + pad)

    if (dotHist == "dotplot") {
        dots <- data.frame(idx = idx) |>
            dplyr::group_by(idx) |>
            dplyr::mutate(y = dplyr::row_number()) |>
            dplyr::ungroup()
        dots$x <- bin_xrange(dots$idx, obs_stat, b$bw, 1, "center")$mid

        # At least 0.6 above the tallest stack so the top dot (semi-height
        # up to 0.45) is never clipped by the y limit
        y_top <- max(max(dots$y) * 1.06, max(dots$y) + 0.6)
        a <- 0.42 * b$bw
        semi_h <- dot_semi_height(a, diff(xlims), y_top)

        p <- ggplot2::ggplot(dots) +
            ggforce::geom_ellipse(
                ggplot2::aes(x0 = x, y0 = y, a = a, b = semi_h, angle = 0),
                fill = "grey35", color = NA, n = 36, show.legend = FALSE) +
            ggplot2::theme_minimal() +
            ggplot2::ylab("count") +
            ggplot2::xlab(xlab) +
            ggplot2::scale_x_continuous(limits = xlims) +
            ggplot2::scale_y_continuous(
                limits = c(0, y_top),
                expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
            ggplot2::theme(text = ggplot2::element_text(size = 14))
    } else {
        p <- ggplot2::ggplot(bars) +
            ggplot2::geom_rect(
                ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = n),
                fill = "grey35", color = "white", linewidth = 0.3,
                show.legend = FALSE) +
            ggplot2::theme_minimal() +
            ggplot2::ylab("count") +
            ggplot2::xlab(xlab) +
            ggplot2::scale_x_continuous(limits = xlims) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0, 0.04))) +
            ggplot2::theme(text = ggplot2::element_text(size = 14))
    }

    if (show_lines)
        p <- p + ggplot2::geom_vline(xintercept = c(cila, ciua), linetype = "dashed", color = "red")
    if (show_caption)
        p <- p + ggplot2::labs(caption = caption) +
            ggplot2::theme(plot.caption = ggplot2::element_text(color = "red", hjust = 0))

    p
}
