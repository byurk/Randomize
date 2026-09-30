#' Bootstrap helper utilities
#'
#' Shared functions for computing bootstrap confidence intervals and
#' plotting bootstrap distributions. Used internally by all bootstrap
#' CI analyses (SingleMeanCI, SinglePropCI, twomeanCI, pairedmeanCI,
#' slopeCI, TwoPropCI).
#'
#' @name bootstrap_utils
#' @import ggplot2
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
#' discreteness of the replicates (see \code{\link{choose_binning}}),
#' chosen so the histogram never shows interior empty bins.  The dashed
#' CI lines are drawn at the exact reported bounds; since no tail is
#' shaded on CI plots, a line falling inside a bar is fine (and
#' preferable to nudging the line away from the true percentile).
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
#' @param show_counts Logical; if \code{TRUE}, print the count above
#'   each bar or dot stack.  Default \code{FALSE}.
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
                           show_caption = TRUE,
                           show_counts = FALSE) {
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

    # Sparse discrete replicates cannot be binned gap-free; draw one bar
    # or dot column per exact value (see plot_null_dist for rationale).
    sparse <- is_sparse_values(boot$stat)

    if (sparse) {
        g <- sparse_groups(boot$stat)
        u <- g$centers
        a <- 0.4 * min(diff(u))
        idx <- g$idx
        dot_x <- u[idx]
        bars <- data.frame(idx = idx) |>
            dplyr::count(idx)
        bars$xmin <- u[bars$idx] - a
        bars$xmax <- u[bars$idx] + a
    } else {
        b <- choose_binning(boot$stat, obs_stat, align = "center")
        idx <- bin_index(boot$stat, obs_stat, b$bw, b$off, 1, "center")
        a <- 0.42 * b$bw
        dot_x <- dot_column_x(boot$stat, idx, b, obs_stat, 1, "center")
        # a bin straddling a domain bound would place its column at an
        # impossible value; clamp the center to the bound instead
        if (!is.null(clamp))
            dot_x <- pmin(pmax(dot_x, clamp[1]), clamp[2])

        bars <- data.frame(idx = idx) |>
            dplyr::count(idx)
        xr <- bin_xrange(bars$idx, obs_stat, b$bw, 1, "center")
        bars$xmin <- xr$xmin
        bars$xmax <- xr$xmax
    }

    # Trim to the statistic's domain (the CI clamp bounds) so no bar
    # implies impossible values, e.g. proportions outside [0, 1].
    # Center-aligned bins put values at bin centers, so a boundary bar
    # trims to an honest half-bar rather than vanishing.
    if (!is.null(clamp)) {
        bars$xmin <- pmax(bars$xmin, clamp[1])
        bars$xmax <- pmin(bars$xmax, clamp[2])
    }

    pad <- 0.03 * (max(bars$xmax, ciu) - min(bars$xmin, cil))
    xlims <- c(min(bars$xmin, cil) - pad, max(bars$xmax, ciu) + pad)

    if (dotHist == "dotplot") {
        dots <- data.frame(idx = idx) |>
            dplyr::group_by(idx) |>
            dplyr::mutate(y = dplyr::row_number()) |>
            dplyr::ungroup()
        dots$x <- dot_x

        # At least 0.6 above the tallest stack so the top dot (semi-height
        # up to 0.45) is never clipped by the y limit; count labels need
        # a little more headroom
        cnt_room <- if (!show_counts) 1.06 else if (counts_vertical(bars$n)) 1.26 else 1.14
        y_top <- max(max(dots$y) * cnt_room,
                     max(dots$y) + 0.6 + (if (show_counts) 0.8 else 0))

        p <- ggplot2::ggplot(dots) +
            stack_dots(dots, a = a, fill = "grey35") +
            # floor-sized dots on very tall stacks may poke past y_top
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::theme_minimal() +
            ggplot2::ylab("count") +
            ggplot2::xlab(xlab) +
            ggplot2::scale_x_continuous(limits = xlims) +
            ggplot2::scale_y_continuous(
                limits = c(0, y_top),
                expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
            ggplot2::theme(text = ggplot2::element_text(size = 14))

        if (show_counts) {
            tops <- dots |>
                dplyr::group_by(idx) |>
                dplyr::summarize(x = x[1], n = max(y), .groups = "drop")
            p <- p + count_labels(tops$x, tops$n + 0.5, tops$n)
        }
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
                expand = ggplot2::expansion(mult = c(0, if (!show_counts) 0.04
                    else if (counts_vertical(bars$n)) 0.24 else 0.12))) +
            ggplot2::theme(text = ggplot2::element_text(size = 14))

        if (show_counts)
            p <- p + count_labels((bars$xmin + bars$xmax) / 2, bars$n, bars$n)
    }

    if (show_lines)
        p <- p + ggplot2::geom_vline(xintercept = c(cil, ciu), linetype = "dashed", color = "red")
    if (show_caption)
        p <- p + ggplot2::labs(caption = caption) +
            ggplot2::theme(plot.caption = ggplot2::element_text(color = "red", hjust = 0))

    p
}

#' Bootstrap distribution of a difference in means (no formulas)
#'
#' Resamples whole rows with replacement, as \code{infer::generate(type =
#' "bootstrap")} does for a two-variable specification, and computes the
#' difference in group means.  See \code{\link{resample_means}} for why
#' \pkg{infer} cannot be used for this in Jamovi 2.7.  Replicates in which
#' a group happens to be empty (possible only for tiny samples) are
#' dropped.
#'
#' @inheritParams resample_means
#' @return A data frame with \code{replicate} and \code{stat} columns.
#' @keywords internal
bootstrap_diff_means <- function(dep, group, levels, reps) {
    group <- as.character(group)
    n <- length(dep)
    stat <- vapply(seq_len(reps), function(i) {
        s <- sample.int(n, n, replace = TRUE)
        d <- dep[s]; g <- group[s]
        mean(d[g == levels[1]]) - mean(d[g == levels[2]])
    }, numeric(1))
    out <- data.frame(replicate = seq_len(reps), stat = stat)
    out[is.finite(out$stat), ]
}
