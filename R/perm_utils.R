#' Null distribution helper utilities
#'
#' Shared functions for seed handling, computing p-values from permutation
#' (or draw-based) null distributions, and plotting those distributions.
#' Used internally by all hypothesis test analyses.
#'
#' @name perm_utils
#' @import ggplot2
#' @import ggforce
#' @import dplyr
NULL

#' Strip infer attributes from a tibble
#'
#' The \pkg{infer} package attaches formula objects with large environment
#' chains to its result tibbles.  When these are stored via
#' \code{jmvcore::setState()}, the serialized size can exceed Jamovi's
#' internal message limits, preventing the frontend from receiving the
#' COMPLETE status.  This helper reduces the tibble to a plain data frame
#' containing only the \code{stat} column.
#'
#' @param x An infer tibble with a \code{stat} column.
#' @return A plain \code{data.frame} with one column \code{stat}.
#' @keywords internal
strip_infer <- function(x) {
    data.frame(stat = x$stat)
}

#' Set the RNG seed conditionally
#'
#' Consolidates the seed-setting pattern used by every analysis.
#' When \code{seedBool} is \code{TRUE}, sets a fixed seed for
#' reproducibility; otherwise resets to random seeding.
#'
#' @param seedBool Logical; whether to use a fixed seed.
#' @param rngSeed Integer seed value (used only when \code{seedBool} is
#'   \code{TRUE}).
#'
#' @keywords internal
set_seed_if <- function(seedBool, rngSeed) {
    if (seedBool) {
        set.seed(rngSeed)
    } else {
        set.seed(NULL)
    }
}

#' Map Jamovi hypothesis option strings to infer direction strings
#'
#' Translates the various hypothesis option names used across the module's
#' analyses into the direction strings expected by
#' \code{\link[infer]{get_p_value}}.
#'
#' @param hypothesis A string from the Jamovi UI: \code{"oneGreater"},
#'   \code{"twoGreater"}, \code{"different"}, \code{"greater"},
#'   \code{"less"}, or \code{"notequal"}.
#'
#' @return One of \code{"greater"}, \code{"less"}, or \code{"two_sided"}.
#'
#' @keywords internal
map_direction <- function(hypothesis) {
    switch(hypothesis,
        "oneGreater" = "greater",
        "greater"    = "greater",
        "twoGreater" = "less",
        "less"       = "less",
        "two_sided"
    )
}

#' Compute a p-value from a null distribution
#'
#' Wrapper around \code{\link[infer]{get_p_value}} that extracts
#' the numeric p-value from the result tibble.
#'
#' @param perms A data frame with a \code{stat} column containing
#'   simulated statistics under the null hypothesis (as produced by
#'   \code{infer::calculate()}).
#' @param obs_stat The observed test statistic.
#' @param direction One of \code{"greater"}, \code{"less"}, or
#'   \code{"two_sided"}.
#'
#' @return A single numeric p-value.
#'
#' @keywords internal
compute_null_pval <- function(perms, obs_stat, direction) {
    perms |>
        infer::get_p_value(obs_stat = obs_stat, direction = direction) |>
        dplyr::pull()
}

#' Plot a null distribution with p-value shading
#'
#' Creates either a dotplot or histogram of simulated null-distribution
#' statistics with the observed statistic marked by a red dashed line.
#' The tail used for the p-value is shaded in red, with a caption
#' explaining the calculation.
#'
#' Bins are anchored so that one bin edge falls exactly at the observed
#' statistic and each simulated value is clamped to the bins on its own
#' side, so no bar (or dot stack) ever mixes values that count toward the
#' p-value with values that do not, and the shaded bars always lie
#' entirely beyond the observed-value line.  The bin width adapts to the
#' discreteness of the statistics (see \code{\link{choose_binning}}) so
#' bars stay contiguous and evenly spaced.
#'
#' @inheritParams compute_null_pval
#' @param dotHist Either \code{"dotplot"} or \code{"histogram"}.
#' @param xlab Label for the x-axis (e.g. \code{"difference (group 1 - group 2)"}).
#' @param obs_label Annotation text placed at the observed statistic line
#'   (e.g. \code{"Observed\\nDifference"}).
#' @param show_line Logical; if \code{FALSE}, omit the dashed vertical line
#'   at the observed statistic. Default \code{TRUE}.
#' @param show_label Logical; if \code{FALSE}, omit the text annotation
#'   labelling the observed statistic. Default \code{TRUE}.
#' @param show_caption Logical; if \code{FALSE}, omit the caption explaining
#'   how the p-value is calculated. Default \code{TRUE}.
#' @param show_tail Logical; if \code{FALSE}, do not shade the tail region
#'   used for the p-value calculation. Default \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
plot_null_dist <- function(perms, obs_stat, direction,
                           dotHist = c("dotplot", "histogram"),
                           xlab = "statistic",
                           obs_label = "Observed\nStatistic",
                           show_line = TRUE,
                           show_label = TRUE,
                           show_caption = TRUE,
                           show_tail = TRUE) {
    dotHist <- match.arg(dotHist)

    if (direction == "less") {
        ptail <- "lt"
        caption <- "One-sided: p-val is proportion of\n results \u2264 observed value"
    } else if (direction == "greater") {
        ptail <- "rt"
        caption <- "One-sided: p-val is proportion of\n results \u2265 observed value"
    } else {
        ptail <- perms |>
            dplyr::summarize(lt = mean(stat < obs_stat), rt = mean(stat > obs_stat)) |>
            dplyr::mutate(pt = dplyr::if_else(rt < lt, "rt", "lt")) |>
            dplyr::pull(pt)
        if (ptail == "lt")
            caption <- "Two-sided: p-val is 2\u00D7 proportion of\n results \u2264 observed value"
        else
            caption <- "Two-sided: p-val is 2\u00D7 proportion of\n results \u2265 observed value"
    }

    sgn <- if (ptail == "lt") -1 else 1

    # Same comparison used for the p-value; this flag (not floating-point
    # bin arithmetic) decides which side of the line a value is drawn on.
    extreme <- if (ptail == "lt") perms$stat <= obs_stat else perms$stat >= obs_stat

    b <- choose_binning(perms$stat, obs_stat, align = "edge", sign = sgn)
    idx <- bin_index(perms$stat, obs_stat, b$bw, b$off, sgn, "edge")
    idx[extreme] <- pmax(idx[extreme], 0L)
    idx[!extreme] <- pmin(idx[!extreme], -1L)

    fill <- extreme & show_tail

    bars <- data.frame(idx = idx, fill = fill) |>
        dplyr::count(idx, fill)
    xr <- bin_xrange(bars$idx, obs_stat, b$bw, sgn, "edge")
    bars$xmin <- xr$xmin
    bars$xmax <- xr$xmax
    bars$mid <- xr$mid

    pad <- 0.03 * (max(bars$xmax, obs_stat) - min(bars$xmin, obs_stat))
    xlims <- c(min(bars$xmin, obs_stat) - pad, max(bars$xmax, obs_stat) + pad)

    fill_scale <- ggplot2::scale_fill_manual(
        values = c("FALSE" = "black", "TRUE" = "#ff8c8c"), guide = "none")

    if (dotHist == "dotplot") {
        dots <- data.frame(idx = idx, fill = fill) |>
            dplyr::group_by(idx) |>
            dplyr::mutate(y = dplyr::row_number()) |>
            dplyr::ungroup()
        dots$x <- bin_xrange(dots$idx, obs_stat, b$bw, sgn, "edge")$mid

        max_stack <- max(dots$y)
        # At least 0.6 above the tallest stack so the top dot (semi-height
        # up to 0.45) is never clipped by the y limit
        y_top <- max(max_stack * (if (show_label) 1.3 else 1.06),
                     max_stack + 0.6)
        a <- 0.42 * b$bw
        semi_h <- dot_semi_height(a, diff(xlims), y_top)

        p <- ggplot2::ggplot(dots) +
            ggforce::geom_ellipse(
                ggplot2::aes(x0 = x, y0 = y, a = a, b = semi_h, angle = 0, fill = fill),
                color = NA, n = 36, show.legend = FALSE) +
            fill_scale +
            ggplot2::theme_minimal() +
            ggplot2::ylab("count") +
            ggplot2::xlab(xlab) +
            ggplot2::scale_x_continuous(limits = xlims) +
            ggplot2::scale_y_continuous(
                limits = c(0, y_top),
                expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
            ggplot2::theme(text = ggplot2::element_text(size = 14))
    } else {
        y_top <- max(bars$n) * (if (show_label) 1.2 else 1.04)

        p <- ggplot2::ggplot(bars) +
            ggplot2::geom_rect(
                ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = n, fill = fill),
                color = "white", linewidth = 0.3, show.legend = FALSE) +
            fill_scale +
            ggplot2::theme_minimal() +
            ggplot2::ylab("count") +
            ggplot2::xlab(xlab) +
            ggplot2::scale_x_continuous(limits = xlims) +
            ggplot2::scale_y_continuous(
                limits = c(0, y_top),
                expand = ggplot2::expansion(mult = c(0, 0.02))) +
            ggplot2::theme(text = ggplot2::element_text(size = 14))
    }

    if (show_line)
        p <- p + ggplot2::geom_vline(xintercept = obs_stat, linetype = "dashed", color = "red")
    if (show_label)
        p <- p + ggplot2::annotate("text", x = obs_stat, y = y_top, vjust = "top",
                                   hjust = inward_hjust(obs_stat, xlims),
                                   label = obs_label, color = "red")
    if (show_caption)
        p <- p + ggplot2::labs(caption = caption) +
            ggplot2::theme(plot.caption = ggplot2::element_text(color = "red", hjust = 0))

    p
}
