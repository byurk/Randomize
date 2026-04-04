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
                           show_caption = TRUE) {
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

    if (dotHist == "dotplot") {
        ndist <- dplyr::n_distinct(perms$stat)
        bw <- perms |>
            dplyr::summarize(min = min(stat), max = max(stat)) |>
            dplyr::mutate(bw = (max - min) / 30) |>
            dplyr::pull(bw)

        # Guard against zero bin width (all values identical)
        if (bw == 0) bw <- abs(obs_stat) * 0.01 + 0.001

        if (ndist <= 30) {
            perms <- perms |> dplyr::mutate(x.bin = stat)
        } else {
            if (ptail == "lt") {
                perms <- perms |>
                    dplyr::mutate(x.bin = obs_stat - ((obs_stat - stat) %/% bw) * bw)
            } else {
                perms <- perms |>
                    dplyr::mutate(x.bin = obs_stat + ((stat - obs_stat) %/% bw) * bw)
            }
        }

        perms <- perms |>
            dplyr::mutate(extreme = (ptail == "lt" & stat <= obs_stat) | (ptail == "rt" & stat >= obs_stat)) |>
            dplyr::group_by(x.bin) |>
            dplyr::mutate(y = seq_along(x.bin)) |>
            dplyr::ungroup()

        lab_ht <- max(max(perms$y) - 1, 3)
        bigx <- max(perms$x.bin)
        littlex <- min(perms$x.bin)

        x_lo <- min(littlex - bw, obs_stat - bw)
        x_hi <- max(bigx + bw, obs_stat + bw)

        p <- ggplot2::ggplot(perms) +
            ggforce::geom_ellipse(ggplot2::aes(x0 = x.bin, y0 = y, a = bw / 3, b = 0.5, angle = 0,
                                               fill = extreme, color = extreme),
                                  show.legend = FALSE) +
            ggplot2::scale_fill_manual(values = c("FALSE" = "black", "TRUE" = "#ff8c8c")) +
            ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#ff8c8c")) +
            ggplot2::theme_minimal() +
            ggplot2::ylab("count") +
            ggplot2::xlab(xlab) +
            ggplot2::ylim(0, lab_ht + 3) +
            ggplot2::xlim(x_lo, x_hi) +
            ggplot2::coord_equal(ratio = bw * 2 / 3) +
            ggplot2::theme(text = ggplot2::element_text(size = 14))

        if (show_line)
            p <- p + ggplot2::geom_vline(xintercept = obs_stat, linetype = "dashed", color = "red")
        if (show_label)
            p <- p + ggplot2::annotate("text", x = obs_stat, y = lab_ht, label = obs_label, color = "red")
        if (show_caption)
            p <- p + ggplot2::labs(caption = caption) +
                ggplot2::theme(plot.caption = ggplot2::element_text(color = "red", hjust = 0))
    } else {
        closed <- ifelse(ptail == "lt", "right", "left")
        perms <- perms |>
            dplyr::mutate(extreme = (ptail == "lt" & stat <= obs_stat) | (ptail == "rt" & stat >= obs_stat))
        p <- ggplot2::ggplot(perms, ggplot2::aes(x = stat, fill = extreme)) +
            ggplot2::geom_histogram(boundary = obs_stat, closed = closed, show.legend = FALSE) +
            ggplot2::scale_fill_manual(values = c("FALSE" = "black", "TRUE" = "#ff8c8c")) +
            ggplot2::theme_minimal() +
            ggplot2::xlab(xlab) +
            ggplot2::ylab("count") +
            ggplot2::theme(text = ggplot2::element_text(size = 14))

        if (show_line)
            p <- p + ggplot2::geom_vline(xintercept = obs_stat, linetype = "dashed", color = "red")
        if (show_caption)
            p <- p + ggplot2::labs(caption = caption) +
                ggplot2::theme(plot.caption = ggplot2::element_text(color = "red", hjust = 0))
        if (show_label) {
            yMax <- ggplot2::layer_scales(p)$y$range$range[2]
            p <- p + ggplot2::annotate("text", x = obs_stat, y = yMax, vjust = "top", label = obs_label, color = "red")
        }
    }
    p
}
