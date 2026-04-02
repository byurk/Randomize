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
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
plot_null_dist <- function(perms, obs_stat, direction,
                           dotHist = c("dotplot", "histogram"),
                           xlab = "statistic",
                           obs_label = "Observed\nStatistic") {
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

    # Determine if observed stat is far from the null distribution
    dist_min <- min(perms$stat)
    dist_max <- max(perms$stat)
    dist_span <- dist_max - dist_min
    if (dist_span == 0) dist_span <- abs(obs_stat) * 0.01 + 0.001
    obs_off_right <- obs_stat > dist_max + 1.5 * dist_span
    obs_off_left  <- obs_stat < dist_min - 1.5 * dist_span
    obs_off_screen <- obs_off_right || obs_off_left

    if (dotHist == "dotplot") {
        ndist <- dplyr::n_distinct(perms$stat)
        bw <- dist_span / 30

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
        dist_range <- bigx - littlex + 2 * bw

        if (obs_off_screen) {
            # Keep x limits tight to the distribution; no misleading line
            x_lo <- littlex - bw
            x_hi <- bigx + bw
            edge_margin <- dist_range * 0.15
            if (obs_off_right) {
                x_hi <- x_hi + edge_margin
                label_x <- bigx + bw
                edge_label <- paste0(obs_label, "\n= ", round(obs_stat, 2), " \u2192")
                label_hjust <- 1
            } else {
                x_lo <- x_lo - edge_margin
                label_x <- littlex - bw
                edge_label <- paste0("\u2190 ", obs_label, "\n= ", round(obs_stat, 2))
                label_hjust <- 0
            }
        } else {
            # obs_stat is near the distribution — include it in range
            x_lo <- min(littlex - bw, obs_stat - bw)
            x_hi <- max(bigx + bw, obs_stat + bw)
        }

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
            ggplot2::labs(caption = caption) +
            ggplot2::theme(
                text = ggplot2::element_text(size = 14),
                plot.caption = ggplot2::element_text(color = "red", hjust = 0)
            )

        p <- p + ggplot2::coord_equal(ratio = bw * 2 / 3)

        if (obs_off_screen) {
            # Annotate at the edge with actual value — no dashed line
            p <- p + ggplot2::annotate("text", x = label_x, y = lab_ht,
                                       label = edge_label, color = "red",
                                       hjust = label_hjust)
        } else {
            p <- p +
                ggplot2::geom_vline(xintercept = obs_stat, linetype = "dashed", color = "red") +
                ggplot2::annotate("text", x = obs_stat, y = lab_ht,
                                  label = obs_label, color = "red")
        }
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
            ggplot2::labs(caption = caption) +
            ggplot2::theme(
                text = ggplot2::element_text(size = 14),
                plot.caption = ggplot2::element_text(color = "red", hjust = 0)
            )

        if (obs_off_screen) {
            # Use coord_cartesian to avoid removing data (no warnings)
            hist_lo <- dist_min - 0.5 * dist_span
            hist_hi <- dist_max + 0.5 * dist_span
            edge_margin <- dist_span * 0.15
            if (obs_off_right) {
                hist_hi <- hist_hi + edge_margin
                label_x <- dist_max + 0.5 * dist_span
                edge_label <- paste0(obs_label, "\n= ", round(obs_stat, 2), " \u2192")
                label_hjust <- 1
            } else {
                hist_lo <- hist_lo - edge_margin
                label_x <- dist_min - 0.5 * dist_span
                edge_label <- paste0("\u2190 ", obs_label, "\n= ", round(obs_stat, 2))
                label_hjust <- 0
            }
            p <- p + ggplot2::coord_cartesian(xlim = c(hist_lo, hist_hi))
            yMax <- ggplot2::layer_scales(p)$y$range$range[2]
            p <- p + ggplot2::annotate("text", x = label_x, y = yMax, vjust = "top",
                                       label = edge_label, color = "red",
                                       hjust = label_hjust)
        } else {
            yMax <- ggplot2::layer_scales(p)$y$range$range[2]
            p <- p +
                ggplot2::geom_vline(xintercept = obs_stat, linetype = "dashed", color = "red") +
                ggplot2::annotate("text", x = obs_stat, y = yMax, vjust = "top",
                                  label = obs_label, color = "red")
        }
    }
    p
}
