#' Null distribution helper utilities
#'
#' Shared functions for seed handling, computing p-values from permutation
#' (or draw-based) null distributions, and plotting those distributions.
#' Used internally by all hypothesis test analyses.
#'
#' @name perm_utils
#' @import ggplot2
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

#' Format a simulation p-value for a results table
#'
#' A simulation p-value of exactly zero means none of the \code{reps}
#' simulated statistics was as extreme as the observed one, so all the
#' simulation can say is that the p-value is below \code{1/reps}.  Jamovi
#' would otherwise print such a value according to the number of decimal
#' places displayed (e.g. \code{0.0000} for 100 reps, or \code{< .001}
#' for 10 reps), which misstates the resolution of the simulation.
#'
#' @param p Numeric p-value from \code{\link{compute_null_pval}}.
#' @param reps Number of simulated samples.
#' @return \code{p} unchanged when positive; otherwise a string such as
#'   \code{"< .01"} (100 reps) or \code{"< .001"} (1000 reps).
#' @keywords internal
format_sim_pval <- function(p, reps) {
    if (is.null(p) || is.na(p) || p > 0)
        return(p)
    thr <- signif(1 / reps, 2)
    txt <- sub("0+$", "", sprintf("%.10f", thr))
    paste0("< ", sub("^0", "", txt))
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
#' @param domain Optional length-2 bounds of the statistic (e.g.
#'   \code{c(0, Inf)} for chi-square or F, \code{c(0, 1)} for a
#'   proportion).  Bars and dot columns are trimmed so nothing is drawn
#'   at impossible values.
#' @param show_counts Logical; if \code{TRUE}, print the count above
#'   each bar or dot stack.  Default \code{FALSE}.
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
                           show_tail = TRUE,
                           domain = NULL,
                           show_counts = FALSE) {
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

    fill <- extreme & show_tail

    # Sparse discrete statistics (e.g. chi-square from a 2x2 table, or a
    # handful of reps) cannot be binned gap-free with equal-width bins:
    # too few achievable values, unevenly spaced.  Draw one bar / dot
    # column per exact value instead -- the space between values is
    # honest, every element sits exactly at its value, and a value equal
    # to the observed statistic lies on the line.
    u <- sort(unique(perms$stat))
    sparse <- length(u) >= 2 && length(u) <= 8

    if (sparse) {
        a <- 0.4 * min(diff(u))
        idx <- match(perms$stat, u)
        dot_x <- u[idx]
        bars <- data.frame(idx = idx, fill = fill) |>
            dplyr::count(idx, fill)
        bars$xmin <- u[bars$idx] - a
        bars$xmax <- u[bars$idx] + a
        if (!is.null(domain)) {
            bars$xmin <- pmax(bars$xmin, domain[1])
            bars$xmax <- pmin(bars$xmax, domain[2])
        }
        bars$mid <- u[bars$idx]
    } else {
        b <- choose_binning(perms$stat, obs_stat, align = "edge", sign = sgn)
        idx <- bin_index(perms$stat, obs_stat, b$bw, b$off, sgn, "edge")
        idx[extreme] <- pmax(idx[extreme], 0L)
        idx[!extreme] <- pmin(idx[!extreme], -1L)
        a <- 0.42 * b$bw
        dot_x <- dot_column_x(perms$stat, idx, b, obs_stat, sgn, "edge")
        # a bin straddling a domain bound would place its column at an
        # impossible value; clamp the center to the bound instead
        if (!is.null(domain))
            dot_x <- pmin(pmax(dot_x, domain[1]), domain[2])

        bars <- data.frame(idx = idx, fill = fill) |>
            dplyr::count(idx, fill)
        xr <- bin_xrange(bars$idx, obs_stat, b$bw, sgn, "edge")
        bars$xmin <- xr$xmin
        bars$xmax <- xr$xmax
        bars$mid <- xr$mid

        # Trim grouped bins to the statistic's domain so no bar implies
        # impossible values (e.g. negative chi-square).  Single-value
        # bins are left alone: each bar starts exactly at its value (the
        # value-at-edge convention), which is already honest.
        single_val <- b$lattice && abs(b$bw - b$res) <= b$res * 1e-9
        if (!is.null(domain) && !single_val) {
            bars$xmin <- pmax(bars$xmin, domain[1])
            bars$xmax <- pmin(bars$xmax, domain[2])
            bars$mid <- (bars$xmin + bars$xmax) / 2
        }
    }

    # Extend the axis to reach the observed value, but never so far that
    # the distribution collapses into a sliver: beyond half the data span
    # past the data, the observed value is indicated with an arrow at the
    # panel edge instead of a line.
    dat_lo <- min(bars$xmin)
    dat_hi <- max(bars$xmax)
    dspan <- dat_hi - dat_lo
    cap_lo <- dat_lo - 0.5 * dspan
    cap_hi <- dat_hi + 0.5 * dspan
    obs_in <- obs_stat >= cap_lo && obs_stat <= cap_hi
    x_lo <- max(min(dat_lo, obs_stat), cap_lo)
    x_hi <- min(max(dat_hi, obs_stat), cap_hi)
    pad <- 0.03 * (x_hi - x_lo)
    xlims <- c(x_lo - pad, x_hi + pad)

    fill_scale <- ggplot2::scale_fill_manual(
        values = c("FALSE" = "black", "TRUE" = "#ff8c8c"), guide = "none")

    if (dotHist == "dotplot") {
        dots <- data.frame(idx = idx, fill = fill) |>
            dplyr::group_by(idx) |>
            dplyr::mutate(y = dplyr::row_number()) |>
            dplyr::ungroup()
        dots$x <- dot_x

        max_stack <- max(dots$y)
        # At least 0.6 above the tallest stack so the top dot (semi-height
        # up to 0.45) is never clipped by the y limit; count labels need
        # a little more headroom
        y_top <- max(max_stack * (if (show_label) 1.3 else 1.06) *
                         (if (show_counts) 1.08 else 1),
                     max_stack + 0.6 + (if (show_counts) 0.8 else 0))

        p <- ggplot2::ggplot(dots) +
            stack_dots(dots, a = a) +
            # floor-sized dots on very tall stacks may poke past y_top
            ggplot2::coord_cartesian(clip = "off") +
            fill_scale +
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
                dplyr::summarize(x = x[1], n = max(y), fill = fill[1],
                                 .groups = "drop")
            p <- p + count_labels(tops$x, tops$n + 0.5, tops$n, tops$fill)
        }
    } else {
        y_top <- max(bars$n) * (if (show_label) 1.2 else 1.04) *
            (if (show_counts) 1.08 else 1)

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

        if (show_counts)
            p <- p + count_labels(bars$mid, bars$n, bars$n, bars$fill)
    }

    if (show_line && obs_in)
        p <- p + ggplot2::geom_vline(xintercept = obs_stat, linetype = "dashed", color = "red")
    W <- diff(xlims)
    if (obs_in) {
        if (show_label) {
            hj <- inward_hjust(obs_stat, xlims)
            # keep the text clear of the dashed line when edge-justified
            lab_x <- obs_stat + (0.5 - hj) * 2 * (0.01 * W)
            p <- p + ggplot2::annotate("text", x = lab_x, y = y_top,
                                       vjust = "top", hjust = hj,
                                       label = obs_label, color = "red")
        }
    } else {
        # Observed value beyond the axis cap: a drawn arrow at the panel
        # edge stands in for the dashed line (drawn geometrically -- text
        # arrows are missing from some graphics-device fonts)
        edge <- if (obs_stat > cap_hi) xlims[2] else xlims[1]
        dir <- if (obs_stat > cap_hi) 1 else -1
        if (show_line)
            p <- p + ggplot2::annotate("segment",
                x = edge - dir * 0.08 * W, xend = edge - dir * 0.005 * W,
                y = y_top * 0.92, yend = y_top * 0.92,
                color = "red", linewidth = 0.6,
                arrow = grid::arrow(length = grid::unit(6, "pt"), type = "closed"))
        if (show_label)
            p <- p + ggplot2::annotate("text",
                x = edge - dir * 0.095 * W, y = y_top, vjust = "top",
                hjust = if (dir > 0) 1 else 0,
                label = obs_label, color = "red")
    }
    if (show_caption)
        p <- p + ggplot2::labs(caption = caption) +
            ggplot2::theme(plot.caption = ggplot2::element_text(color = "red", hjust = 0))

    p
}

#' Count labels above bars or dot stacks
#'
#' @param x,y Label positions (the top of each bar or stack).
#' @param n Counts to print.
#' @param extreme Logical per label; \code{TRUE} labels are drawn in the
#'   tail-shading red so a student can read the p-value numerator
#'   straight off the plot.
#' @return A \code{geom_text} layer.
#' @keywords internal
count_labels <- function(x, y, n, extreme = FALSE) {
    d <- data.frame(x = x, y = y, n = n,
                    col = ifelse(rep_len(extreme, length(x)), "#d9534f", "grey25"))
    ggplot2::geom_text(
        data = d, ggplot2::aes(x = x, y = y, label = n),
        color = d$col, vjust = -0.35, size = 3, inherit.aes = FALSE)
}
