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

#' Round a statistic to 10 significant digits
#'
#' Permutation statistics computed by different arithmetic paths (a
#' difference of two proportions, say) can land one ulp apart from a
#' mathematically equal observed statistic.  Every comparison against
#' the observed value -- the p-value, the tail shading, the tie column --
#' goes through this so such values are ties, as they should be.
#' @param x Numeric vector.
#' @keywords internal
snap_stat <- function(x) signif(x, 10)

#' Replace floating-point dust with exact zero
#'
#' @param x Numeric vector.
#' @param tol Values with absolute value below this become 0.
#' @keywords internal
zap_tiny <- function(x, tol = 1e-10) {
    x[!is.na(x) & abs(x) < tol] <- 0
    x
}

#' Reuse simulations across re-runs that do not change them
#'
#' Jamovi re-runs an analysis from scratch on every option change, so
#' toggling a display-only option (plot type, count labels, descriptives)
#' used to draw a fresh set of simulations -- the plot changed shape for
#' no statistical reason, which confuses students.  The simulations are
#' now kept in the results table's state together with a key describing
#' everything that determines them (the data actually used, the rep
#' count, the seed settings, and any analysis-specific inputs such as
#' the null value).  They are recomputed only when that key changes.
#' The table's own \code{clearWith} list wipes the state for the same
#' options, and the key guards against edits to the data.
#'
#' @param holder A results element (normally the simulation results
#'   table) whose \code{state} carries the cache.
#' @param key A list of everything the simulations depend on; compared
#'   with \code{identical()}.
#' @param compute A function of no arguments returning the simulations
#'   (an \pkg{infer} tibble or a data frame with a \code{stat} column).
#' @param slot Cache slot, for analyses that simulate several things
#'   (one per variable pair).
#' @return A plain data frame with a \code{stat} column.
#' @keywords internal
cached_sims <- function(holder, key, compute, slot = "sims") {
    st <- holder$state
    entry <- if (is.list(st)) st[[slot]] else NULL
    if (!is.null(entry) && identical(entry$key, key) &&
        is.data.frame(entry$sims) && "stat" %in% names(entry$sims) &&
        is.numeric(entry$sims$stat) && nrow(entry$sims) > 0)
        return(entry$sims)
    sims <- strip_infer(compute())
    if (!is.list(st)) st <- list()
    st[[slot]] <- list(key = key, sims = sims)
    holder$setState(st)
    sims
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
    # Snap away floating-point dust first: a simulated difference in
    # proportions that mathematically equals the observed one can differ
    # from it by one ulp (15/25 - 11/25 vs 0.16) and would then be
    # counted on the wrong side of >= / <=
    perms$stat <- snap_stat(perms$stat)
    obs_stat <- snap_stat(obs_stat)
    # infer warns when p = 0; the tables report that case as "< 1/reps"
    # (format_sim_pval), so the warning is noise for R users
    withCallingHandlers(
        perms |>
            infer::get_p_value(obs_stat = obs_stat, direction = direction) |>
            dplyr::pull(),
        warning = function(w) {
            if (grepl("p-value of 0", conditionMessage(w))) invokeRestart("muffleWarning")
        })
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
    txt <- sub("\\.$", "", txt)      # reps = 1: "< 1", not "< 1."
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
#' @param plot_width Logical plot width in pixels, used to size the
#'   count labels (see \code{\link{count_label_size}}).
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
                           show_counts = FALSE,
                           plot_width = 400) {
    dotHist <- match.arg(dotHist)
    # same snapping as the p-value, so the shading and the tie column
    # agree with the number in the table
    perms$stat <- snap_stat(perms$stat)
    obs_stat <- snap_stat(obs_stat)

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
    # to the observed statistic lies on the line.  (Values that all but
    # coincide share a column; see sparse_groups.)
    sparse <- is_sparse_values(perms$stat)
    few <- dotHist == "dotplot" && !sparse && is_few_reps(perms$stat)

    if (sparse) {
        g <- sparse_groups(perms$stat, extreme)
        u <- g$centers
        a <- 0.4 * min(diff(u))
        idx <- g$idx
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
    } else if (few) {
        # A handful of simulations (a step-by-step classroom demo):
        # every dot sits at its own value rather than at a bin centre,
        # so the picture shows exactly where each simulation landed
        # relative to the line.  Exact ties stack; near-ties overlap.
        span <- diff(range(c(perms$stat, obs_stat)))
        if (span <= 0) span <- max(abs(perms$stat[1]), 1)
        a <- 0.02 * span
        # values within a dot diameter of a cluster's first value share a
        # stack (at their mean, never across the line) so neither dots nor
        # count labels pile up; no dot moves more than a diameter
        g <- cluster_near(perms$stat, extreme, tol = 2 * a)
        u <- g$centers
        idx <- g$idx
        dot_x <- u[idx]
        bars <- data.frame(idx = idx, fill = fill) |>
            dplyr::count(idx, fill)
        bars$mid <- u[bars$idx]
        bars$xmin <- bars$mid - a
        bars$xmax <- bars$mid + a
        if (!is.null(domain)) {
            bars$xmin <- pmax(bars$xmin, domain[1])
            bars$xmax <- pmin(bars$xmax, domain[2])
        }
    } else {
        b <- choose_binning(perms$stat, obs_stat, align = "edge", sign = sgn)
        single_val <- b$lattice && abs(b$bw - b$res) <= b$res * 1e-9
        a <- 0.42 * b$bw
        if (single_val) {
            # One achievable value per column: index each value by its
            # lattice offset from the observed value and draw the bar (or
            # dot stack) centred on the value itself, like the dot stacks
            # always were.  A column holds one value, so it is one colour
            # by construction, and the tie column straddles the line
            # (shaded, since ties count toward the p-value) whichever
            # tail is extreme -- an edge-anchored bar would have to jump
            # to the other side of the line with the direction, dragging
            # its neighbour with it.
            idx <- as.integer(round((perms$stat - obs_stat) / b$res))
            vals <- stats::ave(perms$stat, idx, FUN = function(v) v[1])
            dot_x <- vals
            bars <- data.frame(idx = idx, fill = fill, v = vals) |>
                dplyr::count(idx, fill, v)
            bars$mid <- bars$v
            bars$xmin <- bars$v - b$bw / 2
            bars$xmax <- bars$v + b$bw / 2
            bars$v <- NULL
        } else {
            idx <- bin_index(perms$stat, obs_stat, b$bw, b$off, sgn, "edge")
            idx[extreme] <- pmax(idx[extreme], 0L)
            idx[!extreme] <- pmin(idx[!extreme], -1L)
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
        }

        # Trim bars to the statistic's domain so no bar implies
        # impossible values (e.g. a proportion above 1)
        if (!is.null(domain)) {
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
    # With a handful of simulations (a step-by-step classroom demo)
    # there is no distribution to protect from collapsing, so always
    # bring the observed value onto the axis
    if (nrow(bars) <= 2 || nrow(perms) <= 10) { cap_lo <- -Inf; cap_hi <- Inf }
    obs_in <- obs_stat >= cap_lo && obs_stat <= cap_hi
    x_lo <- max(min(dat_lo, obs_stat), cap_lo)
    x_hi <- min(max(dat_hi, obs_stat), cap_hi)
    pad <- 0.03 * (x_hi - x_lo)
    # a dot drawn at a domain bound (few reps) must not be cut by the panel
    if (few) pad <- max(pad, a)
    xlims <- c(x_lo - pad, x_hi + pad)

    fill_scale <- ggplot2::scale_fill_manual(
        values = c("FALSE" = "black", "TRUE" = "#ff8c8c"), guide = "none")
    ts <- text_scale(plot_width)

    if (dotHist == "dotplot") {
        # Dot k of a stack is centred at k - 1/2 so the bottom dot rests on
        # the axis and the stack's top edge is at its count
        dots <- data.frame(idx = idx, fill = fill) |>
            dplyr::group_by(idx) |>
            dplyr::mutate(y = dplyr::row_number() - 0.5) |>
            dplyr::ungroup()
        dots$x <- dot_x

        max_stack <- max(dots$y) + 0.5
        # Headroom above the tallest stack for the label / count labels.
        # The axis never shows fewer than 8 counts: with a handful of dots
        # a unit would otherwise be so tall that the 4 mm radius ceiling
        # leaves the bottom dot hovering above the baseline.
        cnt_room <- if (!show_counts) 1 else if (counts_vertical(bars$n, plot_width)) 1.2 else 1.08
        y_top <- max(max_stack * (if (show_label) 1.3 else 1.06) * cnt_room,
                     max_stack + 0.6 + (if (show_counts) 0.8 else 0),
                     8)

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
            ggplot2::theme(text = ggplot2::element_text(size = 14 * ts))

        cnt <- NULL
        if (show_counts) {
            tops <- dots |>
                dplyr::group_by(idx) |>
                dplyr::summarize(x = x[1], n = max(y) + 0.5, fill = fill[1],
                                 .groups = "drop")
            cnt <- count_labels(tops$x, tops$n, tops$n, tops$fill, plot_width,
                                line_x = if (show_line) obs_stat else NULL)
        }
    } else {
        cnt_room <- if (!show_counts) 1 else if (counts_vertical(bars$n, plot_width)) 1.2 else 1.08
        y_top <- max(bars$n) * (if (show_label) 1.2 else 1.04) * cnt_room

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
            ggplot2::theme(text = ggplot2::element_text(size = 14 * ts))

        cnt <- if (show_counts) count_labels(bars$mid, bars$n, bars$n, bars$fill, plot_width,
                                             line_x = if (show_line) obs_stat else NULL) else NULL
    }

    if (show_line && obs_in)
        p <- p + ggplot2::geom_vline(xintercept = obs_stat, linetype = "dashed", color = "red")
    W <- diff(xlims)
    if (obs_in) {
        if (show_label) {
            hj <- inward_hjust(obs_stat, xlims)
            # keep the text clear of the dashed line when edge-justified
            lab_x <- obs_stat + (0.5 - hj) * 2 * (0.01 * W)
            p <- p + ggplot2::annotate("label", x = lab_x, y = y_top,
                                       vjust = "top", hjust = hj, size = 3.88 * ts,
                                       label = obs_label, color = "red",
                                       fill = "white", label.size = 0,
                                       label.padding = grid::unit(0.12, "lines"))
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
                hjust = if (dir > 0) 1 else 0, size = 3.88 * ts,
                label = obs_label, color = "red")
    }
    if (show_caption)
        p <- p + ggplot2::labs(caption = caption) +
            ggplot2::theme(plot.caption = ggplot2::element_text(color = "red", hjust = 0))
    # count labels last, so their halos sit over the dashed line
    if (!is.null(cnt))
        p <- p + cnt

    p
}

#' Count labels above bars or dot stacks
#'
#' Every bar / stack gets its count.  The label is as large as the axis
#' tick labels (\code{size = 3.9}, about 11 pt) whenever the widest count
#' fits inside its bar at Jamovi's 400-pixel plot width, shrinks with the
#' bar width down to \code{size = 3}, and below that stands up vertically
#' (still size 3 or larger) so neighbours cannot run into each other.
#' Size 3 is a hard floor: Jamovi renders with ragg at 72 ppi, which
#' silently drops smaller text.  If even that is too small to read, the
#' user can switch the counts off.
#'
#' @param x,y Label positions (the top of each bar or stack).
#' @param n Counts to print.
#' @param extreme Logical per label; \code{TRUE} labels are drawn in the
#'   tail-shading red so a student can read the p-value numerator
#'   straight off the plot.
#' @return A \code{geom_text} layer.
#' @keywords internal
count_labels <- function(x, y, n, extreme = FALSE, plot_width = 400, line_x = NULL) {
    d <- data.frame(x = x, y = y, n = n,
                    col = ifelse(rep_len(extreme, length(x)), "#d9534f", "grey25"))
    sz <- count_label_size(n, plot_width)
    common <- list(mapping = ggplot2::aes(x = x, y = y, label = n),
                   size = sz$size, inherit.aes = FALSE)
    if (sz$vertical)
        common <- c(common, list(angle = 90, hjust = -0.15, vjust = 0.5))
    else
        common <- c(common, list(vjust = -0.35))
    layers <- list(do.call(ggplot2::geom_text, c(common, list(data = d, color = d$col))))
    # A white halo (a borderless label box with invisible text) under
    # the labels that a vertical line would otherwise strike through:
    # those within half a column of the observed value / CI limits.
    # Only those -- drawn under every label, the box would erase bits of
    # neighbouring dots in the valleys of a dense dotplot.
    if (!is.null(line_x) && length(line_x) > 0) {
        ux <- sort(unique(x))
        half <- if (length(ux) > 1) 0.5 * min(diff(ux)) else Inf
        near <- vapply(x, function(xi) any(abs(xi - line_x) <= half * (1 + 1e-9)), logical(1))
        if (any(near)) {
            halo <- c(common, list(data = d[near, ], fill = "white", colour = "white",
                                   label.size = 0,
                                   label.padding = grid::unit(0.08, "lines"),
                                   label.r = grid::unit(0, "lines")))
            layers <- c(list(do.call(ggplot2::geom_label, halo)), layers)
        }
    }
    layers
}

#' Is this a step-by-step handful of simulations?
#'
#' With this few replicates a dotplot places every dot at its own value
#' (no binning), so the picture shows exactly where each simulation
#' landed relative to the observed value.  The threshold is shared with
#' the plot invariants.
#' @param stats Numeric vector of simulated statistics.
#' @keywords internal
is_few_reps <- function(stats) length(stats) <= 25

#' Axis label naming the groups (and outcome level) being compared
#'
#' @param what "means" or "proportions".
#' @param levels The two group names, in the order of the subtraction.
#' @param outcome For proportions, the outcome level whose proportion is
#'   compared.
#' @keywords internal
diff_label <- function(what, levels = NULL, outcome = NULL) {
    if (is.null(levels) || length(levels) < 2)
        return(sprintf("difference in %s (group 1 \u2212 group 2)", what))
    lv <- shorten_label(levels)
    if (is.null(outcome)) {
        one <- sprintf("difference in %s (%s \u2212 %s)", what, lv[1], lv[2])
        # long group names: second line, or the title runs off the panel
        if (nchar(one) <= 36) one
        else sprintf("difference in %s\n(%s \u2212 %s)", what, lv[1], lv[2])
    } else {
        sprintf("difference in proportion of %s\n(%s \u2212 %s)", shorten_label(outcome), lv[1], lv[2])
    }
}

#' @rdname diff_label
#' @keywords internal
shorten_label <- function(x, n = 18) {
    x <- as.character(x)
    ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "\u2026"), x)
}

#' Value from a stored plot state, or a default for states written by an
#' older build that lacks it
#' @keywords internal
state_or <- function(x, default) if (is.null(x)) default else x

#' Size and orientation of count labels
#'
#' Fits the widest count into one bar's share of Jamovi's plot width
#' (400 px wide, roughly 340 px of panel): horizontal when that allows at
#' least size 3, capped at the tick-label size (3.9); otherwise vertical,
#' sized by the bar spacing with the same floor and cap.
#'
#' @param n The counts to be printed.
#' @param plot_width Logical width of the plot in pixels (Jamovi's
#'   default is 400; a plot the user has resized in the results panel
#'   reports its new width, so labels grow with it).
#' @return A list with \code{size} (ggplot text size, mm) and
#'   \code{vertical}.
#' @keywords internal
count_label_size <- function(n, plot_width = 400) {
    n_bars <- length(n)
    chars <- max(nchar(as.character(n)))
    cap <- 3.9 * text_scale(plot_width)          # the (scaled) tick-label size
    per_bar <- 0.95 * 0.85 * plot_width / n_bars  # px available per label
    pt_per_size <- 2.845                          # ggplot size (mm) -> pt (= px at 72 ppi)
    size_h <- per_bar / (chars * 0.6 * pt_per_size)   # digits are ~0.6 em wide
    if (size_h >= 3)
        return(list(size = min(cap, size_h), vertical = FALSE))
    size_v <- per_bar / (0.8 * pt_per_size)      # vertical: bar spacing vs glyph height
    list(size = max(3, min(cap, size_v)), vertical = TRUE)
}

#' Text magnification for an enlarged plot
#'
#' When a plot is dragged larger in Jamovi's results panel the bars grow
#' but text drawn at a fixed point size does not, so labels look smaller
#' relative to the plot.  All text in the distribution plots is therefore
#' scaled by the plot width relative to the 400-pixel default (capped at
#' 2x), so enlarging a plot magnifies it as a whole.
#'
#' @param plot_width Logical plot width in pixels.
#' @keywords internal
text_scale <- function(plot_width = 400) {
    min(2, max(1, plot_width / 400))
}

#' Should count labels be drawn vertically?
#' @param n The counts to be printed.
#' @keywords internal
counts_vertical <- function(n, plot_width = 400) count_label_size(n, plot_width)$vertical

#' Permutation null distributions computed without formulas
#'
#' Jamovi 2.7 runs module code under a formula sandbox: its \code{as.formula}
#' rejects any formula that calls a function outside a short allowlist.
#' \code{infer::specify()} builds \code{response_variable(x) ~
#' explanatory_variable(x)} internally (to set theoretical-distribution
#' parameters via \code{t.test()} / \code{lm()}) whenever the response is
#' numeric and the explanatory variable is a factor, so the two-sample and
#' multi-sample mean analyses cannot go through \pkg{infer} there.  These
#' helpers draw the same permutation / bootstrap distributions directly:
#' labels are shuffled (or rows resampled) and the statistic is computed
#' from group means.  The result mimics \code{infer::calculate()}: a data
#' frame with \code{replicate} and \code{stat} columns.
#'
#' @param dep Numeric response.
#' @param group Factor of group labels (same length as \code{dep}).
#' @param levels Length-2 character vector: the statistic is
#'   \code{mean(levels[1]) - mean(levels[2])}.
#' @param reps Number of replicates.
#' @name resample_means
#' @keywords internal
NULL

#' @rdname resample_means
#' @keywords internal
permute_diff_means <- function(dep, group, levels, reps) {
    group <- as.character(group)
    stat <- vapply(seq_len(reps), function(i) {
        g <- sample(group)
        mean(dep[g == levels[1]]) - mean(dep[g == levels[2]])
    }, numeric(1))
    data.frame(replicate = seq_len(reps), stat = stat)
}

#' One-way ANOVA F statistic from group means (no model formula)
#'
#' @inheritParams resample_means
#' @return The F statistic \eqn{(SSB/(k-1)) / (SSW/(N-k))}, identical to
#'   \code{anova(lm(dep ~ group))$F[1]}.
#' @keywords internal
f_stat <- function(dep, group) {
    group <- as.character(group)
    n <- length(dep)
    gm <- tapply(dep, group, mean)
    ng <- tapply(dep, group, length)
    k <- length(gm)
    ssb <- sum(ng * (gm - mean(dep))^2)
    ssw <- sum((dep - gm[group])^2)
    (ssb / (k - 1)) / (ssw / (n - k))
}

#' @rdname resample_means
#' @keywords internal
permute_F <- function(dep, group, reps) {
    group <- as.character(group)
    stat <- vapply(seq_len(reps), function(i) f_stat(dep, sample(group)), numeric(1))
    data.frame(replicate = seq_len(reps), stat = stat)
}
