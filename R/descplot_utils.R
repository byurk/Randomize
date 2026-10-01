#' Descriptive plot helper utilities
#'
#' Shared functions for building mean/median descriptive plots
#' used by twomeanCI, twomeanhtest, pairedmeanCI, pairedmeanhtest,
#' and multimeanhtest.
#'
#' @import ggplot2
NULL

#' Build plot data for grouped mean/median descriptive plots
#'
#' Creates a data frame with mean and median statistics per group,
#' suitable for use with plot_desc_stats().
#'
#' @param dep numeric vector of dependent variable values
#' @param group factor vector of group labels
#' @return data frame with columns: group, stat, cie, type
build_desc_plot_data <- function(dep, group) {
    means <- aggregate(dep, by = list(group),
                       function(x) tryNaN(mean(x)), simplify = FALSE)
    medians <- aggregate(dep, by = list(group),
                         function(x) tryNaN(median(x)), simplify = FALSE)

    meanPlotData <- data.frame(group = means$Group.1)
    meanPlotData <- cbind(meanPlotData, stat = unlist(means$x))
    meanPlotData <- cbind(meanPlotData, cie = NA)
    meanPlotData <- cbind(meanPlotData, type = 'mean')

    medianPlotData <- data.frame(group = medians$Group.1)
    medianPlotData <- cbind(medianPlotData, stat = unlist(medians$x))
    medianPlotData <- cbind(medianPlotData, cie = NA)
    medianPlotData <- cbind(medianPlotData, type = 'median')

    stats <- rbind(meanPlotData, medianPlotData)
    # the observations themselves, drawn behind the markers -- at most
    # DESC_MAX_POINTS per group (see thin_points), or a large data set's
    # plot state would push the analysis past Jamovi's 4 MB message limit
    keep <- unlist(lapply(split(seq_along(dep), group), function(i) thin_points(dep[i], i)))
    out <- rbind(stats, data.frame(group = group[keep], stat = dep[keep], cie = NA, type = 'point'))
    attr(out, "n_total") <- length(dep); attr(out, "n_shown") <- length(keep)
    out
}

# Observations kept for a descriptive plot: all of them up to this many
# per group, beyond that evenly spaced order statistics (the shape of
# the distribution survives; a 100 000-row data set would otherwise
# store ~3 MB of points in the plot state)
DESC_MAX_POINTS <- 2000

#' @rdname build_desc_plot_data
#' @keywords internal
thin_points <- function(x, idx = seq_along(x), max_n = DESC_MAX_POINTS) {
    if (length(x) <= max_n) return(idx)
    o <- order(x)
    idx[o[round(seq(1, length(x), length.out = max_n))]]
}

#' Build plot data for paired mean/median descriptive plots
#'
#' Like build_desc_plot_data() but for paired data where we have
#' two named columns rather than a dep~group structure.
#'
#' @param col1 numeric vector for measure 1
#' @param col2 numeric vector for measure 2
#' @param name1 label for measure 1
#' @param name2 label for measure 2
#' @return data frame with columns: group, stat, cie, type
build_paired_desc_plot_data <- function(col1, col2, name1, name2) {
    means <- c(tryNaN(mean(col1)), tryNaN(mean(col2)))
    medians <- c(tryNaN(median(col1)), tryNaN(median(col2)))

    meanPlotData <- data.frame(group = c(name1, name2))
    meanPlotData <- cbind(meanPlotData, stat = means)
    meanPlotData <- cbind(meanPlotData, cie = NA)
    meanPlotData <- cbind(meanPlotData, type = 'mean')

    medianPlotData <- data.frame(group = c(name1, name2))
    medianPlotData <- cbind(medianPlotData, stat = medians)
    medianPlotData <- cbind(medianPlotData, cie = NA)
    medianPlotData <- cbind(medianPlotData, type = 'median')

    # pairs are thinned together (by the rank of measure 1) so every
    # kept observation still has its partner
    keep <- thin_points(col1)
    c1 <- col1[keep]; c2 <- col2[keep]
    plotData <- rbind(meanPlotData, medianPlotData,
                      data.frame(group = rep(c(name1, name2), c(length(c1), length(c2))),
                                 stat = c(c1, c2), cie = NA, type = 'point'))
    # which observation of measure 1 goes with which of measure 2, so the
    # plot can join each pair with a line
    plotData$pair <- c(rep(NA, 4), seq_along(c1), seq_along(c2))
    plotData$group <- factor(plotData$group, levels = c(name1, name2))
    attr(plotData, "n_total") <- length(col1); attr(plotData, "n_shown") <- length(c1)
    plotData
}

#' Render a grouped mean/median descriptive plot
#'
#' @param plotData data frame from build_desc_plot_data() or build_paired_desc_plot_data()
#' @param xlab x-axis label (group variable name)
#' @param ylab y-axis label (dependent variable name or NULL)
#' @param ggtheme the Jamovi ggtheme
#' @param theme the Jamovi theme (color/fill)
#' @return ggplot object
plot_desc_stats <- function(plotData, xlab, ylab = NULL,
                            ggtheme, theme) {
    pd <- ggplot2::position_dodge(0.2)
    # the observations (type 'point'; absent from states saved by an
    # earlier build) are drawn first, small and translucent, spread
    # sideways by a fixed low-discrepancy sequence rather than random
    # jitter: the picture is identical on every re-render and drawing it
    # never touches the RNG
    pts <- plotData[plotData$type == 'point', , drop = FALSE]
    stats <- plotData[plotData$type != 'point', , drop = FALSE]
    stats$type <- factor(stats$type, levels = c('mean', 'median'))
    lv <- levels(factor(plotData$group))
    stats$xpos <- as.integer(factor(stats$group, levels = lv))

    plot <- ggplot2::ggplot(data = stats, ggplot2::aes(x = xpos, y = stat, shape = type))
    if (nrow(pts) > 0) {
        g <- as.integer(factor(pts$group, levels = lv))
        k <- stats::ave(seq_along(g), g, FUN = seq_along)
        pts$xpos <- g + ((k * 0.6180339887) %% 1 - 0.5) * 0.36
        # paired data: a thin line joins the two members of each pair
        # (both get the same sideways offset, so the lines run cleanly
        # from one measure to the other; with hundreds of pairs they
        # fuse into a band whose tilt still shows the typical change)
        if ("pair" %in% names(pts) && any(!is.na(pts$pair))) {
            a <- pts[g == 1L & !is.na(pts$pair), ]
            b <- pts[g == 2L & !is.na(pts$pair), ]
            seg <- merge(a[, c("pair", "xpos", "stat")], b[, c("pair", "xpos", "stat")], by = "pair", suffixes = c("", "_end"))
            if (nrow(seg) > 0)
                plot <- plot + ggplot2::geom_segment(
                    data = seg, ggplot2::aes(x = xpos, y = stat, xend = xpos_end, yend = stat_end),
                    inherit.aes = FALSE, color = "grey70", alpha = 0.5, linewidth = 0.3)
        }
        plot <- plot + ggplot2::geom_point(
            data = pts, ggplot2::aes(x = xpos, y = stat), inherit.aes = FALSE,
            color = "grey55", alpha = 0.45, size = 1.6)
    }
    n_total <- attr(plotData, "n_total"); n_shown <- attr(plotData, "n_shown")
    thinned <- nrow(pts) > 0 && !is.null(n_total) && !is.null(n_shown) && n_shown < n_total
    # the markers are filled (the scatterplot's blue) with a dark outline
    # so they stand off the grey observations behind them
    plot <- plot +
        ggplot2::geom_point(ggplot2::aes(x = xpos, y = stat, shape = type),
                            color = theme$color[1], fill = theme$color[2],
                            size = 3.2, stroke = 0.8, position = pd) +
        ggplot2::scale_x_continuous(breaks = seq_along(lv), labels = lv,
                                    limits = c(0.5, length(lv) + 0.5)) +
        ggplot2::labs(x = xlab, y = ylab) +
        ggplot2::scale_shape_manual(
            name = '',
            values = c(mean = 21, median = 22),
            labels = c(mean = 'Mean', median = 'Median')
        ) +
        ggtheme +
        ggplot2::theme(
            plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 5.5 * 1.2)),
            plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
        )

    # (after ggtheme, which would otherwise right-align the caption and
    # push a long line off the 400-px image)
    if (thinned)
        plot <- plot + ggplot2::labs(caption = sprintf("showing %s of %s observations\n(evenly spaced by value)",
                                                       formatC(n_shown, format = "d", big.mark = ","),
                                                       formatC(n_total, format = "d", big.mark = ","))) +
            ggplot2::theme(plot.caption = ggplot2::element_text(color = "grey40", hjust = 0, size = 9))

    plot
}
