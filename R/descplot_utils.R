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
    # the observations themselves, drawn behind the markers
    rbind(stats, data.frame(group = group, stat = dep, cie = NA, type = 'point'))
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

    plotData <- rbind(meanPlotData, medianPlotData,
                      data.frame(group = rep(c(name1, name2), c(length(col1), length(col2))),
                                 stat = c(col1, col2), cie = NA, type = 'point'))
    plotData$group <- factor(plotData$group, levels = c(name1, name2))
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
        plot <- plot + ggplot2::geom_point(
            data = pts, ggplot2::aes(x = xpos, y = stat), inherit.aes = FALSE,
            color = "grey55", alpha = 0.45, size = 1.6)
    }
    plot <- plot +
        ggplot2::geom_point(ggplot2::aes(x = xpos, y = stat, shape = type),
                            color = theme$color[1], fill = theme$fill[1],
                            size = 3, position = pd) +
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

    plot
}
