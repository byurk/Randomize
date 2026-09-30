
# This file is a generated template, your changes will not be overwritten

modelBasedClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "modelBasedClass",
    inherit = modelBasedBase,
    private = list(
        #### Init + run functions ----
        .init = function() {

            private$.initAreaTable()
            private$.initMultTable()

        },
        .run = function() {

            if (self$options$areaBool) {

                resultsArea <- private$.computeArea()

                private$.populateAreaTable(resultsArea)
                private$.preparePlot()

            }

            if (self$options$CIBool && private$.symmetric()) {

                resultsMult <- private$.computeMult()

                private$.populateMultTable(resultsMult)

            }

        },

        #### Helpers ----

        # Normal and t are symmetric about zero: they support left / both
        # tails and CI multipliers.  Chi-square and F are right-tail only.
        .symmetric = function() {
            self$options$distro %in% c("ndistro", "tdistro")
        },
        .distName = function() {
            switch(self$options$distro,
                ndistro = "standard normal distribution",
                tdistro = jmvcore::format("t-distribution with df = {df}", df = self$options$dF),
                chisq   = jmvcore::format("chi-square distribution with df = {df}", df = self$options$dF),
                fdistro = jmvcore::format("F-distribution with df₁ = {df1}, df₂ = {df2}",
                                          df1 = self$options$dF, df2 = self$options$dF2))
        },
        # Effective tail: chi-square and F always use the right tail
        .tail = function() {
            if (private$.symmetric()) self$options$tail else "right"
        },

        #### Compute results ----
        .computeArea = function() {

            distro <- self$options$distro
            dF <- self$options$dF
            dF2 <- self$options$dF2
            tail <- private$.tail()
            obs_stat <- self$options$obsStat

            pfun <- switch(distro,
                ndistro = function(q, lower) stats::pnorm(q, lower.tail = lower),
                tdistro = function(q, lower) stats::pt(q, dF, lower.tail = lower),
                chisq   = function(q, lower) stats::pchisq(q, dF, lower.tail = lower),
                fdistro = function(q, lower) stats::pf(q, dF, dF2, lower.tail = lower))

            if (tail == "left") {
                pval <- pfun(obs_stat, TRUE)
            } else if (tail == "right") {
                pval <- pfun(obs_stat, FALSE)
            } else {
                pval <- pfun(obs_stat, TRUE)
                pval <- 2 * min(pval, 1 - pval)
            }

            list(obsStat = obs_stat, df = dF, df1 = dF, df2 = dF2, area = pval)
        },
        .computeMult = function() {

            distro <- self$options$distro
            dF <- self$options$dF
            confLevel <- self$options$confLevel

            if (distro == "ndistro")
                dF <- Inf

            tcrit <- stats::qt(1 - (1 - confLevel/100)/2, dF)

            list(confLev = confLevel, df = self$options$dF, critVal = tcrit)
        },

        #### Init tables functions ----
        .initAreaTable = function() {

            areatable <- self$results$get('areaTable')

            tail <- private$.tail()
            region <- switch(tail,
                right = "Right-tail area (values ≥ observed) under the {d}",
                left  = "Left-tail area (values ≤ observed) under the {d}",
                both  = "Two-tail area (both tails beyond ±|observed|) under the {d}")

            areatable$setNote('area', jmvcore::format(region, d = private$.distName()))

        },
        .initMultTable = function() {

            multtable <- self$results$get('multTable')

            if (private$.symmetric()) {
                multtable$setNote(
                    'mult',
                    jmvcore::format(
                        'Multiplier for {cL}% CI based on the {d}',
                        cL = self$options$confLevel,
                        d = private$.distName()
                    )
                )
            } else {
                multtable$setNote(
                    'mult',
                    'CI multipliers are defined for the standard normal and t distributions only'
                )
            }

        },

        #### Populate tables functions ----
        .populateAreaTable = function(results) {

            table <- self$results$get('areaTable')
            table$deleteRows()

            table$addRow(rowKey=1, values=results)

        },
        .populateMultTable = function(results) {

            table <- self$results$get('multTable')
            table$deleteRows()

            table$addRow(rowKey=1, values=results)

        },

        #### Plot functions ----
        .preparePlot = function() {

            areaPlot <- self$results$Plot

            areaPlot$setState(list(
                distro = self$options$distro,
                dF = self$options$dF,
                dF2 = self$options$dF2,
                tail = private$.tail(),
                obsStat = self$options$obsStat))

        },
        .areaPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            st <- image$state
            plot_model_density(st$distro, st$dF, st$dF2, st$obsStat, st$tail)
        }
        )
)

#' Plot a reference density with the tail area shaded
#'
#' Draws the standard normal, t, chi-square, or F density with the region
#' corresponding to the requested tail shaded and the observed value marked
#' by a dashed line.  The shaded region always starts exactly at the
#' observed value (for a negative Z with a right tail, the shading starts
#' at that negative value, not at its absolute value).
#'
#' @param distro One of \code{"ndistro"}, \code{"tdistro"}, \code{"chisq"},
#'   \code{"fdistro"}.
#' @param dF Degrees of freedom (numerator df for F).
#' @param dF2 Denominator degrees of freedom (F only).
#' @param obs The observed value of the statistic.
#' @param tail \code{"right"}, \code{"left"}, or \code{"both"}.  Chi-square
#'   and F ignore this and shade the right tail.
#' @return A \code{ggplot} object.
#' @keywords internal
plot_model_density <- function(distro, dF, dF2 = NULL, obs = 0, tail = "right") {

    symmetric <- distro %in% c("ndistro", "tdistro")
    if (!symmetric) tail <- "right"

    dfun <- switch(distro,
        ndistro = function(x) stats::dnorm(x),
        tdistro = function(x) stats::dt(x, dF),
        chisq   = function(x) stats::dchisq(x, dF),
        fdistro = function(x) stats::df(x, dF, dF2))

    if (symmetric) {
        xl <- max(3.5, abs(obs)) * 1.5
        x <- seq(-xl, xl, length.out = 601)
    } else {
        qfun <- switch(distro,
            chisq   = function(p) stats::qchisq(p, dF),
            fdistro = function(p) stats::qf(p, dF, dF2))
        xmax <- max(qfun(0.995), obs * 1.15, 1e-6)
        # start a hair above zero: chi-square / F densities with df = 1
        # are unbounded at zero
        x <- seq(xmax / 400, xmax, length.out = 601)
    }
    curve <- data.frame(x = x, y = dfun(x))

    # Shaded tail(s): each region is cut exactly at the observed value so
    # the shading begins at the dashed line
    region <- function(lo, hi) {
        xs <- x[x > lo & x < hi]
        xs <- c(lo, xs, hi)
        xs <- xs[is.finite(xs)]
        data.frame(x = xs, y = dfun(xs))
    }
    shade <- list()
    lines <- numeric()
    if (tail == "right") {
        lo <- if (symmetric) obs else max(obs, min(x))
        if (lo < max(x)) shade[[1]] <- region(lo, max(x))
        lines <- obs
    } else if (tail == "left") {
        if (obs > min(x)) shade[[1]] <- region(min(x), obs)
        lines <- obs
    } else {
        a <- abs(obs)
        if (a < max(x)) {
            shade[[1]] <- region(-max(x), -a)
            shade[[2]] <- region(a, max(x))
        }
        lines <- c(-a, a)
    }

    p <- ggplot2::ggplot(curve, ggplot2::aes(x = x, y = y))
    for (s in shade)
        p <- p + ggplot2::geom_area(data = s, stat = "identity",
                                    fill = "red", alpha = 0.5)
    p <- p +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept = lines, linetype = "dashed", color = "red") +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(NULL, breaks = NULL,
                                    expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::xlab("") +
        ggplot2::theme(text = ggplot2::element_text(size = 18))

    if (!symmetric)
        p <- p + ggplot2::scale_x_continuous(limits = c(0, max(x)),
                                             expand = ggplot2::expansion(mult = c(0, 0.02)))

    p
}
