#' Convenience functions for using Randomize results from R
#'
#' These helpers make it easy to extract plots and tables from
#' analysis results for use in R scripts, Quarto documents,
#' R Markdown, or any other R-based workflow.
#'
#' @import ggplot2
NULL

# ---- Internal: locate the right Image object ----

find_image <- function(x, which) {
    if (which == "sim") {
        # Simulation / bootstrap distribution plot
        # Different analyses use different names
        if (!is.null(x[["simplot"]])) return(x[["simplot"]])
        if (!is.null(x[["Plot"]]))    return(x[["Plot"]])
        return(NULL)
    }

    if (which == "desc") {
        # Descriptive (mean/median) plot — stored in an Array
        arr <- x[["descplot"]]
        if (is.null(arr)) return(NULL)
        # Get the first element, which contains a $desc sub-image
        keys <- arr$itemKeys
        if (length(keys) == 0) return(NULL)
        item <- arr$get(key = keys[[1]])
        if (!is.null(item[["desc"]])) return(item[["desc"]])
        return(item)
    }

    if (which == "line") {
        # Regression scatterplot — stored in an Array
        arr <- x[["linplot"]]
        if (is.null(arr)) return(NULL)
        keys <- arr$itemKeys
        if (length(keys) == 0) return(NULL)
        item <- arr$get(key = keys[[1]])
        if (!is.null(item[["lp"]])) return(item[["lp"]])
        return(item)
    }

    NULL
}

# ---- plot() S3 method for Jamovi result objects ----

#' Extract and display a plot from a Randomize analysis result
#'
#' @param x A result object from any Randomize analysis function
#' @param which Which plot to extract:
#'   \code{"sim"} (default) for the simulation/bootstrap distribution,
#'   \code{"desc"} for the descriptive means/medians plot,
#'   \code{"line"} for the regression scatterplot
#' @param ... Additional arguments (ignored)
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' d <- data.frame(score = rnorm(50), group = factor(rep(c("A","B"), 25)))
#' r <- twomeanhtest(data = d, vars = "score", group = "group",
#'                   hypothesis = "different", reps = 1000,
#'                   dotHist = "histogram", seedBool = TRUE, rngSeed = 42)
#' plot(r)              # permutation distribution
#' plot(r, "desc")      # means/medians by group (requires desc=TRUE, plots=TRUE)
#' }
#'
#' @export
plot.Group <- function(x, which = c("sim", "desc", "line"), ...) {
    which <- match.arg(which)
    img <- find_image(x, which)

    if (is.null(img)) {
        available <- c()
        if (!is.null(find_image(x, "sim")))  available <- c(available, "sim")
        if (!is.null(find_image(x, "desc"))) available <- c(available, "desc")
        if (!is.null(find_image(x, "line"))) available <- c(available, "line")
        if (length(available) == 0) {
            stop("No plots are available in this result object.")
        }
        stop("No '", which, "' plot available. Try: ",
             paste0("'", available, "'", collapse = ", "))
    }

    if (is.null(img$state)) {
        stop("Plot '", which, "' has no data. ",
             "For 'desc' plots, re-run with desc=TRUE and plots=TRUE. ",
             "For 'line' plots, re-run with plots=TRUE.")
    }

    img$plot$fun()
}

# ---- Table extraction helpers ----

#' Extract the main results table as a data frame
#'
#' Returns the primary results table (hypothesis test results or
#' confidence interval bounds) from any Randomize analysis.
#'
#' The columns vary by analysis type:
#'
#' \strong{Hypothesis tests} (all include \code{reps} and \code{p}):
#' \itemize{
#'   \item \code{twomeanhtest}: \code{md} (observed mean difference), \code{reps}, \code{p}
#'   \item \code{pairedmeanhtest}: \code{md} (observed mean difference), \code{reps}, \code{p}
#'   \item \code{multimeanhtest}: \code{oF} (observed F statistic), \code{reps}, \code{p}
#'   \item \code{slopehtest}: \code{b} (observed slope), \code{reps}, \code{p}
#'   \item \code{SinglePropHTest}: \code{obsProp} (observed proportion), \code{reps}, \code{p}
#'   \item \code{TwoPropHTest}: \code{obsDiff} (observed difference in proportions), \code{reps}, \code{p}
#'   \item \code{ContTabHTest}: \code{x2} (observed X\eqn{^2}), \code{reps}, \code{p}
#' }
#'
#' \strong{Model-based} (no simulation):
#' \itemize{
#'   \item \code{modelBased}: \code{obsStat}, \code{area} (tail area/p-value)
#' }
#'
#' \strong{Confidence intervals} (all include \code{reps}, \code{cil}, \code{ciu}):
#' \itemize{
#'   \item \code{SingleMeanCI}: \code{reps}, \code{obsMean}, \code{cil}, \code{ciu}
#'   \item \code{twomeanCI}: \code{reps}, \code{obsDiff}, \code{cil}, \code{ciu}
#'   \item \code{pairedmeanCI}: \code{reps}, \code{obsDiff}, \code{cil}, \code{ciu}
#'   \item \code{slopeCI}: \code{reps}, \code{b} (observed slope), \code{cil}, \code{ciu}
#'   \item \code{SinglePropCI}: \code{reps}, \code{obsProp}, \code{cil}, \code{ciu}
#'   \item \code{TwoPropCI}: \code{reps}, \code{obsDiff}, \code{cil}, \code{ciu}
#' }
#'
#' @param x A result object from any Randomize analysis function
#' @return A data frame with columns specific to the analysis type (see Details)
#'
#' @examples
#' \dontrun{
#' r <- twomeanhtest(data = d, vars = "score", group = "group",
#'                   hypothesis = "different", reps = 1000,
#'                   seedBool = TRUE, rngSeed = 42)
#' rt <- results_table(r)
#' rt$p    # the p-value
#' rt$md   # the observed mean difference
#' }
#'
#' @export
results_table <- function(x) {
    # Try result table names in priority order
    for (nm in c("htest", "CITable", "simtable", "x2tab", "areaTable")) {
        tbl <- x[[nm]]
        if (!is.null(tbl) && inherits(tbl, "Table")) {
            return(tbl$asDF)
        }
    }
    stop("No results table found in this object. Available elements: ",
         paste(names(x), collapse = ", "))
}

#' Extract descriptive statistics as a data frame
#'
#' Returns the descriptive statistics table from analyses that provide
#' group-level summaries.  Bracket-indexed column names from Jamovi
#' (e.g. \code{mean[1]}) are cleaned to R-friendly names (e.g.
#' \code{mean1}).
#'
#' \strong{Important}: For two-sample and paired analyses, descriptive
#' statistics are only computed when \code{desc = TRUE} is passed to
#' the analysis function.  If the table contains all \code{NA} values,
#' re-run the analysis with \code{desc = TRUE}.
#'
#' Columns vary by analysis type:
#' \itemize{
#'   \item \code{twomeanCI}, \code{twomeanhtest}: \code{dep}, \code{group1},
#'     \code{num1}, \code{mean1}, \code{med1}, \code{sd1}, \code{group2},
#'     \code{num2}, \code{mean2}, \code{med2}, \code{sd2}
#'   \item \code{pairedmeanCI}, \code{pairedmeanhtest}: \code{name},
#'     \code{num}, \code{m} (mean), \code{med}, \code{sd}
#'   \item \code{multimeanhtest}: \code{dep}, \code{group}, \code{num},
#'     \code{mean}, \code{median}, \code{sd}
#'   \item \code{SingleMeanCI}: \code{var}, \code{num}, \code{mean},
#'     \code{med}, \code{sd}
#'   \item \code{SinglePropCI}, \code{SinglePropHTest}: \code{var},
#'     \code{level}, \code{count}, \code{total}, \code{prop}
#' }
#'
#' @param x A result object from any Randomize analysis function
#' @return A data frame with cleaned column names (see Details)
#'
#' @examples
#' \dontrun{
#' r <- twomeanhtest(data = d, vars = "score", group = "group",
#'                   hypothesis = "different", reps = 1000,
#'                   seedBool = TRUE, rngSeed = 42, desc = TRUE)
#' dt <- desc_table(r)
#' dt$mean1   # mean of group 1
#' dt$mean2   # mean of group 2
#' dt$num1    # n of group 1
#' }
#'
#' @export
desc_table <- function(x) {
    # Try descriptive table names
    for (nm in c("desc", "summtable", "freqs")) {
        tbl <- x[[nm]]
        if (!is.null(tbl) && inherits(tbl, "Table")) {
            df <- tbl$asDF
            # Check for all-NA content (desc=TRUE was not set)
            data_cols <- setdiff(names(df), c("name", "var", "dep"))
            if (nrow(df) > 0 && all(is.na(df[, data_cols, drop = FALSE]))) {
                stop("Descriptive statistics table is empty. ",
                     "Re-run the analysis with desc=TRUE to populate it.")
            }
            # Clean bracket-indexed column names: mean[1] -> mean1
            names(df) <- gsub("\\[([0-9]+)\\]", "\\1", names(df))
            return(df)
        }
    }
    stop("No descriptive statistics table found. ",
         "This analysis type may not provide descriptive statistics, ",
         "or re-run with desc=TRUE if available.")
}
