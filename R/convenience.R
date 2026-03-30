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
#' @param x A result object from any Randomize analysis function
#' @return A data frame
#'
#' @examples
#' \dontrun{
#' r <- twomeanhtest(data = d, vars = "score", group = "group",
#'                   hypothesis = "different", reps = 1000,
#'                   dotHist = "histogram", seedBool = TRUE, rngSeed = 42)
#' results_table(r)
#' # Returns data frame with columns: var, md, reps, p
#' }
#'
#' @export
results_table <- function(x) {
    # Try result table names in priority order
    for (nm in c("htest", "CITable", "simtable", "x2tab")) {
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
#' Returns the descriptive statistics table (n, mean, median, sd)
#' from analyses that provide group-level summaries.
#'
#' @param x A result object from any Randomize analysis function
#' @return A data frame
#'
#' @examples
#' \dontrun{
#' r <- twomeanhtest(data = d, vars = "score", group = "group",
#'                   hypothesis = "different", reps = 1000,
#'                   dotHist = "histogram", seedBool = TRUE, rngSeed = 42,
#'                   desc = TRUE)
#' desc_table(r)
#' }
#'
#' @export
desc_table <- function(x) {
    # Try descriptive table names
    for (nm in c("desc", "summtable", "freqs")) {
        tbl <- x[[nm]]
        if (!is.null(tbl) && inherits(tbl, "Table")) {
            return(tbl$asDF)
        }
    }
    stop("No descriptive statistics table found. ",
         "Re-run the analysis with desc=TRUE if available.")
}
