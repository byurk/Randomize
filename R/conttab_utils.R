# Contingency table data-cleaning patterns in this file are adapted from
# the jmv package (https://github.com/jamovi/jmv), which is licensed under
# GPL (>= 2). See jmv/R/conttables.b.R for the original implementation.

#' Contingency table helper utilities
#'
#' Shared functions used by twopropCI, twoprophtest, and conttabhtest
#' for preparing contingency table data.
#'
#' @import jmvcore
NULL

#' Prepare a clean data frame for contingency table analysis
#'
#' Converts row and column variables to factors, handles weight/count
#' columns, and returns a clean data frame.
conttab_clean_data <- function(data, rowVarName, colVarName,
                               countsName = NULL, weights = NULL,
                               B64 = FALSE) {
    columns <- list()

    if (!is.null(rowVarName)) {
        columns[[rowVarName]] <- as.factor(data[[rowVarName]])
    }
    if (!is.null(colVarName)) {
        columns[[colVarName]] <- as.factor(data[[colVarName]])
    }

    if (!is.null(countsName)) {
        columns[['.COUNTS']] <- jmvcore::toNumeric(data[[countsName]])
    } else if (!is.null(weights)) {
        columns[['.COUNTS']] <- jmvcore::toNumeric(weights)
    } else {
        columns[['.COUNTS']] <- as.integer(rep(1, nrow(data)))
    }

    if (B64)
        names(columns) <- jmvcore::toB64(names(columns))

    attr(columns, 'row.names') <- paste(seq_len(length(columns[[1]])))
    class(columns) <- 'data.frame'

    columns
}

#' Build contingency table matrices from cleaned data
conttab_matrices <- function(data) {
    list(ftable(xtabs(.COUNTS ~ ., data = data)))
}

#' Build row expansion grid for contingency table display
conttab_grid <- function(data, rowVarName, incRows = FALSE) {
    expand <- list()

    if (incRows) {
        if (is.null(rowVarName))
            expand[['.']] <- c('.', '. ', jmvcore::..('Total'))
        else
            expand[[rowVarName]] <- c(base::levels(data[[rowVarName]]), '.total')
    }

    rev(expand.grid(expand))
}

#' Initialize contingency table columns in a Jamovi freqs table
#'
#' Sets up the column structure for displaying observed counts,
#' expected counts, and row/column/total percentages.
conttab_init_columns <- function(freqs, data, rowVarName, colVarName,
                                  countsName, options) {
    # Row variable column
    title <- if (!is.null(rowVarName)) rowVarName else '.'
    freqs$addColumn(name = title, title = title, type = 'text')

    # Column variable levels
    if (!is.null(colVarName)) {
        superTitle <- colVarName
        levels <- base::levels(data[[colVarName]])
    } else {
        superTitle <- '.'
        levels <- c('.', '.')
    }

    countsType <- `if`(is.integer(data$.COUNTS), 'integer', 'number')

    subNames  <- c('[count]', '[expected]', '[pcRow]', '[pcCol]', '[pcTot]')
    subTitles <- c(
        jmvcore::..('Observed'), jmvcore::..('Expected'),
        jmvcore::..('% within row'), jmvcore::..('% within column'),
        jmvcore::..('% of total')
    )
    visible <- c('(obs)', '(exp)', '(pcRow)', '(pcCol)', '(pcTot)')
    types   <- c(countsType, 'number', 'number', 'number', 'number')
    formats <- c('', '', 'pc', 'pc', 'pc')

    # Type sub-columns
    for (j in seq_along(subNames)) {
        subName <- subNames[[j]]
        v <- if (subName == '[count]')
            '(obs && (exp || pcRow || pcCol || pcTot))'
        else
            visible[j]

        freqs$addColumn(
            name = paste0('type', subName), title = '',
            type = 'text', visible = v)
    }

    # Data columns per level
    for (i in seq_along(levels)) {
        level <- levels[[i]]
        for (j in seq_along(subNames)) {
            subName <- subNames[[j]]
            freqs$addColumn(
                name = paste0(i, subName), title = level,
                superTitle = superTitle, type = types[j],
                format = formats[j], visible = visible[j])
        }
    }

    # Total columns
    totalLabel <- jmvcore::..('Total')
    if (options$obs)
        freqs$addColumn(name = '.total[count]', title = totalLabel, type = countsType)
    if (options$exp)
        freqs$addColumn(name = '.total[exp]', title = totalLabel, type = 'number')
    if (options$pcRow)
        freqs$addColumn(name = '.total[pcRow]', title = totalLabel, type = 'number', format = 'pc')
    if (options$pcCol)
        freqs$addColumn(name = '.total[pcCol]', title = totalLabel, type = 'number', format = 'pc')
    if (options$pcTot)
        freqs$addColumn(name = '.total[pcTot]', title = totalLabel, type = 'number', format = 'pc')

    # Populate row labels
    values <- list()
    for (i in seq_along(subNames))
        values[[paste0('type', subNames[i])]] <- subTitles[i]

    rows <- conttab_grid(data = data, rowVarName = rowVarName, incRows = TRUE)

    nextIsNewGroup <- TRUE
    for (i in seq_len(nrow(rows))) {
        for (name in colnames(rows)) {
            value <- as.character(rows[i, name])
            if (value == '.total')
                value <- totalLabel
            values[[name]] <- value
        }

        key <- paste0(rows[i, ], collapse = '`')
        freqs$addRow(rowKey = key, values = values)

        if (nextIsNewGroup) {
            freqs$addFormat(rowNo = i, 1, Cell.BEGIN_GROUP)
            nextIsNewGroup <- FALSE
        }

        if (as.character(rows[i, name]) == '.total') {
            freqs$addFormat(rowNo = i, 1, Cell.BEGIN_END_GROUP)
            nextIsNewGroup <- TRUE
            if (i > 1)
                freqs$addFormat(rowNo = i - 1, 1, Cell.END_GROUP)
        }
    }
}

#' Populate a contingency table with observed/expected counts and percentages
conttab_populate <- function(freqs, mat, data, rowVarName, colVarName) {
    freqRowNo <- 1

    suppressWarnings({
        test <- try(chisq.test(mat, correct = FALSE))
        n <- sum(mat)
        if (base::inherits(test, 'try-error'))
            exp <- mat
        else
            exp <- test$expected
    })

    nRows  <- base::nlevels(data[[rowVarName]])
    nCols  <- base::nlevels(data[[colVarName]])

    total <- sum(mat)
    colTotals <- apply(mat, 2, sum)
    rowTotals <- apply(mat, 1, sum)

    for (rowNo in seq_len(nRows)) {
        values <- mat[rowNo, ]
        rowTotal <- sum(values)

        pcRow <- values / rowTotal

        values <- as.list(values)
        names(values) <- paste0(1:nCols, '[count]')
        values[['.total[count]']] <- rowTotal

        expValues <- exp[rowNo, ]
        expValues <- as.list(expValues)
        names(expValues) <- paste0(1:nCols, '[expected]')
        expValues[['.total[exp]']] <- sum(exp[rowNo, ])

        pcRow <- as.list(pcRow)
        names(pcRow) <- paste0(1:nCols, '[pcRow]')
        pcRow[['.total[pcRow]']] <- 1

        pcCol <- as.list(mat[rowNo, ] / colTotals)
        names(pcCol) <- paste0(1:nCols, '[pcCol]')
        pcCol[['.total[pcCol]']] <- unname(rowTotals[rowNo] / total)

        pcTot <- as.list(mat[rowNo, ] / total)
        names(pcTot) <- paste0(1:nCols, '[pcTot]')
        pcTot[['.total[pcTot]']] <- sum(mat[rowNo, ] / total)

        values <- c(values, expValues, pcRow, pcCol, pcTot)

        freqs$setRow(rowNo = freqRowNo, values = values)
        freqRowNo <- freqRowNo + 1
    }

    # Total row
    values <- apply(mat, 2, sum)
    rowTotal <- sum(values)
    values <- as.list(values)
    names(values) <- paste0(1:nCols, '[count]')
    values[['.total[count]']] <- rowTotal

    expValues <- apply(mat, 2, sum)
    expValues <- as.list(expValues)
    names(expValues) <- paste0(1:nCols, '[expected]')

    pcRow <- apply(mat, 2, sum) / rowTotal
    pcRow <- as.list(pcRow)
    names(pcRow) <- paste0(1:nCols, '[pcRow]')

    pcCol <- rep(1, nCols)
    pcCol <- as.list(pcCol)
    names(pcCol) <- paste0(1:nCols, '[pcCol]')

    pcTot <- apply(mat, 2, sum) / total
    pcTot <- as.list(pcTot)
    names(pcTot) <- paste0(1:nCols, '[pcTot]')

    expValues[['.total[exp]']] <- total
    pcRow[['.total[pcRow]']] <- 1
    pcCol[['.total[pcCol]']] <- 1
    pcTot[['.total[pcTot]']] <- 1

    values <- c(values, expValues, pcRow, pcCol, pcTot)

    freqs$setRow(rowNo = freqRowNo, values = values)
}

#' Resolve the counts/weights variable name from options and data attributes
conttab_counts_name <- function(countsOption, data) {
    if (!is.null(countsOption)) {
        return(countsOption)
    } else if (!is.null(attr(data, "jmv-weights-name"))) {
        return(attr(data, "jmv-weights-name"))
    }
    NULL
}
