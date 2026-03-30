
#' @importFrom jmvcore .
ContTabHTestClass <- R6::R6Class(
    "ContTabHTestClass",
    inherit=ContTabHTestBase,
    #### Active bindings ----
    active = list(
        countsName = function() {
            if ( ! is.null(self$options$counts)) {
                return(self$options$counts)
            } else if ( ! is.null(attr(self$data, "jmv-weights-name"))) {
                return (attr(self$data, "jmv-weights-name"))
            }
            NULL
        }
    ),
    private=list(
        #### Init + run functions ----
        .init=function() {

            data <- private$.cleanData()

            private$.initConttable(data) # initialize contingency table
            private$.initX2table(data) # initialize difference in proportions table

            private$.initSimTable() # initialize simulation table
            private$.initPlot() # initialize sumulation plot

            countsName <- self$countsName

            if ( ! is.null(countsName)) {
                message <- jmvcore::..('The data is weighted by the variable {}.', countsName)
                type <- NoticeType$WARNING

                weightsNotice <- jmvcore::Notice$new(
                    self$options,
                    name='.weights',
                    type=type,
                    content=message)
                self$results$insert(1, weightsNotice)
            }

        },
        .run=function() {

            rowVarName <- self$options$rows
            colVarName <- self$options$cols
            countsName <- self$options$counts

            if (is.null(rowVarName) || is.null(colVarName))
                return()

            data <- private$.cleanData()

            if (nlevels(data[[rowVarName]]) < 2)
                jmvcore::reject(.("Row variable '{var}' contains fewer than 2 levels"), code='', var=rowVarName)
            if (nlevels(data[[colVarName]]) < 2)
                jmvcore::reject(.("Column variable '{var}' contains fewer than 2 levels"), code='', var=colVarName)

            if ( ! is.null(data$.COUNTS)) {
                if (any(data$.COUNTS < 0, na.rm=TRUE))
                    jmvcore::reject(.('Counts may not be negative'))
                if (any(is.infinite(data$.COUNTS)))
                    jmvcore::reject(.('Counts may not be infinite'))
            }

            mats <- private$.matrices(data) # counts arranged as in a contingency table with stanardized formatting
            mat <- mats[[1]]

            private$.populateContTable(mat) # fill in contingency table

            suppressWarnings({

                x2 <- NULL
                #lor <- NULL

                #if (all(dim(mat) == 2) && all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
                if (all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
                    x2 <- private$.computeX2(mat)
                    #lor <- vcd::loddsratio(mat) # need this? e.g., is it doing some sort of indirect check?
                }

            }) # suppressWarnings

            private$.populateX2table(mat, x2) #, lor) # fill in the difference in proportion table

            #if ( ! is.null(lor)) {

            perms <- private$.computePerms(mat)
            simres <- private$.computePval(perms, x2$x2)

            private$.populateSimTable(simres)
            private$.preparePlot(perms, x2$x2)

            #}

        },

        #### Compute results ----

        .computeX2 = function(mat) {

            if (self$options$compare == "columns")
                mat <- t(mat)

           x2 <- chisq.test(mat)$statistic

            return(list(x2=x2))

        },

        .computePval = function(perms, x2) {
            reps <- self$options$reps
            pval <- compute_null_pval(perms, x2, "greater")
            list(reps = reps, pval = pval)
        },

        .computePerms = function(mat){

            if (self$options$compare == "columns")
                mat <- t(mat)

            dims <- dim(mat)

            # create data from from contingency table for use with infer functions
            df <- tibble::tibble(Group = rep(paste0("G", seq(1, dims[1])), dims[2]),
                                 Outcome = rep(paste0("O", seq(1, dims[2])), each = dims[1]),
                                 Count = as.vector(mat)) %>%
                tidyr::uncount(Count)

            reps <- self$options$reps

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            if(dims[2] == 2){

                perms <- df %>%
                    infer::specify(Outcome ~ Group, success = "O1") %>%
                    infer::hypothesize(null = "independence") %>%
                    infer::generate(reps = reps, type = "permute") %>%
                    infer::calculate(stat = "Chisq")

            } else {

                perms <- df %>%
                    infer::specify(Outcome ~ Group) %>%
                    infer::hypothesize(null = "independence") %>%
                    infer::generate(reps = reps, type = "permute") %>%
                    infer::calculate(stat = "Chisq")

            }


            return(perms)

        },

        #### Init tables/plots functions ----
        .initConttable = function(data) {
            rowVarName <- self$options$rows
            colVarName <- self$options$cols
            countsName <- self$countsName

            freqs <- self$results$freqs

            # add the row column, containing the row variable
            # fill in dots, if no row variable specified

            if ( ! is.null(rowVarName))
                title <- rowVarName
            else
                title <- '.'

            freqs$addColumn(
                name=title,
                title=title,
                type='text')

            # add the column columns (from the column variable)
            # fill in dots, if no column variable specified

            if ( ! is.null(colVarName)) {
                superTitle <- colVarName
                levels <- base::levels(data[[colVarName]])
            }
            else {
                superTitle <- '.'
                levels <- c('.', '.')
            }

            countsType <- `if`(is.integer(data$.COUNTS), 'integer', 'number')

            subNames  <- c('[count]', '[expected]', '[pcRow]', '[pcCol]', '[pcTot]')
            subTitles <- c(.('Observed'), .('Expected'), .('% within row'), .('% within column'), .('% of total'))
            visible   <- c('(obs)', '(exp)', '(pcRow)', '(pcCol)', '(pcTot)')
            types     <- c(countsType, 'number', 'number', 'number', 'number')
            formats   <- c('', '', 'pc', 'pc', 'pc')

            # iterate over the sub rows

            for (j in seq_along(subNames)) {
                subName <- subNames[[j]]
                if (subName == '[count]')
                    v <- '(obs && (exp || pcRow || pcCol || pcTot))'
                else
                    v <- visible[j]

                freqs$addColumn(
                    name=paste0('type', subName),
                    title='',
                    type='text',
                    visible=v)
            }

            for (i in seq_along(levels)) {
                level <- levels[[i]]

                for (j in seq_along(subNames)) {
                    subName <- subNames[[j]]
                    freqs$addColumn(
                        name=paste0(i, subName),
                        title=level,
                        superTitle=superTitle,
                        type=types[j],
                        format=formats[j],
                        visible=visible[j])
                }
            }

            # add the Total column

            if (self$options$obs) {
                freqs$addColumn(
                    name='.total[count]',
                    title=.('Total'),
                    #title='Total',
                    type=countsType)
            }

            if (self$options$exp) {
                freqs$addColumn(
                    name='.total[exp]',
                    title=.('Total'),
                    #title='Total',
                    type='number')
            }

            if (self$options$pcRow) {
                freqs$addColumn(
                    name='.total[pcRow]',
                    title=.('Total'),
                    #title='Total',
                    type='number',
                    format='pc')
            }

            if (self$options$pcCol) {
                freqs$addColumn(
                    name='.total[pcCol]',
                    title=.('Total'),
                    #title='Total',
                    type='number',
                    format='pc')
            }

            if (self$options$pcTot) {
                freqs$addColumn(
                    name='.total[pcTot]',
                    title=.('Total'),
                    #title='Total',
                    type='number',
                    format='pc')
            }

            # populate the first column with levels of the row variable

            values <- list()
            for (i in seq_along(subNames))
                values[[paste0('type', subNames[i])]] <- subTitles[i]

            rows <- private$.grid(data=data, incRows=TRUE)

            nextIsNewGroup <- TRUE

            for (i in seq_len(nrow(rows))) {

                for (name in colnames(rows)) {
                    value <- as.character(rows[i, name])
                    if (value == '.total')
                        value <- .('Total') #value <- 'Total'

                    values[[name]] <- value
                }

                key <- paste0(rows[i,], collapse='`')
                freqs$addRow(rowKey=key, values=values)

                if (nextIsNewGroup) {
                    freqs$addFormat(rowNo=i, 1, Cell.BEGIN_GROUP)
                    nextIsNewGroup <- FALSE
                }

                if (as.character(rows[i, name]) == '.total') {
                    freqs$addFormat(rowNo=i, 1, Cell.BEGIN_END_GROUP)
                    nextIsNewGroup <- TRUE
                    if (i > 1)
                        freqs$addFormat(rowNo=i - 1, 1, Cell.END_GROUP)
                }
            }

        },
        .initX2table = function(data) {
            x2tab <- self$results$x2tab

            rows <- private$.grid(data=data, incRows=FALSE)
            values <- list()

            if (length(rows) == 0) {

                x2tab$addRow(rowKey=1, values=list())


            } else { # dont think this is needed wont happen

                for (i in seq_len(nrow(rows))) {

                    for (name in dimnames(rows)[[2]]) {
                        value <- as.character(rows[i, name])
                        if (value == '.total')
                            value <- .('Total')
                        values[[name]] <- value
                    }

                    x2tab$addRow(rowKey=i, values=values)

                }
            }

        },
        .initPlot = function() {

            permplot <- self$results$Plot

        },
        .initSimTable = function() {

            simtable <- self$results$get('simtable')

        },

        #### Populate tables functions ----

        .populateSimTable = function(simres) {

            table <- self$results$get('simtable')
            table$deleteRows()

            table$addRow(rowKey=1, values=simres)

        },
        .populateContTable = function(mat) {

            freqRowNo <- 1

            suppressWarnings({

                test <- try(chisq.test(mat, correct = FALSE))

                n <- sum(mat)

                if (base::inherits(test, 'try-error'))
                    exp <- mat
                else
                    exp <- test$expected

            })

            freqs <- self$results$freqs

            data <- private$.cleanData()
            rowVarName <- self$options$rows
            colVarName <- self$options$cols

            nRows  <- base::nlevels(data[[rowVarName]])
            nCols  <- base::nlevels(data[[colVarName]])

            total <- sum(mat)
            colTotals <- apply(mat, 2, sum)
            rowTotals <- apply(mat, 1, sum)

            for (rowNo in seq_len(nRows)) {

                values <- mat[rowNo,]
                rowTotal <- sum(values)

                pcRow <- values / rowTotal

                values <- as.list(values)
                names(values) <- paste0(1:nCols, '[count]')
                values[['.total[count]']] <- rowTotal

                expValues <- exp[rowNo,]
                expValues <- as.list(expValues)
                names(expValues) <- paste0(1:nCols, '[expected]')
                expValues[['.total[exp]']] <- sum(exp[rowNo,])

                pcRow <- as.list(pcRow)
                names(pcRow) <- paste0(1:nCols, '[pcRow]')
                pcRow[['.total[pcRow]']] <- 1

                pcCol <- as.list(mat[rowNo,] / colTotals)
                names(pcCol) <- paste0(1:nCols, '[pcCol]')
                pcCol[['.total[pcCol]']] <- unname(rowTotals[rowNo] / total)

                pcTot <- as.list(mat[rowNo,] / total)
                names(pcTot) <- paste0(1:nCols, '[pcTot]')
                pcTot[['.total[pcTot]']] <- sum(mat[rowNo,] / total)

                values <- c(values, expValues, pcRow, pcCol, pcTot)

                freqs$setRow(rowNo=freqRowNo, values=values)
                freqRowNo <- freqRowNo + 1
            }

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

            freqs$setRow(rowNo=freqRowNo, values=values)
            freqRowNo <- freqRowNo + 1

        },
        .populateX2table = function(mat, x2) { #, lor) {

            x2tab <- self$results$x2tab

            x2tab$setRow(rowNo=1, list(`v[x2]`=x2$x2))

            # I think lor is just being used to check if this is a 2 by 2 contingency table
            # shouln't be needed here
            # othRowNo <- 1
            #
            # if ( ! is.null(lor)) {
            #     diffProp$setRow(rowNo=othRowNo, list(
            #         `v[dp]`=dp$dp
            #     ))
            #
            #     footnote <- `if`(self$options$compare == 'rows', .('Rows compared'), .('Columns compared'))
            #     diffProp$addFootnote(rowNo=othRowNo, 'v[dp]', footnote)
            #
            #
            # } else {
            #     diffProp$setRow(rowNo=othRowNo, list(
            #         `v[dp]`=NaN))
            #     diffProp$addFootnote(rowNo=othRowNo, 'v[dp]', .('Available for 2x2 tables only'))
            # }

        },

        #### Plot functions ----

        .preparePlot = function(perms, x2) {
            permplot <- self$results$Plot
            dotHist <- self$options$dotHist
            permplot$setState(list(df=perms, obs_stat=x2, direction="greater", dotHist=dotHist))
        },
        .permPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            st <- image$state
            plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                           xlab = "X\u00B2",
                           obs_label = "Observed\nX\u00B2")
        },

        #### Helper functions ----
        .cleanData = function(B64 = FALSE) {

            data <- self$data

            rowVarName <- self$options$rows
            colVarName <- self$options$cols
            countsName <- self$options$counts

            columns <- list()

            if ( ! is.null(rowVarName)) {
                columns[[rowVarName]] <- as.factor(data[[rowVarName]])
            }
            if ( ! is.null(colVarName)) {
                columns[[colVarName]] <- as.factor(data[[colVarName]])
            }

            if ( ! is.null(countsName)) {
                columns[['.COUNTS']] <- jmvcore::toNumeric(data[[countsName]])
            } else if ( ! is.null(attr(data, "jmv-weights"))) {
                columns[['.COUNTS']] <- jmvcore::toNumeric(attr(data, "jmv-weights"))
            } else {
                columns[['.COUNTS']] <- as.integer(rep(1, nrow(data)))
            }

            if (B64)
                names(columns) <- jmvcore::toB64(names(columns))

            attr(columns, 'row.names') <- paste(seq_len(length(columns[[1]])))
            class(columns) <- 'data.frame'

            columns
        },
        .matrices=function(data) {

            matrices <- list()

            rowVarName <- self$options$rows
            colVarName <- self$options$cols

            matrices <- list(ftable(xtabs(.COUNTS ~ ., data=data)))

            matrices
        },
        .grid=function(data, incRows=FALSE) {

            rowVarName <- self$options$rows

            expand <- list()

            if (incRows) {
                if (is.null(rowVarName))
                    expand[['.']] <- c('.', '. ', .('Total'))
                    #expand[['.']] <- c('.', '. ', 'Total')
                else
                    expand[[rowVarName]] <- c(base::levels(data[[rowVarName]]), '.total')
            }

            rows <- rev(expand.grid(expand))

            rows
        },

        .sourcifyOption = function(option) {
            if (option$name %in% c('rows', 'cols', 'counts'))
                return('')
            super$.sourcifyOption(option)
        }



    )
)
