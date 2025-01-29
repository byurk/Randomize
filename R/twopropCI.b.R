
#' @importFrom jmvcore .
TwoPropCIClass <- R6::R6Class(
    "TwoPropCIClass",
    inherit=TwoPropCIBase,
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
            private$.initDPtable(data) # initialize difference in proportions table

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

                dp <- NULL
                lor <- NULL

                if (all(dim(mat) == 2) && all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
                    dp <- private$.diffProp(mat)
                    lor <- vcd::loddsratio(mat) # need this? e.g., is it doing some sort of indirect check?
                }

            }) # suppressWarnings

            private$.populateDPtable(mat, dp, lor) # fill in the difference in proportion table

            if ( ! is.null(lor)) {

                boots <- private$.computeBoots(mat)
                boots <- tidyr::drop_na(boots)
                simres <- private$.computeCI(boots, dp$dp)

                private$.populateSimTable(simres)
                private$.preparePlot(boots, dp$dp)

            }

        },

        #### Compute results ----

        .diffProp = function(mat) {

            dims <- dim(mat)

            if (dims[1] > 2 || dims[2] > 2)
                return(NULL)

            if (self$options$compare == "columns")
                mat <- t(mat)

            a <- mat[1,1]
            b <- mat[1,2]
            c <- mat[2,1]
            d <- mat[2,2]

            p1 <- a / (a + b)
            p2 <- c / (c + d)

            dp <- p1 - p2

            return(list(dp=dp))

        },

        .computeCI = function(boots, dp) {

            reps <- self$options$reps
            confLevel <- self$options$confLevel/100
            ciType <- self$options$ciType

            obs_stat <- dp

            if(ciType == "bootperc"){

                ci <- boots %>%
                    dplyr::summarize(cil = quantile(stat, (1 - confLevel)/2),
                                     ciu = quantile(stat, 1 - (1 - confLevel)/2))

                se <- NULL
                zcrit <- NULL

            } else {

                zcrit <- qnorm(1 - (1 - confLevel)/2)

                se <- boots %>%
                    dplyr::summarize(se = sd(stat)) %>%
                    dplyr::pull()

                ci <- tibble::tibble(cil = obs_stat - zcrit * se,
                                     ciu = obs_stat + zcrit * se)

            }


            cil <- ci %>% dplyr::pull(cil)
            ciu <- ci %>% dplyr::pull(ciu)

            simres <- list(reps = reps, obsDiff = obs_stat, cil = cil, ciu = ciu, se = se, zcrit = zcrit)
            return(simres)
        },

        .computeBoots = function(mat){

            if (self$options$compare == "columns")
                mat <- t(mat)

            # create data from from contingency table for use with infer functions
            df <- tibble::tibble(Group = c("G1", "G2", "G1", "G2"),
                                 Outcome = c("O1", "O1", "O2", "O2"),
                                 Count = c(mat[1,1], mat[2,1], mat[1,2], mat[2,2])) %>%
                tidyr::uncount(Count)

            reps <- self$options$reps

            if(self$options$seedBool){
                set.seed(self$options$rngSeed)
            } else {
                set.seed(NULL)
            }


            boots <- df %>%
                infer::specify(Outcome ~ Group, success = "O1") %>%
                infer::generate(reps = reps, type = "bootstrap") %>%
                infer::calculate(stat = "diff in props", order = c("G1", "G2"))

            return(boots)

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
        .initDPtable = function(data) {
            diffProp <- self$results$diffProp

            rows <- private$.grid(data=data, incRows=FALSE)
            values <- list()

            if (length(rows) == 0) {

                diffProp$addRow(rowKey=1, values=list())


            } else { # dont think this is needed wont happen

                for (i in seq_len(nrow(rows))) {

                    for (name in dimnames(rows)[[2]]) {
                        value <- as.character(rows[i, name])
                        if (value == '.total')
                            value <- .('Total')
                        values[[name]] <- value
                    }

                    diffProp$addRow(rowKey=i, values=values)

                }
            }

        },
        .initPlot = function() {

            permplot <- self$results$Plot

        },
        .initSimTable = function() {

            simtable <- self$results$get('simtable')

            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            if (ciType == 'bootperc')
                typ <- "bootstrap percentile"
            else
                typ <- 'bootstrap SE'

            simtable$setNote(
                'ci',
                jmvcore::format(
                    '{clev}% {typ} confidence interval',
                    clev=confLevel,
                    typ=typ
                )
            )
        },

        #### Populate tables functions ----

        .populateSimTable = function(simres) {

            table <- self$results$get('simtable')
            table$deleteRows()

            simres <- within(simres, rm(se, zcrit))

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
        .populateDPtable = function(mat, dp, lor) {

            diffProp <- self$results$diffProp

            othRowNo <- 1

            if ( ! is.null(lor)) {
                diffProp$setRow(rowNo=othRowNo, list(
                    `v[dp]`=dp$dp
                ))

                footnote <- `if`(self$options$compare == 'rows', .('Rows compared'), .('Columns compared'))
                diffProp$addFootnote(rowNo=othRowNo, 'v[dp]', footnote)


            } else {
                diffProp$setRow(rowNo=othRowNo, list(
                    `v[dp]`=NaN))
                diffProp$addFootnote(rowNo=othRowNo, 'v[dp]', .('Available for 2x2 tables only'))
            }

        },

        #### Plot functions ----

        .preparePlot = function(boots, dp) {

            bootplot <- self$results$Plot
            dotHist <- self$options$dotHist
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            bootplot$setState(list(df=boots, obs_stat=dp, confLevel = confLevel, ciType = ciType, dotHist=dotHist))

        },
        .bootPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            df <- image$state[["df"]]

            obs_stat = image$state[["obs_stat"]]

            dotHist = image$state[["dotHist"]]

            confLevel <- image$state[["confLevel"]]

            ciType <- image$state[["ciType"]]

            dp <- obs_stat

            boot <- df

            # We've already computed these, but it's fast anyway
            ciList <- private$.computeCI(boot, dp)
            cil <- ciList$cil
            ciu <- ciList$ciu

            if(ciType == "bootperc"){

                caption <- paste0("CI limits (dashed) are ",
                                  round((100-confLevel)/2,1),
                                  "% and ",
                                  round(100 - (100-confLevel)/2,1),
                                  "%",
                                  "\n percentiles of bootstrap differences")

            } else {

                se <- ciList$se
                zcrit <- ciList$zcrit

                caption <- paste0("CI (dashed) is calculated using SE of bootstrap differences",
                                  "\n (SE = ",
                                  round(se, 3),
                                  ", z* = ",
                                  round(zcrit, 3),
                                  ")")

            }


            if(dotHist == "dotplot"){

                ndist <- dplyr::n_distinct(boot$stat)

                # bin width, also used for scaling dots when not binning
                bw <- boot %>%
                    dplyr::summarize(min = min(stat), max = max(stat)) %>%
                    dplyr::mutate(bw = (max - min) / 30 ) %>%
                    dplyr::pull(bw)


                if (ndist <= 30){ # don't bin unless there are more than 30 unique values

                    boot <- boot %>%
                        dplyr::mutate(x.bin = stat)

                    # no adjustment to position of ci on plot
                    cila <- cil
                    ciua <- ciu

                } else {

                    boot <- boot %>%
                        dplyr::mutate(x.bin = dp + ( (stat - bw/2 - dp) %/% bw ) * bw)

                    # adjust limits of CI (round) just like we did with other x positions
                    cila <- dp + ( (cil - bw/2 - dp) %/% bw ) * bw
                    ciua <- dp + ( (ciu - bw/2 - dp) %/% bw ) * bw


                }

                boot <- boot %>%
                    dplyr::group_by(x.bin) %>%
                    dplyr::mutate(y = seq_along(x.bin))

                p <- ggplot2::ggplot(boot) +
                    ggforce::geom_ellipse(aes(x0 = x.bin, y0 = y, a = bw/3, b = 0.5, angle = 0),
                                          show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept = cila, linetype = "dashed", color = "red") +
                    ggplot2::geom_vline(xintercept = ciua, linetype = "dashed", color = "red") +
                    ggplot2::theme_minimal() +
                    ggplot2::ylab("count") +
                    ggplot2::xlab("difference (group 1 - group 2)") +
                    ggplot2::coord_equal(ratio = bw*2/3) +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                                   plot.caption = element_text(color = "red", hjust = 0))



            } else {

                p <- ggplot2::ggplot(data=boot, aes(x=stat)) +
                    ggplot2::geom_histogram(center = dp, show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept = cil, linetype = "dashed", color = "red") +
                    ggplot2::geom_vline(xintercept = ciu, linetype = "dashed", color = "red") +
                    ggplot2::theme_minimal() +
                    ggplot2::xlab("difference (group 1 - group 2)") +
                    ggplot2::ylab("count") +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                                   plot.caption = element_text(color = "red", hjust = 0))

            }


            return(p)
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
