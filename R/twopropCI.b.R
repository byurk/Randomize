
#' @importFrom jmvcore .
TwoPropCIClass <- R6::R6Class(
    "TwoPropCIClass",
    inherit=TwoPropCIBase,
    private=list(
        #### Init + run functions ----
        .init=function() {

            data <- conttab_clean_data(self$data, self$options$rows, self$options$cols,
                                       self$options$counts,
                                       attr(self$data, "jmv-weights"))

            private$.initConttable(data) # initialize contingency table
            private$.initDPtable(data) # initialize difference in proportions table

            private$.initSimTable() # initialize simulation table
            private$.initPlot() # initialize sumulation plot

            countsName <- conttab_counts_name(self$options$counts, self$data)

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

            data <- conttab_clean_data(self$data, self$options$rows, self$options$cols,
                                       self$options$counts,
                                       attr(self$data, "jmv-weights"))

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

            mats <- conttab_matrices(data) # counts arranged as in a contingency table with standardized formatting
            mat <- mats[[1]]

            private$.populateContTable(mat) # fill in contingency table

            dp <- NULL
            is_2x2 <- all(dim(mat) == 2) && all(rowSums(mat) > 0) && all(colSums(mat) > 0)

            if (is_2x2) {
                dp <- private$.diffProp(mat)
            }

            private$.populateDPtable(mat, dp, is_2x2) # fill in the difference in proportion table

            if (is_2x2) {

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
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            ci <- compute_boot_ci(boots, dp, confLevel, ciType)
            simres <- c(list(reps = reps, obsDiff = dp), ci)
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

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            boots <- df %>%
                infer::specify(Outcome ~ Group, success = "O1") %>%
                infer::generate(reps = reps, type = "bootstrap") %>%
                infer::calculate(stat = "diff in props", order = c("G1", "G2"))

            return(boots)

        },

        #### Init tables/plots functions ----
        .initConttable = function(data) {
            countsName <- conttab_counts_name(self$options$counts, self$data)
            conttab_init_columns(self$results$freqs, data, self$options$rows,
                                 self$options$cols, countsName, self$options)
        },
        .initDPtable = function(data) {
            diffProp <- self$results$diffProp

            rows <- conttab_grid(data=data, rowVarName=self$options$rows, incRows=FALSE)
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
            data <- conttab_clean_data(self$data, self$options$rows, self$options$cols,
                                       self$options$counts,
                                       attr(self$data, "jmv-weights"))
            conttab_populate(self$results$freqs, mat, data, self$options$rows, self$options$cols)
        },
        .populateDPtable = function(mat, dp, is_2x2) {

            diffProp <- self$results$diffProp

            othRowNo <- 1

            if (is_2x2) {
                diffProp$setRow(rowNo=othRowNo, list(
                    `v[dp]`=dp$dp
                ))

                footnote <- `if`(self$options$compare == 'rows', .('Rows compared'), .('Columns compared'))
                diffProp$addFootnote(rowNo=othRowNo, 'v[dp]', footnote)


            } else {
                diffProp$setRow(rowNo=othRowNo, list(
                    `v[dp]`=NaN))
                diffProp$addFootnote(rowNo=othRowNo, 'v[dp]', if (all(dim(mat) == 2)) .('Not available: an empty row or column leaves a proportion undefined') else .('Available for 2x2 tables only'))
            }

        },

        #### Plot functions ----

        .preparePlot = function(boots, dp) {

            bootplot <- self$results$Plot
            dotHist <- self$options$dotHist
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            bootplot$setState(list(df=strip_infer(boots), obs_stat=dp, confLevel = confLevel, ciType = ciType, dotHist=dotHist, showCounts=self$options$showCounts,
                                          xlab="difference (group 1 - group 2)", stat_label="bootstrap differences"))

        },
        .bootPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            st <- image$state
            p <- plot_boot_dist(st$df, st$obs_stat, st$confLevel, st$ciType,
                                st$dotHist,
                                xlab = "difference (group 1 - group 2)",
                                stat_label = "bootstrap differences",
                                show_counts = isTRUE(st$showCounts),
                           plot_width = image$width)
            return(p)
        },

        #### Helper functions ----
        .sourcifyOption = function(option) {
            if (option$name %in% c('rows', 'cols', 'counts'))
                return('')
            super$.sourcifyOption(option)
        }



    )
)
