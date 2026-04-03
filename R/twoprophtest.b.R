
#' @importFrom jmvcore .
TwoPropHTestClass <- R6::R6Class(
    "TwoPropHTestClass",
    inherit=TwoPropHTestBase,
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

                perms <- private$.computePerms(mat)
                simres <- private$.computePval(perms, dp$dp)

                private$.populateSimTable(simres)
                private$.preparePlot(perms, dp$dp, simres$direction)

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

        .computePval = function(perms, dp) {
            reps <- self$options$reps
            direction <- map_direction(self$options$hypothesis)
            pval <- compute_null_pval(perms, dp, direction)
            list(reps = reps, p = pval, direction = direction)
        },

        .computePerms = function(mat){

            if (self$options$compare == "columns")
                mat <- t(mat)

            # create data from from contingency table for use with infer functions
            df <- tibble::tibble(Group = c("G1", "G2", "G1", "G2"),
                                 Outcome = c("O1", "O1", "O2", "O2"),
                                 Count = c(mat[1,1], mat[2,1], mat[1,2], mat[2,2])) %>%
                tidyr::uncount(Count)

            reps <- self$options$reps

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            perms <- df %>%
                infer::specify(Outcome ~ Group, success = "O1") %>%
                infer::hypothesize(null = "independence") %>%
                infer::generate(reps = reps, type = "permute") %>%
                infer::calculate(stat = "diff in props", order = c("G1", "G2"))

            return(perms)

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

            alt <- self$options$hypothesis

            if (alt == 'oneGreater')
                op <- '>'
            else if (alt == 'twoGreater')
                op <- '<'
            else
                op <- '\u2260'

            simtable$setNote(
                'hyp',
                jmvcore::format(
                    'H\u2090: p\u2081 - p\u2082 {direction} 0',
                    direction=op
                )
            )
        },

        #### Populate tables functions ----

        .populateSimTable = function(simres) {

            table <- self$results$get('simtable')
            table$deleteRows()

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
                diffProp$addFootnote(rowNo=othRowNo, 'v[dp]', .('Available for 2x2 tables only'))
            }

        },

        #### Plot functions ----

        .preparePlot = function(perms, dp, direction) {
            permplot <- self$results$Plot
            dotHist <- self$options$dotHist
            permplot$setState(list(df=strip_infer(perms), obs_stat=dp, direction=direction, dotHist=dotHist))
        },
        .permPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            st <- image$state
            plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                           xlab = "difference (group 1 - group 2)",
                           obs_label = "Observed\nDifference")
        },

        #### Helper functions ----
        .sourcifyOption = function(option) {
            if (option$name %in% c('rows', 'cols', 'counts'))
                return('')
            super$.sourcifyOption(option)
        }



    )
)
