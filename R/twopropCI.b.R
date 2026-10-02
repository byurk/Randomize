
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
                # a missing count beside present row / column values is an
                # error, not an empty cell (xtabs would carry the NA into
                # chisq.test, which then dies); rows where the row / column
                # value itself is missing are simply omitted
                if (!is.null(countsName) && any(is.na(data$.COUNTS) & !is.na(data[[rowVarName]]) & !is.na(data[[colVarName]])))
                    jmvcore::reject(.("Count variable '{v}' is missing a value in a row where '{rv}' and '{cv}' are present"), code='', v=countsName, rv=rowVarName, cv=colVarName)
            }

            mats <- conttab_matrices(data) # counts arranged as in a contingency table with standardized formatting
            mat <- mats[[1]]

            private$.populateContTable(mat) # fill in contingency table

            # A jamovi filter leaves hidden levels on the factors; they show
            # as empty rows / columns in the table above but play no part
            # in the comparison
            full_dim <- dim(mat)
            # (subsetting an ftable drops its level names; keep them as
            # dimnames so the plot can name the groups and the outcome)
            keep_r <- rowSums(mat) > 0
            keep_c <- colSums(mat) > 0
            lev_r <- attr(mat, "row.vars")[[1]]
            lev_c <- attr(mat, "col.vars")[[1]]
            mat <- mat[keep_r, keep_c, drop = FALSE]
            dimnames(mat) <- list(lev_r[keep_r], lev_c[keep_c])
            attr(mat, "full_dim") <- full_dim

            dp <- NULL
            is_2x2 <- all(dim(mat) == 2)

            if (is_2x2) {
                dp <- private$.diffProp(mat)
            }

            private$.populateDPtable(mat, dp, is_2x2) # fill in the difference in proportion table

            if (is_2x2) {

                boots <- cached_sims(self$results$simtable,
                    list(mat = unname(mat), compare = self$options$compare, reps = self$options$reps, seedBool = self$options$seedBool, rngSeed = self$options$rngSeed),
                    function() tidyr::drop_na(private$.computeBoots(mat)))
                simres <- private$.computeCI(boots, dp$dp)

                private$.populateSimTable(simres)
                private$.preparePlot(boots, dp$dp, mat)

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

            reps <- self$options$reps

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            # Drawn directly rather than through infer: with a handful of
            # rows and one or two reps a resample can contain a single
            # group, and infer then errors ("G2 is not a level of the
            # explanatory variable") instead of dropping that replicate.
            boots <- bootstrap_diff_props(unclass(mat), reps)

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
                diffProp$addFootnote(rowNo=othRowNo, 'v[dp]', if (all(attr(mat, "full_dim") == 2)) .('Not available: an empty row or column leaves a proportion undefined') else .('Available for 2x2 tables only'))
            }

        },

        #### Plot functions ----

        .preparePlot = function(boots, dp, mat = NULL) {
            m <- if (!is.null(mat) && self$options$compare == "columns") t(mat) else mat

            bootplot <- self$results$Plot
            dotHist <- self$options$dotHist
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            bootplot$setState(list(df=strip_infer(boots), obs_stat=dp, confLevel = confLevel, ciType = ciType, dotHist=dotHist, showCounts=self$options$showCounts,
                                          xlab=diff_label("proportions", rownames(m), colnames(m)[1]), stat_label="bootstrap differences", obs_label="Observed\nDifference"))

        },
        .bootPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            st <- image$state
            p <- plot_boot_dist(st$df, st$obs_stat, st$confLevel, st$ciType,
                                st$dotHist,
                                xlab = state_or(st$xlab, diff_label("proportions")),
                                stat_label = "bootstrap differences",
                                obs_label = "Observed\nDifference",
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
