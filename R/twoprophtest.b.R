
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

                perms <- cached_sims(self$results$simtable,
                    list(mat = unname(mat), compare = self$options$compare, reps = self$options$reps, seedBool = self$options$seedBool, rngSeed = self$options$rngSeed),
                    function() private$.computePerms(mat))
                simres <- private$.computePval(perms, dp$dp)

                private$.populateSimTable(simres)
                private$.preparePlot(perms, dp$dp, simres$direction, mat)

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
            list(obsDiff = dp, reps = reps, p = format_sim_pval(pval, reps), direction = direction)
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

            # Only pass table-defined columns (not direction)
            table$addRow(rowKey=1, values=simres[c("obsDiff", "reps", "p")])

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

        .preparePlot = function(perms, dp, direction, mat = NULL) {
            permplot <- self$results$Plot
            dotHist <- self$options$dotHist
            m <- if (!is.null(mat) && self$options$compare == "columns") t(mat) else mat
            permplot$setState(list(df=strip_infer(perms), obs_stat=dp, direction=direction, dotHist=dotHist, showCounts=self$options$showCounts, domain=c(-1, 1),
                                          xlab=diff_label("proportions", rownames(m), colnames(m)[1]), obs_label="Observed\nDifference"))
        },
        .permPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            st <- image$state
            plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                           xlab = state_or(st$xlab, diff_label("proportions")),
                           obs_label = "Observed\nDifference",
                           domain = c(-1, 1),
                           show_counts = isTRUE(st$showCounts),
                           plot_width = image$width)
        },

        #### Helper functions ----
        .sourcifyOption = function(option) {
            if (option$name %in% c('rows', 'cols', 'counts'))
                return('')
            super$.sourcifyOption(option)
        }



    )
)
