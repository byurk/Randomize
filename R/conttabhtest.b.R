
#' @importFrom jmvcore .
ContTabHTestClass <- R6::R6Class(
    "ContTabHTestClass",
    inherit=ContTabHTestBase,
    private=list(
        #### Init + run functions ----
        .init=function() {

            data <- conttab_clean_data(self$data, self$options$rows, self$options$cols,
                                       self$options$counts,
                                       attr(self$data, "jmv-weights"))

            private$.initConttable(data) # initialize contingency table
            private$.initX2table(data) # initialize difference in proportions table

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

            # hidden (filtered-out) levels show as empty rows / columns in
            # the table but play no part in the test
            mat <- mat[rowSums(mat) > 0, colSums(mat) > 0, drop = FALSE]

            suppressWarnings({

                x2 <- NULL

                if (all(dim(mat) >= 2)) {
                    x2 <- private$.computeX2(mat)
                    # chisq.test() returns floating-point dust (e.g.
                    # 1.1e-31) for tables whose statistic is exactly
                    # zero; zap it so the table reads 0
                    x2$x2 <- zap_tiny(x2$x2)
                }

            }) # suppressWarnings

            private$.populateX2table(mat, x2) #, lor) # fill in the difference in proportion table

            if (!is.null(x2)) {
                perms <- cached_sims(self$results$simtable,
                    list(mat = unname(mat), compare = self$options$compare, reps = self$options$reps, seedBool = self$options$seedBool, rngSeed = self$options$rngSeed),
                    function() private$.computePerms(mat))
                # same zap for the simulations, so the >= comparison
                # against the observed value is unaffected
                perms$stat <- zap_tiny(perms$stat)
                simres <- private$.computePval(perms, x2$x2)

                private$.populateSimTable(simres)
                private$.preparePlot(perms, x2$x2)
            }

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
            list(x2 = as.numeric(x2), reps = reps, p = format_sim_pval(pval, reps))
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
            countsName <- conttab_counts_name(self$options$counts, self$data)
            conttab_init_columns(self$results$freqs, data, self$options$rows,
                                 self$options$cols, countsName, self$options)
        },
        .initX2table = function(data) {
            x2tab <- self$results$x2tab

            rows <- conttab_grid(data=data, rowVarName=self$options$rows, incRows=FALSE)
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
            data <- conttab_clean_data(self$data, self$options$rows, self$options$cols,
                                       self$options$counts,
                                       attr(self$data, "jmv-weights"))
            conttab_populate(self$results$freqs, mat, data, self$options$rows, self$options$cols)
        },
        .populateX2table = function(mat, x2) { #, lor) {

            x2tab <- self$results$x2tab

            x2tab$setRow(rowNo=1, list(`v[x2]`=x2$x2))

        },

        #### Plot functions ----

        .preparePlot = function(perms, x2) {
            permplot <- self$results$Plot
            dotHist <- self$options$dotHist
            permplot$setState(list(df=strip_infer(perms), obs_stat=as.numeric(x2), direction="greater", dotHist=dotHist, showCounts=self$options$showCounts, domain=c(0, Inf),
                                          xlab="X\u00B2", obs_label="Observed\nX\u00B2"))
        },
        .permPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            st <- image$state
            plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                           xlab = "X\u00B2",
                           obs_label = "Observed\nX\u00B2",
                           domain = c(0, Inf),
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
