
# This file is a generated template, your changes will not be overwritten

SinglePropHTestClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "SinglePropHTestClass",
    inherit = SinglePropHTestBase,
    private = list(
        #### Init + run functions ----
        .init = function() {

            private$.initSummTable()
            private$.initSimTable()
            private$.initPlot()

        },
        .run = function() {

            private$.errorCheck()

            if (length(self$options$resp) > 0) {

                results <- private$.computeSumm()
                boot <- private$.computeBoots()
                simres <- private$.computePval(boot)

                private$.populateSummTable(results)
                private$.populateSimTable(simres)

                private$.preparePlot(boot, simres$direction)

            }
        },

        #### Compute results ----
        .computeSumm = function() {

            resp <- self$options$resp

            levelList <- list()

            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total
            levels <- results$levels

            for (i in seq_along(counts)) {

                level <- levels[i]
                count <- counts[i]
                prop <- count / total

                levelList[[i]] <- list(var=resp, level=level, count=count, total=total,
                                       prop=prop)
            }

            return(levelList)
        },
        .computeBoots = function() {

            resp <- self$options$resp
            reps <- self$options$reps
            testValue <- self$options$testValue

            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total
            levels <- results$levels

            df <- tibble::tibble(level=levels, count=counts) %>%
                tidyr::uncount(count)

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            boot <- df %>%
                infer::specify(response = level, success = levels[1]) %>%
                infer::hypothesize(null = "point", p = testValue) %>%
                infer::generate(reps = reps, type = "draw") %>%
                infer::calculate(stat = "prop")


            return(boot)
        },
        .computePval = function(boot) {

            resp <- self$options$resp
            reps <- self$options$reps

            direction <- map_direction(self$options$alt)

            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total

            pval <- compute_null_pval(boot, counts[1] / total, direction)

            list(obsProp = counts[1] / total, reps = reps, p = pval, direction = direction)
        },

        #### Init tables/plots functions ----
        .initSummTable = function() {

            summtable <- self$results$get('summtable')

            resp <- self$options$resp

            if(!is.null(resp)){

            varData <- jmvcore::naOmit(self$data[[resp]])

            if (self$options$areCounts) {

                levels <- paste(1:2)

            } else {

                levels <- base::levels(varData)
                if (length(levels) == 0)
                    levels <- c('\u2026', '\u2026 ')
            }

            for (level in levels) {
                key <- paste0(resp, '`', level)
                summtable$addRow(rowKey=key, values=list(var=resp, level=level))
            }

            summtable$addFormat(rowKey=paste0(resp, '`', levels[1]), 'var', jmvcore::Cell.BEGIN_GROUP)
            summtable$addFormat(rowKey=paste0(resp, '`', levels[length(levels)]), 'var', jmvcore::Cell.END_GROUP)
            }

        },
        .initSimTable = function() {

            simtable <- self$results$get('simtable')

            alt <- self$options$alt
            if (alt == 'greater')
                op <- '>'
            else if (alt == 'less')
                op <- '<'
            else
                op <- '\u2260'

            simtable$setNote(
                'hyp',
                jmvcore::format(
                    'H\u2090: p {direction} {testValue}',
                    direction=op,
                    testValue=self$options$testValue
                )
            )
        },
        .initPlot = function() {

            bootplot <- self$results$Plot

        },

        #### Populate tables functions ----
        .populateSummTable = function(results) {

            table <- self$results$get('summtable')
            table$deleteRows()

            resp <- self$options$resp
            varResults <- results

            levels <- character(length(varResults))
            for (i in seq_along(levels))
                levels[i] <- varResults[[i]]$level

            for (i in seq_along(levels)) {

                level <- levels[i]
                key <- paste0(resp, '`', level)
                table$addRow(rowKey=key, values=varResults[[i]])

            }

            table$addFormat(rowKey=paste0(resp, '`', levels[1]), 'var', jmvcore::Cell.BEGIN_GROUP)
            table$addFormat(rowKey=paste0(resp, '`', levels[length(levels)]), 'var', jmvcore::Cell.END_GROUP)

        },
        .populateSimTable = function(simres) {

            table <- self$results$get('simtable')
            table$deleteRows()

            # Only pass table-defined columns (not direction)
            table$addRow(rowKey=1, values=simres[c("obsProp", "reps", "p")])

        },

        #### Plot functions ----
        .preparePlot = function(boot, direction) {
            bootplot <- self$results$Plot
            dotHist <- self$options$dotHist
            resp <- self$options$resp
            results <- private$.counts(resp)
            obs_stat <- results$counts[1] / results$total
            bootplot$setState(list(df=strip_infer(boot), obs_stat=obs_stat, direction=direction, dotHist=dotHist,
                                          xlab="proportion", obs_label="Observed\nProportion"))
        },
        .bootPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            st <- image$state
            plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                           xlab = "proportion",
                           obs_label = "Observed\nProportion")
        },

        #### Helper functions ----
        .errorCheck = function() {

            resp <- self$options$resp
            data <- self$data

            if(!is.null(resp)){

            column <- jmvcore::naOmit(data[[resp]])
            if (length(column) == 0) {
                jmvcore::reject(
                    jmvcore::format("Variable '{resp}' contains no data", resp=resp),
                    code=''
                )
            }
            }

        },
        .counts = function(var) {

            varData <- jmvcore::naOmit(self$data[[var]])

            if (self$options$areCounts) {

                levels <- paste(1:length(varData))

                if (jmvcore::canBeNumeric(varData))
                    counts <- jmvcore::toNumeric(varData)
                else
                    counts <- suppressWarnings(as.numeric(as.character(varData)))

            } else {

                levels <- base::levels(varData)
                if (length(levels) == 0)
                    levels <- paste(sort(unique(varData)))

                counts <- as.vector(table(varData))

            }

            total <- base::sum(counts, na.rm=TRUE)

            list(levels=levels, counts=counts, total=total)
        }
        )
)
