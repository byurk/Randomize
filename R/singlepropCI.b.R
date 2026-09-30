
# This file is a generated template, your changes will not be overwritten

SinglePropCIClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "SinglePropCIClass",
    inherit = SinglePropCIBase,
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
                simres <- private$.computeCI(boot)

                private$.populateSummTable(results)
                private$.populateSimTable(simres)

                private$.preparePlot(boot)

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

            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total
            levels <- results$levels

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            if (any(counts == 0)) {
                # infer::specify() drops unused factor levels, so a level
                # with no observations makes the response look
                # single-level and calculate() refuses to compute a
                # proportion.  Resample directly instead: a bootstrap
                # proportion from a binary sample is binomial with the
                # observed proportion as its success probability.
                return(data.frame(replicate = seq_len(reps),
                                  stat = stats::rbinom(reps, total, counts[1] / total) / total))
            }

            df <- tibble::tibble(level=levels, count=counts) %>%
                tidyr::uncount(count)

            boot <- df %>%
                infer::specify(response = level, success = levels[1]) %>%
                infer::generate(reps = reps, type = "bootstrap") %>%
                infer::calculate(stat = "prop")


            return(boot)
        },
        .computeCI = function(boot) {

            resp <- self$options$resp
            reps <- self$options$reps
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            results <- private$.counts(resp)
            counts <- results$counts
            total  <- results$total

            obs_stat <- counts[1] / total

            ci <- compute_boot_ci(boot, obs_stat, confLevel, ciType,
                                  clamp = c(0, 1))

            simres <- c(list(reps = reps, obsProp = obs_stat), ci)
            return(simres)
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

            simres <- within(simres, rm(se, zcrit))

            table$addRow(rowKey=1, values=simres)

        },

        #### Plot functions ----
        .preparePlot = function(boot) {

            bootplot <- self$results$Plot
            dotHist <- self$options$dotHist
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            resp <- self$options$resp
            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total
            levels <- results$levels

            obs_stat <- counts[1] / total

            bootplot$setState(list(df=strip_infer(boot), obs_stat=obs_stat, confLevel = confLevel, ciType = ciType, dotHist=dotHist, showCounts=self$options$showCounts,
                                          xlab="proportion", stat_label="bootstrap proportions", clamp=c(0, 1)))

        },
        .bootPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            st <- image$state
            p <- plot_boot_dist(st$df, st$obs_stat, st$confLevel, st$ciType,
                                st$dotHist, xlab = "proportion",
                                stat_label = "bootstrap proportions",
                                clamp = c(0, 1),
                                show_counts = isTRUE(st$showCounts),
                           plot_width = image$width)
            return(p)
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

            results <- private$.counts(resp)
            if (length(results$counts) != 2) {
                jmvcore::reject(
                    jmvcore::format(
                        "Variable '{resp}' must have exactly 2 levels (found {n}). A proportion is only defined for a binary variable.",
                        resp=resp, n=length(results$counts)),
                    code=''
                )
            }
            if (any(is.na(results$counts)) || results$total <= 0) {
                jmvcore::reject(
                    jmvcore::format("Counts for '{resp}' must be non-negative numbers with a positive total", resp=resp),
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
