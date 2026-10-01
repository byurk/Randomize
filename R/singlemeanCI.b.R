
# This file is a generated template, your changes will not be overwritten

SingleMeanCIClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "SingleMeanCIClass",
    inherit = SingleMeanCIBase,
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
                xbar <- results$mean

                boot <- cached_sims(self$results$simtable,
                    list(x = jmvcore::naOmit(self$data[[self$options$resp]]), reps = self$options$reps, seedBool = self$options$seedBool, rngSeed = self$options$rngSeed),
                    function() private$.computeBoots())
                simres <- private$.computeCI(boot, xbar)

                private$.populateSummTable(results)
                private$.populateSimTable(simres)

                private$.preparePlot(boot, xbar)

            }
        },

        #### Compute results ----
        .computeSumm = function() {

            resp <- self$options$resp

            varData <- jmvcore::naOmit(self$data[[resp]])

            descriptives <- list(var=resp, num=length(varData), mean=mean(varData), med=median(varData), sd=sd(varData))

            return(descriptives)
        },
        .computeBoots = function() {

            resp <- self$options$resp
            reps <- self$options$reps

            varData <- jmvcore::naOmit(self$data[[resp]])

            df <- tibble::tibble(val = varData)

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            if (length(varData) < 2)
                jmvcore::reject(jmvcore::format("Variable '{resp}' needs at least 2 observations", resp = resp), code = '')

            # infer::specify() runs t.test() internally, which refuses a
            # constant sample; every bootstrap mean of a constant sample
            # is that constant, so draw it directly
            if (stats::sd(varData) == 0)
                return(data.frame(replicate = seq_len(reps), stat = rep(mean(varData), reps)))

            boot <- df %>%
                infer::specify(response = val) %>%
                infer::generate(reps = reps, type = "bootstrap") %>%
                infer::calculate(stat = "mean")


            return(boot)
        },
        .computeCI = function(boot, xbar) {

            reps <- self$options$reps
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            ci <- compute_boot_ci(boot, xbar, confLevel, ciType)
            simres <- c(list(reps = reps, obsMean = xbar), ci)
            return(simres)
        },

        #### Init tables/plots functions ----
        .initSummTable = function() {

            summtable <- self$results$get('summtable')

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


            key <- paste0(resp)

            table$addRow(rowKey=key, values=results)

        },
        .populateSimTable = function(simres) {

            table <- self$results$get('simtable')
            table$deleteRows()

            simres <- within(simres, rm(se, zcrit))

            table$addRow(rowKey=1, values=simres)

        },

        #### Plot functions ----
        .preparePlot = function(boot, xbar) {

            bootplot <- self$results$Plot
            dotHist <- self$options$dotHist
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            resp <- self$options$resp

            obs_stat <- xbar

            bootplot$setState(list(df=strip_infer(boot), obs_stat=obs_stat, confLevel = confLevel, ciType = ciType, dotHist=dotHist, showCounts=self$options$showCounts,
                                          xlab="mean", stat_label="bootstrap means"))

        },
        .bootPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            st <- image$state
            p <- plot_boot_dist(st$df, st$obs_stat, st$confLevel, st$ciType,
                                st$dotHist, xlab = "mean",
                                stat_label = "bootstrap means",
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
            }

        }
        )
)
