
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

                boot <- private$.computeBoots()
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

            if(self$options$seedBool){
                set.seed(self$options$rngSeed)
            } else {
                set.seed(NULL)
            }

            boot <- df %>%
                infer::specify(response = val) %>%
                infer::generate(reps = reps, type = "bootstrap") %>%
                infer::calculate(stat = "mean")


            return(boot)
        },
        .computeCI = function(boot, xbar) {

            resp <- self$options$resp
            reps <- self$options$reps

            confLevel <- self$options$confLevel/100
            ciType <- self$options$ciType

            obs_stat <- xbar

            if(ciType == "bootperc"){

                ci <- boot %>%
                    dplyr::summarize(cil = quantile(stat, (1 - confLevel)/2),
                                     ciu = quantile(stat, 1 - (1 - confLevel)/2))

                se <- NULL
                zcrit <- NULL

            } else {

                zcrit <- qnorm(1 - (1 - confLevel)/2)

                se <- boot %>%
                    dplyr::summarize(se = sd(stat)) %>%
                    dplyr::pull()

                ci <- tibble::tibble(cil = obs_stat - zcrit * se,
                             ciu = obs_stat + zcrit * se)

            }


            cil <- ci %>% dplyr::pull(cil)
            ciu <- ci %>% dplyr::pull(ciu)

            simres <- list(reps = reps, obsMean = obs_stat, cil = cil, ciu = ciu, se = se, zcrit = zcrit)
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

            bootplot$setState(list(df=boot, obs_stat=obs_stat, confLevel = confLevel, ciType = ciType, dotHist=dotHist))

        },
        .bootPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            df <- image$state[["df"]]

            obs_stat <- image$state[["obs_stat"]]

            dotHist = image$state[["dotHist"]]

            confLevel <- image$state[["confLevel"]]

            ciType <- image$state[["ciType"]]

             xbar <- obs_stat

             boot <- df

             # We've already computed these, but it's fast anyway
             ciList <- private$.computeCI(boot, xbar)
             cil <- ciList$cil
             ciu <- ciList$ciu

             if(ciType == "bootperc"){

                 caption <- paste0("CI limits (dashed) are ",
                                   round((100-confLevel)/2,1),
                                   "% and ",
                                   round(100 - (100-confLevel)/2,1),
                                   "%",
                                   "\n percentiles of bootstrap means")

             } else {

                 se <- ciList$se
                 zcrit <- ciList$zcrit

                 caption <- paste0("CI (dashed) is calculated using SE of bootstrap means",
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
                        dplyr::mutate(x.bin = xbar + ( (stat - bw/2 - xbar) %/% bw ) * bw)

                    # adjust limits of CI (round) just like we did with other x positions
                    cila <- xbar + ( (cil - bw/2 - xbar) %/% bw ) * bw
                    ciua <- xbar + ( (ciu - bw/2 - xbar) %/% bw ) * bw

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
                    ggplot2::xlab("mean") +
                    ggplot2::coord_equal(ratio = bw*2/3) +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                          plot.caption = element_text(color = "red", hjust = 0))



            } else {

                p <- ggplot2::ggplot(data=boot, aes(x=stat)) +
                    ggplot2::geom_histogram(center = xbar, show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept = cil, linetype = "dashed", color = "red") +
                    ggplot2::geom_vline(xintercept = ciu, linetype = "dashed", color = "red") +
                    ggplot2::theme_minimal() +
                    ggplot2::xlab("mean") +
                    ggplot2::ylab("count") +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                          plot.caption = element_text(color = "red", hjust = 0))

            }


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
