
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

            df <- tibble::tibble(level=levels, count=counts) %>%
                tidyr::uncount(count)

            if(self$options$seedBool){
                set.seed(self$options$rngSeed)
            } else {
                set.seed(NULL)
            }

            boot <- df %>%
                infer::specify(response = level, success = levels[1]) %>%
                infer::generate(reps = reps, type = "bootstrap") %>%
                infer::calculate(stat = "prop")


            return(boot)
        },
        .computeCI = function(boot) {

            resp <- self$options$resp
            reps <- self$options$reps
            confLevel <- self$options$confLevel/100
            ciType <- self$options$ciType

            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total
            levels <- results$levels

            obs_stat <- counts[1] / total

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

            cil <- max(cil, 0) # for proportion don't want CI to extend beyon [0,1]
            ciu <- min(ciu, 1)

            simres <- list(reps = reps, obsProp = obs_stat, cil = cil, ciu = ciu, se = se, zcrit = zcrit)
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

            bootplot$setState(list(df=boot, obs_stat=obs_stat, confLevel = confLevel, ciType = ciType, dotHist=dotHist))

        },
        .bootPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            df <- image$state[["df"]]

            obs_stat = image$state[["obs_stat"]]

            dotHist = image$state[["dotHist"]]

            confLevel <- image$state[["confLevel"]]

            ciType <- image$state[["ciType"]]

             phat <- obs_stat

             boot <- df

             # We've already computed these, but it's fast anyway
             ciList <- private$.computeCI(boot)
             cil <- ciList$cil
             ciu <- ciList$ciu

             if(ciType == "bootperc"){

                 caption <- paste0("CI limits (dashed) are ",
                                   round((100-confLevel)/2,1),
                                   "% and ",
                                   round(100 - (100-confLevel)/2,1),
                                   "%",
                                   "\n percentiles of bootstrap proportions")

             } else {

                 se <- ciList$se
                 zcrit <- ciList$zcrit

                 caption <- paste0("CI (dashed) is calculated using SE of bootstrap proportions",
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
                        dplyr::mutate(x.bin = phat + ( (stat - bw/2 - phat) %/% bw ) * bw)

                    # adjust limits of CI (round) just like we did with other x positions
                    cila <- phat + ( (cil - bw/2 - phat) %/% bw ) * bw
                    ciua <- phat + ( (ciu - bw/2 - phat) %/% bw ) * bw


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
                    ggplot2::xlab("proportion") +
                    ggplot2::coord_equal(ratio = bw*2/3) +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                          plot.caption = element_text(color = "red", hjust = 0))



            } else {

                p <- ggplot2::ggplot(data=boot, aes(x=stat)) +
                    ggplot2::geom_histogram(center = phat, show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept = cil, linetype = "dashed", color = "red") +
                    ggplot2::geom_vline(xintercept = ciu, linetype = "dashed", color = "red") +
                    ggplot2::theme_minimal() +
                    ggplot2::xlab("proportion") +
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

        },
        .counts = function(var) {

            initing <- nrow(self$data) == 0

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
