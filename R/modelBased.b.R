
# This file is a generated template, your changes will not be overwritten

modelBasedClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "modelBasedClass",
    inherit = modelBasedBase,
    private = list(
        #### Init + run functions ----
        .init = function() {

            private$.initAreaTable()
            private$.initMultTable()
            private$.initPlot()

        },
        .run = function() {

            if (self$options$areaBool) {

                resultsArea <- private$.computeArea()

                private$.populateAreaTable(resultsArea)
                private$.preparePlot()

            }

            if (self$options$CIBool) {

                resultsMult <- private$.computeMult()

                private$.populateMultTable(resultsMult)

            }

        },

        #### Compute results ----
        .computeArea = function() {

            distro <- self$options$distro
            dF <- self$options$dF
            tail <- self$options$tail
            obs_stat <- self$options$obsStat

            if(distro == "ndistro")
                dF = Inf

            if(tail == "left"){
                pval <- pt(obs_stat, dF, lower.tail = TRUE)
            } else if(tail == "right"){
                pval <- pt(obs_stat, dF, lower.tail = FALSE)
            } else {
                pval <- pt(obs_stat, dF, lower.tail = TRUE)
                pval <- 2*min(pval, 1-pval)
            }

            resultList <- list(obsStat = obs_stat, area = pval)

            return(resultList)
        },
        .computeMult = function() {

            distro <- self$options$distro
            dF <- self$options$dF
            confLevel <- self$options$confLevel

            if(distro == "ndistro")
                dF = Inf

            tcrit <- qt(1 - (1 - confLevel/100)/2, dF)

            resultList <- list(confLev = confLevel, critVal = tcrit)

            return(resultList)
        },

        #### Init tables/plots functions ----
        .initAreaTable = function() {

            areatable <- self$results$get('areaTable')

            distro <- self$options$distro
            dF <- self$options$dF

            if(distro == "ndistro"){

                areatable$setNote(
                    'area',
                    'Area of shaded region under standard normal curve'
                )

            } else {

                areatable$setNote(
                    'area',
                    jmvcore::format(
                        'Area under density curve for T-distribution with df = {df}',
                        df=dF
                    )
                )

            }

        },
        .initMultTable = function() {

            multtable <- self$results$get('multTable')

            distro <- self$options$distro
            dF <- self$options$dF
            confLevel <- self$options$confLevel

            if(distro == "ndistro"){

                multtable$setNote(
                    'mult',
                    jmvcore::format(
                        'Multiplier for {cL}% CI based on standard normal distribution',
                        cL=confLevel
                    )
                )

            } else {

                multtable$setNote(
                    'mult',
                    jmvcore::format(
                        'Multiplier for {cL}% CI based on T-distribution with df = {df}',
                        cL=confLevel,
                        df=dF
                    )
                )

            }

        },
        .initPlot = function() {

            areaplot <- self$results$Plot

        },

        #### Populate tables functions ----
        .populateAreaTable = function(results) {

            table <- self$results$get('areaTable')
            table$deleteRows()

            table$addRow(rowKey=1, values=results)

        },
        .populateMultTable = function(results) {

            table <- self$results$get('multTable')
            table$deleteRows()

            table$addRow(rowKey=1, values=results)

        },

        #### Plot functions ----
        .preparePlot = function() {

            areaPlot <- self$results$Plot

            distro <- self$options$distro
            dF <- self$options$dF
            tail <- self$options$tail
            obsStat <- self$options$obsStat

            areaPlot$setState(list(distro = distro, dF = dF, tail = tail, obsStat = obsStat))

        },
        .areaPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            distro <- image$state[["distro"]]
            dF <- image$state[["dF"]]
            tail = image$state[["tail"]]
            obsStat = image$state[["obsStat"]]

            if(distro == "ndistro")
                dF = Inf

            xl <- 3.5

            if(abs(obsStat) > xl)
                xl <- abs(obsStat)

            xl <- xl + xl/2

            p <- ggplot(data.frame(x = c(-xl, xl)), aes(x)) +
                stat_function(fun = dt, args = list(df = dF))

            if(tail != "right"){

                if(tail == "both")
                    obsStat <- -abs(obsStat)

                p <- p +
                    stat_function(fun = dt, args = list(df = dF),
                                  xlim = c(-xl, obsStat),
                                  geom = "area",
                                  fill = "red",
                                  alpha = 0.5)
            }


            if(tail != "left"){

                if(tail == "both")
                    obsStat <- abs(obsStat)

                p <- p +
                    stat_function(fun = dt, args = list(df = dF),
                                  xlim = c(abs(obsStat), xl),
                                  geom = "area",
                                  fill = "red",
                                  alpha = 0.5)

            }


            p <- p +
                theme_classic() +
                scale_y_continuous(NULL, breaks = NULL) +
                xlab("") +
                theme(text = element_text(size = 18))

            return(p)
        }
        )
)
