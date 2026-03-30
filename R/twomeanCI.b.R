
#' @importFrom jmvcore .
twomeanCIClass <- R6::R6Class(
  "twomeanCIClass",
  inherit=twomeanCIBase,
  private=list(
    .run=function() {
      
      groupVarName <- self$options$group
      depVarName <- self$options$vars
      varNames <- c(groupVarName, depVarName)
      
      if (is.null(groupVarName) || is.null(depVarName))
      return()
      
      data <- select(self$data, varNames)
      
      data[[depVarName]] <- jmvcore::toNumeric(data[[depVarName]])
      data[[groupVarName]] <- droplevels(as.factor(data[[groupVarName]]))
      
      CITable <- self$results$CITable
      descTable <- self$results$desc
      
      if (depVarName == groupVarName)
      jmvcore::reject(.("Grouping variable '{a}' must not also be a dependent variable"),
      code="a_is_dependent_variable", a=groupVarName)
      
      # exclude rows with missings in the grouping variable
      data <- data[ ! is.na(data[[groupVarName]]),]
      
      groupLevels <- base::levels(data[[groupVarName]])
      
      if (length(groupLevels) != 2)
      jmvcore::reject(.("Grouping variable '{a}' must have exactly 2 levels"),
      code="grouping_var_must_have_2_levels", a=groupVarName)
      
      dataCI <- data.frame(dep=data[[depVarName]], group=data[[groupVarName]])
      dataCI <- tidyr::drop_na(dataCI)

      groupLevels <- base::levels(dataCI$group)
      v <- tapply(dataCI$dep, dataCI$group, function(x) tryNaN(var(x)))
      n <- tapply(dataCI$dep, dataCI$group, length)
      m <- tapply(dataCI$dep, dataCI$group, function(x) tryNaN(mean(x)))
      med <- tapply(dataCI$dep, dataCI$group, function(x) tryNaN(median(x)))
      sd <- sqrt(v)
      
      n[is.na(n)] <- 0
      m[is.na(m)] <- NaN
      med[is.na(med)] <- NaN
      sd[is.na(sd)] <- NaN
      
      
      if (is.factor(dataCI$dep))
      res <- createError(.('Variable is not numeric'))
      else if (any(is.infinite(dataCI$dep)))
      res <- createError(.('Variable contains infinite values'))
      else
      #res <- try(t.test(dep ~ group, data=dataCI, var.equal=TRUE,
      #                  alternative=Ha, conf.level=confInt), silent=TRUE)
      #res <- try(t.test(dep ~ group, data=dataCI, var.equal=TRUE,
      #  alternative=Ha), silent=TRUE)
      boots <- private$.computeBoots(dataCI)
      boots <- tidyr::drop_na(boots)
      res <- private$.computeCI(boots, m[1]-m[2])
      res <- within(res, rm(se, zcrit))
      private$.preparePlot(boots, m[1]-m[2])
        
        if (isError(res)) {
          
          CITable$setRow(rowKey=depVarName, list(
            "reps"='',
            "obsDiff"='',
            "cil"='',
            "ciu"=''))
            
            message <- extractErrorMessage(res)
            
            CITable$addFootnote(rowKey=depVarName, 'md', message)
            
          } else {
            
            CITable$setRow(rowKey=depVarName, res)
            }
            
            if (self$options$desc) {
              
              descTable$setRow(rowKey=depVarName, list(
                "dep"=depVarName,
                "group[1]"=groupLevels[1],
                "num[1]"=n[1],
                "mean[1]"=m[1],
                "sd[1]"=sd[1],
                "med[1]"=med[1],
                "group[2]"=groupLevels[2],
                "num[2]"=n[2],
                "mean[2]"=m[2],
                "sd[2]"=sd[2],
                "med[2]"=med[2]
              ))
            }
            
            if (self$options$plots) {
              
              image <- self$results$descplot$get(key=depVarName)$desc
              
              if (nrow(dataCI) > 0) {
                
                #ciWidth <- self$options$ciWidth
                #tail <- qnorm(1 - (100 - ciWidth) / 200)
                
                means <- aggregate(dataCI$dep, by=list(dataCI$group),
                function(x) tryNaN(mean(x)), simplify=FALSE)
                #cies  <- aggregate(dataCI$dep, by=list(dataCI$group),
                #                   function(x) { tail * tryNaN(sd(x)) / sqrt(length(x)) }, simplify=FALSE)
                medians <- aggregate(dataCI$dep, by=list(dataCI$group),
                function(x) tryNaN(median(x)), simplify=FALSE)
                
                meanPlotData <- data.frame(group=means$Group.1)
                meanPlotData <- cbind(meanPlotData, stat=unlist(means$x))
                #meanPlotData <- cbind(meanPlotData, cie=unlist(cies$x))
                meanPlotData <- cbind(meanPlotData, cie=NA)
                meanPlotData <- cbind(meanPlotData, type='mean')
                
                medianPlotData <- data.frame(group=medians$Group.1)
                medianPlotData <- cbind(medianPlotData, stat=unlist(medians$x))
                medianPlotData <- cbind(medianPlotData, cie=NA)
                medianPlotData <- cbind(medianPlotData, type='median')
                
                plotData <- rbind(meanPlotData, medianPlotData)
                
                if (all(is.na(plotData$stat)))
                image$setState(NULL)
                else
                image$setState(plotData)
                
              } else {
                
                image$setState(NULL)
              }
            }
          },
          .init=function() {
            
            simtable <- self$results$CITable
            
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
          
          .descplot=function(image, ggtheme, theme, ...) {
            
            if (is.null(image$state))
            return(FALSE)
            
            groupName <- self$options$group
            
            #ciw <- self$options$ciWidth
            
            pd <- position_dodge(0.2)
            
            plot <- ggplot(data=image$state, aes(x=group, y=stat, shape=type)) +
            #geom_errorbar(aes(x=group, ymin=stat-cie, ymax=stat+cie, width=.1),
            #              size=.8, colour=theme$color[2], position=pd) +
            geom_point(aes(x=group, y=stat, shape=type), color=theme$color[1],
            fill=theme$fill[1], size=3, position=pd) +
            labs(x=groupName, y=image$key) +
            scale_shape_manual(
              name='',
              values=c(mean=21, median=22),
              labels=c(
                #mean=jmvcore::format(.('Mean ({ciWidth}% CI)'), ciWidth=ciw),
                mean=.('Mean'),
                median=.('Median')
              )
            ) +
            ggtheme +
            theme(
              plot.title=ggplot2::element_text(margin=ggplot2::margin(b = 5.5 * 1.2)),
              plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
            )
            
            return(plot)
          },
          .sourcifyOption = function(option) {
            if (option$name %in% c('deps', 'group'))
            return('')
            super$.sourcifyOption(option)
          },
          .computeCI = function(boots, dm) {

            reps <- self$options$reps
            confLevel <- self$options$confLevel/100
            ciType <- self$options$ciType

            obs_stat <- dm

            if(ciType == "bootperc"){

                ci <- boots %>%
                    dplyr::summarize(cil = quantile(stat, (1 - confLevel)/2, na.rm=TRUE),
                                     ciu = quantile(stat, 1 - (1 - confLevel)/2, na.rm=TRUE))

                se <- NULL
                zcrit <- NULL

            } else {

                zcrit <- qnorm(1 - (1 - confLevel)/2)

                se <- boots %>%
                    dplyr::summarize(se = sd(stat, na.rm=TRUE)) %>%
                    dplyr::pull()

                ci <- tibble::tibble(cil = obs_stat - zcrit * se,
                                     ciu = obs_stat + zcrit * se)

            }


            cil <- ci %>% dplyr::pull(cil)
            ciu <- ci %>% dplyr::pull(ciu)

            simres <- list(reps = reps, obsDiff = obs_stat, cil = cil, ciu = ciu, se = se, zcrit = zcrit)
            return(simres)
        },

        .computeBoots = function(df){

            reps <- self$options$reps
          
            groupLevels <- base::levels(df$group)

            if(self$options$seedBool){
                set.seed(self$options$rngSeed)
            } else {
                set.seed(NULL)
            }


            boots <- df %>%
                infer::specify(dep ~ group) %>%
                infer::generate(reps = reps, type = "bootstrap") %>%
                infer::calculate(stat = "diff in means", order = c(groupLevels[1], groupLevels[2]))

            return(boots)

        },
          .preparePlot = function(boots, dm) {

            bootplot <- self$results$simplot
            dotHist <- self$options$dotHist
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            bootplot$setState(list(df=boots, obs_stat=dm, confLevel = confLevel, ciType = ciType, dotHist=dotHist))

          },
          .bootPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            df <- image$state[["df"]]

            obs_stat = image$state[["obs_stat"]]

            dotHist = image$state[["dotHist"]]

            confLevel <- image$state[["confLevel"]]

            ciType <- image$state[["ciType"]]

            dm <- obs_stat

            boot <- df

            # We've already computed these, but it's fast anyway
            ciList <- private$.computeCI(boot, dm)
            cil <- ciList$cil
            ciu <- ciList$ciu

            if(ciType == "bootperc"){

                caption <- paste0("CI limits (dashed) are ",
                                  round((100-confLevel)/2,1),
                                  "% and ",
                                  round(100 - (100-confLevel)/2,1),
                                  "%",
                                  "\n percentiles of bootstrap differences")

            } else {

                se <- ciList$se
                zcrit <- ciList$zcrit

                caption <- paste0("CI (dashed) is calculated using SE of bootstrap differences",
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
                        dplyr::mutate(x.bin = dm + ( (stat - bw/2 - dm) %/% bw ) * bw)

                    # adjust limits of CI (round) just like we did with other x positions
                    cila <- dm + ( (cil - bw/2 - dm) %/% bw ) * bw
                    ciua <- dm + ( (ciu - bw/2 - dm) %/% bw ) * bw


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
                    ggplot2::xlab("difference (group 1 - group 2)") +
                    ggplot2::coord_equal(ratio = bw*2/3) +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                                   plot.caption = element_text(color = "red", hjust = 0))



            } else {

                p <- ggplot2::ggplot(data=boot, aes(x=stat)) +
                    ggplot2::geom_histogram(center = dm, show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept = cil, linetype = "dashed", color = "red") +
                    ggplot2::geom_vline(xintercept = ciu, linetype = "dashed", color = "red") +
                    ggplot2::theme_minimal() +
                    ggplot2::xlab("difference (group 1 - group 2)") +
                    ggplot2::ylab("count") +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                                   plot.caption = element_text(color = "red", hjust = 0))

            }


            return(p)
          },
          
          .formula=function() {
            jmvcore:::composeFormula(self$options$vars, self$options$group)
          }
        )
      )