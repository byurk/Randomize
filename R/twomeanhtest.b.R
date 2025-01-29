
#' @importFrom jmvcore .
twomeanhtestClass <- R6::R6Class(
  "twomeanhtestClass",
  inherit=twomeanhtestBase,
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
      
      htestTable <- self$results$htest
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
      
      ## Hypothesis options checking
      if (self$options$hypothesis == 'oneGreater')
      Ha <- "greater"
      else if (self$options$hypothesis == 'twoGreater')
      Ha <- "less"
      else
      Ha <- "two.sided"
      
      dataHTest <- data.frame(dep=data[[depVarName]], group=data[[groupVarName]])
      
      groupLevels <- base::levels(dataHTest$group)
      v <- tapply(dataHTest$dep, dataHTest$group, function(x) tryNaN(var(x)))
      n <- tapply(dataHTest$dep, dataHTest$group, length)
      m <- tapply(dataHTest$dep, dataHTest$group, function(x) tryNaN(mean(x)))
      med <- tapply(dataHTest$dep, dataHTest$group, function(x) tryNaN(median(x)))
      sd <- sqrt(v)
      
      n[is.na(n)] <- 0
      m[is.na(m)] <- NaN
      med[is.na(med)] <- NaN
      sd[is.na(sd)] <- NaN
      
      
      if (is.factor(dataHTest$dep))
      res <- createError(.('Variable is not numeric'))
      else if (any(is.infinite(dataHTest$dep)))
      res <- createError(.('Variable contains infinite values'))
      else
      #res <- try(t.test(dep ~ group, data=dataHTest, var.equal=TRUE,
      #                  alternative=Ha, conf.level=confInt), silent=TRUE)
      #res <- try(t.test(dep ~ group, data=dataHTest, var.equal=TRUE,
      #  alternative=Ha), silent=TRUE)
      perms <- private$.computePerms(dataHTest)
      res <- private$.computePval(perms, m[1]-m[2])
      private$.preparePlot(perms, m[1]-m[2])
        
        if (isError(res)) {
          
          htestTable$setRow(rowKey=depVarName, list(
            "reps"='',
            "md"='',
            "p"=''))
            
            message <- extractErrorMessage(res)
            
            htestTable$addFootnote(rowKey=depVarName, 'md', message)
            
          } else {
            
            htestTable$setRow(rowKey=depVarName, list(
              "reps"=self$options$reps,
              #"md"=res$estimate[1]-res$estimate[2],
              "md"=m[1]-m[2],
              "p"=res$pval))
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
              
              if (nrow(dataHTest) > 0) {
                
                #ciWidth <- self$options$ciWidth
                #tail <- qnorm(1 - (100 - ciWidth) / 200)
                
                means <- aggregate(dataHTest$dep, by=list(dataHTest$group),
                function(x) tryNaN(mean(x)), simplify=FALSE)
                #cies  <- aggregate(dataHTest$dep, by=list(dataHTest$group),
                #                   function(x) { tail * tryNaN(sd(x)) / sqrt(length(x)) }, simplify=FALSE)
                medians <- aggregate(dataHTest$dep, by=list(dataHTest$group),
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
            
            hypothesis <- self$options$hypothesis
            groupName <- self$options$group
            
            groups <- NULL
            if ( ! is.null(groupName))
            groups <- base::levels(self$data[[groupName]])
            if (length(groups) != 2)
            groups <- c('Group 1', 'Group 2')
            
            table <- self$results$htest
            
            if (hypothesis == 'oneGreater')
            table$setNote("hyp", jmvcore::format("H\u2090 \u03BC\u2009<sub>{}</sub> > \u03BC\u2009<sub>{}</sub>", groups[1], groups[2]))
            else if (hypothesis == 'twoGreater')
            table$setNote("hyp", jmvcore::format("H\u2090 \u03BC\u2009<sub>{}</sub> < \u03BC\u2009<sub>{}</sub>", groups[1], groups[2]))
            else
            table$setNote("hyp", jmvcore::format("H\u2090 \u03BC\u2009<sub>{}</sub> \u2260 \u03BC\u2009<sub>{}</sub>", groups[1], groups[2]))
            
          },
          .computePval = function(perms, dm) {

            reps <- self$options$reps
            alt <- self$options$hypothesis

            if (alt == "oneGreater")
                direction <- "greater"
            else if (alt == "twoGreater")
                direction <- "less"
            else
                direction <- "two_sided"

            pval <- perms %>%
                infer::get_p_value(obs_stat = dm, direction = direction) %>%
                dplyr::pull()

            simres <- list(reps = reps, pval = pval)
            return(simres)
          },
          .computePerms = function(dataHTest){

            groupLevels <- base::levels(dataHTest$group)
            
            reps <- self$options$reps

            if(self$options$seedBool){
                set.seed(self$options$rngSeed)
            } else {
                set.seed(NULL)
            }

            perms <- dataHTest %>%
                infer::specify(dep ~ group) %>%
                infer::hypothesize(null = "independence") %>%
                infer::generate(reps = reps, type = "permute") %>%
                infer::calculate(stat = "diff in means", order = c(groupLevels[1], groupLevels[2]))

            return(perms)

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
          .preparePlot = function(perms, dm) {

            permplot <- self$results$simplot
            dotHist <- self$options$dotHist
            alt <- self$options$hypothesis

            permplot$setState(list(df=perms, obs_stat=dm, alt = alt, dotHist=dotHist))

        },
          .permPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            df <- image$state[["df"]]

            obs_stat = image$state[["obs_stat"]]

            dotHist = image$state[["dotHist"]]

            alt <- image$state[["alt"]]

            dm <- obs_stat

            # which tail of the distro should be used to calculate a p-val?
            if (alt == "twoGreater"){
                ptail <- "lt"
                caption <- "One-sided: p-val is proportion of\n results \u2264 observed value"
            } else if (alt == "oneGreater"){
                ptail <- "rt"
                caption <- "One-sided: p-val is proportion of\n results \u2265 observed value"
            } else {
                ptail <- df %>%
                    dplyr::summarize(lt = mean(stat < dm), rt = mean(stat > dm)) %>%
                    dplyr::mutate(pt = dplyr::if_else(rt < lt, "rt", "lt")) %>%
                    dplyr::pull(pt)
                if (ptail == "lt")
                    caption <- "Two-sided: p-val is 2\u00D7 proportion of\n results \u2264 observed value"
                else
                    caption <- "Two-sided: p-val is 2\u00D7 proportion of\n results \u2265 observed value"
            }

            if(dotHist == "dotplot"){

                ndist <- dplyr::n_distinct(df$stat)

                # bin width, also used for scaling dots when not binning
                bw <- df %>%
                    dplyr::summarize(min = min(stat), max = max(stat)) %>%
                    #mutate(bw = (max - min) / 30 * (1 + 3*.Machine$double.eps)) %>%
                    dplyr::mutate(bw = (max - min) / 30 ) %>%
                    dplyr::pull(bw)


                if (ndist <= 30){ # don't bin unless there are more than 30 unique values

                    df <- df %>%
                        dplyr::mutate(x.bin = stat)

                } else {

                    if(ptail == "lt"){

                        df <- df %>%
                            dplyr::mutate(x.bin = dm - ( (dm - stat) %/% bw ) * bw)

                    } else {

                        df <- df %>%
                            dplyr::mutate(x.bin = dm + ( (stat - dm) %/% bw ) * bw)

                    }


                }


                # color the dots in the tail that will be used to calculate the p-val
                # put vertical line at value of dm

                df <- df %>%
                    dplyr::mutate(extreme = (ptail == "lt" & stat <= dm) | (ptail == "rt" & stat >= dm))


                # dotplot based on tjebo answer from:
                # https://stackoverflow.com/questions/53697235/ggplot-dotplot-what-is-the-proper-use-of-geom-dotplot


                df <- df %>%
                    dplyr::group_by(x.bin) %>%
                    dplyr::mutate(y = seq_along(x.bin))

                labht <- df %>%
                    dplyr::ungroup() %>%
                    dplyr::summarize(labht = max(y)) %>%
                    dplyr::pull()

                lab_ht <- max(labht - 1, 3)

                bigx <- df %>%
                    dplyr::ungroup() %>%
                    dplyr::summarize(bigx = max(x.bin)) %>%
                    dplyr::pull()

                littlex <- df %>%
                    dplyr::ungroup() %>%
                    dplyr::summarize(littlex = min(x.bin)) %>%
                    dplyr::pull()

                # do not show legend
                p <- ggplot2::ggplot(df) +
                    ggforce::geom_ellipse(aes(x0 = x.bin, y0 = y, a = bw/3, b = 0.5, angle = 0, fill = extreme, color = extreme),
                                          show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept = dm, linetype = "dashed", color = "red") +
                    ggplot2::scale_fill_manual(values = c("FALSE" = "black", "TRUE" = "#ff8c8c")) +
                    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#ff8c8c")) +
                    ggplot2::annotate("text", x = dm, y = lab_ht, label = "Observed\nDifference", color = "red") +
                    ggplot2::theme_minimal() +
                    ggplot2::ylab("count") +
                    ggplot2::xlab("difference (group 1 - group 2)") +
                    ggplot2::ylim(0, lab_ht + 3) +
                    ggplot2::xlim(min(littlex - bw, dm - bw), max(bigx + bw, dm + bw)) +
                    ggplot2::coord_equal(ratio = bw*2/3) +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                                   plot.caption = element_text(color = "red", hjust = 0))



            } else {

                closed <- dplyr::if_else(ptail == "lt", "right", "left")

                df <- df %>%
                    dplyr::mutate(extreme = (ptail == "lt" & stat <= dm) | (ptail == "rt" & stat >= dm))

                p <- ggplot2::ggplot(data=df, aes(x=stat, fill = extreme)) +
                    ggplot2::geom_histogram(boundary = dm, closed = closed, show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept=obs_stat, linetype='dashed', color = "red") +
                    ggplot2::scale_fill_manual(values = c("FALSE" = "black", "TRUE" = "#ff8c8c")) +
                    ggplot2::theme_minimal() +
                    ggplot2::xlab("difference (group 1 - group 2)") +
                    ggplot2::ylab("count") +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                                   plot.caption = element_text(color = "red", hjust = 0))

                yMax <- ggplot2::layer_scales(p)$y$range$range[2]  # upper y-limit
                p <- p + ggplot2::annotate("text", x = dm, y = yMax, vjust = "top", label = "Observed\nDifference", color = "red")
            }


            return(p)
        },
          .formula=function() {
            jmvcore:::composeFormula(self$options$vars, self$options$group)
          }
        )
      )