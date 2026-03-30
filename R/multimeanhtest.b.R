
#' @importFrom jmvcore .
multimeanhtestClass <- R6::R6Class(
  "multimeanhtestClass",
  inherit=multimeanhtestBase,
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
      
      if (length(groupLevels) < 3)
      jmvcore::reject(.("Grouping variable '{a}' must have at least 3 levels"),
      code="grouping_var_must_have_3_levels", a=groupVarName)
      
      dataHTest <- data.frame(dep=data[[depVarName]], group=data[[groupVarName]])
      dataHTest <- tidyr::drop_na(dataHTest)

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

      Fobs <- dataHTest %>%
        infer::specify(dep ~ group) %>%
        infer::hypothesize(null = "independence") %>%
        infer::calculate("F") %>%
        pull()

      perms <- private$.computePerms(dataHTest)
      res <- private$.computePval(perms, Fobs)
      private$.preparePlot(perms, Fobs)
        
        if (isError(res)) {
          
          htestTable$setRow(rowKey=depVarName, list(
            "reps"='',
            "oF"='',
            "p"=''))
            
            message <- extractErrorMessage(res)
            
            htestTable$addFootnote(rowKey=depVarName, 'oF', message)
            
          } else {
            
            htestTable$setRow(rowKey=depVarName, list(
              "reps"=self$options$reps,
              #"md"=res$estimate[1]-res$estimate[2],
              "oF"=Fobs,
              "p"=res$pval))
            }
            
            if (self$options$desc) {

              desc <- tapply(dataHTest$dep, dataHTest$group, function (x) {
                n <- length(x)
                mean <- mean(x)
                median <- median(x)
                sd <- sd(x)
                
                return(c(n=n, mean=mean, median = median, sd=sd))
            })

            for( level in groupLevels){

              row <- list(
                "num" = as.numeric(desc[[level]]['n']),
                "mean" = as.numeric(desc[[level]]['mean']),
                "median" = as.numeric(desc[[level]]['median']),
                "sd" = as.numeric(desc[[level]]['sd'])
            )

            descTable$setRow(rowKey=paste0(depVarName,level), row)

            }

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
            
            groupName <- self$options$group
            
            #groups <- NULL
            #if ( ! is.null(groupName))
            #groups <- base::levels(self$data[[groupName]])
            #if (length(groups) != 2)
            #groups <- c('Group 1', 'Group 2')
            
            #table <- self$results$htest

            table <- self$results$desc

            group <- self$options$group

            if (is.null(group))
                return()

            levels <- levels(self$data[[group]])

            table$getColumn('group')$setTitle(group)

            depVarName <- self$options$vars

                for (i in seq_along(levels)) {
                    table$addRow(paste0(depVarName,levels[i]), list(dep=depVarName, group=levels[i]))

                    if (i == 1)
                        table$addFormat(rowKey=paste0(depVarName,levels[i]), col=1, jmvcore::Cell.BEGIN_GROUP)
                }
            
            
          },
          .computePval = function(perms, oF) {
            reps <- self$options$reps
            pval <- compute_null_pval(perms, oF, "greater")
            list(reps = reps, pval = pval)
          },
          .computePerms = function(dataHTest){
            
            reps <- self$options$reps

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            perms <- dataHTest %>%
                infer::specify(dep ~ group) %>%
                infer::hypothesize(null = "independence") %>%
                infer::generate(reps = reps, type = "permute") %>%
                infer::calculate(stat = "F")

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
          .preparePlot = function(perms, oF) {
            permplot <- self$results$simplot
            dotHist <- self$options$dotHist
            permplot$setState(list(df=perms, obs_stat=oF, direction="greater", dotHist=dotHist))
        },
          .permPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            st <- image$state
            plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                           xlab = "F",
                           obs_label = "Observed\nF")
        },
          .formula=function() {
            jmvcore:::composeFormula(self$options$vars, self$options$group)
          }
        )
      )