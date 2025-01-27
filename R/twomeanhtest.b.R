
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
      res <- try(t.test(dep ~ group, data=dataHTest, var.equal=TRUE,
        alternative=Ha), silent=TRUE)
        
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
              "p"=res$p.value))
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
          .computePerms = function(dataHTest){

            groupLevels <- base::levels(dataHTest$group) # Can use [1] , [2] to get order?
            
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
                infer::calculate(stat = "diff in means", order = c("G1", "G2"))

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
          .formula=function() {
            jmvcore:::composeFormula(self$options$vars, self$options$group)
          }
        )
      )