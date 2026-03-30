
#' @importFrom jmvcore .
pairedmeanhtestClass <- R6::R6Class(
  "pairedmeanhtestClass",
  inherit=pairedmeanhtestBase,
  private=list(
    .run=function() {
      
      data <- self$data
      
      htestTable <- self$results$htest
      descTable <- self$results$desc
      
      pairs <- self$options$get('pairs')
      
      if(length(pairs) > 0){
        
        
        pair <- pairs[[1]]
        
        if(!any(sapply(pair,length)== 0)){
          
          ## Hypothesis options checking
          if (self$options$hypothesis == 'oneGreater')
          Ha <- "greater"
          else if (self$options$hypothesis == 'twoGreater')
          Ha <- "less"
          else
          Ha <- "two.sided"
          
          name1 <- pair$i1
          name2 <- pair$i2
          
          
          
          data[[name1]] <- jmvcore::toNumeric(data[[name1]])
          data[[name2]] <- jmvcore::toNumeric(data[[name2]])

          # Remove rows where either variable has NA before computing stats
          complete <- complete.cases(data[[name1]], data[[name2]])
          column1 <- data[[name1]][complete]
          column2 <- data[[name2]][complete]

          n    <- length(column1)
          m1   <- tryNaN(mean(column1))
          m2   <- tryNaN(mean(column2))
          med1 <- tryNaN(median(column1))
          med2 <- tryNaN(median(column2))
          sd1  <- tryNaN(stats::sd(column1))
          sd2  <- tryNaN(stats::sd(column2))

          dataHTest <- data.frame(dif = column1 - column2)

          if (is.factor(column1) | is.factor(column2))
          res <- createError(.('One or both variables are not numeric'))
          else if (any(is.infinite(column1)) | any(is.infinite(column2)))
          res <- createError(.('One or both variables contain infinite values'))
          else
          #res <- try(t.test(dep ~ group, data=dataHTest, var.equal=TRUE,
          #                  alternative=Ha, conf.level=confInt), silent=TRUE)
          #res <- try(t.test(dep ~ group, data=dataHTest, var.equal=TRUE,
          #  alternative=Ha), silent=TRUE)
          perms <- private$.computePerms(dataHTest)
          res <- private$.computePval(perms, m1-m2)
          private$.preparePlot(perms, m1-m2, res$direction)
          
          if (isError(res)) {
            
            htestTable$setRow(rowKey=pair, list(
              "reps"='',
              "md"='',
              "p"=''))
              
              message <- extractErrorMessage(res)
              
              htestTable$addFootnote(rowKey=pair, 'md', message)
              
            } else {
              
              htestTable$setRow(rowKey=pair, list(
                "reps"=self$options$reps,
                #"md"=res$estimate[1]-res$estimate[2],
                "md"=m1-m2,
                "p"=res$pval))
              }
              
              if (self$options$desc) {
                
                row1Key <- paste0(pair$i1, 1)
                row2Key <- paste0(pair$i2, 1)
                
                descTable$setRow(rowKey=row1Key, list(
                  "name"=name1,
                  "num"=n,
                  "m"=m1,
                  "med"=med1,
                  "sd"=sd1))
                  
                  descTable$setRow(rowKey=row2Key, list(
                    "name"=name2,
                    "num"=n,
                    "m"=m2,
                    "med"=med2,
                    "sd"=sd2))
                    
                    descTable$addFormat(col='name', rowKey=row1Key, Cell.BEGIN_GROUP)
                    descTable$addFormat(col='name', rowKey=row2Key, Cell.END_GROUP)
                    
                  }
                  
                  if (self$options$plots) {
                    
                    image <- self$results$descplot$get(key=pair)$desc
                    
                    if (nrow(dataHTest) > 0) {
                      
                      #ciWidth <- self$options$ciWidth
                      #tail <- qnorm(1 - (100 - ciWidth) / 200)
                      
                      means <- c(tryNaN(mean(column1)), tryNaN(mean(column2)))
                      #cies  <- aggregate(dataHTest$dep, by=list(dataHTest$group),
                      #                   function(x) { tail * tryNaN(sd(x)) / sqrt(length(x)) }, simplify=FALSE)
                      medians <- c(tryNaN(median(column1)), tryNaN(median(column2)))
                      
                      meanPlotData <- data.frame(group=c(name1, name2))
                      meanPlotData <- cbind(meanPlotData, stat=means)
                      #meanPlotData <- cbind(meanPlotData, cie=unlist(cies$x))
                      meanPlotData <- cbind(meanPlotData, cie=NA)
                      meanPlotData <- cbind(meanPlotData, type='mean')
                      
                      medianPlotData <- data.frame(group=c(name1, name2))
                      medianPlotData <- cbind(medianPlotData, stat=medians)
                      medianPlotData <- cbind(medianPlotData, cie=NA)
                      medianPlotData <- cbind(medianPlotData, type='median')
                      
                      plotData <- rbind(meanPlotData, medianPlotData)
                      plotData$group <- factor(plotData$group, levels=unique(plotData$group))
                      
                      if (all(is.na(plotData$stat)))
                      image$setState(NULL)
                      else
                      image$setState(plotData)
                      
                    } else {
                      
                      image$setState(NULL)
                    }
                  }
                }
              }
              
            },
            .init=function() {
              
              hypothesis <- self$options$hypothesis
              
              table <- self$results$htest
              
              if (hypothesis == 'oneGreater') {
                table$setNote(
                  "hyp",
                  jmvcore::format("H\u2090 \u03BC\u2009<sub>{}</sub> > 0", .("Measure 1 - Measure 2"))
                )
              } else if (hypothesis == 'twoGreater') {
                table$setNote(
                  "hyp",
                  jmvcore::format("H\u2090 \u03BC\u2009<sub>{}</sub> < 0", .("Measure 1 - Measure 2"))
                )
              } else {
                table$setNote(
                  "hyp",
                  jmvcore::format("H\u2090 \u03BC\u2009<sub>{}</sub> \u2260 0", .("Measure 1 - Measure 2"))
                )
              }
              
              descTable <- self$results$desc
              
              plots <- self$results$descplot
              
              pairs <- self$options$get('pairs')
              
              if(length(pairs) > 0){
                
                pair <- pairs[[1]]
                
                table$setRow(rowKey=pair, list(
                  `var1`=pair$i1,
                  `var2`=pair$i2))
                  
                  row1Key <- paste0(pair$i1, 1)
                  row2Key <- paste0(pair$i2, 1)
                  descTable$addRow(row1Key)
                  descTable$addRow(row2Key)
                  
                  plots$get(pair)$setTitle(paste0(pair, collapse=' - '))
                  
                }
                
              },
              .computePval = function(perms, dm) {

                reps <- self$options$reps
                direction <- map_direction(self$options$hypothesis)
                pval <- compute_null_pval(perms, dm, direction)

                simres <- list(reps = reps, pval = pval, direction = direction)
                return(simres)
              },
              .computePerms = function(dataHTest){
                
                reps <- self$options$reps
                
                set_seed_if(self$options$seedBool, self$options$rngSeed)
                
                perms <- dataHTest %>%
                infer::specify(response = dif) %>%
                infer::hypothesize(null = "paired independence") %>%
                infer::generate(reps = reps, type = "permute") %>%
                infer::calculate(stat = "mean")
                
                return(perms)
                
              },
              .descplot=function(image, ggtheme, theme, ...) {
                
                if (is.null(image$state))
                return(FALSE)
                
                groupName <- self$options$get('group')
                
                #ciw <- self$options$ciWidth
                
                pd <- position_dodge(0.2)
                
                plot <- ggplot(data=image$state, aes(x=group, y=stat, shape=type)) +
                #geom_errorbar(aes(x=group, ymin=stat-cie, ymax=stat+cie, width=.1),
                #              size=.8, colour=theme$color[2], position=pd) +
                geom_point(aes(x=group, y=stat, shape=type), color=theme$color[1],
                fill=theme$fill[1], size=3, position=pd) +
                labs(x=groupName, y=NULL) +
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
              .preparePlot = function(perms, dm, direction) {
                permplot <- self$results$simplot
                dotHist <- self$options$dotHist
                permplot$setState(list(df=perms, obs_stat=dm, direction=direction, dotHist=dotHist))
              },
              .permPlot = function(image, ggtheme, theme, ...) {
                if (is.null(image$state))
                    return(FALSE)
                st <- image$state
                plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                               xlab = "mean difference",
                               obs_label = "Observed\nDifference")
              }
            )
          )