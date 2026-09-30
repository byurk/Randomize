
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

          if (is.factor(column1) | is.factor(column2)) {
            res <- createError(.('One or both variables are not numeric'))
          } else if (any(is.infinite(column1)) | any(is.infinite(column2))) {
            res <- createError(.('One or both variables contain infinite values'))
          } else {
            perms <- private$.computePerms(dataHTest)
            res <- private$.computePval(perms, m1-m2)
            private$.preparePlot(perms, m1-m2, res$direction)
          }
          
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
                "md"=m1-m2,
                "p"=format_sim_pval(res$pval, self$options$reps)))
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
                      plotData <- build_paired_desc_plot_data(column1, column2, name1, name2)
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
                plot_desc_stats(image$state, xlab = NULL, ylab = NULL,
                                ggtheme = ggtheme, theme = theme)
              },
              .preparePlot = function(perms, dm, direction) {
                permplot <- self$results$simplot
                dotHist <- self$options$dotHist
                permplot$setState(list(df=strip_infer(perms), obs_stat=dm, direction=direction, dotHist=dotHist, showCounts=self$options$showCounts,
                                          xlab="mean difference", obs_label="Observed\nDifference"))
              },
              .permPlot = function(image, ggtheme, theme, ...) {
                if (is.null(image$state))
                    return(FALSE)
                st <- image$state
                plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                               xlab = "mean difference",
                               obs_label = "Observed\nDifference",
                               show_counts = isTRUE(st$showCounts),
                           plot_width = image$width)
              }
            )
          )