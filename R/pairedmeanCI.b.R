
#' @importFrom jmvcore .
pairedmeanCIClass <- R6::R6Class(
  "pairedmeanCIClass",
  inherit=pairedmeanCIBase,
  private=list(
    .run=function() {
      
      data <- self$data
      
      CITable <- self$results$CITable
      descTable <- self$results$desc
      
      pairs <- self$options$get('pairs')
      
      if(length(pairs) > 0){
        
        
        for (pair in pairs) {
        
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

          dataCI <- data.frame(dif = column1 - column2)

          if (is.factor(column1) | is.factor(column2)) {
            res <- createError(.('One or both variables are not numeric'))
          } else if (any(is.infinite(column1)) | any(is.infinite(column2))) {
            res <- createError(.('One or both variables contain infinite values'))
          } else if (n < 2) {
            res <- createError(.('At least 2 complete pairs are needed (rows with a missing value are dropped)'))
          } else {
            boots <- cached_sims(CITable,
                list(dif = dataCI$dif, reps = self$options$reps, seedBool = self$options$seedBool, rngSeed = self$options$rngSeed),
                function() tidyr::drop_na(private$.computeBoots(dataCI)), slot = paste(name1, name2))
            res <- private$.computeCI(boots, m1-m2)
            res <- within(res, rm(se, zcrit))
            private$.preparePlot(pair, boots, m1-m2)
          }
          
          if (isError(res)) {
            
            CITable$setRow(rowKey=pair, list(
              "reps"='',
              "obsDiff"='',
              "cil"='',
              "ciu"=''))
              
              message <- extractErrorMessage(res)
              
              CITable$addFootnote(rowKey=pair, 'obsDiff', message)
              
            } else {
              
              CITable$setRow(rowKey=pair, list(
                "reps"=self$options$reps,
                "obsDiff"=m1-m2,
                "cil"=res$cil,
                "ciu"=res$ciu))
              }
              
              if (self$options$desc) {
                
                row1Key <- paste(pair$i1, pair$i2, 1)
                row2Key <- paste(pair$i1, pair$i2, 2)
                
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
                    if (nrow(dataCI) > 0) {
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
                  }  # for (pair in pairs)
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
              
              descTable <- self$results$desc
              
              plots <- self$results$descplot
              
              pairs <- self$options$get('pairs')
              
              if(length(pairs) > 0){
                
                for (pair in pairs) {
                
                simtable$setRow(rowKey=pair, list(
                  `var1`=pair$i1,
                  `var2`=pair$i2))
                  
                  row1Key <- paste(pair$i1, pair$i2, 1)
                  row2Key <- paste(pair$i1, pair$i2, 2)
                  descTable$addRow(row1Key)
                  descTable$addRow(row2Key)
                  
                  plots$get(pair)$setTitle(paste0(pair, collapse=' - '))
                  
                    }  # for (pair in pairs)
                }
                
              },

              .computeCI = function(boots, dm) {

                reps <- self$options$reps
                confLevel <- self$options$confLevel
                ciType <- self$options$ciType

                ci <- compute_boot_ci(boots, dm, confLevel, ciType)
                simres <- c(list(reps = reps, obsDiff = dm), ci)
                return(simres)
            },
    
            .computeBoots = function(df){
    
                reps <- self$options$reps
    
                set_seed_if(self$options$seedBool, self$options$rngSeed)

                # infer::specify() runs t.test() internally, which refuses
                # constant differences; every bootstrap mean of a constant
                # sample is that constant, so draw it directly
                if (stats::sd(df$dif) == 0)
                    return(data.frame(replicate = seq_len(reps), stat = rep(df$dif[1], reps)))

                boots <- df %>%
                    infer::specify(response = dif) %>%
                    infer::generate(reps = reps, type = "bootstrap") %>%
                    infer::calculate(stat = "mean")
    
                return(boots)
    
            },

              .descplot=function(image, ggtheme, theme, ...) {
                if (is.null(image$state))
                    return(FALSE)
                plot_desc_stats(image$state, xlab = NULL, ylab = NULL,
                                ggtheme = ggtheme, theme = theme)
              },
              .preparePlot = function(pair, boots, dm) {

                bootplot <- self$results$simplot$get(key=pair)
                bootplot$setTitle(paste(pair$i1, "\u2212", pair$i2))
                dotHist <- self$options$dotHist
                confLevel <- self$options$confLevel
                ciType <- self$options$ciType
    
                bootplot$setState(list(df=strip_infer(boots), obs_stat=dm, confLevel = confLevel, ciType = ciType, dotHist=dotHist, showCounts=self$options$showCounts,
                                          xlab="mean difference", stat_label="bootstrap differences"))
                
              },

              .bootPlot = function(image, ggtheme, theme, ...) {

                if (is.null(image$state))
                    return(FALSE)

                st <- image$state
                p <- plot_boot_dist(st$df, st$obs_stat, st$confLevel, st$ciType,
                                    st$dotHist,
                                    xlab = "mean difference",
                                    stat_label = "bootstrap differences",
                                    show_counts = isTRUE(st$showCounts),
                           plot_width = image$width)
                return(p)
              }
            )
          )