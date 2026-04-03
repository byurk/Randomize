
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

          dataCI <- data.frame(dif = column1 - column2)

          if (is.factor(column1) | is.factor(column2)) {
            res <- createError(.('One or both variables are not numeric'))
          } else if (any(is.infinite(column1)) | any(is.infinite(column2))) {
            res <- createError(.('One or both variables contain infinite values'))
          } else {
            boots <- private$.computeBoots(dataCI)
            boots <- tidyr::drop_na(boots)
            res <- private$.computeCI(boots, m1-m2)
            res <- within(res, rm(se, zcrit))
            private$.preparePlot(boots, m1-m2)
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
                
                pair <- pairs[[1]]
                
                simtable$setRow(rowKey=pair, list(
                  `var1`=pair$i1,
                  `var2`=pair$i2))
                  
                  row1Key <- paste0(pair$i1, 1)
                  row2Key <- paste0(pair$i2, 1)
                  descTable$addRow(row1Key)
                  descTable$addRow(row2Key)
                  
                  plots$get(pair)$setTitle(paste0(pair, collapse=' - '))
                  
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
              .preparePlot = function(boots, dm) {
                
                bootplot <- self$results$simplot
                dotHist <- self$options$dotHist
                confLevel <- self$options$confLevel
                ciType <- self$options$ciType
    
                bootplot$setState(list(df=strip_infer(boots), obs_stat=dm, confLevel = confLevel, ciType = ciType, dotHist=dotHist))
                
              },

              .bootPlot = function(image, ggtheme, theme, ...) {

                if (is.null(image$state))
                    return(FALSE)

                st <- image$state
                p <- plot_boot_dist(st$df, st$obs_stat, st$confLevel, st$ciType,
                                    st$dotHist,
                                    xlab = "mean difference",
                                    stat_label = "bootstrap differences")
                return(p)
              }
            )
          )