
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
          
          column1 <- data[[name1]]
          column2 <- data[[name2]]
          
          var1 <- tryNaN(var(column1))
          var2 <- tryNaN(var(column2))
          n    <- length(column1)
          m1   <- tryNaN(mean(column1))
          m2   <- tryNaN(mean(column2))
          med1 <- tryNaN(median(column1))
          med2 <- tryNaN(median(column2))
          sd1  <- tryNaN(stats::sd(column1))
          sd2  <- tryNaN(stats::sd(column2))
          
          dataCI <- data.frame(dif = column1 - column2)
          dataCI <- tidyr::drop_na(dataCI)
          
          if (is.factor(column1) | is.factor(column2))
          res <- createError(.('One or both variables are not numeric'))
          else if (any(is.infinite(column1)) | any(is.infinite(column2)))
          res <- createError(.('ne or both variables contain infinite values'))
          else
          #res <- try(t.test(dep ~ group, data=dataCI, var.equal=TRUE,
          #                  alternative=Ha, conf.level=confInt), silent=TRUE)
          #res <- try(t.test(dep ~ group, data=dataCI, var.equal=TRUE,
          #  alternative=Ha), silent=TRUE)
          boots <- private$.computeBoots(dataCI)
          boots <- tidyr::drop_na(boots)
          res <- private$.computeCI(boots, m1-m2)
          res <- within(res, rm(se, zcrit))
          private$.preparePlot(boots, m1-m2)
          
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
                #"md"=res$estimate[1]-res$estimate[2],
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
                      
                      #ciWidth <- self$options$ciWidth
                      #tail <- qnorm(1 - (100 - ciWidth) / 200)
                      
                      means <- c(tryNaN(mean(column1)), tryNaN(mean(column2)))
                      #cies  <- aggregate(dataCI$dep, by=list(dataCI$group),
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
    
                if(self$options$seedBool){
                    set.seed(self$options$rngSeed)
                } else {
                    set.seed(NULL)
                }
    
    
                boots <- df %>%
                    infer::specify(response = dif) %>%
                    infer::generate(reps = reps, type = "bootstrap") %>%
                    infer::calculate(stat = "mean")
    
                return(boots)
    
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

                st <- image$state
                p <- plot_boot_dist(st$df, st$obs_stat, st$confLevel, st$ciType,
                                    st$dotHist,
                                    xlab = "mean difference",
                                    stat_label = "bootstrap differences")
                return(p)
              }
            )
          )