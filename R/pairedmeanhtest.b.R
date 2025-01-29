
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
          
          dataHTest <- data.frame(dif = column1 - column2)
          dataHTest <- tidyr::drop_na(dataHTest)
          
          if (is.factor(column1) | is.factor(column2))
          res <- createError(.('One or both variables are not numeric'))
          else if (any(is.infinite(column1)) | any(is.infinite(column2)))
          res <- createError(.('ne or both variables contain infinite values'))
          else
          #res <- try(t.test(dep ~ group, data=dataHTest, var.equal=TRUE,
          #                  alternative=Ha, conf.level=confInt), silent=TRUE)
          #res <- try(t.test(dep ~ group, data=dataHTest, var.equal=TRUE,
          #  alternative=Ha), silent=TRUE)
          perms <- private$.computePerms(dataHTest)
          res <- private$.computePval(perms, m1-m2)
          private$.preparePlot(perms, m1-m2)
          
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
                
                reps <- self$options$reps
                
                if(self$options$seedBool){
                  set.seed(self$options$rngSeed)
                } else {
                  set.seed(NULL)
                }
                
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
                  ggplot2::xlab("mean difference") +
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
                  ggplot2::xlab("mean difference") +
                  ggplot2::ylab("count") +
                  ggplot2::labs(caption = caption) +
                  ggplot2::theme(text = element_text(size = 14),
                  plot.caption = element_text(color = "red", hjust = 0))
                  
                  yMax <- ggplot2::layer_scales(p)$y$range$range[2]  # upper y-limit
                  p <- p + ggplot2::annotate("text", x = dm, y = yMax, vjust = "top", label = "Observed\nDifference", color = "red")
                }
                
                
                return(p)
              }
            )
          )