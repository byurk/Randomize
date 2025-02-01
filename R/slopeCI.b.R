
#' @importFrom jmvcore .
slopeCIClass <- R6::R6Class(
  "slopeCIClass",
  inherit=slopeCIBase,
  private=list(
    .run=function() {
      
      indepVarName <- self$options$indep
      depVarName <- self$options$dep
      
      varNames <- c(indepVarName, depVarName)
      
      if (is.null(indepVarName) || is.null(depVarName))
      return()
      
      data <- select(self$data, varNames)
      
      data[[depVarName]] <- jmvcore::toNumeric(data[[depVarName]])
      data[[indepVarName]] <- jmvcore::toNumeric(data[[indepVarName]])
      
      CITable <- self$results$CITable
      coefTable <- self$results$coef
      modelfitTable <- self$results$modelfit
      
      if (depVarName == indepVarName)
      jmvcore::reject(.("Independent variable '{a}' must not also be a dependent variable"),
      code="a_is_dependent_variable", a=indepVarName)
      
      # exclude rows with missings
      data <- data %>% tidyr::drop_na()
      
      dataCI <- data.frame(dep=data[[depVarName]], indep=data[[indepVarName]])

      if (is.factor(dataCI$dep))
      res <- createError(.('Dependent variable is not numeric'))
      else if (any(is.infinite(dataCI$dep)))
      res <- createError(.('Dependent variable contains infinite values'))
      else if (is.factor(dataCI$indep))
      res <- createError(.('Independent variable is not numeric'))
      else if (any(is.infinite(dataCI$indep)))
      res <- createError(.('Independent variable contains infinite values'))
      else
        
      lm1 <- lm(dep ~ indep, data=dataCI)
      coef <- lm1$coefficients
      b <- as.numeric(coef[2])
      r2 <- summary(lm1)$r.squared
      r <- cor(dataCI$indep, dataCI$dep)
      
      boots <- private$.computeBoots(dataCI)
      boots <- tidyr::drop_na(boots)
      res <- private$.computeCI(boots, b)
      res <- within(res, rm(se, zcrit))
      private$.preparePlot(boots, b)
        
        if (isError(res)) {
          
          CITable$setRow(rowKey=depVarName, list(
            "reps"='',
            "b"='',
            "cil"='',
            "ciu"=''))
            
            message <- extractErrorMessage(res)
            
            CITable$addFootnote(rowKey=depVarName, 'b', message)
            
          } else {
            
            CITable$setRow(rowKey=depVarName, res)
            }
            
            if (self$options$coef) {
              
              #coefTable$addRow(rowKey="Intercept", list(
              coefTable$setRow(rowKey="1", list(
                "term"="Intercept",
                "est"=as.numeric(coef[1])
              ))

              #coefTable$addRow(rowKey="Slope", list(
              coefTable$setRow(rowKey="2", list(
                "term"=indepVarName,
                "est"=as.numeric(coef[2])
              ))
            }
      
            if (self$options$modelfit) {
              
              modelfitTable$setRow(rowKey=depVarName, list(
                "r"=r,
                "r2"=r2
              ))
            }
            
            if (self$options$plots) {
              
              image <- self$results$linplot$get(key=depVarName)$lp
              
              if (nrow(dataCI) > 0) {
                      
                image$setState(dataCI)
                
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

          .computeCI = function(boots, b) {

            reps <- self$options$reps
            confLevel <- self$options$confLevel/100
            ciType <- self$options$ciType

            obs_stat <- b

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

            simres <- list(reps = reps, b = b, cil = cil, ciu = ciu, se = se, zcrit = zcrit)
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
                infer::specify(dep ~ indep) %>%
                infer::generate(reps = reps, type = "bootstrap") %>%
                infer::calculate(stat = "slope")

            return(boots)

        },
         
          .linplot=function(image, ggtheme, theme, ...) {
            
            if (is.null(image$state))
            return(FALSE)
            
            depName <- self$options$dep
            indepName <- self$options$indep
            
            plot <- ggplot2::ggplot(data=image$state, aes(x=indep, y=dep)) +
            ggplot2::geom_point(alpha = 0.8, size=2.5, shape = 21, color = theme$color[1], fill = theme$color[2]) +
            ggplot2::geom_smooth(method = "lm", se = FALSE, color = theme$color[1]) +
            ggplot2::labs(x=indepName, y=depName) +
            ggtheme
            
            return(plot)
          },
          .sourcifyOption = function(option) {
            if (option$name %in% c('deps', 'group'))
            return('')
            super$.sourcifyOption(option)
          },
          .preparePlot = function(boots, b) {

            bootplot <- self$results$simplot
            dotHist <- self$options$dotHist
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            bootplot$setState(list(df=boots, obs_stat=b, confLevel = confLevel, ciType = ciType, dotHist=dotHist))

        },

        .bootPlot = function(image, ggtheme, theme, ...) {

          if (is.null(image$state))
              return(FALSE)

          df <- image$state[["df"]]

          obs_stat = image$state[["obs_stat"]]

          dotHist = image$state[["dotHist"]]

          confLevel <- image$state[["confLevel"]]

          ciType <- image$state[["ciType"]]

          b <- obs_stat

          boot <- df

          # We've already computed these, but it's fast anyway
          ciList <- private$.computeCI(boot, b)
          cil <- ciList$cil
          ciu <- ciList$ciu

          if(ciType == "bootperc"){

              caption <- paste0("CI limits (dashed) are ",
                                round((100-confLevel)/2,1),
                                "% and ",
                                round(100 - (100-confLevel)/2,1),
                                "%",
                                "\n percentiles of bootstrap slopes")

          } else {

              se <- ciList$se
              zcrit <- ciList$zcrit

              caption <- paste0("CI (dashed) is calculated using SE of bootstrap slopes",
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
                      dplyr::mutate(x.bin = b + ( (stat - bw/2 - b) %/% bw ) * bw)

                  # adjust limits of CI (round) just like we did with other x positions
                  cila <- b + ( (cil - bw/2 - b) %/% bw ) * bw
                  ciua <- b + ( (ciu - bw/2 - b) %/% bw ) * bw


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
                  ggplot2::xlab("slope") +
                  ggplot2::coord_equal(ratio = bw*2/3) +
                  ggplot2::labs(caption = caption) +
                  ggplot2::theme(text = element_text(size = 14),
                                 plot.caption = element_text(color = "red", hjust = 0))



          } else {

              p <- ggplot2::ggplot(data=boot, aes(x=stat)) +
                  ggplot2::geom_histogram(center = b, show.legend = FALSE) +
                  ggplot2::geom_vline(xintercept = cil, linetype = "dashed", color = "red") +
                  ggplot2::geom_vline(xintercept = ciu, linetype = "dashed", color = "red") +
                  ggplot2::theme_minimal() +
                  ggplot2::xlab("slope") +
                  ggplot2::ylab("count") +
                  ggplot2::labs(caption = caption) +
                  ggplot2::theme(text = element_text(size = 14),
                                 plot.caption = element_text(color = "red", hjust = 0))

          }


          return(p)
        },   
          
        .formula=function() {
            jmvcore:::composeFormula(self$options$dep, self$options$indep)
          }
        )
      )