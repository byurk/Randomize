
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
            confLevel <- self$options$confLevel
            ciType <- self$options$ciType

            ci <- compute_boot_ci(boots, b, confLevel, ciType)
            simres <- c(list(reps = reps, b = b), ci)
            return(simres)
        },

          .computeBoots = function(df){

            reps <- self$options$reps

            set_seed_if(self$options$seedBool, self$options$rngSeed)


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

          st <- image$state
          p <- plot_boot_dist(st$df, st$obs_stat, st$confLevel, st$ciType,
                              st$dotHist,
                              xlab = "slope",
                              stat_label = "bootstrap slopes")
          return(p)
        },
          
        .formula=function() {
            jmvcore:::composeFormula(self$options$dep, self$options$indep)
          }
        )
      )