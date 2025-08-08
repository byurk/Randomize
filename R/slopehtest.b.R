
#' @importFrom jmvcore .
slopehtestClass <- R6::R6Class(
  "slopehtestClass",
  inherit=slopehtestBase,
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
      
      htestTable <- self$results$htest
      coefTable <- self$results$coef
      modelfitTable <- self$results$modelfit
      
      if (depVarName == indepVarName)
      jmvcore::reject(.("Independent variable '{a}' must not also be a dependent variable"),
      code="a_is_dependent_variable", a=indepVarName)
      
      # exclude rows with missings
      data <- data %>% tidyr::drop_na()
      
      ## Hypothesis options checking
      if (self$options$hypothesis == 'greater')
      Ha <- "greater"
      else if (self$options$hypothesis == 'less')
      Ha <- "less"
      else
      Ha <- "two.sided"
      
      dataHTest <- data.frame(dep=data[[depVarName]], indep=data[[indepVarName]])

      if (is.factor(dataHTest$dep))
      res <- createError(.('Dependent variable is not numeric'))
      else if (any(is.infinite(dataHTest$dep)))
      res <- createError(.('Dependent variable contains infinite values'))
      else if (is.factor(dataHTest$indep))
      res <- createError(.('Independent variable is not numeric'))
      else if (any(is.infinite(dataHTest$indep)))
      res <- createError(.('Independent variable contains infinite values'))
      else
        
      lm1 <- lm(dep ~ indep, data=dataHTest)
      coef <- lm1$coefficients
      b <- as.numeric(coef[2])
      r2 <- summary(lm1)$r.squared
      r <- cor(dataHTest$indep, dataHTest$dep)
      
      perms <- private$.computePerms(dataHTest)
      res <- private$.computePval(perms, b)
      private$.preparePlot(perms, b, res$direction)
        
        if (isError(res)) {
          
          htestTable$setRow(rowKey=depVarName, list(
            "reps"='',
            "b"='',
            "p"=''))
            
            message <- extractErrorMessage(res)
            
            htestTable$addFootnote(rowKey=depVarName, 'b', message)
            
          } else {
            
            htestTable$setRow(rowKey=depVarName, list(
              "reps"=self$options$reps,
              #"md"=res$estimate[1]-res$estimate[2],
              "b"=b,
              "p"=res$pval))
            }
            
            if (self$options$coef) {
              
              coefTable$setRow(rowKey="1", list(
                "term"="Intercept",
                "est"=as.numeric(coef[1])
              ))

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
              
              if (nrow(dataHTest) > 0) {
                      
                image$setState(dataHTest)
                
              } else {
                
                image$setState(NULL)
              }
            }
    },
          .init=function() {
            
            hypothesis <- self$options$hypothesis
            
            table <- self$results$htest
            
            if (hypothesis == 'greater')
            table$setNote("hyp","H\u2090: \u03B2\u2009 > 0")
            else if (hypothesis == 'less')
            table$setNote("hyp", "H\u2090: \u03B2\u2009 < 0")
            else
            table$setNote("hyp","H\u2090: \u03B2\u2009 \u2260 0")

          },
          .computePval = function(perms, b) {

            alt <- self$options$hypothesis

            if (alt == "greater")
                direction <- "greater"
            else if (alt == "less")
                direction <- "less"
            else
                direction <- "two_sided"

            pval <- compute_perm_pval(perms, b, direction)

            list(reps = self$options$reps, pval = pval, direction = direction)
          },
          .computePerms = function(dataHTest){

            reps <- self$options$reps

            if(self$options$seedBool){
                set.seed(self$options$rngSeed)
            } else {
                set.seed(NULL)
            }

            perms <- dataHTest %>%
                infer::specify(dep ~ indep) %>%
                infer::hypothesize(null = "independence") %>%
                infer::generate(reps = reps, type = "permute") %>%
                infer::calculate(stat = "slope")

            return(perms)

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
          .preparePlot = function(perms, b, direction) {

            permplot <- self$results$simplot
            dotHist <- self$options$dotHist

            permplot$setState(list(df=perms, obs_stat=b, direction=direction, dotHist=dotHist))

        },
                    .permPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            df <- image$state[["df"]]
            obs_stat <- image$state[["obs_stat"]]
            dotHist <- image$state[["dotHist"]]
            direction <- image$state[["direction"]]

            plot_perm_dist(df, obs_stat, direction, dotHist,
                           xlab = "slope", obs_label = "Observed\nSlope")
        }
          .formula=function() {
            jmvcore:::composeFormula(self$options$dep, self$options$indep)
          }
        )
      )
