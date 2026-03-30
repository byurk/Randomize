
#' @importFrom jmvcore .
twomeanCIClass <- R6::R6Class(
  "twomeanCIClass",
  inherit=twomeanCIBase,
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
      
      CITable <- self$results$CITable
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
      
      dataCI <- data.frame(dep=data[[depVarName]], group=data[[groupVarName]])
      dataCI <- tidyr::drop_na(dataCI)

      groupLevels <- base::levels(dataCI$group)
      v <- tapply(dataCI$dep, dataCI$group, function(x) tryNaN(var(x)))
      n <- tapply(dataCI$dep, dataCI$group, length)
      m <- tapply(dataCI$dep, dataCI$group, function(x) tryNaN(mean(x)))
      med <- tapply(dataCI$dep, dataCI$group, function(x) tryNaN(median(x)))
      sd <- sqrt(v)
      
      n[is.na(n)] <- 0
      m[is.na(m)] <- NaN
      med[is.na(med)] <- NaN
      sd[is.na(sd)] <- NaN
      
      
      if (is.factor(dataCI$dep))
      res <- createError(.('Variable is not numeric'))
      else if (any(is.infinite(dataCI$dep)))
      res <- createError(.('Variable contains infinite values'))
      else
      #res <- try(t.test(dep ~ group, data=dataCI, var.equal=TRUE,
      #                  alternative=Ha, conf.level=confInt), silent=TRUE)
      #res <- try(t.test(dep ~ group, data=dataCI, var.equal=TRUE,
      #  alternative=Ha), silent=TRUE)
      boots <- private$.computeBoots(dataCI)
      boots <- tidyr::drop_na(boots)
      res <- private$.computeCI(boots, m[1]-m[2])
      res <- within(res, rm(se, zcrit))
      private$.preparePlot(boots, m[1]-m[2])
        
        if (isError(res)) {
          
          CITable$setRow(rowKey=depVarName, list(
            "reps"='',
            "obsDiff"='',
            "cil"='',
            "ciu"=''))
            
            message <- extractErrorMessage(res)
            
            CITable$addFootnote(rowKey=depVarName, 'md', message)
            
          } else {
            
            CITable$setRow(rowKey=depVarName, res)
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
              if (nrow(dataCI) > 0) {
                plotData <- build_desc_plot_data(dataCI$dep, dataCI$group)
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
          
          .descplot=function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            groupName <- self$options$group
            plot_desc_stats(image$state, xlab = groupName, ylab = image$key,
                            ggtheme = ggtheme, theme = theme)
          },
          .sourcifyOption = function(option) {
            if (option$name %in% c('deps', 'group'))
            return('')
            super$.sourcifyOption(option)
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
          
            groupLevels <- base::levels(df$group)

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            boots <- df %>%
                infer::specify(dep ~ group) %>%
                infer::generate(reps = reps, type = "bootstrap") %>%
                infer::calculate(stat = "diff in means", order = c(groupLevels[1], groupLevels[2]))

            return(boots)

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
                                xlab = "difference (group 1 - group 2)",
                                stat_label = "bootstrap differences")
            return(p)
          },
          
          .formula=function() {
            jmvcore:::composeFormula(self$options$vars, self$options$group)
          }
        )
      )