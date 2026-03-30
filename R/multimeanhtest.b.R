
#' @importFrom jmvcore .
multimeanhtestClass <- R6::R6Class(
  "multimeanhtestClass",
  inherit=multimeanhtestBase,
  private=list(
    .run=function() {
      
      groupVarName <- self$options$group
      depVarName <- self$options$vars
      varNames <- c(groupVarName, depVarName)
      
      if (is.null(groupVarName) || is.null(depVarName))
      return()
      
      data <- dplyr::select(self$data, dplyr::all_of(varNames))
      
      data[[depVarName]] <- jmvcore::toNumeric(data[[depVarName]])
      data[[groupVarName]] <- droplevels(as.factor(data[[groupVarName]]))
      
      htestTable <- self$results$htest
      descTable <- self$results$desc
      
      if (depVarName == groupVarName)
      jmvcore::reject(.("Grouping variable '{a}' must not also be a dependent variable"),
      code="a_is_dependent_variable", a=groupVarName)
      
      # exclude rows with missings in the grouping variable
      data <- data[ ! is.na(data[[groupVarName]]),]
      
      groupLevels <- base::levels(data[[groupVarName]])
      
      if (length(groupLevels) < 3)
      jmvcore::reject(.("Grouping variable '{a}' must have at least 3 levels"),
      code="grouping_var_must_have_3_levels", a=groupVarName)
      
      dataHTest <- data.frame(dep=data[[depVarName]], group=data[[groupVarName]])
      dataHTest <- tidyr::drop_na(dataHTest)

      groupLevels <- base::levels(dataHTest$group)
      v <- tapply(dataHTest$dep, dataHTest$group, function(x) tryNaN(var(x)))
      n <- tapply(dataHTest$dep, dataHTest$group, length)
      m <- tapply(dataHTest$dep, dataHTest$group, function(x) tryNaN(mean(x)))
      med <- tapply(dataHTest$dep, dataHTest$group, function(x) tryNaN(median(x)))
      sd <- sqrt(v)
      
      n[is.na(n)] <- 0
      m[is.na(m)] <- NaN
      med[is.na(med)] <- NaN
      sd[is.na(sd)] <- NaN
      
      
      if (is.factor(dataHTest$dep))
      res <- createError(.('Variable is not numeric'))
      else if (any(is.infinite(dataHTest$dep)))
      res <- createError(.('Variable contains infinite values'))
      else
      Fobs <- dataHTest %>%
        infer::specify(dep ~ group) %>%
        infer::hypothesize(null = "independence") %>%
        infer::calculate("F") %>%
        pull()

      perms <- private$.computePerms(dataHTest)
      res <- private$.computePval(perms, Fobs)
      private$.preparePlot(perms, Fobs)
        
        if (isError(res)) {
          
          htestTable$setRow(rowKey=depVarName, list(
            "reps"='',
            "oF"='',
            "p"=''))
            
            message <- extractErrorMessage(res)
            
            htestTable$addFootnote(rowKey=depVarName, 'oF', message)
            
          } else {
            
            htestTable$setRow(rowKey=depVarName, list(
              "reps"=self$options$reps,
              "oF"=Fobs,
              "p"=res$pval))
            }
            
            if (self$options$desc) {

              desc <- tapply(dataHTest$dep, dataHTest$group, function (x) {
                n <- length(x)
                mean <- mean(x)
                median <- median(x)
                sd <- sd(x)
                
                return(c(n=n, mean=mean, median = median, sd=sd))
            })

            for( level in groupLevels){

              row <- list(
                "num" = as.numeric(desc[[level]]['n']),
                "mean" = as.numeric(desc[[level]]['mean']),
                "median" = as.numeric(desc[[level]]['median']),
                "sd" = as.numeric(desc[[level]]['sd'])
            )

            descTable$setRow(rowKey=paste0(depVarName,level), row)

            }

            }
            
            if (self$options$plots) {
              image <- self$results$descplot$get(key=depVarName)$desc
              if (nrow(dataHTest) > 0) {
                plotData <- build_desc_plot_data(dataHTest$dep, dataHTest$group)
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
            
            groupName <- self$options$group
            
            table <- self$results$desc

            group <- self$options$group

            if (is.null(group))
                return()

            levels <- levels(self$data[[group]])

            table$getColumn('group')$setTitle(group)

            depVarName <- self$options$vars

                for (i in seq_along(levels)) {
                    table$addRow(paste0(depVarName,levels[i]), list(dep=depVarName, group=levels[i]))

                    if (i == 1)
                        table$addFormat(rowKey=paste0(depVarName,levels[i]), col=1, jmvcore::Cell.BEGIN_GROUP)
                }
            
            
          },
          .computePval = function(perms, oF) {
            reps <- self$options$reps
            pval <- compute_null_pval(perms, oF, "greater")
            list(reps = reps, pval = pval)
          },
          .computePerms = function(dataHTest){
            
            reps <- self$options$reps

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            perms <- dataHTest %>%
                infer::specify(dep ~ group) %>%
                infer::hypothesize(null = "independence") %>%
                infer::generate(reps = reps, type = "permute") %>%
                infer::calculate(stat = "F")

            return(perms)

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
          .preparePlot = function(perms, oF) {
            permplot <- self$results$simplot
            dotHist <- self$options$dotHist
            permplot$setState(list(df=perms, obs_stat=oF, direction="greater", dotHist=dotHist))
        },
          .permPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            st <- image$state
            plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                           xlab = "F",
                           obs_label = "Observed\nF")
        },
          .formula=function() {
            jmvcore:::composeFormula(self$options$vars, self$options$group)
          }
        )
      )