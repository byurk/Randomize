
#' @importFrom jmvcore .
twomeanhtestClass <- R6::R6Class(
  "twomeanhtestClass",
  inherit=twomeanhtestBase,
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
      
      if (length(groupLevels) != 2)
      jmvcore::reject(.("Grouping variable '{a}' must have exactly 2 levels"),
      code="grouping_var_must_have_2_levels", a=groupVarName)
      
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
      
      
      if (is.factor(dataHTest$dep)) {
        res <- createError(.('Variable is not numeric'))
      } else if (any(is.infinite(dataHTest$dep))) {
        res <- createError(.('Variable contains infinite values'))
      } else if (any(n == 0)) {
        # every value of one group missing: the level survives droplevels()
        # (it is applied before the NA rows go) and every resampled
        # statistic would be NaN
        res <- createError(jmvcore::format(.('Group \'{g}\' has no non-missing observations'), g = groupLevels[n == 0][1]))
      } else {
        perms <- cached_sims(htestTable,
            list(dep = dataHTest$dep, group = as.character(dataHTest$group), reps = self$options$reps, seedBool = self$options$seedBool, rngSeed = self$options$rngSeed),
            function() private$.computePerms(dataHTest))
        res <- private$.computePval(perms, m[1]-m[2])
        private$.preparePlot(perms, m[1]-m[2], res$direction)
      }
        
        if (isError(res)) {
          
          htestTable$setRow(rowKey=depVarName, list(
            "reps"='',
            "md"='',
            "p"=''))
            
            message <- extractErrorMessage(res)
            
            htestTable$addFootnote(rowKey=depVarName, 'md', message)
            
          } else {
            
            htestTable$setRow(rowKey=depVarName, list(
              "reps"=self$options$reps,
              "md"=m[1]-m[2],
              "p"=format_sim_pval(res$pval, self$options$reps)))
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
            
            hypothesis <- self$options$hypothesis
            groupName <- self$options$group
            
            groups <- NULL
            if ( ! is.null(groupName))
            groups <- base::levels(self$data[[groupName]])
            if (length(groups) != 2)
            groups <- c('Group 1', 'Group 2')
            
            table <- self$results$htest
            
            if (hypothesis == 'oneGreater')
            table$setNote("hyp", jmvcore::format("H\u2090 \u03BC\u2009<sub>{}</sub> > \u03BC\u2009<sub>{}</sub>", groups[1], groups[2]))
            else if (hypothesis == 'twoGreater')
            table$setNote("hyp", jmvcore::format("H\u2090 \u03BC\u2009<sub>{}</sub> < \u03BC\u2009<sub>{}</sub>", groups[1], groups[2]))
            else
            table$setNote("hyp", jmvcore::format("H\u2090 \u03BC\u2009<sub>{}</sub> \u2260 \u03BC\u2009<sub>{}</sub>", groups[1], groups[2]))
            
          },
          .computePval = function(perms, dm) {

            reps <- self$options$reps
            direction <- map_direction(self$options$hypothesis)

            pval <- compute_null_pval(perms, dm, direction)

            list(reps = reps, pval = pval, direction = direction)
          },
          .computePerms = function(dataHTest){

            groupLevels <- base::levels(dataHTest$group)
            
            reps <- self$options$reps

            set_seed_if(self$options$seedBool, self$options$rngSeed)

            # Not infer::specify(dep ~ group): under Jamovi 2.7's formula
            # sandbox that call is rejected (see ?resample_means).  Same
            # permutation distribution, drawn directly.
            perms <- permute_diff_means(dataHTest$dep, dataHTest$group, groupLevels, reps)

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
          .preparePlot = function(perms, dm, direction) {

            permplot <- self$results$simplot
            dotHist <- self$options$dotHist

            permplot$setState(list(df=strip_infer(perms), obs_stat=dm, direction=direction, dotHist=dotHist, showCounts=self$options$showCounts,
                                          xlab="difference (group 1 - group 2)", obs_label="Observed\nDifference"))

        },
          .permPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            st <- image$state
            plot_null_dist(st$df, st$obs_stat, st$direction, st$dotHist,
                           xlab = "difference (group 1 - group 2)",
                           obs_label = "Observed\nDifference",
                           show_counts = isTRUE(st$showCounts),
                           plot_width = image$width)
        },
          .formula=function() {
            jmvcore:::composeFormula(self$options$vars, self$options$group)
          }
        )
      )