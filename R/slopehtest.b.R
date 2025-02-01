
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
      private$.preparePlot(perms, b)
        
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

            reps <- self$options$reps
            alt <- self$options$hypothesis

            if (alt == "greater")
                direction <- "greater"
            else if (alt == "less")
                direction <- "less"
            else
                direction <- "two_sided"

            pval <- perms %>%
                infer::get_p_value(obs_stat = b, direction = direction) %>%
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
          .preparePlot = function(perms, b) {

            permplot <- self$results$simplot
            dotHist <- self$options$dotHist
            alt <- self$options$hypothesis

            permplot$setState(list(df=perms, obs_stat=b, alt = alt, dotHist=dotHist))

        },
          .permPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            df <- image$state[["df"]]

            obs_stat = image$state[["obs_stat"]]

            dotHist = image$state[["dotHist"]]

            alt <- image$state[["alt"]]

            b <- obs_stat

            # which tail of the distro should be used to calculate a p-val?
            if (alt == "less"){
                ptail <- "lt"
                caption <- "One-sided: p-val is proportion of\n results \u2264 observed value"
            } else if (alt == "greater"){
                ptail <- "rt"
                caption <- "One-sided: p-val is proportion of\n results \u2265 observed value"
            } else {
                ptail <- df %>%
                    dplyr::summarize(lt = mean(stat < b), rt = mean(stat > b)) %>%
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
                            dplyr::mutate(x.bin = b - ( (b - stat) %/% bw ) * bw)

                    } else {

                        df <- df %>%
                            dplyr::mutate(x.bin = b + ( (stat - b) %/% bw ) * bw)

                    }


                }


                # color the dots in the tail that will be used to calculate the p-val
                # put vertical line at value of b

                df <- df %>%
                    dplyr::mutate(extreme = (ptail == "lt" & stat <= b) | (ptail == "rt" & stat >= b))


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
                    ggplot2::geom_vline(xintercept = b, linetype = "dashed", color = "red") +
                    ggplot2::scale_fill_manual(values = c("FALSE" = "black", "TRUE" = "#ff8c8c")) +
                    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#ff8c8c")) +
                    ggplot2::annotate("text", x = b, y = lab_ht, label = "Observed\nSlope", color = "red") +
                    ggplot2::theme_minimal() +
                    ggplot2::ylab("count") +
                    ggplot2::xlab("slope") +
                    ggplot2::ylim(0, lab_ht + 3) +
                    ggplot2::xlim(min(littlex - bw, b - bw), max(bigx + bw, b + bw)) +
                    ggplot2::coord_equal(ratio = bw*2/3) +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                                   plot.caption = element_text(color = "red", hjust = 0))



            } else {

                closed <- dplyr::if_else(ptail == "lt", "right", "left")

                df <- df %>%
                    dplyr::mutate(extreme = (ptail == "lt" & stat <= b) | (ptail == "rt" & stat >= b))

                p <- ggplot2::ggplot(data=df, aes(x=stat, fill = extreme)) +
                    ggplot2::geom_histogram(boundary = b, closed = closed, show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept=obs_stat, linetype='dashed', color = "red") +
                    ggplot2::scale_fill_manual(values = c("FALSE" = "black", "TRUE" = "#ff8c8c")) +
                    ggplot2::theme_minimal() +
                    ggplot2::xlab("slope") +
                    ggplot2::ylab("count") +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                                   plot.caption = element_text(color = "red", hjust = 0))

                yMax <- ggplot2::layer_scales(p)$y$range$range[2]  # upper y-limit
                p <- p + ggplot2::annotate("text", x = b, y = yMax, vjust = "top", label = "Observed\nSlope", color = "red")
            }


            return(p)
        },
          .formula=function() {
            jmvcore:::composeFormula(self$options$dep, self$options$indep)
          }
        )
      )