
# This file is a generated template, your changes will not be overwritten

SinglePropHTestClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "SinglePropHTestClass",
    inherit = SinglePropHTestBase,
    private = list(
        #### Init + run functions ----
        .init = function() {

            private$.initSummTable()
            private$.initSimTable()
            private$.initPlot()

        },
        .run = function() {

            private$.errorCheck()

            if (length(self$options$resp) > 0) {

                results <- private$.computeSumm()
                boot <- private$.computeBoots()
                simres <- private$.computePval(boot)

                private$.populateSummTable(results)
                private$.populateSimTable(simres)

                private$.preparePlot(boot)

            }
        },

        #### Compute results ----
        .computeSumm = function() {

            resp <- self$options$resp

            levelList <- list()

            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total
            levels <- results$levels

            for (i in seq_along(counts)) {

                level <- levels[i]
                count <- counts[i]
                prop <- count / total

                levelList[[i]] <- list(var=resp, level=level, count=count, total=total,
                                       prop=prop)
            }

            return(levelList)
        },
        .computeBoots = function() {

            resp <- self$options$resp
            reps <- self$options$reps
            testValue <- self$options$testValue

            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total
            levels <- results$levels

            df <- tibble::tibble(level=levels, count=counts) %>%
                tidyr::uncount(count)

            if(self$options$seedBool){
                set.seed(self$options$rngSeed)
            } else {
                set.seed(NULL)
            }

            boot <- df %>%
                infer::specify(response = level, success = levels[1]) %>%
                infer::hypothesize(null = "point", p = testValue) %>%
                infer::generate(reps = reps, type = "draw") %>%
                infer::calculate(stat = "prop")


            return(boot)
        },
        .computePval = function(boot) {

            resp <- self$options$resp
            reps <- self$options$reps
            alt <- self$options$alt

            if (alt == "greater")
                direction <- "greater"
            else if (alt == "less")
                direction <- "less"
            else
                direction <- "two_sided"

            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total
            levels <- results$levels

            pval <- boot %>%
                infer::get_p_value(obs_stat = counts[1] / total, direction = direction) %>%
                dplyr::pull()


            simres <- list(reps = reps, pval = pval)
            return(simres)
        },

        #### Init tables/plots functions ----
        .initSummTable = function() {

            summtable <- self$results$get('summtable')

            resp <- self$options$resp

            if(!is.null(resp)){

            varData <- jmvcore::naOmit(self$data[[resp]])

            if (self$options$areCounts) {

                levels <- paste(1:2)

            } else {

                levels <- base::levels(varData)
                if (length(levels) == 0)
                    levels <- c('\u2026', '\u2026 ')
            }

            for (level in levels) {
                key <- paste0(resp, '`', level)
                summtable$addRow(rowKey=key, values=list(var=resp, level=level))
            }

            summtable$addFormat(rowKey=paste0(resp, '`', levels[1]), 'var', jmvcore::Cell.BEGIN_GROUP)
            summtable$addFormat(rowKey=paste0(resp, '`', levels[length(levels)]), 'var', jmvcore::Cell.END_GROUP)
            }

        },
        .initSimTable = function() {

            simtable <- self$results$get('simtable')

            alt <- self$options$alt
            if (alt == 'greater')
                op <- '>'
            else if (alt == 'less')
                op <- '<'
            else
                op <- '\u2260'

            simtable$setNote(
                'hyp',
                jmvcore::format(
                    'H\u2090: p {direction} {testValue}',
                    direction=op,
                    testValue=self$options$testValue
                )
            )
        },
        .initPlot = function() {

            bootplot <- self$results$Plot

        },

        #### Populate tables functions ----
        .populateSummTable = function(results) {

            table <- self$results$get('summtable')
            table$deleteRows()

            resp <- self$options$resp
            varResults <- results

            levels <- character(length(varResults))
            for (i in seq_along(levels))
                levels[i] <- varResults[[i]]$level

            for (i in seq_along(levels)) {

                level <- levels[i]
                key <- paste0(resp, '`', level)
                table$addRow(rowKey=key, values=varResults[[i]])

            }

            table$addFormat(rowKey=paste0(resp, '`', levels[1]), 'var', jmvcore::Cell.BEGIN_GROUP)
            table$addFormat(rowKey=paste0(resp, '`', levels[length(levels)]), 'var', jmvcore::Cell.END_GROUP)

        },
        .populateSimTable = function(simres) {

            table <- self$results$get('simtable')
            table$deleteRows()

            table$addRow(rowKey=1, values=simres)

        },

        #### Plot functions ----
        .preparePlot = function(boot) {

            bootplot <- self$results$Plot
            dotHist <- self$options$dotHist
            alt <- self$options$alt

            resp <- self$options$resp
            results <- private$.counts(resp)

            counts <- results$counts
            total  <- results$total
            levels <- results$levels

            obs_stat <- counts[1] / total

            bootplot$setState(list(df=boot, obs_stat=obs_stat, alt = alt, dotHist=dotHist))

        },
        .bootPlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            df <- image$state[["df"]]

            obs_stat = image$state[["obs_stat"]]

            dotHist = image$state[["dotHist"]]

            alt <- image$state[["alt"]]

             phat <- obs_stat

             boot <- df

             # which tail of the distro should be used to calculate a p-val?
             if (alt == "less"){
                 ptail <- "lt"
                 caption <- "One-sided: p-val is proportion of\n results \u2264 observed value"
             } else if (alt == "greater"){
                 ptail <- "rt"
                 caption <- "One-sided: p-val is proportion of\n results \u2265 observed value"
             } else {
                 ptail <- boot %>%
                     dplyr::summarize(lt = mean(stat < phat), rt = mean(stat > phat)) %>%
                     dplyr::mutate(pt = dplyr::if_else(rt < lt, "rt", "lt")) %>%
                     dplyr::pull(pt)
                 if (ptail == "lt")
                     caption <- "Two-sided: p-val is 2\u00D7 proportion of\n results \u2264 observed value"
                 else
                     caption <- "Two-sided: p-val is 2\u00D7 proportion of\n results \u2265 observed value"
             }

             if(dotHist == "dotplot"){

                ndist <- dplyr::n_distinct(boot$stat)

                # bin width, also used for scaling dots when not binning
                bw <- boot %>%
                    dplyr::summarize(min = min(stat), max = max(stat)) %>%
                    #mutate(bw = (max - min) / 30 * (1 + 3*.Machine$double.eps)) %>%
                    dplyr::mutate(bw = (max - min) / 30 ) %>%
                    dplyr::pull(bw)


                if (ndist <= 30){ # don't bin unless there are more than 30 unique values

                    boot <- boot %>%
                        dplyr::mutate(x.bin = stat)

                } else {

                    if(ptail == "lt"){

                        boot <- boot %>%
                            dplyr::mutate(x.bin = phat - ( (phat - stat) %/% bw ) * bw)

                    } else {

                        boot <- boot %>%
                            dplyr::mutate(x.bin = phat + ( (stat - phat) %/% bw ) * bw)

                    }


                }


                # color the dots in the tail that will be used to calculate the p-val
                # put vertical line at value of phat

                boot <- boot %>%
                    dplyr::mutate(extreme = (ptail == "lt" & stat <= phat) | (ptail == "rt" & stat >= phat))


                # dotplot based on tjebo answer from:
                # https://stackoverflow.com/questions/53697235/ggplot-dotplot-what-is-the-proper-use-of-geom-dotplot


                boot <- boot %>%
                    dplyr::group_by(x.bin) %>%
                    dplyr::mutate(y = seq_along(x.bin))

                labht <- boot %>%
                    dplyr::ungroup() %>%
                    dplyr::summarize(labht = max(y)) %>%
                    dplyr::pull()

                lab_ht <- max(labht - 1, 3)

                bigx <- boot %>%
                    dplyr::ungroup() %>%
                    dplyr::summarize(bigx = max(x.bin)) %>%
                    dplyr::pull()

                littlex <- boot %>%
                    dplyr::ungroup() %>%
                    dplyr::summarize(littlex = min(x.bin)) %>%
                    dplyr::pull()

                # do not show legend
                p <- ggplot2::ggplot(boot) +
                    ggforce::geom_ellipse(aes(x0 = x.bin, y0 = y, a = bw/3, b = 0.5, angle = 0, fill = extreme, color = extreme),
                                          show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept = phat, linetype = "dashed", color = "red") +
                    ggplot2::scale_fill_manual(values = c("black", "#ff8c8c")) +
                    ggplot2::scale_color_manual(values = c("black", "#ff8c8c")) +
                    ggplot2::annotate("text", x = phat, y = lab_ht, label = "Observed\nProportion", color = "red") +
                    ggplot2::theme_minimal() +
                    ggplot2::ylab("count") +
                    ggplot2::xlab("proportion") +
                    ggplot2::ylim(0, lab_ht + 3) +
                    ggplot2::xlim(min(littlex - bw, phat - bw), max(bigx + bw, phat + bw)) +
                    ggplot2::coord_equal(ratio = bw*2/3) +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                          plot.caption = element_text(color = "red", hjust = 0))



            } else {

                closed <- dplyr::if_else(ptail == "lt", "right", "left")

                boot <- boot %>%
                    dplyr::mutate(extreme = (ptail == "lt" & stat <= phat) | (ptail == "rt" & stat >= phat))

                p <- ggplot2::ggplot(data=boot, aes(x=stat, fill = extreme)) +
                    ggplot2::geom_histogram(boundary = phat, closed = closed, show.legend = FALSE) +
                    ggplot2::geom_vline(xintercept=obs_stat, linetype='dashed', color = "red") +
                    ggplot2::scale_fill_manual(values = c("black", "#ff8c8c")) +
                    ggplot2::theme_minimal() +
                    ggplot2::xlab("proportion") +
                    ggplot2::ylab("count") +
                    ggplot2::labs(caption = caption) +
                    ggplot2::theme(text = element_text(size = 14),
                          plot.caption = element_text(color = "red", hjust = 0))

                yMax <- ggplot2::layer_scales(p)$y$range$range[2]  # upper y-limit
                p <- p + ggplot2::annotate("text", x = phat, y = yMax, vjust = "top", label = "Observed\nProportion", color = "red")
            }


            return(p)
        },

        #### Helper functions ----
        .errorCheck = function() {

            resp <- self$options$resp
            data <- self$data

            if(!is.null(resp)){

            column <- jmvcore::naOmit(data[[resp]])
            if (length(column) == 0) {
                jmvcore::reject(
                    jmvcore::format("Variable '{resp}' contains no data", resp=resp),
                    code=''
                )
            }
            }

        },
        .counts = function(var) {

            initing <- nrow(self$data) == 0

            varData <- jmvcore::naOmit(self$data[[var]])

            if (self$options$areCounts) {

                levels <- paste(1:length(varData))

                if (jmvcore::canBeNumeric(varData))
                    counts <- jmvcore::toNumeric(varData)
                else
                    counts <- suppressWarnings(as.numeric(as.character(varData)))

            } else {

                levels <- base::levels(varData)
                if (length(levels) == 0)
                    levels <- paste(sort(unique(varData)))

                counts <- as.vector(table(varData))

            }

            total <- base::sum(counts, na.rm=TRUE)

            list(levels=levels, counts=counts, total=total)
        }
        )
)
