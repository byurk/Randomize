# Render end-to-end plots from the real analysis functions for eyeball
# inspection (complements generate_plots.R, which feeds synthetic
# statistics straight into the plot utilities).
#
# Usage: Rscript tests/visual/generate_analysis_plots.R <tag>
# Output: tests/visual/out/<tag>/analysis_sheet_NN.png (2 x 2 per sheet)
#
# Every scenario is rendered at Jamovi's plot size (400 x 350) in dotplot
# and histogram mode, with and without count labels.

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "analyses"

suppressMessages(suppressWarnings(devtools::load_all(".", quiet = TRUE)))
suppressMessages({library(grid); library(ggplot2)})

out_dir <- file.path("tests/visual/out", tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
W <- 400; H <- 350; SCALE <- 2

set.seed(2026)
mk_two <- function(n1, n2, d = 1, sd = 2) data.frame(
    score = c(rnorm(n1, 10, sd), rnorm(n2, 10 + d, sd)),
    group = factor(rep(c("A", "B"), c(n1, n2))))
mk_pair <- function(n, d = 1) { x <- rnorm(n, 50, 10); data.frame(pre = x, post = x + rnorm(n, d, 4)) }
mk_multi <- function(n, k = 3, d = 1) data.frame(
    score = unlist(lapply(seq_len(k), function(i) rnorm(n, 10 + d * (i - 1), 2))),
    group = factor(rep(LETTERS[1:k], each = n)))
mk_reg <- function(n, b = 0.5) { x <- rnorm(n, 5, 2); data.frame(x = x, y = 3 + b * x + rnorm(n, 0, 2)) }
mk_cat <- function(n1, p1, n2, p2) data.frame(
    outcome = factor(c(rep(c("Yes", "No"), round(c(p1, 1 - p1) * n1)),
                       rep(c("Yes", "No"), round(c(p2, 1 - p2) * n2)))),
    group = factor(rep(c("T", "C"), c(n1, n2))))
mk_tab <- function(mat) {
    d <- as.data.frame(as.table(mat)); names(d) <- c("r", "c", "n")
    d[rep(seq_len(nrow(d)), d$n), c("r", "c")]
}
skew <- function(n) data.frame(v = rexp(n, 1 / 20))

S <- list()
# Each scenario's data set is generated once so the four renders of a
# scenario (dotplot / histogram, with / without counts) show the same
# simulation.
add <- function(name, make, fun) {
    d <- make()
    S[[name]] <<- function(dh, sc) fun(dh, sc, d)
}

# ---- hypothesis tests ----
add("twomean_n25_1000", function() mk_two(25, 25), function(dh, sc, d) twomeanhtest(data = d, vars = "score", group = "group",
    hypothesis = "different", reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 1, showCounts = sc))
add("twomean_n8_100", function() mk_two(8, 8, d = 2), function(dh, sc, d) twomeanhtest(data = d, vars = "score", group = "group",
    hypothesis = "oneGreater", reps = 100, dotHist = dh, seedBool = TRUE, rngSeed = 2, showCounts = sc))
add("twomean_n60_5000_p0", function() mk_two(60, 60, d = 3), function(dh, sc, d) twomeanhtest(data = d, vars = "score", group = "group",
    hypothesis = "different", reps = 5000, dotHist = dh, seedBool = TRUE, rngSeed = 3, showCounts = sc))
add("twomean_reps10", function() mk_two(20, 20), function(dh, sc, d) twomeanhtest(data = d, vars = "score", group = "group",
    hypothesis = "twoGreater", reps = 10, dotHist = dh, seedBool = TRUE, rngSeed = 4, showCounts = sc))
add("paired_n30_1000", function() mk_pair(30, 2), function(dh, sc, d) pairedmeanhtest(data = d, pairs = list(list(i1 = "pre", i2 = "post")),
    hypothesis = "different", reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 5, showCounts = sc))
add("paired_n6_200", function() mk_pair(6, 1), function(dh, sc, d) pairedmeanhtest(data = d, pairs = list(list(i1 = "pre", i2 = "post")),
    hypothesis = "twoGreater", reps = 200, dotHist = dh, seedBool = TRUE, rngSeed = 6, showCounts = sc))
add("multi_3g_1000", function() mk_multi(15, 3, 1.2), function(dh, sc, d) multimeanhtest(data = d, vars = "score", group = "group",
    reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 7, showCounts = sc))
add("multi_5g_500", function() mk_multi(10, 5, 0.5), function(dh, sc, d) multimeanhtest(data = d, vars = "score", group = "group",
    reps = 500, dotHist = dh, seedBool = TRUE, rngSeed = 8, showCounts = sc))
add("slope_n40_1000", function() mk_reg(40, 0.4), function(dh, sc, d) slopehtest(data = d, dep = "y", indep = "x",
    hypothesis = "notequal", reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 9, showCounts = sc))
add("slope_n12_100_greater", function() mk_reg(12, 0.8), function(dh, sc, d) slopehtest(data = d, dep = "y", indep = "x",
    hypothesis = "greater", reps = 100, dotHist = dh, seedBool = TRUE, rngSeed = 10, showCounts = sc))
add("sprop_n20_100", function() data.frame(x = c(14, 6)), function(dh, sc, d) SinglePropHTest(data = d, resp = "x", areCounts = TRUE,
    testValue = 0.5, alt = "greater", reps = 100, dotHist = dh, seedBool = TRUE, rngSeed = 11, showCounts = sc))
add("sprop_n30_1000_two", function() data.frame(x = c(20, 10)), function(dh, sc, d) SinglePropHTest(data = d, resp = "x", areCounts = TRUE,
    testValue = 0.5, alt = "notequal", reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 12, showCounts = sc))
add("sprop_n100_2000_less", function() data.frame(x = c(41, 59)), function(dh, sc, d) SinglePropHTest(data = d, resp = "x", areCounts = TRUE,
    testValue = 0.5, alt = "less", reps = 2000, dotHist = dh, seedBool = TRUE, rngSeed = 13, showCounts = sc))
add("sprop_n10_10000", function() data.frame(x = c(8, 2)), function(dh, sc, d) SinglePropHTest(data = d, resp = "x", areCounts = TRUE,
    testValue = 0.5, alt = "greater", reps = 10000, dotHist = dh, seedBool = TRUE, rngSeed = 14, showCounts = sc))
add("sprop_n1000_1000", function() data.frame(x = c(530, 470)), function(dh, sc, d) SinglePropHTest(data = d, resp = "x", areCounts = TRUE,
    testValue = 0.5, alt = "greater", reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 15, showCounts = sc))
add("sprop_5_0", function() data.frame(x = c(5, 0)), function(dh, sc, d) SinglePropHTest(data = d, resp = "x", areCounts = TRUE,
    testValue = 0.5, alt = "greater", reps = 100, dotHist = dh, seedBool = TRUE, rngSeed = 16, showCounts = sc))
add("sprop_p1_testval01", function() data.frame(x = c(12, 8)), function(dh, sc, d) SinglePropHTest(data = d, resp = "x", areCounts = TRUE,
    testValue = 0.1, alt = "less", reps = 500, dotHist = dh, seedBool = TRUE, rngSeed = 17, showCounts = sc))
add("twoprop_25_25", function() mk_cat(25, 0.72, 25, 0.48), function(dh, sc, d) TwoPropHTest(data = d, rows = "group", cols = "outcome",
    hypothesis = "different", reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 18, showCounts = sc, compare = "rows"))
add("twoprop_30_20", function() mk_cat(30, 0.6, 20, 0.35), function(dh, sc, d) TwoPropHTest(data = d, rows = "group", cols = "outcome",
    hypothesis = "oneGreater", reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 19, showCounts = sc, compare = "rows"))
add("twoprop_12_8_100", function() mk_cat(12, 0.75, 8, 0.5), function(dh, sc, d) TwoPropHTest(data = d, rows = "group", cols = "outcome",
    hypothesis = "different", reps = 100, dotHist = dh, seedBool = TRUE, rngSeed = 20, showCounts = sc, compare = "rows"))
add("twoprop_200_150", function() mk_cat(200, 0.55, 150, 0.45), function(dh, sc, d) TwoPropHTest(data = d, rows = "group", cols = "outcome",
    hypothesis = "different", reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 21, showCounts = sc, compare = "rows"))
add("chisq_2x2", function() mk_tab(matrix(c(18, 7, 12, 13), 2)), function(dh, sc, d) ContTabHTest(data = d, rows = "r", cols = "c",
    reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 22, showCounts = sc, compare = "rows"))
add("chisq_3x3", function() mk_tab(matrix(c(10, 5, 5, 6, 12, 4, 4, 5, 9), 3)), function(dh, sc, d) ContTabHTest(data = d, rows = "r", cols = "c",
    reps = 1000, dotHist = dh, seedBool = TRUE, rngSeed = 23, showCounts = sc, compare = "rows"))
add("chisq_4x2_200", function() mk_tab(matrix(c(20, 15, 10, 5, 10, 15, 20, 25), 4)), function(dh, sc, d) ContTabHTest(data = d, rows = "r", cols = "c",
    reps = 200, dotHist = dh, seedBool = TRUE, rngSeed = 24, showCounts = sc, compare = "rows"))

# ---- confidence intervals ----
add("smean_n30_1000", function() skew(30), function(dh, sc, d) SingleMeanCI(data = d, resp = "v", reps = 1000, confLevel = 95,
    ciType = "bootperc", dotHist = dh, seedBool = TRUE, rngSeed = 31, showCounts = sc))
add("smean_n8_100_se", function() skew(8), function(dh, sc, d) SingleMeanCI(data = d, resp = "v", reps = 100, confLevel = 90,
    ciType = "bootse", dotHist = dh, seedBool = TRUE, rngSeed = 32, showCounts = sc))
add("twomeanCI_n25", function() mk_two(25, 25), function(dh, sc, d) twomeanCI(data = d, vars = "score", group = "group", reps = 1000,
    confLevel = 95, ciType = "bootperc", dotHist = dh, seedBool = TRUE, rngSeed = 33, showCounts = sc))
add("pairedCI_n20_99", function() mk_pair(20, 2), function(dh, sc, d) pairedmeanCI(data = d, pairs = list(list(i1 = "pre", i2 = "post")),
    reps = 1000, confLevel = 99, ciType = "bootperc", dotHist = dh, seedBool = TRUE, rngSeed = 34, showCounts = sc))
add("slopeCI_n40", function() mk_reg(40, 0.4), function(dh, sc, d) slopeCI(data = d, dep = "y", indep = "x", reps = 1000,
    confLevel = 95, ciType = "bootperc", dotHist = dh, seedBool = TRUE, rngSeed = 35, showCounts = sc))
add("spropCI_n20", function() data.frame(x = c(14, 6)), function(dh, sc, d) SinglePropCI(data = d, resp = "x", areCounts = TRUE, reps = 1000,
    confLevel = 95, ciType = "bootperc", dotHist = dh, seedBool = TRUE, rngSeed = 36, showCounts = sc))
add("spropCI_n30_edge", function() data.frame(x = c(29, 1)), function(dh, sc, d) SinglePropCI(data = d, resp = "x", areCounts = TRUE, reps = 1000,
    confLevel = 95, ciType = "bootperc", dotHist = dh, seedBool = TRUE, rngSeed = 37, showCounts = sc))
add("spropCI_n100_se", function() data.frame(x = c(57, 43)), function(dh, sc, d) SinglePropCI(data = d, resp = "x", areCounts = TRUE, reps = 2000,
    confLevel = 95, ciType = "bootse", dotHist = dh, seedBool = TRUE, rngSeed = 38, showCounts = sc))
add("spropCI_5_0", function() data.frame(x = c(5, 0)), function(dh, sc, d) SinglePropCI(data = d, resp = "x", areCounts = TRUE, reps = 100,
    confLevel = 95, ciType = "bootperc", dotHist = dh, seedBool = TRUE, rngSeed = 39, showCounts = sc))
add("twopropCI_25_25", function() mk_cat(25, 0.72, 25, 0.48), function(dh, sc, d) TwoPropCI(data = d, rows = "group", cols = "outcome", reps = 1000,
    confLevel = 95, ciType = "bootperc", dotHist = dh, seedBool = TRUE, rngSeed = 40, showCounts = sc, compare = "rows"))
add("twopropCI_30_20_se", function() mk_cat(30, 0.6, 20, 0.35), function(dh, sc, d) TwoPropCI(data = d, rows = "group", cols = "outcome", reps = 1000,
    confLevel = 90, ciType = "bootse", dotHist = dh, seedBool = TRUE, rngSeed = 41, showCounts = sc, compare = "rows"))

plots <- list(); errors <- character(); timings <- character()
for (nm in names(S)) {
    for (dh in c("dotplot", "histogram")) {
        for (sc in c(FALSE, TRUE)) {
            key <- sprintf("%s_%s%s", nm, substr(dh, 1, 4), if (sc) "_cnt" else "")
            t0 <- Sys.time()
            p <- tryCatch({
                r <- S[[nm]](dh, sc)
                plot(r) + ggtitle(key) + theme(plot.title = element_text(size = 9))
            }, error = function(e) { errors <<- c(errors, paste(key, "-", conditionMessage(e))); NULL })
            if (!is.null(p)) {
                plots[[key]] <- p
                png(file.path(out_dir, paste0(key, ".png")), width = W * SCALE, height = H * SCALE, res = 72 * SCALE)
                ok <- tryCatch({ print(p); TRUE }, error = function(e) {
                    errors <<- c(errors, paste(key, "- render:", conditionMessage(e))); FALSE })
                dev.off()
            }
            timings <- c(timings, sprintf("%-32s %5.1fs", key, as.numeric(Sys.time() - t0, units = "secs")))
        }
    }
}

keys <- names(plots); n_per <- 4
n_sheet <- ceiling(length(keys) / n_per)
for (i in seq_len(n_sheet)) {
    sk <- keys[((i - 1) * n_per + 1):min(i * n_per, length(keys))]
    png(file.path(out_dir, sprintf("analysis_sheet_%02d.png", i)), width = 2 * W * SCALE, height = 2 * H * SCALE, res = 72 * SCALE)
    grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 2)))
    for (j in seq_along(sk)) {
        vp <- viewport(layout.pos.row = ceiling(j / 2), layout.pos.col = ((j - 1) %% 2) + 1)
        tryCatch(print(plots[[sk[j]]], vp = vp), error = function(e) {
            pushViewport(vp); grid.text(paste(sk[j], "ERROR:", conditionMessage(e))); popViewport() })
    }
    dev.off()
}
cat("Rendered", length(plots), "plots on", n_sheet, "sheets to", out_dir, "\n")
slow <- timings[as.numeric(sub(".* ([0-9.]+)s$", "\\1", timings)) > 5]
if (length(slow)) { cat("Slow (>5s):\n"); cat(paste(" ", slow), sep = "\n") }
if (length(errors)) { cat("ERRORS:\n"); cat(paste(" -", errors), sep = "\n") }
