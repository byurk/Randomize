# Test script for Randomize module
#
# This script exercises each analysis with clean data and with NA data
# to verify the refactored code produces correct results.
#
# Run with: Rscript tests/test_analyses.R
# Requires the package to be installed first via jmvtools::install()

suppressWarnings(devtools::load_all('.', quiet = TRUE))

passed <- 0
failed <- 0

test <- function(name, expr) {
    result <- tryCatch({
        expr
        cat("  PASS:", name, "\n")
        passed <<- passed + 1
        TRUE
    }, error = function(e) {
        cat("  FAIL:", name, "-", e$message, "\n")
        failed <<- failed + 1
        FALSE
    })
}

assert_close <- function(actual, expected, tol = 0.01) {
    if (is.na(actual) || is.na(expected))
        stop(paste("NA result: actual =", actual, ", expected =", expected))
    if (abs(actual - expected) > tol)
        stop(paste("Values differ: actual =", round(actual, 4), ", expected =", round(expected, 4)))
}

assert_not_na <- function(val, label = "value") {
    if (is.null(val) || (length(val) == 1 && is.na(val)))
        stop(paste(label, "is NA or NULL"))
}

# The simulation plot: a single Image for most analyses, but the paired
# analyses draw one per pair (an Array keyed by the pairs); take the first
sim_img <- function(r) {
    img <- if (!is.null(r[["simplot"]])) r$simplot else r$Plot
    if (inherits(img, "Array")) img <- img$get(key = img$itemKeys[[1]])
    img
}

# A simulation p-value of 0 is reported as the string "< 1/reps"
# (e.g. "< .001"); recover a numeric bound for range checks
p_num <- function(p) {
    if (is.character(p)) as.numeric(paste0("0", sub("^<\\s*", "", p))) else p
}

# ============================================================
# Test data
# ============================================================

set.seed(42)
n <- 50

# Continuous data with two groups
clean_data <- data.frame(
    score = c(rnorm(25, mean = 10, sd = 2), rnorm(25, mean = 12, sd = 2)),
    group = factor(rep(c("A", "B"), each = 25)),
    measure1 = rnorm(n, 50, 10),
    measure2 = rnorm(n, 55, 10)
)

# Same data but with some NAs
na_data <- clean_data
na_data$score[c(3, 7, 28, 45)] <- NA
na_data$measure1[c(5, 15)] <- NA
na_data$measure2[c(10, 20)] <- NA

# Categorical data for proportions
cat_data <- data.frame(
    outcome = factor(sample(c("Yes", "No"), n, replace = TRUE, prob = c(0.6, 0.4))),
    group = factor(rep(c("Treatment", "Control"), each = 25))
)

# Multi-group data for ANOVA
multi_data <- data.frame(
    score = c(rnorm(20, 10, 2), rnorm(20, 12, 2), rnorm(20, 11, 2)),
    group = factor(rep(c("A", "B", "C"), each = 20))
)

# Regression data
reg_data <- data.frame(
    x = rnorm(n, 5, 2),
    y = NA_real_
)
reg_data$y <- 3 + 2 * reg_data$x + rnorm(n, 0, 1)

cat("\n=== Testing SingleMeanCI ===\n")
test("clean data", {
    r <- SingleMeanCI(data = clean_data, resp = "score", reps = 500,
                      confLevel = 95, ciType = "bootperc", dotHist = "histogram",
                      seedBool = TRUE, rngSeed = 123)
    tbl <- r$simtable$asDF
    assert_not_na(tbl$cil, "cil")
    assert_not_na(tbl$ciu, "ciu")
    stopifnot(tbl$cil < tbl$ciu)
})
test("data with NAs", {
    r <- SingleMeanCI(data = na_data, resp = "score", reps = 500,
                      confLevel = 95, ciType = "bootperc", dotHist = "histogram",
                      seedBool = TRUE, rngSeed = 123)
    tbl <- r$simtable$asDF
    assert_not_na(tbl$cil, "cil")
    assert_not_na(tbl$ciu, "ciu")
})
test("bootstrap SE method", {
    r <- SingleMeanCI(data = clean_data, resp = "score", reps = 500,
                      confLevel = 95, ciType = "bootse", dotHist = "dotplot",
                      seedBool = TRUE, rngSeed = 123)
    tbl <- r$simtable$asDF
    assert_not_na(tbl$cil, "cil")
})

cat("\n=== Testing TwoMeanCI ===\n")
test("clean data", {
    r <- twomeanCI(data = clean_data, vars = "score", group = "group",
                   reps = 500, confLevel = 95, ciType = "bootperc",
                   dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$CITable$asDF
    assert_not_na(tbl$cil, "cil")
    assert_not_na(tbl$ciu, "ciu")
})
test("data with NAs", {
    r <- twomeanCI(data = na_data, vars = "score", group = "group",
                   reps = 500, confLevel = 95, ciType = "bootperc",
                   dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$CITable$asDF
    assert_not_na(tbl$cil, "cil")
})

cat("\n=== Testing TwoMeanHTest ===\n")
test("two-sided", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                      hypothesis = "different", reps = 500,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$htest$asDF
    assert_not_na(tbl$p, "p-value")
    stopifnot(tbl$p >= 0 && tbl$p <= 1)
})
test("one-sided (oneGreater)", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                      hypothesis = "oneGreater", reps = 500,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$htest$asDF
    assert_not_na(tbl$p, "p-value")
})
test("data with NAs", {
    r <- twomeanhtest(data = na_data, vars = "score", group = "group",
                      hypothesis = "different", reps = 500,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$htest$asDF
    assert_not_na(tbl$p, "p-value")
})

cat("\n=== Testing PairedMeanCI ===\n")
test("clean data", {
    r <- pairedmeanCI(data = clean_data, pairs = list(list(i1 = "measure1", i2 = "measure2")),
                      reps = 500, confLevel = 95, ciType = "bootperc",
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$CITable$asDF
    assert_not_na(tbl$cil, "cil")
    assert_not_na(tbl$ciu, "ciu")
})
test("data with NAs", {
    r <- pairedmeanCI(data = na_data, pairs = list(list(i1 = "measure1", i2 = "measure2")),
                      reps = 500, confLevel = 95, ciType = "bootperc",
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$CITable$asDF
    assert_not_na(tbl$cil, "cil")
})

cat("\n=== Testing PairedMeanHTest ===\n")
test("two-sided", {
    r <- pairedmeanhtest(data = clean_data, pairs = list(list(i1 = "measure1", i2 = "measure2")),
                         hypothesis = "different", reps = 500,
                         dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$htest$asDF
    assert_not_na(tbl$p, "p-value")
})
test("data with NAs", {
    r <- pairedmeanhtest(data = na_data, pairs = list(list(i1 = "measure1", i2 = "measure2")),
                         hypothesis = "different", reps = 500,
                         dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$htest$asDF
    assert_not_na(tbl$p, "p-value")
})

cat("\n=== Testing SlopeCI ===\n")
test("clean data", {
    r <- slopeCI(data = reg_data, dep = "y", indep = "x",
                 reps = 500, confLevel = 95, ciType = "bootperc",
                 dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$CITable$asDF
    assert_not_na(tbl$cil, "cil")
    assert_not_na(tbl$ciu, "ciu")
    # Slope should be near 2
    assert_close(tbl$b, 2, tol = 0.5)
})

cat("\n=== Testing SlopeHTest ===\n")
test("two-sided", {
    r <- slopehtest(data = reg_data, dep = "y", indep = "x",
                    hypothesis = "notequal", reps = 500,
                    dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$htest$asDF
    assert_not_na(tbl$p, "p-value")
    # Should be significant (true slope = 2); with 500 reps and no
    # permutation as extreme the p-value is reported as "< .002"
    stopifnot(p_num(tbl$p) < 0.05)
})

cat("\n=== Testing SinglePropCI ===\n")
test("clean data", {
    r <- SinglePropCI(data = cat_data, resp = "outcome",
                      reps = 500, confLevel = 95, ciType = "bootperc",
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$simtable$asDF
    assert_not_na(tbl$cil, "cil")
    assert_not_na(tbl$ciu, "ciu")
    # CI should be in [0, 1]
    stopifnot(tbl$cil >= 0 && tbl$ciu <= 1)
})

cat("\n=== Testing SinglePropHTest ===\n")
test("two-sided", {
    r <- SinglePropHTest(data = cat_data, resp = "outcome",
                         testValue = 0.5, alt = "notequal", reps = 500,
                         dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    tbl <- r$simtable$asDF
    assert_not_na(tbl$p, "p-value")
    stopifnot(tbl$p >= 0 && tbl$p <= 1)
})

cat("\n=== Testing MultiMeanHTest ===\n")
test("clean data", {
    r <- multimeanhtest(data = multi_data, vars = "score", group = "group",
                        reps = 500, dotHist = "histogram",
                        seedBool = TRUE, rngSeed = 123)
    tbl <- r$htest$asDF
    assert_not_na(tbl$p, "p-value")
    stopifnot(p_num(tbl$p) >= 0 && p_num(tbl$p) <= 1)
})

cat("\n=== Testing TwoPropCI ===\n")
test("clean data", {
    r <- TwoPropCI(data = cat_data, rows = "group", cols = "outcome",
                   reps = 500, confLevel = 95, ciType = "bootperc",
                   dotHist = "histogram", seedBool = TRUE, rngSeed = 123,
                   compare = "rows")
    tbl <- r$simtable$asDF
    assert_not_na(tbl$cil, "cil")
    assert_not_na(tbl$ciu, "ciu")
})

cat("\n=== Testing TwoPropHTest ===\n")
test("two-sided", {
    r <- TwoPropHTest(data = cat_data, rows = "group", cols = "outcome",
                      hypothesis = "different", reps = 500,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123,
                      compare = "rows")
    tbl <- r$simtable$asDF
    assert_not_na(tbl$p, "p-value")
    stopifnot(tbl$p >= 0 && tbl$p <= 1)
})

cat("\n=== Testing ContTabHTest ===\n")
test("clean data", {
    r <- ContTabHTest(data = cat_data, rows = "group", cols = "outcome",
                      reps = 500, dotHist = "histogram",
                      seedBool = TRUE, rngSeed = 123, compare = "rows")
    tbl <- r$simtable$asDF
    assert_not_na(tbl$p, "p-value")
    stopifnot(p_num(tbl$p) >= 0 && p_num(tbl$p) <= 1)
})

# ============================================================
# Plot display toggle tests
# ============================================================

assert_is_ggplot <- function(p) {
    if (!inherits(p, "ggplot"))
        stop(paste("Expected ggplot, got", paste(class(p), collapse = ", ")))
}

count_layers <- function(p, geom_class) {
    sum(vapply(p$layers, function(l) inherits(l$geom, geom_class), logical(1)))
}

has_caption <- function(p) {
    cap <- p$labels$caption
    !is.null(cap) && nzchar(cap)
}

cat("\n=== Testing plot() display toggles (hypothesis test) ===\n")
test("default plot returns ggplot", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                      hypothesis = "different", reps = 500,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    p <- plot(r)
    assert_is_ggplot(p)
    stopifnot(count_layers(p, "GeomVline") >= 1)
    stopifnot(has_caption(p))
})
test("show_text = FALSE keeps lines and shading, removes text", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                      hypothesis = "different", reps = 500,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    p <- plot(r, show_text = FALSE)
    assert_is_ggplot(p)
    stopifnot(count_layers(p, "GeomVline") >= 1)
    stopifnot(count_layers(p, "GeomText") == 0)
    stopifnot(!has_caption(p))
})
test("show_lines = FALSE gives bare distribution", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                      hypothesis = "different", reps = 500,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    p <- plot(r, show_lines = FALSE)
    assert_is_ggplot(p)
    stopifnot(count_layers(p, "GeomVline") == 0)
    stopifnot(count_layers(p, "GeomText") == 0)
    stopifnot(!has_caption(p))
})
test("dotplot variant with toggles", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                      hypothesis = "different", reps = 500,
                      dotHist = "dotplot", seedBool = TRUE, rngSeed = 123)
    p_full <- plot(r)
    p_bare <- plot(r, show_lines = FALSE)
    assert_is_ggplot(p_full)
    assert_is_ggplot(p_bare)
    stopifnot(count_layers(p_full, "GeomVline") >= 1)
    stopifnot(count_layers(p_bare, "GeomVline") == 0)
})

cat("\n=== Testing plot() display toggles (bootstrap CI) ===\n")
test("default plot returns ggplot", {
    r <- twomeanCI(data = clean_data, vars = "score", group = "group",
                   reps = 500, confLevel = 95, ciType = "bootperc",
                   dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    p <- plot(r)
    assert_is_ggplot(p)
    stopifnot(count_layers(p, "GeomVline") >= 1)
    stopifnot(has_caption(p))
})
test("show_lines = FALSE gives bare distribution", {
    r <- twomeanCI(data = clean_data, vars = "score", group = "group",
                   reps = 500, confLevel = 95, ciType = "bootperc",
                   dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    p <- plot(r, show_lines = FALSE)
    assert_is_ggplot(p)
    stopifnot(count_layers(p, "GeomVline") == 0)
    stopifnot(!has_caption(p))
})
test("show_text = FALSE keeps lines, removes caption", {
    r <- twomeanCI(data = clean_data, vars = "score", group = "group",
                   reps = 500, confLevel = 95, ciType = "bootperc",
                   dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    p <- plot(r, show_text = FALSE)
    assert_is_ggplot(p)
    stopifnot(!has_caption(p))
    stopifnot(count_layers(p, "GeomVline") >= 1)
})
test("dotplot variant with toggles", {
    r <- twomeanCI(data = clean_data, vars = "score", group = "group",
                   reps = 500, confLevel = 95, ciType = "bootperc",
                   dotHist = "dotplot", seedBool = TRUE, rngSeed = 123)
    p_full <- plot(r)
    p_bare <- plot(r, show_lines = FALSE)
    assert_is_ggplot(p_full)
    assert_is_ggplot(p_bare)
    stopifnot(count_layers(p_full, "GeomVline") >= 1)
    stopifnot(count_layers(p_bare, "GeomVline") == 0)
})

# ============================================================
# Issue regressions
# ============================================================

cat("\n=== #12: a level with zero count ===\n")
test("SinglePropHTest, counts 5 / 0", {
    r <- SinglePropHTest(data = data.frame(x = c(5, 0)), resp = "x", areCounts = TRUE,
                         testValue = 0.5, alt = "greater", reps = 200,
                         dotHist = "dotplot", seedBool = TRUE, rngSeed = 1)
    tbl <- r$simtable$asDF
    stopifnot(tbl$obsProp == 1)
    # P(5 heads in 5 fair flips) = 1/32: the simulated p must land near it
    stopifnot(p_num(tbl$p) > 0, p_num(tbl$p) < 0.15)
    stopifnot(all(r$Plot$state$df$stat >= 0 & r$Plot$state$df$stat <= 1))
    assert_is_ggplot(plot(r))
})
test("SinglePropHTest, counts 0 / 5", {
    r <- SinglePropHTest(data = data.frame(x = c(0, 5)), resp = "x", areCounts = TRUE,
                         testValue = 0.5, alt = "less", reps = 200,
                         dotHist = "histogram", seedBool = TRUE, rngSeed = 1)
    stopifnot(r$simtable$asDF$obsProp == 0)
    assert_is_ggplot(plot(r))
})
test("SinglePropHTest, factor level never observed", {
    d <- data.frame(y = factor(rep("H", 5), levels = c("H", "T")))
    r <- SinglePropHTest(data = d, resp = "y", testValue = 0.5, alt = "greater",
                         reps = 200, seedBool = TRUE, rngSeed = 1)
    stopifnot(r$simtable$asDF$obsProp == 1)
    stopifnot(r$summtable$asDF$count == c(5, 0))
})
test("SinglePropCI, counts 5 / 0", {
    r <- SinglePropCI(data = data.frame(x = c(5, 0)), resp = "x", areCounts = TRUE,
                      reps = 200, confLevel = 95, ciType = "bootperc",
                      dotHist = "dotplot", seedBool = TRUE, rngSeed = 1)
    tbl <- r$simtable$asDF
    stopifnot(tbl$obsProp == 1, tbl$cil == 1, tbl$ciu == 1)
    assert_is_ggplot(plot(r))
})
test("SinglePropHTest rejects a variable with 3 levels", {
    d <- data.frame(y = factor(c("a", "b", "c", "a")))
    msg <- tryCatch({
        SinglePropHTest(data = d, resp = "y", testValue = 0.5, reps = 50)
        ""
    }, error = function(e) conditionMessage(e))
    stopifnot(grepl("exactly 2 levels", msg))
})

cat("\n=== Chi-square: exact-zero statistic reads 0, not 1e-31 ===\n")
test("ContTabHTest zaps floating-point dust in X2", {
    d <- data.frame(group = factor(rep(c("A", "B"), c(13, 12))),
                    outcome = factor(c(rep("No", 4), rep("Yes", 9), rep("No", 3), rep("Yes", 9))))
    r <- ContTabHTest(data = d, rows = "group", cols = "outcome", reps = 200,
                      dotHist = "dotplot", seedBool = TRUE, rngSeed = 3, compare = "rows")
    stopifnot(as.numeric(r$x2tab$asDF[["v[x2]"]]) == 0)
    stopifnot(as.numeric(results_table(r)$x2) == 0)
    st <- r$Plot$state
    stopifnot(all(st$df$stat >= 0), sum(st$df$stat == 0) > 0)
    stopifnot(p_num(results_table(r)$p) == 1)
})

cat("\n=== Floating-point ties count as ties ===\n")
test("a simulated difference one ulp below the observed one is a tie", {
    s <- c(15/25 - 11/25, 0.16, 4/25, 0.60 - 0.44)   # all mathematically 0.16
    stopifnot(any(s < 0.16))                           # the arithmetic really does differ
    stopifnot(compute_null_pval(data.frame(stat = s), 0.16, "greater") == 1)
    stopifnot(compute_null_pval(data.frame(stat = s), 0.16, "less") == 1)
    p <- plot_null_dist(data.frame(stat = s), 0.16, "greater", "histogram")
    b <- ggplot2::ggplot_build(p)$data[[1]]
    stopifnot(nrow(b) == 1, b$fill == "#ff8c8c")       # one bar, shaded
})

cat("\n=== Filtered-out levels (unused factor levels) ===\n")
test("proportion analyses drop unused levels beyond two, keep an empty second level", {
    sp <- data.frame(ans = factor(c(rep("Yes", 13), rep("No", 7)), levels = c("No", "Yes", "Maybe")))
    r <- SinglePropHTest(data = sp, resp = "ans", testValue = 0.5, reps = 100, seedBool = TRUE, rngSeed = 1)
    stopifnot(nrow(r$summtable$asDF) == 2, abs(results_table(r)$obsProp - 0.35) < 1e-9)
    r <- SinglePropCI(data = sp, resp = "ans", reps = 100, seedBool = TRUE, rngSeed = 1)
    stopifnot(nrow(r$summtable$asDF) == 2)
    one <- data.frame(ans = factor(rep("H", 5), levels = c("H", "T")))
    r <- SinglePropHTest(data = one, resp = "ans", testValue = 0.5, reps = 100, seedBool = TRUE, rngSeed = 1)
    stopifnot(nrow(r$summtable$asDF) == 2, results_table(r)$obsProp == 1)
})
test("two-proportion and chi-square analyses ignore empty rows/columns from hidden levels", {
    tab <- data.frame(g = factor(rep(c("A","B","C"), each = 10)), o = factor(rep(c("Y","N"), 15)))
    tab <- tab[tab$g != "C", ]
    r <- TwoPropHTest(data = tab, rows = "g", cols = "o", hypothesis = "different", reps = 100, seedBool = TRUE, rngSeed = 1, compare = "rows")
    stopifnot(nrow(r$simtable$asDF) == 1, nrow(r$freqs$asDF) == 3 + 1)   # table still shows C (plus total)
    r <- TwoPropCI(data = tab, rows = "g", cols = "o", reps = 100, seedBool = TRUE, rngSeed = 1, compare = "rows")
    stopifnot(nrow(r$simtable$asDF) == 1)
    r <- ContTabHTest(data = tab, rows = "g", cols = "o", reps = 100, seedBool = TRUE, rngSeed = 1, compare = "rows")
    stopifnot(nrow(r$simtable$asDF) == 1, is.finite(results_table(r)$x2))
    # a genuinely empty level of a 2-level factor is still reported, not silently dropped
    z <- data.frame(g = factor(c("A","A","A"), levels = c("A","B")), o = factor(c("Y","N","Y")))
    r <- TwoPropHTest(data = z, rows = "g", cols = "o", hypothesis = "different", reps = 100, compare = "rows")
    stopifnot(nrow(r$simtable$asDF) == 0)
    txt <- gsub("\\s+", " ", paste(capture.output(print(r$diffProp)), collapse = " "))
    stopifnot(grepl("empty row or column", txt))
})
test("cache ignores a corrupt cached object", {
    holder <- jmvcore::Table$new(options = jmvcore::Options$new(), name = "t", title = "t")
    holder$setState(list(sims = list(key = list(k = 1), sims = "not a data frame")))
    r <- cached_sims(holder, list(k = 1), function() data.frame(stat = 1:3))
    stopifnot(is.data.frame(r), nrow(r) == 3)
})

cat("\n=== jamovi weights variable ===\n")
test("a jamovi weights column weights the table analyses once, not twice, and does not crash", {
    d <- data.frame(g = factor(c("A","A","B","B")), o = factor(c("Y","N","Y","N")), w = c(30, 10, 15, 25))
    attr(d, "jmv-weights") <- d$w; attr(d, "jmv-weights-name") <- "w"
    r <- ContTabHTest(data = d, rows = "g", cols = "o", reps = 100, seedBool = TRUE, rngSeed = 1, compare = "rows")
    ft <- r$freqs$asDF
    stopifnot(abs(results_table(r)$x2 - chisq.test(matrix(c(10, 25, 30, 15), 2))$statistic) < 1e-9)
    r2 <- TwoPropHTest(data = d, rows = "g", cols = "o", hypothesis = "different", reps = 100, seedBool = TRUE, rngSeed = 1, compare = "rows")
    stopifnot(abs(results_table(r2)$obsDiff - (10/40 - 25/40)) < 1e-9)
    r3 <- TwoPropCI(data = d, rows = "g", cols = "o", reps = 100, seedBool = TRUE, rngSeed = 1, compare = "rows")
    stopifnot(abs(results_table(r3)$obsDiff - (10/40 - 25/40)) < 1e-9)
    dm <- data.frame(score = c(1, 2, 10, 11), group = factor(c("A","A","B","B")), w = c(3, 1, 1, 3))
    attr(dm, "jmv-weights") <- dm$w; attr(dm, "jmv-weights-name") <- "w"
    r4 <- twomeanhtest(data = dm, vars = "score", group = "group", hypothesis = "different", reps = 50, seedBool = TRUE, rngSeed = 1, desc = TRUE)
    stopifnot(r4$desc$asDF[["num[1]"]] == 4, abs(results_table(r4)$md - (mean(c(1,1,1,2)) - mean(c(10,11,11,11)))) < 1e-9)
})

cat("\n=== #10: p = 0 reported as < 1/reps ===\n")
test("format_sim_pval", {
    stopifnot(identical(format_sim_pval(0, 100), "< .01"))
    stopifnot(identical(format_sim_pval(0, 1000), "< .001"))
    stopifnot(identical(format_sim_pval(0, 5000), "< .0002"))
    stopifnot(identical(format_sim_pval(0, 300), "< .0033"))
    stopifnot(identical(format_sim_pval(0.034, 100), 0.034))
})
test("SinglePropHTest reports < .01 with 100 reps and no extreme sims", {
    r <- SinglePropHTest(data = data.frame(x = c(19, 1)), resp = "x", areCounts = TRUE,
                         testValue = 0.5, alt = "greater", reps = 100,
                         seedBool = TRUE, rngSeed = 1)
    stopifnot(identical(r$simtable$asDF$p, "< .01"))
    stopifnot(identical(results_table(r)$p, "< .01"))
})
test("twomeanhtest reports < .001 with 1000 reps and no extreme sims", {
    d <- data.frame(score = c(rnorm(30, 0, 1), rnorm(30, 8, 1)),
                    group = factor(rep(c("A", "B"), each = 30)))
    r <- twomeanhtest(data = d, vars = "score", group = "group",
                      hypothesis = "different", reps = 1000,
                      seedBool = TRUE, rngSeed = 1)
    stopifnot(identical(r$htest$asDF$p, "< .001"))
})

cat("\n=== #13 / #8: model-based calculator ===\n")
test("normal, Z = -1.5, right tail: area and shading both start at -1.5", {
    r <- modelBased(distro = "ndistro", areaBool = TRUE, obsStat = -1.5, tail = "right")
    assert_close(r$areaTable$asDF$area, pnorm(-1.5, lower.tail = FALSE), 1e-6)
    b <- ggplot2::ggplot_build(plot(r))
    # the shaded area layer comes first; its left edge must be the observed value
    shade <- b$data[[1]]
    assert_close(min(shade$x), -1.5, 1e-9)
    stopifnot(max(shade$x) > 3)
})
test("normal, Z = -1.5, left tail", {
    r <- modelBased(distro = "ndistro", areaBool = TRUE, obsStat = -1.5, tail = "left")
    assert_close(r$areaTable$asDF$area, pnorm(-1.5), 1e-6)
    shade <- ggplot2::ggplot_build(plot(r))$data[[1]]
    assert_close(max(shade$x), -1.5, 1e-9)
})
test("t, both tails", {
    r <- modelBased(distro = "tdistro", dF = 7, areaBool = TRUE, obsStat = 2.2, tail = "both")
    assert_close(r$areaTable$asDF$area, 2 * pt(-2.2, 7), 1e-6)
    stopifnot(r$areaTable$asDF$df == 7)
})
test("chi-square: right tail regardless of tail option", {
    r <- modelBased(distro = "chisq", dF = 4, areaBool = TRUE, obsStat = 9.5, tail = "left")
    assert_close(r$areaTable$asDF$area, pchisq(9.5, 4, lower.tail = FALSE), 1e-6)
    stopifnot(r$areaTable$asDF$df == 4)
    b <- ggplot2::ggplot_build(plot(r))
    shade <- b$data[[1]]
    assert_close(min(shade$x), 9.5, 1e-9)
    stopifnot(min(b$data[[2]]$x) >= 0)   # density curve never below zero
    assert_is_ggplot(plot(r))
})
test("F: area, df columns, no CI multiplier", {
    r <- modelBased(distro = "fdistro", dF = 3, dF2 = 24, areaBool = TRUE, obsStat = 3.1,
                    CIBool = TRUE, confLevel = 95)
    tbl <- r$areaTable$asDF
    assert_close(tbl$area, pf(3.1, 3, 24, lower.tail = FALSE), 1e-6)
    stopifnot(tbl$df1 == 3, tbl$df2 == 24)
    stopifnot(nrow(r$multTable$asDF) == 0)
})
test("CI multiplier for normal and t", {
    r <- modelBased(distro = "ndistro", CIBool = TRUE, confLevel = 90)
    assert_close(r$multTable$asDF$critVal, qnorm(0.95), 1e-6)
    r <- modelBased(distro = "tdistro", dF = 12, CIBool = TRUE, confLevel = 95)
    assert_close(r$multTable$asDF$critVal, qt(0.975, 12), 1e-6)
    stopifnot(r$multTable$asDF$df == 12)
})

cat("\n=== #7: counts on bars / dot stacks ===\n")
test("showCounts option adds one label per bar (histogram)", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                      hypothesis = "different", reps = 300,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123,
                      showCounts = TRUE)
    p <- plot(r)
    # the observed-value label is itself a text layer; counts add one more
    r0 <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                       hypothesis = "different", reps = 300,
                       dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    n0 <- count_layers(plot(r0), "GeomText")
    stopifnot(count_layers(p, "GeomText") == n0 + 1)
    b <- ggplot2::ggplot_build(p)
    bars <- b$data[[1]]
    labs <- b$data[[which(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1)))[1]]]
    stopifnot(nrow(labs) == nrow(bars))
    stopifnot(sum(as.numeric(labs$label)) == 300)
    # the option can also be forced from R on an analysis run without it
    stopifnot(count_layers(plot(r0, show_counts = TRUE), "GeomText") == n0 + 1)
})
test("showCounts labels every dot stack (bootstrap dotplot)", {
    r <- SinglePropCI(data = cat_data, resp = "outcome", reps = 300,
                      confLevel = 95, ciType = "bootperc", dotHist = "dotplot",
                      seedBool = TRUE, rngSeed = 123, showCounts = TRUE)
    p <- plot(r)
    b <- ggplot2::ggplot_build(p)
    labs <- b$data[[which(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1)))[1]]]
    stopifnot(sum(as.numeric(labs$label)) == 300)
    stopifnot(length(unique(labs$x)) == length(unique(b$data[[1]]$x)))
})
test("an enlarged plot magnifies all text (plot_width)", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                      hypothesis = "different", reps = 300,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123)
    st <- sim_img(r)$state
    p400 <- plot_null_dist(st$df, st$obs_stat, st$direction, "histogram", show_counts = TRUE)
    p800 <- plot_null_dist(st$df, st$obs_stat, st$direction, "histogram", show_counts = TRUE,
                           plot_width = 800)
    stopifnot(p400$theme$text$size == 14, p800$theme$text$size == 28)
    txt <- function(p) {
        b <- ggplot2::ggplot_build(p)
        b$data[[which(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1)))[1]]]
    }
    stopifnot(unique(txt(p800)$size) > unique(txt(p400)$size))
    # count labels never exceed the (scaled) tick-label size
    stopifnot(unique(txt(p800)$size) <= 3.9 * 2 + 1e-9)
    stopifnot(text_scale(400) == 1, text_scale(600) == 1.5, text_scale(2000) == 2)
})
test("bare mode still draws no counts", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group",
                      hypothesis = "different", reps = 300,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 123,
                      showCounts = TRUE)
    stopifnot(count_layers(plot(r, show_lines = FALSE, show_counts = FALSE), "GeomText") == 0)
})

cat("\n=== Degenerate input gives a clear message, not a stack trace ===\n")
expect_reject <- function(expr, pattern) {
    msg <- tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))
    if (!grepl(pattern, msg)) stop(paste0("expected message matching '", pattern, "', got: ", substr(msg, 1, 120)))
}
# The mean/slope analyses report data problems as a table footnote (the
# analysis itself succeeds), so look for the message in the printed table
expect_footnote <- function(tbl, pattern) {
    txt <- gsub("\\s+", " ", paste(capture.output(print(tbl)), collapse = " "))
    if (!grepl(pattern, txt)) stop(paste0("expected footnote matching '", pattern, "'"))
}
test("F test: constant response / one obs per group", {
    r <- multimeanhtest(data = data.frame(score = rep(3, 30), group = factor(rep(c("A","B","C"), 10))),
                        vars = "score", group = "group", reps = 50)
    expect_footnote(r$htest, "does not vary")
    r <- multimeanhtest(data = data.frame(score = c(1, 2, 4), group = factor(c("A","B","C"))),
                        vars = "score", group = "group", reps = 50)
    expect_footnote(r$htest, "at least 2 observations")
})
test("slope: constant predictor / too few rows", {
    r <- slopehtest(data = data.frame(x = rep(2, 20), y = rnorm(20)), dep = "y", indep = "x", reps = 50)
    expect_footnote(r$htest, "constant")
    r <- slopeCI(data = data.frame(x = rep(2, 20), y = rnorm(20)), dep = "y", indep = "x", reps = 50)
    expect_footnote(r$CITable, "constant")
    r <- slopehtest(data = data.frame(x = c(1, 2), y = c(3, 5)), dep = "y", indep = "x", reps = 50)
    expect_footnote(r$htest, "At least 3")
})
test("single mean CI: n = 1 rejects, constant sample gives a point CI", {
    expect_reject(SingleMeanCI(data = data.frame(v = 7), resp = "v", reps = 50), "at least 2")
    r <- SingleMeanCI(data = data.frame(v = rep(7, 15)), resp = "v", reps = 100, seedBool = TRUE, rngSeed = 1)
    stopifnot(r$simtable$asDF$cil == 7, r$simtable$asDF$ciu == 7)
    assert_is_ggplot(plot(r))
})
test("proportion counts: NA / negative / fractional / wrong row count", {
    expect_reject(SinglePropHTest(data = data.frame(x = c(3, NA)), resp = "x", areCounts = TRUE, reps = 50), "missing")
    expect_reject(SinglePropHTest(data = data.frame(x = c(3, -2)), resp = "x", areCounts = TRUE, reps = 50), "negative")
    expect_reject(SinglePropHTest(data = data.frame(x = c(3.5, 2.5)), resp = "x", areCounts = TRUE, reps = 50), "whole numbers")
    expect_reject(SinglePropHTest(data = data.frame(x = c(5, 3, 2)), resp = "x", areCounts = TRUE, reps = 50), "exactly 2 counts")
    expect_reject(SinglePropCI(data = data.frame(x = c(0, 0)), resp = "x", areCounts = TRUE, reps = 50), "both zero")
})
test("reps must be a positive integer; confidence level in (0, 100)", {
    expect_reject(twomeanhtest(data = clean_data, vars = "score", group = "group", reps = 0), "reps|Reps|greater|minimum|min")
    expect_reject(twomeanhtest(data = clean_data, vars = "score", group = "group", reps = -5), "reps|Reps|greater|minimum|min")
    expect_reject(SingleMeanCI(data = clean_data, resp = "score", reps = 50, confLevel = 100), "confLevel|maximum|max")
    expect_reject(modelBased(distro = "tdistro", dF = 10, CIBool = TRUE, confLevel = 0), "confLevel|minimum|min")
})
test("p = 0 no longer emits infer's warning", {
    w <- NULL
    withCallingHandlers(
        SinglePropHTest(data = data.frame(x = c(19, 1)), resp = "x", areCounts = TRUE, testValue = 0.5,
                        alt = "greater", reps = 100, seedBool = TRUE, rngSeed = 1),
        warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    stopifnot(!any(grepl("p-value of 0", w)))
})
test("dot stacks rest on the axis", {
    for (r in list(
        twomeanhtest(data = clean_data, vars = "score", group = "group", hypothesis = "different",
                     reps = 3, dotHist = "dotplot", seedBool = TRUE, rngSeed = 1),
        SingleMeanCI(data = clean_data, resp = "score", reps = 3, dotHist = "dotplot",
                     seedBool = TRUE, rngSeed = 1))) {
        b <- ggplot2::ggplot_build(plot(r))
        d <- b$data[[1]]
        stopifnot(min(d$y) == 0.5)                       # bottom dot centred half a unit up
        stopifnot(b$layout$panel_params[[1]]$y.range[2] >= 8)  # axis floor
    }
})
test("a single simulation keeps the observed value on the axis", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group", hypothesis = "different",
                      reps = 1, dotHist = "dotplot", seedBool = TRUE, rngSeed = 1)
    p <- plot(r); b <- ggplot2::ggplot_build(p)
    xr <- b$layout$panel_params[[1]]$x.range
    obs <- sim_img(r)$state$obs_stat
    stopifnot(obs >= xr[1], obs <= xr[2])
    stopifnot(count_layers(p, "GeomVline") == 1)   # dashed line, not an edge arrow
})

cat("\n=== Simulations are reused when only display options change ===\n")
test("cached_sims: same key reuses, changed key recomputes", {
    holder <- jmvcore::Table$new(options = jmvcore::Options$new(), name = "t", title = "t")
    n <- 0
    f <- function() { n <<- n + 1; data.frame(stat = rnorm(5)) }
    a <- cached_sims(holder, list(x = 1:3, reps = 5), f)
    b <- cached_sims(holder, list(x = 1:3, reps = 5), f)
    stopifnot(identical(a, b), n == 1)
    c3 <- cached_sims(holder, list(x = 1:4, reps = 5), f)
    stopifnot(!identical(a, c3), n == 2)
    cached_sims(holder, list(x = 1:4, reps = 5), f, slot = "other"); stopifnot(n == 3)
    cached_sims(holder, list(x = 1:4, reps = 5), f, slot = "other"); stopifnot(n == 3)
})
test("re-running an analysis with unchanged inputs keeps the simulations (unseeded)", {
    # jamovi calls .run() again on every option change; the cache lives in
    # the results table's state, which jamovi preserves unless a clearWith
    # option changed
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group", hypothesis = "different", reps = 200,
                      dotHist = "histogram", seedBool = FALSE)
    s1 <- sim_img(r)$state$df$stat
    stopifnot(identical(r$htest$state$sims$sims$stat, s1))
    r$analysis$.__enclos_env__$private$.run()
    stopifnot(identical(sim_img(r)$state$df$stat, s1))
    # and the cache key covers the data: a different data set gets new draws
    r2 <- twomeanhtest(data = na_data, vars = "score", group = "group", hypothesis = "different", reps = 200,
                       dotHist = "histogram", seedBool = FALSE)
    stopifnot(!identical(sim_img(r2)$state$df$stat, s1))
})
test("edited data invalidates the cache even when no option changed", {
    holder <- jmvcore::Table$new(options = jmvcore::Options$new(), name = "t", title = "t")
    k1 <- list(dep = c(1, 2, 3), reps = 10); k2 <- list(dep = c(1, 2, 4), reps = 10)
    a <- cached_sims(holder, k1, function() data.frame(stat = 1)); b <- cached_sims(holder, k2, function() data.frame(stat = 2))
    stopifnot(a$stat == 1, b$stat == 2)
})

cat("\n=== Option-change matrix: only simulation inputs redraw the simulations ===\n")
# For every analysis, change each option in place and re-run (as jamovi does on
# every option change). Options that do not determine the simulations must leave
# them untouched; options that do must produce new draws (unseeded).
sim_state <- function(r) sim_img(r)$state$df$stat
set_opt <- function(r, name, value) { o <- r$analysis$options$option(name); o$value <- value }
rerun <- function(r) r$analysis$.__enclos_env__$private$.run()
check_matrix <- function(label, r, keep, change) {
    test(paste("matrix:", label), {
        base <- sim_state(r)
        for (nm in names(keep)) {
            set_opt(r, nm, keep[[nm]]); rerun(r)
            if (!identical(sim_state(r), base)) stop(paste("changing", nm, "redrew the simulations"))
        }
        for (nm in names(change)) {
            before <- sim_state(r); set_opt(r, nm, change[[nm]]); rerun(r)
            if (identical(sim_state(r), before)) stop(paste("changing", nm, "did not redraw the simulations"))
        }
    })
}
cat_tab <- data.frame(group = factor(rep(c("A","B"), each = 20)),
                      outcome = factor(c(rep("Yes", 15), rep("No", 5), rep("Yes", 8), rep("No", 12))))
check_matrix("twomeanhtest",
    twomeanhtest(data = clean_data, vars = "score", group = "group", hypothesis = "different", reps = 150, seedBool = FALSE),
    keep = list(hypothesis = "oneGreater", dotHist = "histogram", showCounts = TRUE, desc = TRUE, plots = TRUE),
    change = list(reps = 160, seedBool = TRUE, rngSeed = 5))
check_matrix("twomeanCI",
    twomeanCI(data = clean_data, vars = "score", group = "group", reps = 150, seedBool = FALSE),
    keep = list(confLevel = 90, ciType = "bootse", dotHist = "histogram", showCounts = TRUE, desc = TRUE, plots = TRUE),
    change = list(reps = 160))
check_matrix("pairedmeanhtest",
    pairedmeanhtest(data = clean_data, pairs = list(list(i1 = "measure1", i2 = "measure2")), hypothesis = "different", reps = 150, seedBool = FALSE),
    keep = list(hypothesis = "twoGreater", dotHist = "histogram", showCounts = TRUE, desc = TRUE, plots = TRUE),
    change = list(reps = 160))
check_matrix("pairedmeanCI",
    pairedmeanCI(data = clean_data, pairs = list(list(i1 = "measure1", i2 = "measure2")), reps = 150, seedBool = FALSE),
    keep = list(confLevel = 99, ciType = "bootse", dotHist = "histogram", showCounts = TRUE),
    change = list(reps = 160))
check_matrix("multimeanhtest",
    multimeanhtest(data = multi_data, vars = "score", group = "group", reps = 150, seedBool = FALSE),
    keep = list(dotHist = "histogram", showCounts = TRUE, desc = TRUE, plots = TRUE),
    change = list(reps = 160))
check_matrix("slopehtest",
    slopehtest(data = reg_data, dep = "y", indep = "x", hypothesis = "notequal", reps = 150, seedBool = FALSE),
    keep = list(hypothesis = "greater", dotHist = "histogram", showCounts = TRUE, coef = TRUE, modelfit = TRUE, plots = TRUE),
    change = list(reps = 160))
check_matrix("slopeCI",
    slopeCI(data = reg_data, dep = "y", indep = "x", reps = 150, seedBool = FALSE),
    keep = list(confLevel = 90, ciType = "bootse", dotHist = "histogram", showCounts = TRUE, coef = TRUE),
    change = list(reps = 160))
check_matrix("SingleMeanCI",
    SingleMeanCI(data = clean_data, resp = "score", reps = 150, seedBool = FALSE),
    keep = list(confLevel = 90, ciType = "bootse", dotHist = "histogram", showCounts = TRUE),
    change = list(reps = 160))
check_matrix("SinglePropHTest",
    SinglePropHTest(data = cat_data, resp = "outcome", testValue = 0.5, alt = "notequal", reps = 150, seedBool = FALSE),
    keep = list(alt = "greater", dotHist = "histogram", showCounts = TRUE),
    change = list(testValue = 0.4, reps = 160))
check_matrix("SinglePropCI",
    SinglePropCI(data = cat_data, resp = "outcome", reps = 150, seedBool = FALSE),
    keep = list(confLevel = 90, ciType = "bootse", dotHist = "histogram", showCounts = TRUE),
    change = list(reps = 160))
check_matrix("TwoPropHTest",
    TwoPropHTest(data = cat_tab, rows = "group", cols = "outcome", hypothesis = "different", reps = 150, seedBool = FALSE, compare = "rows"),
    keep = list(hypothesis = "oneGreater", dotHist = "histogram", showCounts = TRUE, obs = TRUE, exp = TRUE, pcRow = TRUE),
    change = list(compare = "columns", reps = 160))
check_matrix("TwoPropCI",
    TwoPropCI(data = cat_tab, rows = "group", cols = "outcome", reps = 150, seedBool = FALSE, compare = "rows"),
    keep = list(confLevel = 90, ciType = "bootse", dotHist = "histogram", showCounts = TRUE, exp = TRUE),
    change = list(compare = "columns", reps = 160))
check_matrix("ContTabHTest",
    ContTabHTest(data = cat_tab, rows = "group", cols = "outcome", reps = 150, seedBool = FALSE, compare = "rows"),
    keep = list(dotHist = "histogram", showCounts = TRUE, obs = TRUE, exp = TRUE, pcRow = TRUE, pcCol = TRUE, pcTot = TRUE),
    change = list(reps = 160))

test("cache holders' clearWith lists contain only simulation inputs (jamovi wipes state on these)", {
    holders <- c(twomeanhtest = "htest", pairedmeanhtest = "htest", multimeanhtest = "htest", slopehtest = "htest",
                 twomeanCI = "CITable", pairedmeanCI = "CITable", slopeCI = "CITable", singlemeanCI = "simtable",
                 singleprophtest = "simtable", singlepropCI = "simtable", twoprophtest = "simtable",
                 twopropCI = "simtable", conttabhtest = "simtable")
    never <- c("hypothesis", "alt", "confLevel", "ciType", "dotHist", "showCounts", "desc", "plots",
               "coef", "modelfit", "obs", "exp", "pcRow", "pcCol", "pcTot")
    for (a in names(holders)) {
        txt <- readLines(sprintf("jamovi/%s.r.yaml", a), warn = FALSE, encoding = "UTF-8")
        start <- grep(sprintf("^\\s+- name:\\s+%s\\s*$", holders[[a]]), txt)[1]
        nxt <- grep("^\\s+- name:", txt); nxt <- nxt[nxt > start]; end <- if (length(nxt)) nxt[1] - 1 else length(txt)
        blk <- txt[start:end]; cw <- grep("^\\s+clearWith:", blk)
        items <- sub("^\\s+- ", "", grep("^\\s+- \\w+\\s*$", blk[(cw + 1):length(blk)], value = TRUE))
        items <- items[seq_len(match(FALSE, grepl("^\\w+$", items), nomatch = length(items) + 1) - 1)]
        bad <- intersect(items, never)
        if (length(bad)) stop(paste0(a, "/", holders[[a]], " clearWith has display/inference options: ", paste(bad, collapse = ", ")))
    }
})

cat("\n=== Formula-free resampling matches the model-based statistics ===\n")
test("f_stat equals anova(lm()) F", {
    f_ref <- anova(lm(score ~ group, data = multi_data))$F[1]
    assert_close(f_stat(multi_data$score, multi_data$group), f_ref, 1e-9)
    # unbalanced groups too
    d <- multi_data[-c(1:7), ]
    f_ref <- anova(lm(score ~ group, data = d))$F[1]
    assert_close(f_stat(d$score, d$group), f_ref, 1e-9)
})
test("permutation and bootstrap helpers return well-formed distributions", {
    set.seed(5)
    lv <- levels(clean_data$group)
    pm <- permute_diff_means(clean_data$score, clean_data$group, lv, 400)
    stopifnot(nrow(pm) == 400, all(c("replicate", "stat") %in% names(pm)))
    # under the null the permutation distribution is centred near zero
    stopifnot(abs(mean(pm$stat)) < 0.3)
    bt <- bootstrap_diff_means(clean_data$score, clean_data$group, lv, 400)
    stopifnot(nrow(bt) == 400)
    obs <- mean(clean_data$score[clean_data$group == lv[1]]) -
        mean(clean_data$score[clean_data$group == lv[2]])
    stopifnot(abs(mean(bt$stat) - obs) < 0.4)
    pf <- permute_F(multi_data$score, multi_data$group, 400)
    stopifnot(nrow(pf) == 400, all(pf$stat >= 0))
})

# ============================================================
# Missing values (bug-hunt round 7)
# ============================================================
cat("\n=== Missing values ===\n")

# jamovi pads shorter columns with blanks, so a 2-row count column beside
# longer data arrives with trailing NAs; the "values are counts" check
# must look at the non-missing values (it used to reject "found 5 rows")
test("counts column padded with blanks is accepted", {
    d <- data.frame(x = c(12, 8, NA, NA, NA), other = 1:5)
    r <- SinglePropHTest(data = d, resp = "x", areCounts = TRUE, testValue = 0.5, reps = 200, seedBool = TRUE, rngSeed = 1)
    assert_close(results_table(r)$obsProp, 0.6, 1e-9)
    r <- SinglePropCI(data = d, resp = "x", areCounts = TRUE, reps = 200, seedBool = TRUE, rngSeed = 1)
    assert_close(results_table(r)$obsProp, 0.6, 1e-9)
    expect_reject(SinglePropHTest(data = data.frame(x = c(3, NA)), resp = "x", areCounts = TRUE, reps = 50), "exactly 2 counts")
})

# padded rows (NA in every column) are dropped; a missing count beside
# present row / column values is an error, not an empty cell (it used
# to reach chisq.test as NA: "all entries of 'x' must be nonnegative")
test("table analyses: padded rows dropped, a missing count rejected", {
    padded <- data.frame(g = c("A","A","B","B", NA, NA), o = c("y","n","y","n", NA, NA), n = c(10, 5, 4, 11, NA, NA))
    holed  <- data.frame(g = c("A","A","B","B"), o = c("y","n","y","n"), n = c(10, NA, 4, 11))
    for (f in list(TwoPropHTest, TwoPropCI, ContTabHTest)) {
        r <- f(data = padded, rows = "g", cols = "o", counts = "n", reps = 100, seedBool = TRUE, rngSeed = 1)
        if (r$freqs$asDF[[".total[count]"]][3] != 30) stop("padded table total should be 30")
        expect_reject(f(data = holed, rows = "g", cols = "o", counts = "n", reps = 100), "missing a value in a row where")
    }
})

# a group whose values are all missing keeps its level (droplevels runs
# before the NA rows go); every resampled statistic was NaN and infer's
# get_p_value() threw "All calculated statistics were NaN"
test("two-means analyses: a group with every value missing gets a footnote", {
    d <- data.frame(score = c(rnorm(10), rep(NA, 10)), group = factor(rep(c("A", "B"), each = 10)))
    r <- twomeanhtest(data = d, vars = "score", group = "group", reps = 100, desc = TRUE)
    expect_footnote(r$htest, "no non-missing")
    r <- twomeanCI(data = d, vars = "score", group = "group", reps = 100)
    expect_footnote(r$CITable, "no non-missing")
})

# a single complete pair reached infer's specify() -> t.test() and died
# with "not enough 'x' observations"
test("paired analyses: fewer than 2 complete pairs gets a footnote", {
    set.seed(5)
    d <- data.frame(pre = c(1.2, NA, NA, NA), post = c(0.8, 1.1, NA, 2))
    r <- pairedmeanhtest(data = d, pairs = list(list(i1 = "pre", i2 = "post")), reps = 100)
    expect_footnote(r$htest, "At least 2 complete pairs")
    r <- pairedmeanCI(data = d, pairs = list(list(i1 = "pre", i2 = "post")), reps = 100)
    expect_footnote(r$CITable, "At least 2 complete pairs")
    d$pre[2] <- 0.5
    r <- pairedmeanhtest(data = d, pairs = list(list(i1 = "pre", i2 = "post")), reps = 100, seedBool = TRUE, rngSeed = 1, desc = TRUE)
    if (desc_table(r)$num[1] != 2) stop("two complete pairs expected")
    assert_not_na(p_num(results_table(r)$p))
})

# the descriptives, the observed statistic and the test all use the same
# complete cases (missing in the response OR the group)
test("two-means descriptives and statistic use the complete cases", {
    set.seed(3)
    d <- data.frame(score = rnorm(40, 10, 2), group = factor(rep(c("A", "B"), each = 20)))
    d$score[c(1, 25)] <- NA; d$group[c(2, 30)] <- NA
    cc <- d[complete.cases(d), ]
    r <- twomeanhtest(data = d, vars = "score", group = "group", reps = 200, seedBool = TRUE, rngSeed = 1, desc = TRUE)
    dt <- desc_table(r)
    if (!all(c(dt$num1, dt$num2) == as.vector(table(cc$group)))) stop("group sizes differ from complete cases")
    assert_close(results_table(r)$md, diff(rev(tapply(cc$score, cc$group, mean))), 1e-9)
    assert_not_na(p_num(results_table(r)$p))
})

# constant differences (every post = pre + c, or pre == post): infer's
# specify() runs t.test(), which refuses them ("data are essentially
# constant"); the sign-flip permutation and the bootstrap are drawn directly
test("paired analyses: constant differences give the sign-flip null / point CI", {
    d <- data.frame(pre = 1:10 + 0.5, post = 1:10 + 2.5)
    r <- pairedmeanhtest(data = d, pairs = list(list(i1 = "pre", i2 = "post")), reps = 2000, seedBool = TRUE, rngSeed = 1, hypothesis = "twoGreater")
    assert_close(results_table(r)$md, -2, 1e-9)
    p <- p_num(results_table(r)$p)            # exact sign-test p = 0.5^10
    if (p > 0.01) stop(paste("p should be about 0.001, got", p))
    st <- sim_img(r)$state$df$stat
    if (length(unique(st)) > 11 || any(abs(st) > 2 + 1e-9)) stop("sign-flip null should take at most 11 values in [-2, 2]")
    r <- pairedmeanCI(data = d, pairs = list(list(i1 = "pre", i2 = "post")), reps = 200, seedBool = TRUE, rngSeed = 1)
    rt <- results_table(r); assert_close(rt$cil, -2, 1e-9); assert_close(rt$ciu, -2, 1e-9)
    d$post <- d$pre
    r <- pairedmeanhtest(data = d, pairs = list(list(i1 = "pre", i2 = "post")), reps = 100)
    assert_close(p_num(results_table(r)$p), 1, 1e-9)
})

# every pair gets its own row, simulation plot and descriptive plot (only
# the first pair used to be analysed; the other rows stayed blank)
test("paired analyses: every pair is analysed, one simulation plot per pair", {
    set.seed(2); d <- data.frame(a = rnorm(20), b = rnorm(20, 0.5), c = rnorm(20, 1))
    pairs <- list(list(i1 = "a", i2 = "b"), list(i1 = "a", i2 = "c"), list(i1 = "b", i2 = "a"))
    r <- pairedmeanhtest(data = d, pairs = pairs, reps = 300, seedBool = TRUE, rngSeed = 1, desc = TRUE, plots = TRUE, hypothesis = "oneGreater")
    tb <- r$htest$asDF
    if (any(is.na(tb$md)) || nrow(tb) != 3) stop("every pair needs a filled row")
    assert_close(tb$md[3], -tb$md[1], 1e-12)
    if (length(r$simplot$itemKeys) != 3) stop("one simulation plot per pair")
    for (k in r$simplot$itemKeys) if (is.null(r$simplot$get(key = k)$state)) stop("plot state missing for a pair")
    if (nrow(desc_table(r)) != 6) stop("two descriptive rows per pair")
    if (!inherits(plot(r), "ggplot")) stop("plot(r) must return the first pair's plot")
    r1 <- pairedmeanhtest(data = d, pairs = pairs[1], reps = 300, seedBool = TRUE, rngSeed = 1, hypothesis = "oneGreater")
    if (!identical(r1$htest$asDF$p[1], tb$p[1])) stop("single-pair result changed")
    r2 <- pairedmeanCI(data = d, pairs = pairs, reps = 300, seedBool = TRUE, rngSeed = 1, desc = TRUE, plots = TRUE)
    tb2 <- r2$CITable$asDF
    if (any(is.na(tb2$cil)) || length(r2$simplot$itemKeys) != 3) stop("CI: every pair needs a row and a plot")
    assert_close(tb2$cil[3], -tb2$ciu[1], 1e-12)
})

# a group whose values are all missing (level kept, rows gone) made the
# descriptives table choke on a zero-length summary
test("F test: a group with every value missing gets n = 0 and a footnote", {
    d <- data.frame(y = c(rnorm(6), rep(NA, 3)), g = factor(rep(c("A", "B", "C"), each = 3)))
    r <- multimeanhtest(data = d, vars = "y", group = "g", reps = 50, desc = TRUE, plots = TRUE)
    expect_footnote(r$htest, "at least 2 observations")
    dt <- desc_table(r)
    if (dt$num[dt$group == "C" | grepl("C$", rownames(dt))][1] != 0) stop("group C should report n = 0")
})

# a bootstrap resample of a handful of rows can contain one group only;
# infer errored ("G2 is not a level of the explanatory variable") when
# that was the only replicate, instead of dropping it
test("two-proportion CI: one rep on a tiny table runs; bootstrap matches infer's", {
    d <- data.frame(g = c("A","A","B","B"), o = c("y","n","y","n"), n = c(1, 2, 2, 0))
    for (reps in c(1, 2, 5)) {
        r <- TwoPropCI(data = d, rows = "g", cols = "o", counts = "n", reps = reps, seedBool = TRUE, rngSeed = 25)
        assert_not_na(results_table(r)$obsDiff)
    }
    d2 <- data.frame(g = c("A","A","B","B"), o = c("y","n","y","n"), n = c(18, 12, 9, 21))
    r <- TwoPropCI(data = d2, rows = "g", cols = "o", counts = "n", reps = 3000, seedBool = TRUE, rngSeed = 1)
    st <- r$Plot$state$df$stat
    df <- tibble::tibble(Group = c("G1","G2","G1","G2"), Outcome = c("O1","O1","O2","O2"), Count = c(12, 21, 18, 9)) |> tidyr::uncount(Count)
    set.seed(2)
    ib <- df |> infer::specify(Outcome ~ Group, success = "O1") |> infer::generate(reps = 3000, type = "bootstrap") |> infer::calculate(stat = "diff in props", order = c("G1", "G2"))
    if (suppressWarnings(ks.test(st, ib$stat)$p.value) < 0.001) stop("bootstrap distribution differs from infer's")
    assert_close(mean(st), -0.3, 0.02)
})

# a jamovi column of integers with a nominal measure type arrives as a
# factor with a "values" attribute; the single-mean CI fed it to median()
# ("need numeric data") while every other mean analysis converts it
test("single mean CI: nominal-integer column (jamovi factor with values) is converted", {
    v <- c(12, 8, 5, 7, 9, 11, 6, 10)
    f <- factor(v, levels = sort(unique(v))); attr(f, "values") <- sort(unique(v))
    d <- data.frame(y = f)
    r <- SingleMeanCI(data = d, resp = "y", reps = 100, seedBool = TRUE, rngSeed = 1)
    assert_close(results_table(r)$obsMean, mean(v), 1e-9)
    assert_not_na(results_table(r)$cil)
    ft <- factor(c("a", "b", "a")); attr(ft, "values") <- NULL
    expect_reject(SingleMeanCI(data = data.frame(y = ft), resp = "y", reps = 50), "not numeric|requires a numeric")
})

# plot states saved by an earlier build lack xlab / obs_label / area;
# the render callbacks must fall back to the old wording, not fail
test("render callbacks cope with a plot state from an earlier build", {
    r <- twomeanhtest(data = clean_data, vars = "score", group = "group", reps = 100, seedBool = TRUE, rngSeed = 1)
    img <- sim_img(r); st <- img$state; st$xlab <- NULL; st$obs_label <- NULL; img$setState(st)
    p <- img$plot$fun(); b <- ggplot2::ggplot_build(p)
    if (!grepl("group 1", p$labels$x)) stop("expected the generic axis label, got: ", p$labels$x)
    r <- SinglePropCI(data = cat_data, resp = "outcome", reps = 100, seedBool = TRUE, rngSeed = 1)
    img <- r$Plot; st <- img$state; st$xlab <- NULL; st$obs_label <- NULL; img$setState(st)
    p <- img$plot$fun(); ggplot2::ggplot_build(p)
    if (p$labels$x != "proportion") stop("expected 'proportion', got: ", p$labels$x)
    r <- modelBased(distro = "tdistro", dF = 9, areaBool = TRUE, obsStat = 2)
    img <- r$Plot; st <- img$state; st$area <- NULL; img$setState(st)
    p <- img$plot$fun(); ggplot2::ggplot_build(p)
    r2 <- TwoPropHTest(data = cat_data, rows = "group", cols = "outcome", reps = 100, seedBool = TRUE, rngSeed = 1)
    img <- r2$Plot; st <- img$state; st$xlab <- NULL; img$setState(st)
    p <- img$plot$fun(); ggplot2::ggplot_build(p)
    if (!grepl("group 1", p$labels$x)) stop("two-prop fallback label missing")
})

# ============================================================
# Jamovi 2.7 formula sandbox
# ============================================================
# Jamovi 2.7 replaces stats::as.formula with jmvcore::asFormula, which
# rejects any formula calling a function outside a short allowlist.
# infer::specify() builds response_variable(x) ~ explanatory_variable(x)
# internally for numeric-response / factor-explanatory specifications,
# which broke the two-means and multiple-means analyses in the app.
# Replicate the sandbox here (same allowlist as jmvcore 2.7.35) and run
# every analysis under it.  Runs last because the override is global.

cat("\n=== All analyses under a replica of Jamovi 2.7's formula sandbox ===\n")
sandbox_allowed <- c("log", "log2", "log10", "log1p", "exp", "expm1", "sqrt", "abs",
    "sign", "sin", "cos", "tan", "asin", "acos", "atan", "atan2", "floor", "ceiling",
    "round", "trunc", "mean", "sd", "var", "median", "min", "max", "sum", "length",
    "rank", "scale", "poly", "ns", "bs", "I", "cbind", "rbind", "c", "as.numeric",
    "as.integer", "as.factor", "as.character", "factor", "ordered", "cut", "offset",
    "relevel", "interaction", "pmin", "pmax", "ifelse", "Error", "Surv", "strata",
    "cluster", "frailty", "tt", "pspline", "Hist", "Event", "s", "te", "ti", "t2",
    "rcs", "lsp", "pdSymm", "pdDiag", "item")
orig_as_formula <- stats::as.formula
sandbox_as_formula <- function(object, env = parent.frame()) {
    fmla <- orig_as_formula(object, env)
    txt <- paste(deparse(fmla, width.cutoff = 500), collapse = " ")
    txt <- gsub("`[^`]+`", "", txt)
    calls <- regmatches(txt, gregexpr("\\.?[a-zA-Z_][a-zA-Z0-9_.]*(?=\\s*\\()", txt, perl = TRUE))[[1]]
    bad <- setdiff(calls, sandbox_allowed)
    if (length(bad))
        stop("Security violation: Unsafe function '", bad[1], "'.", call. = FALSE)
    fmla
}
install_sandbox <- function(fun) {
    ns <- asNamespace("stats")
    unlockBinding("as.formula", ns)
    assign("as.formula", fun, envir = ns)
    lockBinding("as.formula", ns)
}
install_sandbox(sandbox_as_formula)
test("sandbox replica rejects infer-style formulas", {
    msg <- tryCatch({ stats::as.formula("response_variable(x) ~ explanatory_variable(x)"); "" },
                    error = function(e) conditionMessage(e))
    stopifnot(grepl("Security violation", msg))
    stats::as.formula("y ~ x + log(z)")  # plain formulas still work
})
sandbox_run <- function(name, expr) test(paste("sandbox:", name), { force(expr); invisible() })
sandbox_run("twomeanhtest", twomeanhtest(data = clean_data, vars = "score", group = "group",
    hypothesis = "different", reps = 100, seedBool = TRUE, rngSeed = 1, desc = TRUE, plots = TRUE))
sandbox_run("twomeanCI", twomeanCI(data = clean_data, vars = "score", group = "group",
    reps = 100, seedBool = TRUE, rngSeed = 1, desc = TRUE, plots = TRUE))
sandbox_run("pairedmeanhtest", pairedmeanhtest(data = clean_data, pairs = list(list(i1 = "measure1", i2 = "measure2")),
    hypothesis = "different", reps = 100, seedBool = TRUE, rngSeed = 1))
sandbox_run("pairedmeanCI", pairedmeanCI(data = clean_data, pairs = list(list(i1 = "measure1", i2 = "measure2")),
    reps = 100, seedBool = TRUE, rngSeed = 1))
sandbox_run("multimeanhtest", multimeanhtest(data = multi_data, vars = "score", group = "group",
    reps = 100, seedBool = TRUE, rngSeed = 1, desc = TRUE, plots = TRUE))
sandbox_run("slopehtest", slopehtest(data = reg_data, dep = "y", indep = "x",
    hypothesis = "notequal", reps = 100, seedBool = TRUE, rngSeed = 1, plots = TRUE))
sandbox_run("slopeCI", slopeCI(data = reg_data, dep = "y", indep = "x",
    reps = 100, seedBool = TRUE, rngSeed = 1, plots = TRUE))
sandbox_run("SingleMeanCI", SingleMeanCI(data = clean_data, resp = "score", reps = 100, seedBool = TRUE, rngSeed = 1))
sandbox_run("SinglePropHTest", SinglePropHTest(data = cat_data, resp = "outcome", testValue = 0.5,
    reps = 100, seedBool = TRUE, rngSeed = 1))
sandbox_run("SinglePropCI", SinglePropCI(data = cat_data, resp = "outcome", reps = 100, seedBool = TRUE, rngSeed = 1))
sandbox_run("TwoPropHTest", TwoPropHTest(data = cat_data, rows = "group", cols = "outcome",
    hypothesis = "different", reps = 100, seedBool = TRUE, rngSeed = 1, compare = "rows"))
sandbox_run("TwoPropCI", TwoPropCI(data = cat_data, rows = "group", cols = "outcome",
    reps = 100, seedBool = TRUE, rngSeed = 1, compare = "rows"))
sandbox_run("ContTabHTest", ContTabHTest(data = cat_data, rows = "group", cols = "outcome",
    reps = 100, seedBool = TRUE, rngSeed = 1, compare = "rows"))
sandbox_run("modelBased", modelBased(distro = "fdistro", dF = 2, dF2 = 20, areaBool = TRUE, obsStat = 3))
install_sandbox(orig_as_formula)

cat("\n\n============================\n")
cat("Results:", passed, "passed,", failed, "failed\n")
cat("============================\n\n")

if (failed > 0) quit(status = 1)
