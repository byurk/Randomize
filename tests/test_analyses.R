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
    # Should be significant (true slope = 2)
    stopifnot(tbl$p < 0.05)
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
    assert_not_na(tbl$pval, "p-value")
    stopifnot(tbl$pval >= 0 && tbl$pval <= 1)
})

cat("\n=== Testing MultiMeanHTest ===\n")
test("clean data", {
    r <- multimeanhtest(data = multi_data, vars = "score", group = "group",
                        reps = 500, dotHist = "histogram",
                        seedBool = TRUE, rngSeed = 123)
    tbl <- r$htest$asDF
    assert_not_na(tbl$p, "p-value")
    stopifnot(tbl$p >= 0 && tbl$p <= 1)
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
    assert_not_na(tbl$pval, "p-value")
    stopifnot(tbl$pval >= 0 && tbl$pval <= 1)
})

cat("\n=== Testing ContTabHTest ===\n")
test("clean data", {
    r <- ContTabHTest(data = cat_data, rows = "group", cols = "outcome",
                      reps = 500, dotHist = "histogram",
                      seedBool = TRUE, rngSeed = 123, compare = "rows")
    tbl <- r$simtable$asDF
    assert_not_na(tbl$pval, "p-value")
    stopifnot(tbl$pval >= 0 && tbl$pval <= 1)
})

cat("\n\n============================\n")
cat("Results:", passed, "passed,", failed, "failed\n")
cat("============================\n\n")

if (failed > 0) quit(status = 1)
