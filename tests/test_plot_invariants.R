# Pedagogical invariant tests for the distribution plots
#
# These tests inspect the GEOMETRY ggplot actually draws (via
# ggplot_build), not just the input data, and assert the properties the
# plots must never violate:
#
#   1. Side purity: no histogram bar or dot stack ever mixes simulations
#      that count toward the p-value ("extreme") with ones that do not.
#      Every red element lies entirely on the extreme side of the
#      observed-value line; every black element entirely on the other.
#   2. Count conservation: bar heights / dot counts sum to the number of
#      simulations, and the red elements sum to exactly the p-value
#      numerator (so a student counting red dots recovers the p-value).
#   3. Toggle stability: bare/no-text modes draw bars at identical
#      positions to the full plot.
#   4. Bootstrap lattice CI lines never slice through a bar.
#
# Run with: Rscript tests/test_plot_invariants.R

suppressWarnings(suppressMessages(devtools::load_all('.', quiet = TRUE)))
suppressMessages(library(ggplot2))

source("tests/visual/scenarios.R")

passed <- 0
failed <- 0

test <- function(name, expr) {
    tryCatch({
        expr
        cat("  PASS:", name, "\n")
        passed <<- passed + 1
    }, error = function(e) {
        cat("  FAIL:", name, "-", e$message, "\n")
        failed <<- failed + 1
    })
}

RED <- "#ff8c8c"

# Build a plot, failing if any drawn geometry was dropped (e.g. a dot
# clipped by the axis limits); returns the built plot for inspection
build_checked <- function(p) {
    warns <- character()
    built <- withCallingHandlers(
        ggplot_build(p),
        warning = function(w) {
            warns <<- c(warns, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    bad <- grepl("[Rr]emoved", warns)
    if (any(bad))
        stop(paste("geometry clipped:", warns[bad][1]))
    built
}

# The drawn distribution must fill most of the panel: at least 55% of
# the x range (axis may legitimately extend toward a distant observed
# value or CI line, but never so far the data become a sliver) and at
# least 70% of the y range
assert_axis_economy <- function(built, el) {
    xr <- built$layout$panel_params[[1]]$x.range
    yr <- built$layout$panel_params[[1]]$y.range
    x_use <- (max(el$xmax) - min(el$xmin)) / diff(xr)
    if (x_use < 0.55)
        stop(paste0("distribution occupies only ",
                    round(100 * x_use), "% of the x axis"))
    if ("ytop" %in% names(el)) {
        y_use <- max(el$ytop) / yr[2]
        if (y_use < 0.7)
            stop(paste0("distribution occupies only ",
                        round(100 * y_use), "% of the y axis"))
    }
}

# Extract per-element (bar or ellipse) x-extent and fill from the first
# layer of a built plot.  For geom_rect each row is a bar; for
# geom_ellipse the rows are polygon vertices grouped per dot.
layer_elements <- function(built) {
    d <- built$data[[1]]
    if (all(c("xmin", "xmax") %in% names(d)) && !all(is.na(d$xmin))) {
        data.frame(xmin = d$xmin, xmax = d$xmax, fill = d$fill,
                   height = d$ymax - d$ymin, ytop = d$ymax)
    } else {
        agg <- aggregate(d$x, by = list(group = d$group), FUN = min)
        names(agg) <- c("group", "xmin")
        agg$xmax <- aggregate(d$x, by = list(d$group), FUN = max)$x
        agg$fill <- aggregate(d$fill, by = list(d$group), FUN = function(f) f[1])$x
        agg$height <- 1
        agg$ytop <- aggregate(d$y, by = list(d$group), FUN = max)$x
        agg
    }
}

vline_positions <- function(built) {
    for (d in built$data) {
        if ("xintercept" %in% names(d)) return(d$xintercept)
    }
    numeric(0)
}

check_null_invariants <- function(s, mode) {
    df <- data.frame(stat = s$stats)
    p <- plot_null_dist(df, s$obs, s$direction, mode,
                        xlab = s$xlab, obs_label = "Observed\nValue")
    built <- build_checked(p)
    el <- layer_elements(built)
    assert_axis_economy(built, el)

    # Determine which tail is shaded, mirroring plot_null_dist
    if (s$direction == "less") {
        ptail <- "lt"
    } else if (s$direction == "greater") {
        ptail <- "rt"
    } else {
        ptail <- if (mean(s$stats > s$obs) < mean(s$stats < s$obs)) "rt" else "lt"
    }
    extreme <- if (ptail == "lt") s$stats <= s$obs else s$stats >= s$obs

    span <- diff(range(c(s$stats, s$obs)))
    eps <- span * 1e-9 + 1e-12

    red <- el[el$fill == RED, ]
    black <- el[el$fill != RED, ]

    # 1. Side purity of everything actually drawn
    if (ptail == "rt") {
        if (nrow(red) > 0 && min(red$xmin) < s$obs - eps)
            stop("red element extends left of the observed value")
        if (nrow(black) > 0 && max(black$xmax) > s$obs + eps)
            stop("black element extends right of the observed value")
    } else {
        if (nrow(red) > 0 && max(red$xmax) > s$obs + eps)
            stop("red element extends right of the observed value")
        if (nrow(black) > 0 && min(black$xmin) < s$obs - eps)
            stop("black element extends left of the observed value")
    }

    # 2. Count conservation: total and per-color
    if (abs(sum(el$height) - length(s$stats)) > 1e-6)
        stop(paste("element heights sum to", sum(el$height),
                   "but there are", length(s$stats), "simulations"))
    if (abs(sum(red$height) - sum(extreme)) > 1e-6)
        stop(paste("red elements sum to", sum(red$height),
                   "but", sum(extreme), "simulations are extreme"))

    # 3. Bars (or dot columns) must not overlap each other.  Dots stacked
    # in the same column share identical x-extents, so collapse to unique
    # columns first.
    cols <- unique(el[, c("xmin", "xmax")])
    if (nrow(cols) > 1) {
        ord <- cols[order(cols$xmin), ]
        if (any(ord$xmin[-1] - ord$xmax[-nrow(ord)] < -eps - 1e-9 * span))
            stop("elements overlap horizontally")
    }
}

cat("\n=== Null distribution invariants (histogram) ===\n")
null_sc <- make_null_scenarios()
for (nm in names(null_sc)) {
    test(paste0(nm, " [hist]"), check_null_invariants(null_sc[[nm]], "histogram"))
}

cat("\n=== Null distribution invariants (dotplot) ===\n")
for (nm in names(null_sc)) {
    test(paste0(nm, " [dot]"), check_null_invariants(null_sc[[nm]], "dotplot"))
}

cat("\n=== Toggle stability: bare mode draws identical bars ===\n")
for (nm in c("prop_n20_greater", "cont_two_sided", "chisq_2x2")) {
    test(nm, {
        s <- null_sc[[nm]]
        df <- data.frame(stat = s$stats)
        full <- plot_null_dist(df, s$obs, s$direction, "histogram", xlab = s$xlab)
        bare <- plot_null_dist(df, s$obs, s$direction, "histogram", xlab = s$xlab,
                               show_line = FALSE, show_label = FALSE,
                               show_caption = FALSE, show_tail = FALSE)
        ef <- layer_elements(build_checked(full))
        eb <- layer_elements(build_checked(bare))
        if (!isTRUE(all.equal(sort(ef$xmin), sort(eb$xmin))) ||
            !isTRUE(all.equal(sort(ef$height), sort(eb$height))))
            stop("bare-mode bars differ from full-mode bars")
        if (any(eb$fill == RED))
            stop("bare mode still shades the tail")
    })
}

cat("\n=== Bootstrap invariants ===\n")
boot_sc <- make_boot_scenarios()
for (nm in names(boot_sc)) {
    for (mode in c("histogram", "dotplot")) {
        test(paste0(nm, " [", substr(mode, 1, 4), "]"), {
            s <- boot_sc[[nm]]
            df <- data.frame(stat = s$stats)
            p <- plot_boot_dist(df, s$obs, s$conf, s$ci_type, mode,
                                xlab = s$xlab, stat_label = s$stat_label,
                                clamp = s$clamp)
            built <- build_checked(p)
            el <- layer_elements(built)
            assert_axis_economy(built, el)
            span <- diff(range(s$stats))
            eps <- span * 1e-9 + 1e-12

            # Count conservation
            if (abs(sum(el$height) - length(s$stats)) > 1e-6)
                stop("element heights do not sum to the number of replicates")

            # CI lines must not slice a bar when stats are on a lattice
            b <- choose_binning(df$stat, s$obs, align = "center")
            if (b$lattice) {
                vl <- vline_positions(built)
                for (v in vl) {
                    inside <- el$xmin + eps < v & v < el$xmax - eps
                    if (any(inside))
                        stop(paste("CI line at", v, "slices through a bar"))
                }
            }
        })
    }
}

cat("\n=== End-to-end through a real analysis ===\n")
test("TwoPropHTest histogram plot maintains invariants", {
    set.seed(7)
    n <- 60
    cat_data <- data.frame(
        outcome = factor(sample(c("Yes", "No"), n, replace = TRUE, prob = c(0.65, 0.35))),
        group = factor(rep(c("Treatment", "Control"), each = 30))
    )
    r <- TwoPropHTest(data = cat_data, rows = "group", cols = "outcome",
                      hypothesis = "different", reps = 500,
                      dotHist = "histogram", seedBool = TRUE, rngSeed = 99,
                      compare = "rows")
    img <- r$Plot
    st <- img$state
    s <- list(stats = st$df$stat, obs = st$obs_stat, direction = st$direction,
              xlab = "difference in proportions")
    check_null_invariants(s, "histogram")
    check_null_invariants(s, "dotplot")
})

test("SinglePropHTest dotplot plot maintains invariants", {
    set.seed(8)
    cat_data <- data.frame(
        outcome = factor(sample(c("Yes", "No"), 40, replace = TRUE, prob = c(0.7, 0.3)))
    )
    r <- SinglePropHTest(data = cat_data, resp = "outcome",
                         testValue = 0.5, alt = "notequal", reps = 500,
                         dotHist = "dotplot", seedBool = TRUE, rngSeed = 99)
    img <- r$Plot
    st <- img$state
    s <- list(stats = st$df$stat, obs = st$obs_stat, direction = st$direction,
              xlab = "proportion")
    check_null_invariants(s, "histogram")
    check_null_invariants(s, "dotplot")
})

cat("\n\n============================\n")
cat("Results:", passed, "passed,", failed, "failed\n")
cat("============================\n\n")

if (failed > 0) quit(status = 1)
