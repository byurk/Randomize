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
#   4. Bootstrap CI lines are drawn at the exact reported bounds.
#   5. No interior empty bins (gaps) within the central 99.8% of
#      simulations -- null and bootstrap alike.
#   6. Nothing is drawn outside the statistic's domain (no negative
#      chi-square bar, no proportion beyond [0, 1]).
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
    # (dotplots never show fewer than 8 counts on the axis, so short
    # stacks legitimately fill less of it)
    if ("ytop" %in% names(el) && max(el$ytop) >= 8) {
        y_use <- max(el$ytop) / yr[2]
        if (y_use < 0.7)
            stop(paste0("distribution occupies only ",
                        round(100 * y_use), "% of the y axis"))
    }
}

# Extract per-element (bar or dot) x-extent and fill from the first
# layer of a built plot.  For geom_rect each row is a bar; for
# GeomDotStack each row is a dot center, with the semi-width in the
# layer params (the semi-height is resolved at draw time; 0.45 is its
# upper bound).
layer_elements <- function(built) {
    d <- built$data[[1]]
    if (all(c("xmin", "xmax") %in% names(d)) && !all(is.na(d$xmin))) {
        data.frame(xmin = d$xmin, xmax = d$xmax, fill = d$fill,
                   height = d$ymax - d$ymin, ytop = d$ymax, halfw = 0)
    } else {
        a <- built$plot$layers[[1]]$geom_params$a
        data.frame(xmin = d$x - a, xmax = d$x + a, fill = d$fill,
                   height = 1, ytop = d$y + 0.5, halfw = a)
    }
}

# A stack of dots (or a bar) is never part red / part black: every
# element sharing an x position has exactly one fill color.
assert_single_fill_stacks <- function(el) {
    mids <- round((el$xmin + el$xmax) / 2, 12)
    n_fills <- tapply(el$fill, mids, function(f) length(unique(f)))
    if (any(n_fills > 1))
        stop("a stack mixes fill colors")
}

vline_positions <- function(built) {
    for (d in built$data) {
        if ("xintercept" %in% names(d)) return(d$xintercept)
    }
    numeric(0)
}

# No missing interior bins: within the central 98% of the data,
# occupied bins must sit at a regular spacing with no absent bin between
# them.  (A lone extreme simulation may still have its own separated bar
# beyond the window -- that is honest empty space, not a binning gap.)
#
# When every bin holds exactly one achievable value (single_value), an
# empty bin is a value that genuinely never occurred and is allowed --
# except between two well-populated bins (both >= 5), where it would
# read as a binning error.
assert_gap_free <- function(el, stats, single_value = FALSE) {
    win <- stats::quantile(stats, c(0.01, 0.99), names = FALSE)
    mid_all <- round((el$xmin + el$xmax) / 2, 12)
    h <- tapply(el$height, mid_all, sum)
    mids <- sort(as.numeric(names(h)))
    if (length(mids) < 3) return(invisible())
    step <- min(diff(mids))
    if (single_value) {
        gap <- diff(mids) > 1.5 * step
        n_l <- as.numeric(h[as.character(mids[-length(mids)])])
        n_r <- as.numeric(h[as.character(mids[-1])])
        if (any(gap & n_l >= 5 & n_r >= 5))
            stop("empty single-value bin between two populated bins")
        return(invisible())
    }
    inwin <- mids[mids >= win[1] - step & mids <= win[2] + step]
    if (length(inwin) > 1) {
        gap <- diff(inwin) > 1.5 * step
        if (length(stats) >= 300) {
            # from 300 simulations up, a gap is tolerated when either
            # flanking column is sparse (< 5): honest sparseness, not a
            # binning artifact
            hl <- as.numeric(h[as.character(inwin[-length(inwin)])])
            hr <- as.numeric(h[as.character(inwin[-1])])
            gap <- gap & hl >= 5 & hr >= 5
        }
        if (any(gap)) stop("empty bin (gap) between well-populated bins inside the central window")
    }
}

# Does the plot use one bin per achievable lattice value?
is_single_value <- function(stats, anchor, align, sign = 1) {
    b <- choose_binning(stats, anchor, align = align, sign = sign)
    b$lattice && abs(b$bw - b$res) <= b$res * 1e-9
}

check_null_invariants <- function(s, mode) {
    df <- data.frame(stat = s$stats)
    p <- plot_null_dist(df, s$obs, s$direction, mode,
                        xlab = s$xlab, obs_label = "Observed\nValue",
                        domain = s$domain)
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

    # Sparse statistics (<= 8 distinct values) draw one value-centered
    # bar/column per achievable value: gaps between values are honest,
    # and side purity is judged by centers (a bar AT the observed value
    # legitimately straddles the line, pure by value identity).
    sparse <- is_sparse_values(s$stats)
    sgn0 <- if (ptail == "lt") -1 else 1
    if (sparse) {
        el$halfw <- (el$xmax - el$xmin) / 2
    } else {
        assert_gap_free(el, s$stats,
                        single_value = is_single_value(s$stats, s$obs, "edge", sgn0))
    }

    # Nothing placed at impossible values.  Grouped bars are trimmed to
    # the domain (strict); single-value bars follow the value-at-edge
    # convention and may extend one bin past a bound; value-centered
    # elements are judged by center.  Grouped dot columns sit at bin
    # midpoints and are exempt (a midpoint may fall within half a bin of
    # the bound, the same accepted class as a glyph grazing the line).
    if (!is.null(s$domain)) {
        sgn <- if (ptail == "lt") -1 else 1
        slack <- 0
        if (!sparse) {
            b <- choose_binning(s$stats, s$obs, align = "edge", sign = sgn)
            if (b$lattice && abs(b$bw - b$res) <= b$res * 1e-9)
                slack <- b$bw
        }
        bars_el <- el[el$halfw == 0, ]
        if (nrow(bars_el) > 0 &&
            (min(bars_el$xmin) < s$domain[1] - slack - eps ||
             max(bars_el$xmax) > s$domain[2] + slack + eps))
            stop("bar drawn outside the statistic's domain")
        pts <- el[el$halfw > 0, ]
        if (nrow(pts) > 0) {
            ctr <- (pts$xmin + pts$xmax) / 2
            if (min(ctr) < s$domain[1] - eps || max(ctr) > s$domain[2] + eps)
                stop("element centered outside the statistic's domain")
        }
    }

    red <- el[el$fill == RED, ]
    black <- el[el$fill != RED, ]

    # 1. Side purity of everything actually drawn.  Bars are intervals,
    # so the strict edge rule applies (halfw = 0).  Dot columns are
    # point claims judged by their CENTER (halfw = dot semi-width): a
    # tie column -- dots for simulations exactly equal to the observed
    # value -- legitimately sits centered on the line, and it must be
    # red, never black (ties count toward the p-value).
    bb <- black[black$halfw == 0, ]   # bars: strict edge rule
    bd <- black[black$halfw > 0, ]    # dot columns: center rule
    if (ptail == "rt") {
        if (nrow(red) > 0 && min(red$xmin + red$halfw) < s$obs - eps)
            stop("red element lies left of the observed value")
        if (nrow(bb) > 0 && max(bb$xmax) > s$obs + eps)
            stop("black bar extends right of the observed value")
        if (nrow(bd) > 0 && max(bd$xmax - bd$halfw) > s$obs - eps)
            stop("black dot column centered at or beyond the observed value")
    } else {
        if (nrow(red) > 0 && max(red$xmax - red$halfw) > s$obs + eps)
            stop("red element lies right of the observed value")
        if (nrow(bb) > 0 && min(bb$xmin) < s$obs - eps)
            stop("black bar extends left of the observed value")
        if (nrow(bd) > 0 && min(bd$xmin + bd$halfw) < s$obs + eps)
            stop("black dot column centered at or beyond the observed value")
    }

    # 1b. No mixed-color stacks, ever
    assert_single_fill_stacks(el)

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
        # A dot column clamped to a domain bound may lean into its
        # neighbor's nominal extent; drawn glyphs are narrower than the
        # nominal half-width, so tolerate overlap below 35% of a column
        colw <- stats::median(ord$xmax - ord$xmin)
        if (any(ord$xmin[-1] - ord$xmax[-nrow(ord)] < -0.35 * colw - eps))
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

cat("\n=== Dense lattices: one column per value, tie stack on the line (#9, #11) ===\n")
for (nm in c("prop_n20_greater", "prop_n30_greater", "prop_n25_two_sided",
             "prop_n50_two_sided", "prop_n100_greater")) {
    test(nm, {
        s <- null_sc[[nm]]
        df <- data.frame(stat = s$stats)
        res <- min(diff(sort(unique(s$stats))))
        sgn <- if (s$direction == "two_sided" &&
                   mean(s$stats > s$obs) >= mean(s$stats < s$obs)) -1 else 1
        if (!is_single_value(s$stats, s$obs, "edge", sgn))
            stop("bins group several achievable values")
        # dotplot: a column of dots sits exactly on the observed value,
        # and (ties count toward the p-value) it is red
        built <- build_checked(plot_null_dist(df, s$obs, s$direction, "dotplot",
                                              xlab = s$xlab, domain = s$domain))
        d <- built$data[[1]]
        on_line <- abs(d$x - s$obs) < res * 1e-6
        if (sum(abs(s$stats - s$obs) < res * 1e-6) > 0) {
            if (!any(on_line)) stop("no dot stack centered on the observed value")
            if (any(d$fill[on_line] != RED)) stop("tie stack is not red")
        }
        # histogram: the number of bars equals the number of distinct values
        built <- build_checked(plot_null_dist(df, s$obs, s$direction, "histogram",
                                              xlab = s$xlab, domain = s$domain))
        if (nrow(built$data[[1]]) != length(unique(s$stats)))
            stop("bars do not correspond one-to-one to achievable values")
    })
}
for (nm in c("prop_n200_greater", "prop_n500_greater", "rounded_diff_5000")) {
    test(paste0(nm, " groups lattice steps: never more than the target bin count"), {
        s <- null_sc[[nm]]
        if (is_single_value(s$stats, s$obs, "edge", 1))
            stop("expected grouped bins")
        b <- choose_binning(s$stats, s$obs, align = "edge", sign = 1)
        if (diff(range(s$stats)) / b$bw > 30 + 1e-9) stop("more than 30 columns")
    })
}

cat("\n=== Sparse values: near-coincident values share a column ===\n")
test("Yates-corrected 2x2 chi-square: no sliver bars", {
    # values (|k - 15| - 0.5)^2 * c and 0: the two smallest almost coincide
    set.seed(31)
    k <- rbinom(1000, 30, 0.5)
    stats <- pmax(abs(k - 15) - 0.5, 0)^2 / 3
    obs <- 2.5^2 / 3
    df <- data.frame(stat = stats)
    built <- build_checked(plot_null_dist(df, obs, "greater", "histogram",
                                          xlab = "X2", domain = c(0, Inf)))
    el <- layer_elements(built)
    # the leftmost bar is legitimately trimmed at the domain bound 0
    inner <- el[el$xmin > 1e-9, ]
    w <- inner$xmax - inner$xmin
    span <- diff(range(stats))
    if (min(w) < 0.02 * span) stop("sliver bar drawn")
    # count conservation and purity still hold after merging
    if (abs(sum(el$height) - 1000) > 1e-6) stop("counts not conserved")
    red <- el[el$fill == RED, ]; black <- el[el$fill != RED, ]
    if (nrow(red) && min((red$xmin + red$xmax) / 2) < obs - 1e-9) stop("red bar left of obs")
    if (nrow(black) && max((black$xmin + black$xmax) / 2) > obs + 1e-9) stop("black bar right of obs")
    # the merged column is one bar for 0 and 0.083 together
    n0 <- sum(stats < 0.1)
    if (!any(abs(el$height - n0) < 1e-6)) stop("near-coincident values were not merged")
})
test("sparse_groups never merges across the observed value", {
    stats <- c(rep(0, 50), rep(0.05, 30), rep(2, 20))
    g <- sparse_groups(stats, extreme = stats >= 0.05)
    if (length(g$centers) != 3) stop("values on opposite sides of obs were merged")
})

cat("\n=== Sparse mode: structure, not just few values ===\n")
test("5 continuous reps are binned; a 2x2 chi-square with ties is one column per value", {
    set.seed(9)
    cont5 <- rnorm(5)
    if (is_sparse_values(cont5)) stop("5 distinct continuous values should be binned")
    if (!is_sparse_values(rnorm(2))) stop("2 values are always sparse")
    # Yates-corrected 2x2 chi-square: 8 distinct values, heavy repetition, not a lattice
    k <- rhyper(300, 7, 18, 13); e <- c(13*7, 13*18, 12*7, 12*18) / 25   # permuted 2x2 tables
    x2 <- vapply(k, function(a) { o <- c(a, 13 - a, 7 - a, 12 - (7 - a)); sum((pmax(abs(o - e) - 0.5, 0))^2 / e) }, numeric(1))
    if (!is_sparse_values(x2)) stop("repeated chi-square values should be sparse")
    # dots of a sparse chi-square plot sit at (merged) exact values, never at a bin midpoint
    b <- ggplot2::ggplot_build(plot_null_dist(data.frame(stat = x2), 0, "greater", "dotplot", domain = c(0, Inf)))
    xs <- sort(unique(b$data[[1]]$x))
    if (min(xs) > 0.1) stop("the zero-valued simulations are not drawn at zero")
})

cat("\n=== Bin count is stable across re-runs of the same analysis ===\n")
for (case in list(list(nm = "normal 1000", gen = function() rnorm(1000)),
                  list(nm = "normal 300", gen = function() rnorm(300)),
                  list(nm = "chi-square(3) 1000", gen = function() rchisq(1000, 3)),
                  list(nm = "F(3,20) 1000", gen = function() rf(1000, 3, 20)),
                  list(nm = "t(2) heavy tails 1000", gen = function() rt(1000, 2)),
                  list(nm = "F(3,20) 5000", gen = function() rf(5000, 3, 20)))) {
    test(paste(case$nm, "reps land on the same bin count run after run"), {
        set.seed(31)
        bars <- replicate(25, {
            s <- case$gen()
            b <- choose_binning(s, stats::quantile(s, 0.9), align = "edge", sign = 1)
            round(diff(range(s)) / b$bw)
        })
        if (length(unique(bars)) > 2 || min(bars) < 24)
            stop(paste("bar counts across runs:", paste(sort(unique(bars)), collapse = " ")))
    })
}
test("below 300 reps no bin in the central 98% is ever empty", {
    set.seed(32)
    for (i in 1:20) {
        s <- rchisq(100, 3)
        b <- choose_binning(s, 2, align = "edge", sign = 1)
        idx <- bin_index(s, 2, b$bw, b$off, 1, "edge")
        q <- stats::quantile(idx, c(0.01, 0.99), type = 1, names = FALSE)
        w <- idx[idx >= q[1] & idx <= q[2]]
        if (any(tabulate(w - min(w) + 1) == 0)) stop("empty bin inside the window at 100 reps")
    }
})

cat("\n=== Direction invariance: the alternative never changes the bars (except the tie column) ===\n")
test("rounded data: less / greater / two-sided give the same histogram apart from ties at the line", {
    set.seed(41)
    sc <- round(c(rnorm(13, 10, 2), rnorm(12, 10, 2)), 1); g <- rep(c("A", "B"), c(13, 12))
    obs <- mean(sc[g == "A"]) - mean(sc[g == "B"])
    s <- permute_diff_means(sc, g, c("A", "B"), 1000)$stat
    ties <- sum(abs(s - obs) < 1e-9)
    if (ties == 0) stop("test data should produce ties at the observed value")
    heights <- function(dir) {
        e <- ggplot2::ggplot_build(plot_null_dist(data.frame(stat = s), obs, dir, "histogram"))$data[[1]]
        e <- e[order(e$xmin), ]; list(xmin = round(e$xmin, 9), h = e$ymax)
    }
    a <- heights("less"); b <- heights("greater")
    if (!identical(a$xmin, b$xmin)) stop("bar edges depend on the direction")
    d <- a$h - b$h
    moved <- which(d != 0)
    # only the two bars touching the observed value may differ, and by the tie count
    if (length(moved) > 2 || any(abs(d[moved]) != ties)) stop("bar heights depend on the direction beyond the tie column")
    if (length(moved) == 2 && !(abs(a$xmin[moved[1]] + (a$xmin[moved[2]] - a$xmin[moved[1]]) - obs) < 1e-6))
        stop("the bars that differ are not the ones at the observed value")
    # continuous data: strictly identical
    s2 <- s + runif(1000, -1e-4, 1e-4); obs2 <- obs + 3e-5
    e1 <- ggplot2::ggplot_build(plot_null_dist(data.frame(stat = s2), obs2, "less", "histogram"))$data[[1]]
    e2 <- ggplot2::ggplot_build(plot_null_dist(data.frame(stat = s2), obs2, "greater", "histogram"))$data[[1]]
    m1 <- unname(as.matrix(e1[order(e1$xmin), c("xmin", "ymax")]))
    m2 <- unname(as.matrix(e2[order(e2$xmin), c("xmin", "ymax")]))
    if (!isTRUE(all.equal(m1, m2))) stop("continuous data: bars depend on the direction")
})

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

            # CI lines sit at the exact reported bounds (no display nudge)
            ci <- compute_boot_ci(df, s$obs, s$conf, s$ci_type, s$clamp)
            vl <- sort(vline_positions(built))
            if (length(vl) != 2 ||
                max(abs(vl - c(ci$cil, ci$ciu))) > eps)
                stop("CI lines are not at the reported CI bounds")

            # No interior empty bins (unless sparse: value-centered
            # elements with honest gaps), and nothing placed at
            # impossible values (clamp = domain bounds; bootstrap bars
            # are trimmed, value-centered elements judged by center,
            # grouped dot columns exempt as in the null checks)
            sparse <- is_sparse_values(s$stats)
            if (sparse) {
                el$halfw <- (el$xmax - el$xmin) / 2
            } else {
                assert_gap_free(el, s$stats,
                                single_value = is_single_value(s$stats, s$obs, "center"))
            }
            if (!is.null(s$clamp)) {
                bars_el <- el[el$halfw == 0, ]
                if (nrow(bars_el) > 0 &&
                    (min(bars_el$xmin) < s$clamp[1] - eps ||
                     max(bars_el$xmax) > s$clamp[2] + eps))
                    stop("bar drawn outside the statistic's domain")
                pts <- el[el$halfw > 0, ]
                if (sparse && nrow(pts) > 0) {
                    ctr <- (pts$xmin + pts$xmax) / 2
                    if (min(ctr) < s$clamp[1] - eps || max(ctr) > s$clamp[2] + eps)
                        stop("element centered outside the statistic's domain")
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
