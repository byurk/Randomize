# Visual test scenarios for Randomize plotting
#
# Defines a battery of null-distribution and bootstrap-distribution
# scenarios spanning the kinds of statistics the module produces:
# continuous (mean diffs, slopes), lattice-discrete (proportions),
# mixed-lattice (two-prop diffs), irregular-discrete (chi-square),
# and edge cases (p = 0, p = 1, near-constant).
#
# Sourced by generate_plots.R. Each scenario is a list with:
#   $stats     numeric vector of simulated statistics
#   $obs       observed statistic
#   $direction "greater" / "less" / "two_sided"   (null scenarios)
#   $conf, $ci_type, $clamp                        (boot scenarios)

make_null_scenarios <- function() {
    sc <- list()

    set.seed(101)
    sc$cont_two_sided <- list(
        stats = rnorm(1000, 0, 1), obs = 1.8, direction = "two_sided",
        xlab = "difference in means")

    set.seed(102)
    sc$cont_greater <- list(
        stats = rnorm(1000, 0, 1), obs = 1.2, direction = "greater",
        xlab = "difference in means")

    set.seed(103)
    sc$cont_less <- list(
        stats = rnorm(1000, 0, 1), obs = -1.2, direction = "less",
        xlab = "difference in means")

    set.seed(104)
    sc$cont_two_sided_left <- list(
        # two-sided where the LEFT tail is the smaller one (lt branch)
        stats = rnorm(1000, 0, 1), obs = -1.5, direction = "two_sided",
        xlab = "difference in means")

    set.seed(105)
    sc$prop_n20_greater <- list(
        # single proportion, n = 20: lattice spacing 0.05
        stats = rbinom(1000, 20, 0.5) / 20, obs = 0.65, direction = "greater",
        xlab = "proportion", domain = c(0, 1))

    set.seed(106)
    sc$prop_n20_two_sided <- list(
        stats = rbinom(1000, 20, 0.5) / 20, obs = 0.65, direction = "two_sided",
        xlab = "proportion", domain = c(0, 1))

    set.seed(107)
    sc$prop_n50_two_sided <- list(
        # lattice spacing 0.02
        stats = rbinom(1000, 50, 0.5) / 50, obs = 0.6, direction = "two_sided",
        xlab = "proportion", domain = c(0, 1))

    set.seed(108)
    sc$prop_n100_greater <- list(
        # lattice spacing 0.01 -- finer than default bin width
        stats = rbinom(1000, 100, 0.5) / 100, obs = 0.57, direction = "greater",
        xlab = "proportion", domain = c(0, 1))

    set.seed(109)
    sc$prop_n10_greater <- list(
        # very coarse: at most 11 distinct values
        stats = rbinom(1000, 10, 0.5) / 10, obs = 0.7, direction = "greater",
        xlab = "proportion", domain = c(0, 1))

    set.seed(126)
    sc$prop_n30_greater <- list(
        # n = 30, 1000 reps: 18 distinct values -- must get one column per
        # value (a tail gap is honest, not a reason to double the width)
        stats = rbinom(1000, 30, 0.5) / 30, obs = 20 / 30, direction = "greater",
        xlab = "proportion", domain = c(0, 1))

    set.seed(127)
    sc$prop_n25_two_sided <- list(
        stats = rbinom(500, 25, 0.5) / 25, obs = 18 / 25, direction = "two_sided",
        xlab = "proportion", domain = c(0, 1))

    set.seed(128)
    sc$prop_n200_greater <- list(
        # ~47 lattice steps: still one column per value
        stats = rbinom(1000, 200, 0.5) / 200, obs = 115 / 200, direction = "greater",
        xlab = "proportion", domain = c(0, 1))

    set.seed(129)
    sc$prop_n500_greater <- list(
        # too many lattice steps for one column each: grouped bins
        stats = rbinom(1000, 500, 0.5) / 500, obs = 270 / 500, direction = "greater",
        xlab = "proportion", domain = c(0, 1))

    set.seed(110)
    sc$diffprop_25_25 <- list(
        # two-prop diff, equal n: lattice spacing 0.04
        stats = rbinom(1000, 25, 0.5) / 25 - rbinom(1000, 25, 0.5) / 25,
        obs = 0.16, direction = "two_sided",
        xlab = "difference in proportions", domain = c(-1, 1))

    set.seed(111)
    sc$diffprop_30_20 <- list(
        # unequal n: values on 1/60 lattice but unevenly occupied
        stats = rbinom(1000, 30, 0.4) / 30 - rbinom(1000, 20, 0.4) / 20,
        obs = 0.21, direction = "greater",
        xlab = "difference in proportions", domain = c(-1, 1))

    set.seed(112)
    sc$chisq_2x2 <- list(
        # chi-square from permuted 2x2 tables: irregular discrete values
        stats = local({
            rows <- c(30, 30); cols <- c(25, 35)
            tabs <- r2dtable(1000, rows, cols)
            vapply(tabs, function(tt) {
                e <- outer(rows, cols) / 60
                sum((tt - e)^2 / e)
            }, numeric(1))
        }),
        obs = 5.4, direction = "greater",
        xlab = "X-squared", domain = c(0, Inf))

    set.seed(130)
    sc$chisq_2x2_uneven <- list(
        # 2x2 table whose two smallest chi-square values nearly coincide:
        # they must share one bar rather than sizing every bar by that gap
        stats = local({
            rows <- c(25, 25); cols <- c(30, 20)
            tabs <- r2dtable(1000, rows, cols)
            vapply(tabs, function(tt) {
                e <- outer(rows, cols) / 50
                sum((tt - e)^2 / e)
            }, numeric(1))
        }),
        obs = 2.1, direction = "greater",
        xlab = "X-squared", domain = c(0, Inf))

    set.seed(113)
    sc$chisq_3x3 <- list(
        stats = local({
            rows <- c(20, 20, 20); cols <- c(18, 22, 20)
            tabs <- r2dtable(1000, rows, cols)
            vapply(tabs, function(tt) {
                e <- outer(rows, cols) / 60
                sum((tt - e)^2 / e)
            }, numeric(1))
        }),
        obs = 9.2, direction = "greater",
        xlab = "X-squared", domain = c(0, Inf))

    set.seed(114)
    sc$slope_two_sided <- list(
        stats = rnorm(1000, 0, 0.35), obs = 0.8, direction = "two_sided",
        xlab = "slope")

    set.seed(115)
    sc$pval_zero <- list(
        # observed beyond every simulated value
        stats = rnorm(1000, 0, 1), obs = 4.5, direction = "greater",
        xlab = "difference in means")

    set.seed(116)
    sc$pval_one <- list(
        # observed below every simulated value, direction greater: p = 1
        stats = rnorm(1000, 0, 1), obs = -4.5, direction = "greater",
        xlab = "difference in means")

    set.seed(117)
    sc$few_distinct <- list(
        # only a handful of distinct values
        stats = sample(c(-0.2, -0.1, 0, 0.1, 0.2), 1000, replace = TRUE,
                       prob = c(0.1, 0.2, 0.4, 0.2, 0.1)),
        obs = 0.1, direction = "greater",
        xlab = "difference")

    set.seed(118)
    sc$skewed_right <- list(
        # skewed null (like an F or X2 statistic, but continuous)
        stats = rchisq(1000, df = 3), obs = 7.5, direction = "greater",
        xlab = "F statistic")

    set.seed(119)
    sc$reps_100 <- list(
        stats = rbinom(100, 20, 0.5) / 20, obs = 0.65, direction = "greater",
        xlab = "proportion")

    set.seed(120)
    sc$reps_5000 <- list(
        stats = rbinom(5000, 20, 0.5) / 20, obs = 0.65, direction = "greater",
        xlab = "proportion")

    set.seed(121)
    sc$reps_10 <- list(
        # classroom demo scale
        stats = rbinom(10, 20, 0.5) / 20, obs = 0.65, direction = "greater",
        xlab = "proportion")

    set.seed(122)
    sc$reps_25_cont <- list(
        stats = rnorm(25), obs = 1.1, direction = "two_sided",
        xlab = "difference")

    set.seed(123)
    sc$far_obs <- list(
        # observed value far beyond the null: axis must not collapse the
        # distribution into a sliver (arrow at panel edge instead)
        stats = rnorm(1000), obs = 15, direction = "greater",
        xlab = "slope")

    set.seed(124)
    sc$far_obs_left <- list(
        stats = rnorm(1000), obs = -15, direction = "less",
        xlab = "difference")

    set.seed(125)
    sc$outlier <- list(
        # a lone extreme simulation must not force chunky bins or crush
        # the axis
        stats = c(rnorm(999), 8), obs = 1.8, direction = "two_sided",
        xlab = "difference")

    sc
}

make_boot_scenarios <- function() {
    sc <- list()

    set.seed(201)
    sc$boot_mean_perc <- list(
        stats = rnorm(1000, 10.5, 0.4), obs = 10.5,
        conf = 95, ci_type = "bootperc", clamp = NULL,
        xlab = "mean", stat_label = "bootstrap means")

    set.seed(202)
    sc$boot_mean_se <- list(
        stats = rnorm(1000, 10.5, 0.4), obs = 10.5,
        conf = 95, ci_type = "bootse", clamp = NULL,
        xlab = "mean", stat_label = "bootstrap means")

    set.seed(203)
    sc$boot_prop_n20_perc <- list(
        # lattice spacing 0.05
        stats = rbinom(1000, 20, 0.6) / 20, obs = 0.6,
        conf = 95, ci_type = "bootperc", clamp = c(0, 1),
        xlab = "proportion", stat_label = "bootstrap proportions")

    set.seed(204)
    sc$boot_prop_n50_perc <- list(
        stats = rbinom(1000, 50, 0.6) / 50, obs = 0.6,
        conf = 95, ci_type = "bootperc", clamp = c(0, 1),
        xlab = "proportion", stat_label = "bootstrap proportions")

    set.seed(205)
    sc$boot_prop_edge <- list(
        # observed proportion near 1; CI clamped at 1
        stats = rbinom(1000, 30, 0.95) / 30, obs = 29 / 30,
        conf = 95, ci_type = "bootperc", clamp = c(0, 1),
        xlab = "proportion", stat_label = "bootstrap proportions")

    set.seed(206)
    sc$boot_diffprop <- list(
        stats = rbinom(1000, 25, 0.7) / 25 - rbinom(1000, 25, 0.5) / 25,
        obs = 0.2,
        conf = 95, ci_type = "bootperc", clamp = c(-1, 1),
        xlab = "difference in proportions",
        stat_label = "bootstrap differences")

    set.seed(207)
    sc$boot_slope <- list(
        stats = rnorm(1000, 2, 0.15), obs = 2,
        conf = 90, ci_type = "bootperc", clamp = NULL,
        xlab = "slope", stat_label = "bootstrap slopes")

    set.seed(208)
    sc$boot_reps_20 <- list(
        # classroom demo scale
        stats = rbinom(20, 25, 0.6) / 25, obs = 0.6,
        conf = 95, ci_type = "bootperc", clamp = c(0, 1),
        xlab = "proportion", stat_label = "bootstrap proportions")

    set.seed(209)
    sc$boot_outlier <- list(
        # a lone extreme replicate must not force chunky bins
        stats = c(rnorm(999, 2, 0.15), 4), obs = 2,
        conf = 95, ci_type = "bootperc", clamp = NULL,
        xlab = "slope", stat_label = "bootstrap slopes")

    sc
}
