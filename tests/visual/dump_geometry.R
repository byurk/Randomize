# Dump the drawn geometry and results of both visual batteries for a
# before/after regression comparison that survives cosmetic changes
# (byte-for-byte `cmp` of the PNGs cannot).
#
#   Rscript tests/visual/dump_geometry.R <package dir> <out.rds>
#
# Typical use: `git worktree add /tmp/prev <commit>`, dump from /tmp/prev
# and from ".", then `Rscript tests/visual/compare_geometry.R a.rds b.rds`.
# Dumped per scenario: the first layer of the built plot (bar extents and
# heights, or dot centres, with fills) for every synthetic null / bootstrap
# scenario in both modes, and for every end-to-end analysis scenario its
# results table and plot geometry.
args <- commandArgs(TRUE)
pkg <- args[1]; out <- args[2]
suppressMessages(devtools::load_all(pkg, quiet = TRUE))
suppressMessages(library(ggplot2))
source(file.path(pkg, "tests/visual/scenarios.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
geom <- function(p) {
    d <- ggplot_build(p)$data[[1]]
    d <- d[order(d$x %||% d$xmin, d$y %||% d$ymax),
           intersect(c("x", "y", "xmin", "xmax", "ymin", "ymax", "fill"), names(d))]
    rownames(d) <- NULL
    d
}
res <- list()
for (nm in names(make_null_scenarios())) {
    s <- make_null_scenarios()[[nm]]
    for (m in c("dotplot", "histogram"))
        res[[paste("null", nm, m)]] <- geom(plot_null_dist(data.frame(stat = s$stats), s$obs, s$direction, m, domain = s$domain))
}
for (nm in names(make_boot_scenarios())) {
    s <- make_boot_scenarios()[[nm]]
    for (m in c("dotplot", "histogram"))
        res[[paste("boot", nm, m)]] <- geom(plot_boot_dist(data.frame(stat = s$stats), s$obs, s$conf, s$ci_type, m, clamp = s$clamp))
}
# the analysis battery's scenario definitions: every top-level expression
# of its script except the ones that load a package or render
env <- new.env()
for (e in parse(file.path(pkg, "tests/visual/generate_analysis_plots.R"))) {
    txt <- paste(deparse(e), collapse = "")
    if (!grepl("png\\(|out_dir|dev.off|commandArgs|load_all|library\\(|^cat\\(|^for \\(", txt))
        eval(e, envir = env)
}
Sl <- get("S", envir = env)
for (nm in names(Sl)) {
    r <- Sl[[nm]]("dotplot", FALSE)
    res[[paste("tab", nm)]] <- results_table(r)
    res[[paste("plot", nm)]] <- geom(plot(r))
}
saveRDS(res, out)
cat("dumped", length(res), "objects to", out, "\n")
