# Compare two dumps from dump_geometry.R and list every object whose
# drawn geometry or results table differs.
#   Rscript tests/visual/compare_geometry.R before.rds after.rds
args <- commandArgs(TRUE)
a <- readRDS(args[1]); b <- readRDS(args[2])
stopifnot(identical(names(a), names(b)))
same <- mapply(function(x, y) isTRUE(all.equal(x, y, check.attributes = FALSE)), a, b)
cat("objects compared:", length(a), " (tables:", sum(grepl("^tab ", names(a))), ")\n")
cat("objects that differ:", sum(!same), "\n")
for (d in names(a)[!same]) cat("  ", d, "\n")
