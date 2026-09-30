# Render the model-based calculator's density plots across distributions,
# degrees of freedom, tails, and observed values (including edge cases).
#
# Usage: Rscript tests/visual/generate_model_plots.R <tag>
# Output: tests/visual/out/<tag>/model_sheet_NN.png (3 x 2 per sheet)

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "model"
suppressMessages(suppressWarnings(devtools::load_all(".", quiet = TRUE)))
suppressMessages({library(grid); library(ggplot2)})
out_dir <- file.path("tests/visual/out", tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
W <- 450; H <- 250; SCALE <- 2

cases <- list(
    list(distro = "ndistro", obs = 1.96, tail = "right"),
    list(distro = "ndistro", obs = -1.96, tail = "right"),
    list(distro = "ndistro", obs = -0.5, tail = "left"),
    list(distro = "ndistro", obs = 0, tail = "both"),
    list(distro = "ndistro", obs = -2.5, tail = "both"),
    list(distro = "ndistro", obs = 6, tail = "right"),
    list(distro = "ndistro", obs = -6, tail = "left"),
    list(distro = "tdistro", dF = 1, obs = 2, tail = "right"),
    list(distro = "tdistro", dF = 3, obs = -2.2, tail = "both"),
    list(distro = "tdistro", dF = 30, obs = 1.5, tail = "left"),
    list(distro = "tdistro", dF = 200, obs = 2.6, tail = "both"),
    list(distro = "chisq", dF = 1, obs = 3.84),
    list(distro = "chisq", dF = 2, obs = 0.5),
    list(distro = "chisq", dF = 4, obs = 9.49),
    list(distro = "chisq", dF = 10, obs = 25),
    list(distro = "chisq", dF = 50, obs = 40),
    list(distro = "chisq", dF = 3, obs = 0),
    list(distro = "chisq", dF = 3, obs = -1),
    list(distro = "chisq", dF = 3, obs = 40),
    list(distro = "fdistro", dF = 1, dF2 = 10, obs = 4.96),
    list(distro = "fdistro", dF = 2, dF2 = 27, obs = 3.35),
    list(distro = "fdistro", dF = 4, dF2 = 40, obs = 1.2),
    list(distro = "fdistro", dF = 10, dF2 = 200, obs = 2.0),
    list(distro = "fdistro", dF = 3, dF2 = 3, obs = 9.28),
    list(distro = "fdistro", dF = 5, dF2 = 20, obs = 15))
plots <- list(); errors <- character()
for (cs in cases) {
    dF <- if (is.null(cs$dF)) 30 else cs$dF; dF2 <- if (is.null(cs$dF2)) 30 else cs$dF2
    tail <- if (is.null(cs$tail)) "right" else cs$tail
    key <- sprintf("%s_df%d_%d_obs%g_%s", cs$distro, dF, dF2, cs$obs, tail)
    p <- tryCatch({
        r <- modelBased(distro = cs$distro, dF = dF, dF2 = dF2, areaBool = TRUE, obsStat = cs$obs, tail = tail)
        a <- r$areaTable$asDF$area
        plot(r) + ggtitle(sprintf("%s   area = %.4g", key, a)) + theme(plot.title = element_text(size = 9))
    }, error = function(e) { errors <<- c(errors, paste(key, "-", conditionMessage(e))); NULL })
    if (!is.null(p)) plots[[key]] <- p
}
keys <- names(plots); n_per <- 6
n_sheet <- ceiling(length(keys) / n_per)
for (i in seq_len(n_sheet)) {
    sk <- keys[((i - 1) * n_per + 1):min(i * n_per, length(keys))]
    png(file.path(out_dir, sprintf("model_sheet_%02d.png", i)), width = 2 * W * SCALE, height = 3 * H * SCALE, res = 72 * SCALE)
    grid.newpage(); pushViewport(viewport(layout = grid.layout(3, 2)))
    for (j in seq_along(sk)) {
        vp <- viewport(layout.pos.row = ceiling(j / 2), layout.pos.col = ((j - 1) %% 2) + 1)
        tryCatch(print(plots[[sk[j]]], vp = vp), error = function(e) {
            pushViewport(vp); grid.text(paste(sk[j], "ERROR:", conditionMessage(e))); popViewport() })
    }
    dev.off()
}
cat("Rendered", length(plots), "model plots on", n_sheet, "sheets\n")
if (length(errors)) { cat("ERRORS:\n"); cat(paste(" -", errors), sep = "\n") }
