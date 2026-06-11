# Render all visual test scenarios to PNG for inspection.
#
# Usage: Rscript tests/visual/generate_plots.R <tag>
# Output: tests/visual/out/<tag>/  (individual PNGs + contact sheets)
#
# Renders each scenario in both histogram and dotplot mode at the
# same logical size Jamovi uses (400 x 350).

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1) args[1] else "current"

suppressMessages(suppressWarnings(devtools::load_all(".", quiet = TRUE)))
suppressMessages(library(grid))

source("tests/visual/scenarios.R")

out_dir <- file.path("tests/visual/out", tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

W <- 400; H <- 350; SCALE <- 2  # render at 2x for legibility

save_plot <- function(p, file) {
    png(file, width = W * SCALE, height = H * SCALE, res = 72 * SCALE)
    ok <- tryCatch({ print(p); TRUE },
                   error = function(e) { grid.newpage(); grid.text(paste("ERROR:", e$message)); FALSE })
    dev.off()
    ok
}

# Build every plot, collect for contact sheets
plots <- list()   # name -> ggplot
errors <- character()

null_sc <- make_null_scenarios()
for (nm in names(null_sc)) {
    s <- null_sc[[nm]]
    df <- data.frame(stat = s$stats)
    for (mode in c("histogram", "dotplot")) {
        key <- paste0("null_", nm, "_", substr(mode, 1, 4))
        p <- tryCatch(
            plot_null_dist(df, s$obs, s$direction, mode,
                           xlab = s$xlab, obs_label = "Observed\nValue") +
                ggplot2::ggtitle(key) +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 9)),
            error = function(e) { errors <<- c(errors, paste(key, "-", e$message)); NULL })
        if (!is.null(p)) {
            plots[[key]] <- p
            save_plot(p, file.path(out_dir, paste0(key, ".png")))
        }
    }
}

boot_sc <- make_boot_scenarios()
for (nm in names(boot_sc)) {
    s <- boot_sc[[nm]]
    df <- data.frame(stat = s$stats)
    for (mode in c("histogram", "dotplot")) {
        key <- paste0(nm, "_", substr(mode, 1, 4))
        p <- tryCatch(
            plot_boot_dist(df, s$obs, s$conf, s$ci_type, mode,
                           xlab = s$xlab, stat_label = s$stat_label,
                           clamp = s$clamp) +
                ggplot2::ggtitle(key) +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 9)),
            error = function(e) { errors <<- c(errors, paste(key, "-", e$message)); NULL })
        if (!is.null(p)) {
            plots[[key]] <- p
            save_plot(p, file.path(out_dir, paste0(key, ".png")))
        }
    }
}

# Contact sheets: 2 x 2 per sheet
keys <- names(plots)
n_per <- 4
n_sheet <- ceiling(length(keys) / n_per)
for (i in seq_len(n_sheet)) {
    sheet_keys <- keys[((i - 1) * n_per + 1):min(i * n_per, length(keys))]
    file <- file.path(out_dir, sprintf("sheet_%02d.png", i))
    png(file, width = 2 * W * SCALE, height = 2 * H * SCALE, res = 72 * SCALE)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2)))
    for (j in seq_along(sheet_keys)) {
        row <- ceiling(j / 2); col <- ((j - 1) %% 2) + 1
        vp <- viewport(layout.pos.row = row, layout.pos.col = col)
        tryCatch(print(plots[[sheet_keys[j]]], vp = vp),
                 error = function(e) {
                     pushViewport(vp)
                     grid.text(paste(sheet_keys[j], "ERROR:", e$message))
                     popViewport()
                 })
    }
    dev.off()
}

cat("Rendered", length(plots), "plots to", out_dir, "\n")
cat("Contact sheets:", n_sheet, "\n")
if (length(errors)) {
    cat("ERRORS:\n")
    for (e in errors) cat(" -", e, "\n")
}
