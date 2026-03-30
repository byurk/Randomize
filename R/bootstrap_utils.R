#' Bootstrap helper utilities
#'
#' Functions used across analyses for computing bootstrap confidence intervals
#' and for plotting bootstrap distributions.
#'
#' @import ggplot2
#' @import ggforce
#' @import dplyr
#' @import tibble
NULL

compute_boot_ci <- function(boot, obs_stat, conf_level, ci_type,
                            clamp = NULL) {
    if (conf_level > 1)
        conf_level <- conf_level / 100

    if (ci_type == "bootperc") {
        ci <- boot |>
            dplyr::summarise(
                cil = quantile(stat, (1 - conf_level) / 2, na.rm = TRUE),
                ciu = quantile(stat, 1 - (1 - conf_level) / 2, na.rm = TRUE)
            )
        se <- NULL
        zcrit <- NULL
    } else {
        zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
        se <- boot |>
            dplyr::summarise(se = sd(stat, na.rm = TRUE)) |>
            dplyr::pull()
        ci <- tibble::tibble(
            cil = obs_stat - zcrit * se,
            ciu = obs_stat + zcrit * se
        )
    }

    cil <- ci |> dplyr::pull(cil)
    ciu <- ci |> dplyr::pull(ciu)

    if (!is.null(clamp)) {
        cil <- max(cil, clamp[1])
        ciu <- min(ciu, clamp[2])
    }

    list(
        cil = cil,
        ciu = ciu,
        se = se,
        zcrit = zcrit
    )
}

plot_boot_dist <- function(boot, obs_stat, conf_level, ci_type,
                           dotHist = c("dotplot", "histogram"),
                           xlab = "statistic",
                           stat_label = "bootstrap statistics",
                           clamp = NULL) {
    dotHist <- match.arg(dotHist)
    ci <- compute_boot_ci(boot, obs_stat, conf_level, ci_type, clamp)
    cil <- ci$cil
    ciu <- ci$ciu

    if (conf_level > 1) conf_level <- conf_level / 100

    if (ci_type == "bootperc") {
        caption <- paste0(
            "CI limits (dashed) are ",
            round((1 - conf_level) / 2 * 100, 1),
            "% and ",
            round((1 - (1 - conf_level) / 2) * 100, 1),
            "%",
            "\n percentiles of ", stat_label
        )
    } else {
        caption <- paste0(
            "CI (dashed) is calculated using SE of ", stat_label,
            "\n (SE = ", round(ci$se, 3), ", z* = ", round(ci$zcrit, 3), ")"
        )
    }

    if (dotHist == "dotplot") {
        ndist <- dplyr::n_distinct(boot$stat)
        bw <- boot |>
            dplyr::summarise(min = min(stat), max = max(stat)) |>
            dplyr::mutate(bw = (max - min) / 30) |>
            dplyr::pull(bw)

        # Guard against zero bin width (all values identical)
        if (bw == 0) bw <- abs(obs_stat) * 0.01 + 0.001

        if (ndist <= 30) {
            boot <- boot |> dplyr::mutate(x.bin = stat)
            cila <- cil
            ciua <- ciu
        } else {
            boot <- boot |>
                dplyr::mutate(x.bin = obs_stat + ((stat - bw / 2 - obs_stat) %/% bw) * bw)
            cila <- obs_stat + ((cil - bw / 2 - obs_stat) %/% bw) * bw
            ciua <- obs_stat + ((ciu - bw / 2 - obs_stat) %/% bw) * bw
        }

        boot <- boot |>
            dplyr::group_by(x.bin) |>
            dplyr::mutate(y = seq_along(x.bin))

        p <- ggplot2::ggplot(boot) +
            ggforce::geom_ellipse(ggplot2::aes(x0 = x.bin, y0 = y, a = bw / 3, b = 0.5, angle = 0),
                                  show.legend = FALSE) +
            ggplot2::geom_vline(xintercept = c(cila, ciua), linetype = "dashed", color = "red") +
            ggplot2::theme_minimal() +
            ggplot2::ylab("count") +
            ggplot2::xlab(xlab) +
            ggplot2::coord_equal(ratio = bw * 2 / 3) +
            ggplot2::labs(caption = caption) +
            ggplot2::theme(
                text = ggplot2::element_text(size = 14),
                plot.caption = ggplot2::element_text(color = "red", hjust = 0)
            )
    } else {
        p <- ggplot2::ggplot(boot, ggplot2::aes(x = stat)) +
            ggplot2::geom_histogram(center = obs_stat, show.legend = FALSE) +
            ggplot2::geom_vline(xintercept = c(cil, ciu), linetype = "dashed", color = "red") +
            ggplot2::theme_minimal() +
            ggplot2::xlab(xlab) +
            ggplot2::ylab("count") +
            ggplot2::labs(caption = caption) +
            ggplot2::theme(
                text = ggplot2::element_text(size = 14),
                plot.caption = ggplot2::element_text(color = "red", hjust = 0)
            )
    }

    p
}
