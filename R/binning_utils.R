#' Discreteness-aware binning utilities
#'
#' Shared binning logic for the null-distribution and bootstrap-distribution
#' plots.  Simulated statistics are often discrete (proportions are multiples
#' of 1/n, chi-square statistics take few distinct values), and naive
#' fixed-count binning produces histograms with irregular gaps between bars.
#' These helpers pick a bin width adapted to the discreteness of the data
#' while keeping bins aligned to an anchor point (the observed statistic),
#' so that no bar ever mixes values from both sides of the anchor.
#'
#' @name binning_utils
NULL

#' Choose a bin width adapted to the data's discreteness
#'
#' Detects whether the simulated statistics lie on a regular lattice
#' (e.g. multiples of 1/n for proportions).  If so, the bin width is
#' constrained to an integer multiple of the lattice spacing so that every
#' bin covers the same number of possible values -- this removes the
#' irregular gaps and uneven bar widths that aliasing produces.  For
#' sparse or irregular discrete data (e.g. chi-square statistics), the
#' number of bins is reduced until interior empty bins are rare.
#'
#' @param stats Numeric vector of simulated statistics.
#' @param anchor Alignment point for the bin grid (the observed statistic).
#' @param align \code{"edge"} places a bin edge at the anchor (used for
#'   null distributions, so bars never straddle the observed value);
#'   \code{"center"} centers a bin on the anchor (used for bootstrap
#'   distributions).
#' @param sign +1 bins rightward from the anchor, -1 mirrors (used when
#'   the extreme tail is the lower one).
#' @param target_bins Preferred number of bins for continuous data.
#'
#' @return A list with \code{bw} (bin width), \code{off} (half-resolution
#'   offset that keeps lattice values away from bin edges), \code{lattice}
#'   (logical), and \code{res} (lattice resolution).
#'
#' @keywords internal
choose_binning <- function(stats, anchor, align = c("edge", "center"),
                           sign = 1, target_bins = 30) {
    align <- match.arg(align)
    u <- sort(unique(stats))
    span <- u[length(u)] - u[1]

    # Degenerate: all simulated values identical (up to floating-point noise)
    if (span <= max(abs(u[1]), abs(u[length(u)])) * 1e-10) {
        bw <- max(abs(anchor - u[1]) / 5, abs(u[1]) * 0.01, 0.001)
        return(list(bw = bw, off = bw / 2, lattice = TRUE, res = bw))
    }

    # The same lattice point computed from different inputs (e.g.
    # a/30 - b/20) can differ by one ulp, creating phantom near-duplicate
    # "unique" values.  Gaps that are floating-point noise must be ignored
    # or the detected resolution collapses and binning aliases badly.
    d <- diff(u)
    d <- d[d > span * 1e-8]
    res <- min(d)
    mult <- d / res
    lattice <- all(abs(mult - round(mult)) < 0.01)

    bin_counts <- function(bw) {
        idx <- bin_index(stats, anchor, bw, min(res, bw) / 2, sign, align)
        tabulate(idx - min(idx) + 1L)
    }

    frac_empty <- function(counts) {
        1 - sum(counts > 0) / length(counts)
    }

    # Fraction of interior bins that drop to half (or less) of both
    # substantial neighbors.  Detects sawtooth patterns, e.g. mixed
    # lattices like p1 - p2 with unequal n, where adjacent multiples of
    # the combined lattice spacing have systematically different
    # probabilities.
    frac_zigzag <- function(counts) {
        if (length(counts) < 5) return(0)
        i <- 2:(length(counts) - 1)
        dip <- counts[i - 1] >= 3 & counts[i + 1] >= 3 &
            2 * counts[i] <= counts[i - 1] & 2 * counts[i] <= counts[i + 1]
        mean(dip)
    }

    # Few simulations cannot fill many bins; scale the target with sample
    # size (~2*sqrt(n), capped) so classroom-scale runs get chunky bars
    target_bins <- min(target_bins, max(8, ceiling(2 * sqrt(length(stats)))))

    # Try successively fewer bins until interior empty bins are rare and
    # there is no strong sawtooth.  For lattice data the candidate width
    # snaps to a multiple of the lattice spacing, which removes aliasing
    # immediately; for irregular discrete data (e.g. chi-square) wider
    # bins absorb the gaps.
    bw <- span / target_bins
    ks <- unique(c(target_bins, 24, 20, 16, 13, 10, 8))
    for (k in ks[ks <= target_bins]) {
        cand <- if (lattice) res * max(1, round((span / k) / res)) else span / k
        bw <- cand
        counts <- bin_counts(cand)
        if (frac_empty(counts) <= 0.25 && frac_zigzag(counts) <= 0.2) break
    }

    list(bw = bw, off = min(res, bw) / 2, lattice = lattice, res = res)
}

#' Assign values to bins on an anchored grid
#'
#' @param stats Numeric vector to bin.
#' @param anchor Alignment point of the grid.
#' @param bw Bin width.
#' @param off Offset (typically half the lattice resolution) added before
#'   flooring so lattice values sit strictly inside bins rather than on
#'   edges, making assignment robust to floating-point noise.
#' @param sign +1 or -1; -1 mirrors the grid (bins extend leftward).
#' @param align \code{"edge"} or \code{"center"} (see
#'   \code{\link{choose_binning}}).
#'
#' @return Integer bin indices.  With \code{align = "edge"}, bin \code{j}
#'   covers \code{[anchor + j*bw, anchor + (j+1)*bw)} for \code{sign = 1}
#'   (mirrored for \code{sign = -1}), so non-negative indices lie on the
#'   extreme side of the anchor.
#'
#' @keywords internal
bin_index <- function(stats, anchor, bw, off, sign = 1, align = "edge") {
    s <- sign * (stats - anchor) + if (align == "center") bw / 2 else 0
    as.integer(floor((s + off) / bw))
}

#' Convert bin indices back to x-axis intervals
#'
#' @inheritParams bin_index
#' @param idx Integer bin indices from \code{\link{bin_index}}.
#'
#' @return A list with \code{xmin}, \code{xmax}, and \code{mid} vectors.
#'
#' @keywords internal
bin_xrange <- function(idx, anchor, bw, sign = 1, align = "edge") {
    if (align == "center") {
        xmin <- anchor + (idx - 0.5) * bw
    } else if (sign > 0) {
        xmin <- anchor + idx * bw
    } else {
        xmin <- anchor - (idx + 1) * bw
    }
    list(xmin = xmin, xmax = xmin + bw, mid = xmin + bw / 2)
}

#' Horizontal justification that keeps a label inside the panel
#'
#' @param x Label position.
#' @param lims Length-2 panel x-limits.
#' @return An hjust value: labels near the right edge are right-justified,
#'   near the left edge left-justified, otherwise centered.
#' @keywords internal
inward_hjust <- function(x, lims) {
    pos <- (x - lims[1]) / (lims[2] - lims[1])
    if (pos > 0.82) 1 else if (pos < 0.18) 0 else 0.5
}

#' Vertical semi-axis for dotplot dots
#'
#' Chooses the dot half-height so dots render roughly circular at the
#' default Jamovi plot size (400 x 350), capped so stacked dots never
#' overlap.  With tall stacks the dots flatten into ovals rather than
#' compressing the x-axis (which is what a forced equal-coordinate
#' aspect ratio used to do).
#'
#' @param a Dot semi-width in x units.
#' @param x_span,y_span Current axis spans.
#' @return Dot semi-height in y units.
#' @keywords internal
dot_semi_height <- function(a, x_span, y_span) {
    panel_px_x <- 345
    panel_px_y <- 255
    min(0.45, a * (panel_px_x / x_span) * (y_span / panel_px_y))
}
