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
#' @param max_empty Highest tolerated fraction of interior empty bins
#'   for non-lattice data.  The default 0 means histograms and dotplots
#'   of continuous statistics never show gaps; a genuinely detached
#'   outlier still gets its own bar with honest empty space via the
#'   windowed second pass.  Lattice data are handled differently: see
#'   Details.
#'
#' @details
#' For densely occupied lattice data (proportions \eqn{k/n}, where every
#' lattice point in the body of the distribution is achievable) the bin
#' width is always a whole number of lattice steps, so an empty bin is
#' never a binning artifact -- it is a value that genuinely never
#' occurred in the simulation.  Such honest gaps (typically in the
#' tails) therefore do not force coarser bins.  Binning starts from one
#' bin per achievable value whenever the central 99\% of
#' the simulations spans no more than about 1.6 times the target bin
#' count, and coarsens only when the result shows a sawtooth (mixed
#' lattices such as a difference in proportions with unequal group sizes)
#' or an empty bin in the body of the distribution (both nearest occupied
#' neighbours holding at least 5 simulations).  One bin per value is what
#' lets a dot stack sit exactly on the observed-value line.
#'
#' @return A list with \code{bw} (bin width), \code{off} (half-resolution
#'   offset that keeps lattice values away from bin edges), \code{lattice}
#'   (logical), and \code{res} (lattice resolution).
#'
#' @keywords internal
choose_binning <- function(stats, anchor, align = c("edge", "center"),
                           sign = 1, target_bins = 30, max_empty = 0) {
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

    bin_counts <- function(bw, lo, hi) {
        idx <- bin_index(stats, anchor, bw, min(res, bw) / 2, sign, align)
        qs <- stats::quantile(idx, c(lo, hi), type = 1, names = FALSE)
        win <- idx[idx >= qs[1] & idx <= qs[2]]
        tabulate(win - min(win) + 1L)
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

    # An empty bin between two well-populated bins (each holding at least
    # min_n simulations).  For lattice data this is the only kind of gap
    # worth coarsening for: a hole in the body of the distribution looks
    # like a binning error even when it is real, whereas an empty value in
    # a sparse tail is exactly what a student should see.
    body_gap <- function(counts, min_n = 5) {
        occ <- which(counts > 0)
        if (length(occ) < 2) return(FALSE)
        left <- occ[-length(occ)]
        right <- occ[-1]
        any(right - left > 1 & counts[left] >= min_n & counts[right] >= min_n)
    }

    # Few simulations cannot fill many bins; scale the target with sample
    # size (~2*sqrt(n), capped) so classroom-scale runs get chunky bars
    target_bins <- min(target_bins, max(8, ceiling(2 * sqrt(length(stats)))))

    ks <- unique(c(target_bins, 24, 20, 16, 13, 10, 8, 6, 5))
    ks <- ks[ks <= target_bins]

    # A lattice is "dense" when the achievable values actually fill it:
    # proportions k/n occupy every lattice point in the body of the
    # distribution.  Chi-square statistics also lie on a (very fine)
    # rational lattice but occupy only a sparse, irregular subset of it,
    # so their empty bins are binning artifacts, not honest zero counts;
    # they take the non-lattice ladder below.
    qs <- stats::quantile(stats, c(0.005, 0.995), names = FALSE)
    steps <- (qs[2] - qs[1]) / res
    if (lattice) {
        inwin <- u[u >= qs[1] & u <= qs[2]]
        dense <- length(inwin) >= 0.6 * (steps + 1)
    } else {
        dense <- FALSE
    }

    if (dense) {
        # One bin per achievable value if the central 99% of the
        # simulations fits in ~1.6x the target bin count; otherwise the
        # smallest whole number of lattice steps that does.  Coarsen from
        # there only for a sawtooth or a body gap, and never beyond the
        # coarsest width the non-lattice ladder would use.
        max_cols <- ceiling(1.6 * target_bins)
        m <- max(1L, as.integer(ceiling(steps / max_cols)))
        m_max <- max(m, as.integer(round((span / min(ks)) / res)))
        while (m < m_max) {
            counts <- bin_counts(m * res, 0, 1)
            if (m == 1L) {
                # single-value bins: gaps are honest zero counts
                ok <- !body_gap(counts)
            } else {
                # grouped bins: same gap-free rule as the ladder below
                # (full range, else the central 99.8%)
                ok <- frac_empty(counts) <= max_empty ||
                    frac_empty(bin_counts(m * res, 0.001, 0.999)) <= max_empty
            }
            if (ok && frac_zigzag(counts) <= 0.2) break
            m <- m + 1L
        }
        bw <- m * res
        return(list(bw = bw, off = min(res, bw) / 2, lattice = TRUE, res = res))
    }

    # Everything else: try successively fewer bins until interior empty
    # bins are rare and there is no strong sawtooth.  For sparse lattice
    # data (e.g. chi-square) the candidate width still snaps to a
    # multiple of the lattice spacing, which removes aliasing
    # immediately; wider bins absorb the gaps.
    #
    # Two passes: the first evaluates occupancy over the FULL range, so
    # ordinary distributions come out with no gaps anywhere.  Only when
    # no width can manage that (a genuinely detached outlier) does the
    # second pass exclude the extreme 0.1% per side -- the outlier keeps
    # its own bar with honest empty space rather than forcing chunky
    # bins on the whole distribution.
    try_ladder <- function(lo, hi) {
        for (k in ks) {
            cand <- if (lattice) res * max(1, round((span / k) / res)) else span / k
            counts <- bin_counts(cand, lo, hi)
            if (frac_empty(counts) <= max_empty && frac_zigzag(counts) <= 0.2)
                return(cand)
        }
        NULL
    }
    bw <- try_ladder(0, 1)
    if (is.null(bw)) bw <- try_ladder(0.001, 0.999)
    if (is.null(bw))
        bw <- if (lattice) res * max(1, round((span / min(ks)) / res)) else span / min(ks)

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
    # With center alignment, lattice values already sit at bin centers
    # when bw equals the lattice spacing; adding the half-resolution
    # nudge there sums to exactly one bin width (bw/2 + res/2 == bw) and
    # used to shift every value into the next bin -- drawing the whole
    # distribution one lattice step to the right of the data.  The nudge
    # is only needed (and only safe) when it is strictly less than half
    # a bin.
    if (align == "center" && 2 * off >= bw * (1 - 1e-9)) off <- 0
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

#' X-positions for dotplot columns
#'
#' When each bin holds exactly one achievable value (bin width equals
#' the lattice spacing), every dot column sits at its exact value -- so
#' simulations equal to the observed statistic stack directly on the
#' observed-value line instead of at a bin midpoint half a step away
#' (which students can misread as "no simulation matched the observed
#' value").  A single-value column can never mix extreme and non-extreme
#' simulations, so stacks stay one color.  When bins group a range of
#' values, columns sit at bin midpoints as before -- centering the grid
#' on the observed value there would put a mixed-color stack under the
#' line.
#'
#' @param stats Simulated statistics (same order as \code{idx}).
#' @param idx Bin indices from \code{\link{bin_index}}.
#' @param b Binning description from \code{\link{choose_binning}}.
#' @inheritParams bin_index
#' @return X position for each element of \code{stats}.
#' @keywords internal
dot_column_x <- function(stats, idx, b, anchor, sign, align) {
    if (b$lattice && abs(b$bw - b$res) <= b$res * 1e-9)
        return(stats::ave(stats, idx, FUN = function(v) v[1]))
    bin_xrange(idx, anchor, b$bw, sign, align)$mid
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

#' Dotplot dots that stay circular at any device size
#'
#' Draws each simulation as a true circle.  The radius is resolved at
#' draw time from the actual panel geometry as the smaller of the bin
#' semi-width (\code{a}, in x units, so dots never exceed their bin) and
#' the vertical cap (\code{max_b} y-units, so stacked dots never
#' overlap); when the cap binds, the whole dot shrinks -- columns gain a
#' little horizontal breathing room but dots are never squashed into
#' ovals.  Deciding this at draw time (via \code{grid::makeContent}, so
#' it re-runs on every resize) is what keeps dots circular at every
#' device size -- a fixed formula tuned to one panel size distorts
#' everywhere else.
#'
#' @keywords internal
GeomDotStack <- ggplot2::ggproto("GeomDotStack", ggplot2::Geom,
    required_aes = c("x", "y"),
    default_aes = ggplot2::aes(fill = "grey35"),
    draw_key = ggplot2::draw_key_point,
    draw_panel = function(data, panel_params, coord, a = 1, max_b = 0.45) {
        coords <- coord$transform(data, panel_params)
        shifted <- transform(data, x = x + a, y = y + max_b)
        sc <- coord$transform(shifted, panel_params)
        grid::gTree(
            cx = coords$x, cy = coords$y,
            a_npc = sc$x[1] - coords$x[1],
            bcap_npc = sc$y[1] - coords$y[1],
            fill = coords$fill,
            cl = "randomize_dotstack")
    }
)

#' @rdname GeomDotStack
#' @keywords internal
makeContent.randomize_dotstack <- function(x) {
    a_mm <- grid::convertWidth(grid::unit(x$a_npc, "npc"), "mm", valueOnly = TRUE)
    bcap_mm <- grid::convertHeight(grid::unit(x$bcap_npc, "npc"), "mm", valueOnly = TRUE)
    # Radius floor: for very tall stacks the no-overlap cap would shrink
    # dots below visibility, so from there dots keep a legible size and
    # overlap vertically like stacked coins instead (the plot itself
    # disables panel clipping so the top dot survives intact).
    r_mm <- min(a_mm, max(bcap_mm, 1.2))
    dots <- grid::circleGrob(
        x = grid::unit(x$cx, "npc"), y = grid::unit(x$cy, "npc"),
        r = grid::unit(r_mm, "mm"),
        gp = grid::gpar(fill = x$fill, col = NA))
    grid::setChildren(x, grid::gList(dots))
}

#' Layer constructor for \code{GeomDotStack}
#'
#' @param dots Data frame with \code{x}, \code{y} (and optionally
#'   \code{fill}) columns, one row per dot.
#' @param a Dot semi-width in x units.
#' @param max_b Vertical semi-height cap in y units; stacked dots sit 1
#'   y-unit apart, so any value below 0.5 prevents overlap.
#' @param fill Constant fill color; when \code{NULL}, \code{dots$fill}
#'   is mapped through the plot's fill scale instead.
#' @keywords internal
stack_dots <- function(dots, a, max_b = 0.45, fill = NULL) {
    if (is.null(fill)) {
        mapping <- ggplot2::aes(x = x, y = y, fill = fill)
        params <- list(a = a, max_b = max_b)
    } else {
        mapping <- ggplot2::aes(x = x, y = y)
        params <- list(a = a, max_b = max_b, fill = fill)
    }
    ggplot2::layer(
        geom = GeomDotStack, data = dots, mapping = mapping,
        stat = "identity", position = "identity",
        params = params, inherit.aes = FALSE, show.legend = FALSE)
}
