# Randomize

A [Jamovi](https://www.jamovi.org/) module and R package for teaching
randomization-based statistical inference in introductory statistics courses.

## What it does

Randomize provides **permutation tests** and **bootstrap confidence intervals**
for the most common introductory scenarios, plus a model-based calculator
(tail areas and CI multipliers for the normal, t, chi-square, and F
distributions):

| Analysis | CI | Hypothesis test |
|---|---|---|
| Single mean | `SingleMeanCI` | — |
| Difference in means (independent) | `twomeanCI` | `twomeanhtest` |
| Difference in means (paired) | `pairedmeanCI` | `pairedmeanhtest` |
| Multiple means (ANOVA) | — | `multimeanhtest` |
| Slope (simple linear regression) | `slopeCI` | `slopehtest` |
| Single proportion | `SinglePropCI` | `SinglePropHTest` |
| Difference in proportions | `TwoPropCI` | `TwoPropHTest` |
| Contingency table (Chi-square) | — | `ContTabHTest` |
| Model-based (normal / t / chi-square / F) | — | `modelBased` |

Every analysis produces:

- **Simulation plots** — dotplot or histogram of the bootstrap/permutation
  distribution with CI bounds or p-value tail shading; axes name the groups
  and the level being compared, bootstrap plots mark the observed statistic,
  and up to 25 simulations are drawn one dot per value
- **Results tables** — observed statistics, confidence intervals or p-values
- **Descriptive statistics** — sample sizes, means, medians, SDs by group

## Installation

### In Jamovi

Install the `.jmo` file from the Jamovi module library, or build from source:

```r
# Requires jmvtools
install.packages("jmvtools", repos = c("https://repo.jamovi.org", "https://cran.r-project.org"))
jmvtools::install()
```

### As an R package

```r
# Install from GitHub (requires devtools)
devtools::install_github("byurk/Randomize")
```

## Usage from R

Every analysis can be called directly from R — useful for generating output
in Quarto slides, R Markdown labs, or any other R-based workflow.

### Quick example

```r
library(Randomize)

# Create some data
d <- data.frame(
    score = c(rnorm(25, 10, 2), rnorm(25, 12, 2)),
    group = factor(rep(c("Control", "Treatment"), each = 25))
)

# Two-sample permutation test (desc/plots = TRUE also fill the
# descriptive table and the means plot used further down)
r <- twomeanhtest(
    data = d,
    vars = "score",
    group = "group",
    hypothesis = "different",
    reps = 1000,
    dotHist = "dotplot",
    seedBool = TRUE,
    rngSeed = 42,
    desc = TRUE,
    plots = TRUE
)

# Print the results table
r

# Show the permutation distribution
plot(r)

# Extract results as a data frame
results_table(r)
```

### Extracting plots

Use `plot()` to extract ggplot objects from any analysis result:

```r
plot(r)              # simulation/bootstrap distribution (all analyses)
plot(r, show_counts = TRUE)  # label each bar / dot stack with its count
plot(r, show_text = FALSE)   # no labels or caption (bootstrap: keeps the lines)
plot(r, "desc")      # observations with mean/median markers (mean-comparison analyses)
plot(r, "line")      # regression scatterplot with the fitted equation (slope analyses)
```

Since these return standard ggplot objects, you can modify them:

```r
plot(r) + ggplot2::ggtitle("Permutation Test Results")
```

### Extracting tables

```r
results_table(r)     # main results (p-value, CI bounds, etc.)
desc_table(r)        # descriptive statistics (n, mean, median, sd)
```

These return data frames, so they work with `knitr::kable()`,
`gt::gt()`, or any other table-formatting package.

### Common options

All simulation-based analyses share these options:

| Option | Description | Default |
|--------|-------------|---------|
| `reps` | Number of bootstrap/permutation replicates | 1000 |
| `dotHist` | Plot type: `"dotplot"` or `"histogram"` | `"dotplot"` |
| `showCounts` | Print the count above each bar / dot stack | `FALSE` |
| `seedBool` | Use a fixed random seed? | `FALSE` |
| `rngSeed` | The seed value (when `seedBool = TRUE`) | 8675309 |

A simulation p-value of exactly 0 (no simulated statistic as extreme as the
observed one) is reported as `"< 1/reps"`, e.g. `< .001` for 1000 reps.

Bootstrap CI analyses also accept:

| Option | Description | Default |
|--------|-------------|---------|
| `confLevel` | Confidence level (percentage) | 95 |
| `ciType` | `"bootperc"` (percentile) or `"bootse"` (SE method) | `"bootperc"` |

Hypothesis test analyses accept:

| Option | Description |
|--------|-------------|
| `hypothesis` | `"different"`, `"oneGreater"`, `"twoGreater"` (two-sample); `"notequal"`, `"greater"`, `"less"` (slope/proportion) |

Proportion analyses report the proportion of the first level of the response
variable (alphabetical for text data), matching jamovi's own binomial test.

### Analysis examples

#### Bootstrap CI for a single mean

```r
r <- SingleMeanCI(data = d, resp = "score",
                  reps = 1000, confLevel = 95, ciType = "bootperc",
                  dotHist = "histogram", seedBool = TRUE, rngSeed = 42)
plot(r)
results_table(r)
```

#### Paired means permutation test

```r
r <- pairedmeanhtest(data = d,
                     pairs = list(list(i1 = "pretest", i2 = "posttest")),
                     hypothesis = "different", reps = 1000,
                     dotHist = "dotplot", seedBool = TRUE, rngSeed = 42)
plot(r)
```

Several pairs can be given; each gets its own table row and simulation
plot (`r$simplot` is then an array keyed by the pairs, and `plot(r)` shows
the first pair's).

#### Slope hypothesis test

```r
r <- slopehtest(data = d, dep = "y", indep = "x",
                hypothesis = "notequal", reps = 1000,
                dotHist = "dotplot", seedBool = TRUE, rngSeed = 42,
                plots = TRUE)
plot(r)          # permutation distribution for slope
plot(r, "line")  # scatterplot with regression line
```

#### Model-based calculator

```r
r <- modelBased(distro = "chisq", dF = 3, areaBool = TRUE, obsStat = 7.8)
results_table(r)   # observed value, df, right-tail area
plot(r)            # density with the tail shaded from the observed value

r <- modelBased(distro = "tdistro", dF = 24, CIBool = TRUE, confLevel = 95)
r$multTable$asDF   # t* multiplier
```

#### Chi-square contingency table test

```r
r <- ContTabHTest(data = d, rows = "treatment", cols = "outcome",
                  reps = 1000, dotHist = "histogram",
                  seedBool = TRUE, rngSeed = 42, compare = "rows")
plot(r)
```

## Upgrading from an earlier version

Jamovi files (`.omv`) saved with an older Randomize keep their stored
results when reopened; they are not recomputed automatically. After
installing a new version, click an analysis and change any option (or
re-add a variable) to re-run it with the new code.

## Development

```bash
# Clone the repo
git clone git@github.com:byurk/Randomize.git
cd Randomize

# Load for development (without installing)
Rscript -e "devtools::load_all('.')"

# Build and install the Jamovi module
Rscript -e "jmvtools::install()"

# Run tests
Rscript tests/test_analyses.R
```

## Acknowledgments

This package is built on the [Jamovi](https://www.jamovi.org/) module
framework ([jmvcore](https://github.com/jamovi/jmvcore)). Contingency table
data-handling code is adapted from the
[jmv](https://github.com/jamovi/jmv) package. Both are licensed under GPL.

## License

GPL (>= 3)
