# gsDose

gsDose is an R package that provides functions to calculate combination test statistics along with the closed testing procedure
for a survival endpoint in a multi-arm multi-stage group sequential phase 2/3 trial with dose selection.

It implements the design proposed in:

> Zhang, E. (P.) and Jin, M. (2025). A Multi-Arm Multi-Stage Group Sequential Phase 2/3 Design with Dose Selection for Oncology Trials. *Statistics in Biopharmaceutical Research*. [doi:10.1080/19466315.2025.2539831](https://doi.org/10.1080/19466315.2025.2539831)

## Installation

```r
remotes::install_github("eric-zhang16/gsDose", build_vignettes = TRUE, force = TRUE)
```

## Overview

The package provides three functions:

| Function | Purpose |
|----------|---------|
| `gsBoundary()` | Group sequential efficacy (or futility) boundaries for the inverse normal combination test, using an alpha spending function (O'Brien-Fleming or Hwang-Shih-DeCani) |
| `combTest()` | Inverse normal combination test statistics for all intersection hypotheses containing the selected dose, with the closed testing decision |
| `combDunnett()` | Inverse normal combination z statistic based on Dunnett's test for one intersection hypothesis (called by `combTest()`) |

## How to use the package

The example below follows the package vignette: a trial starts with three dose arms and a common
control arm (1:1:1:1). At IA1, one dose is selected; enrollment then continues 1:1 between the selected
dose arm and control until 700 patients are enrolled across the two arms, of which 60 come from cohort 1
(patients enrolled before dose selection). OS is tested at IA2, IA3 and FA, with the final analysis planned
at 520 OS events and a one-sided alpha of 0.025.

```r
library(gsDose)

# weights combining cohort 1 (enrolled before dose selection) and cohort 2
# (enrolled after), proportional to the sample size of each cohort
h <- c(sqrt(60/700), sqrt(1 - 60/700))
```

### Step 1: boundary at the first efficacy analysis (IA2)

Suppose dose 2 is selected at IA1, and IA2 occurs at 270 OS events, 231 of which come from cohort 2.
At the first testing stage all the cohort-1 weight is on the data accrued up to IA2, so `w = 1`:

```r
gsBoundary(w = 1, h = h, d = 270, d2 = 231, s = 3, planD = 520,
           alpha = 0.025, sf = 'OF', bt = 'upper')
#> $Z
#>   stage   Zbound
#> 1     1 2.899754
```

### Step 2: closed testing at IA2

Provide the cohort-1 standardized log-rank z statistics and event counts for each dose arm versus
control, and the cohort-2 log-rank z statistic for the selected dose. `combTest()` evaluates every
intersection hypothesis containing the selected dose (`sid = 2`):

```r
z1.ia2 <- matrix(c(0.8, 1, 0.7), nrow = 1)   # dose 1, 2, 3 vs control (cohort 1)
d1.ia2 <- matrix(c(40, 36, 38), nrow = 1)    # cohort-1 events per comparison

combTest(w = 1, h = h, z1 = z1.ia2, d1 = d1.ia2, z2 = 2.01, sid = 2, zbound = 2.9)
#>    dose 1 dose 2 dose 3 stage DunnettZ1 LogrankZ2   ComboZ Reject
#> H1      1      1      1     1 0.4614985      2.01 2.057041     No
#> H2      1      1      0     1 0.6594718      2.01 2.115001     No
#> H3      0      1      1     1 0.6594718      2.01 2.115001     No
#> H4      0      1      0     1 1.0000000      2.01 2.214698     No
```

No intersection hypothesis is rejected, so the trial continues to IA3.

### Step 3: update the boundary at IA3

IA3 is triggered at 430 events (385 from cohort 2). The cohort-1 weights across stages are
pre-specified proportional to the expected incremental events (here 85% / 15%):

```r
w <- matrix(c(1,          0,
              sqrt(0.85), sqrt(0.15)), nrow = 2, byrow = TRUE)
h.m <- matrix(rep(h, 2), nrow = 2, byrow = TRUE)

gsBoundary(w = w, h = h.m, d = c(270, 430), d2 = c(231, 385), s = 3,
           planD = 520, alpha = 0.025, sf = 'OF', bt = 'upper')
#> $Z
#>   stage   Zbound
#> 1     1 2.899754
#> 2     2 2.223579
```

### Step 4: closed testing at IA3

```r
z1.ia3 <- matrix(c(0.8, 1,   0.7,
                   0.9, 1.2, 0.81), nrow = 2, byrow = TRUE)
d1.ia3 <- matrix(c(40, 36, 38,
                   51, 42, 49), nrow = 2, byrow = TRUE)

combTest(w = c(sqrt(0.85), sqrt(0.15)), h = h, z1 = z1.ia3, d1 = d1.ia3,
         z2 = 2.35, sid = 2, zbound = 2.22)
#>    dose 1 dose 2 dose 3 stage DunnettZ1 LogrankZ2   ComboZ Reject
#> H1      1      1      1     2 0.4841621      2.35 2.388778    Yes
#> H2      1      1      0     2 0.7480401      2.35 2.466033    Yes
#> H3      0      1      1     2 0.7480401      2.35 2.466033    Yes
#> H4      0      1      0     2 1.2029052      2.35 2.599204    Yes
```

All intersection hypotheses containing H02 are rejected, so H02 is rejected overall: the selected
dose demonstrates statistically significant improvement in OS over control.

### Unequal randomization

By default the Dunnett test assumes equal allocation across arms (correlation 0.5). For unequal
dose-to-control randomization, pass the allocation ratio (n_dose / n_control); the correlation is
derived as corr(i, j) = sqrt(lambda_i * lambda_j) with lambda = ratio / (1 + ratio):

```r
# dose arms randomized 2:1 against control
combTest(w = c(sqrt(0.85), sqrt(0.15)), h = h, z1 = z1.ia3, d1 = d1.ia3,
         z2 = 2.35, sid = 2, zbound = 2.22, ratio = 2)
```

### Futility bounds

Lower (futility) bounds are obtained by spending the complement of alpha with `bt = 'lower'`,
e.g. `gsBoundary(..., alpha = 1 - 0.025, sf = 'HSD', sfpar = -6, bt = 'lower')`.

## Full worked example

The vignette walks through the complete flow, including deriving the pre-specified weights from the
expected number of events using gsDesign2:

```r
vignette("Introduction", package = "gsDose")
```
