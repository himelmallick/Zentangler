# Zentangler

Zentangler is an R package for multi-view parallel mediation analysis using `MultiAssayExperiment` input.

The package is designed for studies where an exposure is related to many candidate mediators across multiple assay views, and the goal is to rank mediator signals that connect the exposure to an outcome.

```text
Exposure X -> mediator feature in assay view k -> Outcome Y
```

Each assay in the `MultiAssayExperiment` is treated as one view. Examples of views include microbiome, metabolome, transcriptome, proteome, methylation, pathway abundance, or any numeric feature block stored as an assay.

## Installation

```r
# install.packages("remotes")
remotes::install_github("YOUR_USERNAME/Zentangler")
```

## Main Function

```r
fit_multiview_parallel_zentangler()
```

The function takes a `MultiAssayExperiment`, extracts assay views, aligns samples through MAE sample identifiers, estimates mediator paths, fits a multiview outcome model, and returns ranked mediator tables.

## Quick Start

```r
library(Zentangler)

fit <- fit_multiview_parallel_zentangler(
  mae = mae,
  x_var = "X",
  y_var = "Y",
  sis_n = 50,
  screen_method = "sis",
  fusion_mode = "early",
  b_inference = "debiased_lasso",
  bootstrap_repeats = 100
)

head(fit$combined_mediators)
```

`x_var` and `y_var` must be columns in `colData(mae)`.

## Input Structure

Zentangler expects a `MultiAssayExperiment` with:

- experiments containing numeric assay features
- samples represented in MAE sample identifiers
- exposure, outcome, and optional covariates in `colData(mae)`

For `SummarizedExperiment` assays, the standard Bioconductor orientation is used:

```text
features x samples
```

Zentangler internally converts each assay to:

```text
samples x features
```

for model fitting.

## Method Overview

For each assay view `k` and feature `j`, Zentangler estimates an a-path:

```text
M_jk ~ X + covariates
```

Then it either screens mediators within each view using sure independence screening or keeps all mediators if `screen_method = "none"`.

Selected mediators are passed to a multiview outcome model:

```text
Y ~ X + selected mediators across views + covariates
```

Each mediator receives:

```text
score_jk = a_jk * b_jk
```

where:

- `a_jk` is the exposure-to-mediator coefficient
- `b_jk` is the mediator-to-outcome coefficient from the multiview outcome model
- `score_jk` is the mediation ranking score

## Fusion Modes

### Early Fusion

All selected mediators from all assay views are concatenated and fit in one penalized model.

```r
fusion_mode = "early"
```

### Intermediate Fusion

A cooperative-learning style model fits view-specific sparse predictors and encourages agreement across views.

```r
fusion_mode = "intermediate"
```

This mode currently supports Gaussian outcomes.

### Late Fusion

Each view is modeled separately first, then view-level predictions are combined in a second-stage model.

```r
fusion_mode = "late"
```

## Inference Options

### De-Biased Lasso B-Path Inference

```r
b_inference = "debiased_lasso"
```

Uses a de-biased lasso approximation for B-path p-values in Gaussian outcome models.

### Active-Set Refit

```r
b_inference = "refit"
```

Refits an ordinary model on the selected active mediator set.

### Bootstrap Uncertainty

```r
bootstrap_repeats = 100
```

Adds bootstrap summaries for mediator scores, including:

- `score_boot_mean`
- `score_boot_sd`
- `score_boot_low`
- `score_boot_high`
- `score_boot_selection_freq`

For repeated-measures or clustered data, use a subject or cluster column in `colData(mae)`:

```r
fit <- fit_multiview_parallel_zentangler(
  mae = mae,
  x_var = "X",
  y_var = "Y",
  bootstrap_repeats = 100,
  bootstrap_id = "SubjectID"
)
```

## Simulation Example

Zentangler includes simulation utilities for method development.

```r
library(Zentangler)

sim <- gen_simmba(
  nsample = 100,
  nrep = 1,
  outcome.type = "continuous",
  seed = 1
)

mae_train <- sim$trainMae[[1]]
truth <- sim$truthDat[[1]]

fit <- fit_multiview_parallel_zentangler(
  mae = mae_train,
  x_var = "A",
  y_var = "Y",
  sis_n = 30,
  screen_method = "sis",
  fusion_mode = "early",
  b_inference = "debiased_lasso"
)

head(fit$combined_mediators)
```

A longer simulation runner is installed with the package:

```r
system.file("scripts/run_simmba_zentangler_lasso_fusion.R", package = "Zentangler")
```

## Output Object

The fit object is a list containing:

- `settings`: model settings and selected options
- `sample_ids`: samples used after alignment and filtering
- `views`: per-view mediator results
- `combined_mediators`: stacked mediator ranking table across views
- `x_to_y_coef`: direct exposure coefficient from the selected mediator model
- `effect_decomposition`: direct, indirect, and total-effect summaries
- `bootstrap`: bootstrap matrices and uncertainty summaries when enabled
- `fits`: fitted model objects when `return_fits = TRUE`

The main table is:

```r
fit$combined_mediators
```

Common columns include:

- `omics`: assay view name
- `mediator`: feature name
- `a`: exposure-to-mediator estimate
- `p_a`: a-path p-value
- `q_a`: BH-adjusted a-path p-value
- `b`: mediator-to-outcome estimate
- `p_b`: b-path p-value
- `selected_by_screen`: whether the mediator was retained for the B-stage
- `score`: `a * b`
- `abs_score`: absolute mediation score
- `p_primary`: primary joint evidence p-value
- `q_primary`: BH-adjusted primary p-value

## Package Functions

- `fit_multiview_parallel_zentangler()`: fit the multiview parallel mediation model
- `gen_simmba()`: generate MAE simulation data with mediation truth
- `run_intersim_zentangler()`: run InterSIM/SIMMBA simulation benchmarks across fusion modes
- `trigger_InterSIM()`: generate the null InterSIM data object used by the simulator
