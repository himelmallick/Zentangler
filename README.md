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
remotes::install_github("himelmallick/Zentangler")
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
  method_preset = "fast_lasso",
  sis_n = 50,
  fdr_method = "BH",
  fdr_scope = "global",
  bootstrap_repeats = 100
)

zentangler_top_mediators(fit, n = 10)
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

## Method Presets

Named presets provide coherent starting configurations while keeping every low-level option available.

```r
zentangler_method_presets()
```

Available presets are:

- `custom`: preserve individually supplied options
- `fast_lasso`: LM A-stage, SIS screening, early-fusion lasso
- `elastic_net`: LM A-stage, SIS screening, early-fusion elastic net
- `longitudinal_maaslin2`: MaAsLin2 A-stage, SIS screening, early-fusion lasso
- `full_exploratory`: LM A-stage, no hard screening, early-fusion elastic net

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

The final mediator table contains both global and view-specific q-values so users can decide whether the multiple-testing family should be all mediators together or mediators within each assay view.

## A-Stage Options

The default A-stage is a HIMA-like univariate linear model:

```r
a_stage_model = "lm"
```

For longitudinal or repeated-measures data, Zentangler can use MaAsLin2 for the exposure-to-mediator leg:

```r
fit <- fit_multiview_parallel_zentangler(
  mae = mae,
  x_var = "X",
  y_var = "Y",
  a_stage_model = "maaslin2",
  maaslin2_random_effect = "SubjectID",
  maaslin2_normalization = "NONE",
  maaslin2_transform = "NONE",
  maaslin2_analysis_method = "LM"
)
```

In MaAsLin2 mode, the A-stage model is fit as:

```text
mediator feature ~ X + covariates + optional random effects
```

The B/Y-stage fusion model remains Zentangler's multiview model.

## Fusion Modes

### Early Fusion

All selected mediators from all assay views are concatenated and fit in one penalized `glmnet` model.

```r
fusion_mode = "early"
```

Use `glmnet_alpha` to move from lasso to elastic net:

```r
glmnet_alpha = 1    # lasso
glmnet_alpha = 0.5  # elastic net
glmnet_alpha = 0    # ridge
```

### Intermediate Fusion

A cooperative-learning style model fits view-specific sparse predictors and encourages agreement across views.

```r
fusion_mode = "intermediate"
```

This mode currently supports Gaussian outcomes and uses the cooperative lasso-style update.

### Late Fusion

Each view is modeled separately first using `glmnet`, then view-level predictions are combined in a second-stage model.

```r
fusion_mode = "late"
```

`glmnet_alpha` is also used in the first-stage view-specific `glmnet` models for late fusion.

## Inference Options

### Multiple-Testing Scope

Zentangler supports two final FDR scopes:

```r
fdr_scope = "global"
```

computes the final `q_primary` across all candidate mediators from all assay views.

```r
fdr_scope = "within_view"
```

computes the final `q_primary` separately within each assay view.

The output table always keeps both columns:

- `q_primary_global`: FDR correction across all views together
- `q_primary_within_view`: FDR correction separately within each view
- `q_primary`: the q-value family selected by `fdr_scope`

This lets the same fitted model support two interpretations:

- global FDR asks, "Among all active mediators across the entire multiview study, what proportion may be false discoveries?"
- within-view FDR asks, "Within each assay view, what proportion of active mediators may be false discoveries?"

The correction method is controlled separately:

```r
fdr_method = "BH"  # Benjamini-Hochberg
fdr_method = "BY"  # Benjamini-Yekutieli
```

`BH` is the usual default. `BY` is more conservative and can be useful when dependence across tests is a major concern.

### De-Biased Lasso B-Path Inference

```r
b_inference = "debiased_lasso"
```

Uses a de-biased lasso approximation for B-path p-values in Gaussian outcome models.

If `glmnet_alpha` is set below 1, early/late coefficient fitting uses elastic net, while this p-value routine remains the current HIMA-style lasso de-biasing approximation.

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
  glmnet_alpha = 0.75,
  fdr_method = "BH",
  fdr_scope = "global",
  b_inference = "debiased_lasso"
)

head(fit$combined_mediators)
```

The simulation benchmark runner can evaluate multiple q-thresholds from one fitted model. This is useful on HPC because changing only the q cutoff should not require refitting the same model repeatedly.

```r
res <- run_intersim_zentangler(
  nrep = 10,
  nsample = 1000,
  fusion_modes = c("early", "intermediate", "late"),
  a_stage_model = "maaslin2",
  screen_method = "sis",
  sis_n = 200,
  lambda_choice = "lambda.min",
  glmnet_alpha = 0.5,
  fdr_method = "BH",
  fdr_scope = "both",
  q_threshold = seq(0.05, 0.25, by = 0.05)
)

res$summary
```

With `fdr_scope = "both"`, the summary reports both global and within-view q-value evaluations. With a threshold vector, the summary reports one row per threshold rather than launching a separate fit per threshold.

## Output Object

The fit object is a list containing:

- `settings`: model settings and selected options
- `sample_ids`: samples used after alignment and filtering
- `diagnostics`: runtime, sample counts, feature counts, screening counts, and active B-path counts
- `views`: per-view mediator results
- `combined_mediators`: stacked mediator ranking table across views
- `mediators_all`: alias for the full mediator table
- `mediators_active`: active mediator table using the default q <= 0.25 rule
- `mediators_top`: top 20 mediators by absolute score
- `view_summary`: per-view tested/screened/active mediator counts
- `model_summary`: one-row summary of the run configuration and diagnostics
- `x_to_y_coef`: direct exposure coefficient from the selected mediator model
- `settings$glmnet_alpha`: Y-stage `glmnet` mixing value used for early/late fusion
- `settings$fdr_method`: multiple-testing correction used for q-values (`"BH"` or `"BY"`)
- `settings$fdr_scope`: final q-value scope used by `q_primary` (`"global"` or `"within_view"`)
- `effect_decomposition`: direct, indirect, and total-effect summaries
- `bootstrap`: bootstrap matrices and uncertainty summaries when enabled
- `fits`: fitted model objects when `return_fits = TRUE`

The main table is:

```r
fit$combined_mediators
```

Convenience helpers:

```r
zentangler_model_summary(fit)
zentangler_view_summary(fit)
zentangler_top_mediators(fit, n = 20)
zentangler_active_mediators(fit, q_threshold = 0.25)
summarize_zentangler(fit)
summarize_zentangler(fit, q_threshold = seq(0.05, 0.25, by = 0.05))$threshold_summary
```

Common columns include:

- `omics`: assay view name
- `mediator`: feature name
- `a`: exposure-to-mediator estimate
- `p_a`: a-path p-value
- `q_a`: FDR-adjusted a-path p-value using `settings$fdr_method`
- `b`: mediator-to-outcome estimate
- `p_b`: b-path p-value
- `selected_by_screen`: whether the mediator was retained for the B-stage
- `score`: `a * b`
- `abs_score`: absolute mediation score
- `p_primary`: primary joint evidence p-value
- `q_primary_global`: FDR-adjusted primary p-value across all views
- `q_primary_within_view`: FDR-adjusted primary p-value separately within each assay view
- `q_primary`: selected final q-value according to `settings$fdr_scope`

## Package Functions

- `fit_multiview_parallel_zentangler()`: fit the multiview parallel mediation model
- `zentangler_method_presets()`: inspect named analysis presets
- `zentangler_top_mediators()`: inspect top-ranked mediators
- `zentangler_active_mediators()`: inspect active mediators under a q-value threshold
- `zentangler_view_summary()`: summarize tested/screened/active mediators by view
- `zentangler_model_summary()`: summarize model settings and diagnostics
- `summarize_zentangler()`: return model, view, threshold, top, and active summaries
- `gen_simmba()`: generate MAE simulation data with mediation truth
- `run_intersim_zentangler()`: run InterSIM/SIMMBA simulation benchmarks across fusion modes
- `trigger_InterSIM()`: generate the null InterSIM data object used by the simulator
