# Zentangler

Zentangler is an R package for multiview parallel mediation analysis with
`MultiAssayExperiment` input.

It is designed for studies where one exposure is related to many possible
mediators measured across multiple omics or feature views, and the goal is to
rank mediator signals connecting the exposure to an outcome.

```text
Exposure X -> mediator feature M in view k -> outcome Y
```

Each experiment in the `MultiAssayExperiment` is treated as one view. Views can
be microbiome species, KOs, metabolites, transcripts, proteins, methylation
features, pathway scores, or any numeric feature block.

## Installation

```r
# install.packages("remotes")
remotes::install_github("himelmallick/Zentangler")
```

## Quick Start

```r
library(Zentangler)

fit <- fit_multiview_parallel_zentangler(
  mae = mae,
  x_var = "X",
  y_var = "Y",
  method_preset = "fast_lasso",
  sis_n = 50,
  fusion_mode = "early",
  y_family = "gaussian",
  fdr_method = "BH",
  fdr_scope = "global"
)

zentangler_top_mediators(fit, n = 10)
zentangler_view_summary(fit, q_threshold = 0.25)
```

The main fitting function is:

```r
fit_multiview_parallel_zentangler()
```

It aligns samples across MAE views, estimates exposure-to-mediator paths,
screens candidate mediators, fits a multiview outcome model, computes mediator
scores, applies FDR correction, and returns ranked mediator tables.

## Input Data

Zentangler expects a `MultiAssayExperiment` with:

- numeric assays in one or more experiments
- sample identifiers that can be aligned through MAE sample mapping
- outcome, exposure/design variables, and optional covariates in `colData(mae)`

For `SummarizedExperiment` assays, the usual Bioconductor orientation is:

```text
features x samples
```

Zentangler converts each view internally to:

```text
samples x features
```

## What Zentangler Fits

For each feature in each view, Zentangler first estimates the A-path:

```text
M_jk ~ X + covariates
```

Candidate mediators are then screened within each view, unless
`screen_method = "none"`.

The selected mediators are passed to the multiview B/Y-stage model:

```text
Y ~ X + selected mediators across views + covariates
```

Each mediator receives:

```text
score_jk = a_jk * b_jk
```

where:

- `a_jk` is the exposure-to-mediator estimate
- `b_jk` is the mediator-to-outcome estimate from the multiview model
- `score_jk` is the mediation ranking score

The output keeps global and within-view FDR columns so the same fitted model can
be interpreted under different multiple-testing families.

## Supported Study Designs

The default `study_design = "standard"` uses a numeric exposure column directly.

```r
fit <- fit_multiview_parallel_zentangler(
  mae = mae,
  x_var = "X",
  y_var = "Y",
  study_design = "standard"
)
```

For case-control studies, Zentangler can create the 0/1 exposure internally:

```r
fit_case <- fit_multiview_parallel_zentangler(
  mae = mae,
  y_var = "Y",
  study_design = "case_control",
  case_var = "group",
  control_level = "control",
  case_level = "case"
)
```

For time-point contrasts:

```r
fit_time <- fit_multiview_parallel_zentangler(
  mae = mae,
  y_var = "Y",
  study_design = "time",
  time_var = "visit",
  time_ref = "T0",
  time_compare = "T1"
)
```

For combined case-control and time designs, choose whether the modeled exposure
is the case contrast, the time contrast, or the interaction:

```r
fit_case_time <- fit_multiview_parallel_zentangler(
  mae = mae,
  y_var = "Y",
  study_design = "case_control_time",
  case_var = "group",
  control_level = "control",
  case_level = "case",
  time_var = "visit",
  time_ref = "T0",
  time_compare = "T1",
  exposure_role = "interaction",
  add_design_covariates = TRUE
)
```

## Supported Outcomes

Zentangler supports Gaussian, binomial, and survival outcomes.

| Outcome | `y_family` | Required columns | Main model |
| --- | --- | --- | --- |
| Continuous | `"gaussian"` | `y_var` | linear / Gaussian `glmnet` |
| Binary | `"binomial"` | `y_var` coded numerically | logistic / binomial `glmnet` |
| Time-to-event | `"survival"` | `survival_time_var`, `survival_event_var` | Cox / Cox `glmnet` |

Continuous outcome:

```r
fit_gaussian <- fit_multiview_parallel_zentangler(
  mae,
  x_var = "X",
  y_var = "Y",
  y_family = "gaussian"
)
```

Binary outcome:

```r
fit_binary <- fit_multiview_parallel_zentangler(
  mae,
  x_var = "X",
  y_var = "Y",
  y_family = "binomial",
  b_inference = "debiased_logistic_lasso"
)
```

Survival outcome:

```r
fit_survival <- fit_multiview_parallel_zentangler(
  mae,
  x_var = "X",
  y_family = "survival",
  survival_time_var = "time",
  survival_event_var = "status",
  fusion_mode = "early",
  b_inference = "refit"
)
```

For survival, `survival_event_var` should be coded as 0/1. Early and late
fusion are currently supported for survival outcomes. Intermediate cooperative
fusion is currently Gaussian-only.

## A-Stage Models

The A-stage estimates exposure-to-mediator associations inside each view.

| Option | Meaning |
| --- | --- |
| `a_stage_model = "lm"` | HIMA-like univariate linear models |
| `a_stage_model = "maaslin2"` | MaAsLin2 fixed/random-effect models |

Default A-stage:

```r
a_stage_model = "lm"
```

Repeated-measures or longitudinal A-stage with MaAsLin2:

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

In MaAsLin2 mode, each mediator is modeled approximately as:

```text
mediator ~ X + covariates + optional random effects
```

## Screening

Screening controls which mediators enter the multiview B/Y-stage.

| Option | Meaning |
| --- | --- |
| `screen_method = "sis"` | Sure independence screening within each view |
| `screen_method = "none"` | Keep all usable mediators after numeric/variance filtering |

Use `sis_n` to choose how many mediators are retained per view:

```r
screen_method = "sis"
sis_n = 50
sis_rank = "abs_a"   # or "pvalue"
```

## Fusion Paradigms

Fusion controls how selected mediators from multiple views enter the outcome
model.

| Fusion mode | What it does | Outcome support |
| --- | --- | --- |
| `early` | concatenates selected mediators from all views into one penalized model | Gaussian, binomial, survival |
| `intermediate` | cooperative learning-style view-specific sparse models with agreement penalty | Gaussian |
| `late` | fits view-specific penalized models, then combines view predictions in a meta-model | Gaussian, binomial, survival |

### Early Fusion

```r
fit_early <- fit_multiview_parallel_zentangler(
  mae,
  x_var = "X",
  y_var = "Y",
  fusion_mode = "early"
)
```

Early fusion is the most direct multiview model: all selected mediators are
stacked together.

### Intermediate Fusion

```r
fit_intermediate <- fit_multiview_parallel_zentangler(
  mae,
  x_var = "X",
  y_var = "Y",
  fusion_mode = "intermediate",
  coop_rho = 0.2
)
```

Intermediate fusion uses a cooperative learning-style update. The `coop_rho`
parameter controls how strongly view-specific fitted signals are encouraged to
agree.

### Late Fusion

```r
fit_late <- fit_multiview_parallel_zentangler(
  mae,
  x_var = "X",
  y_var = "Y",
  fusion_mode = "late"
)
```

Late fusion fits each view separately, then combines view-level predictions in a
second-stage model.

## Penalized Outcome Models

Early and late fusion use `glmnet`. The `glmnet_alpha` parameter controls the
penalty for early and late fusion:

```r
glmnet_alpha = 1    # lasso
glmnet_alpha = 0.5  # elastic net
glmnet_alpha = 0    # ridge
```

Choose the cross-validated lambda with:

```r
lambda_choice = "lambda.1se"  # more conservative
lambda_choice = "lambda.min"  # less conservative
```

Intermediate fusion currently uses its cooperative lasso-style update and treats
`glmnet_alpha` as lasso-like.

## B-Stage Inference

The B-stage inference option controls how Zentangler assigns p-values to the
mediator-to-outcome coefficient `b`.

| Option | Main use |
| --- | --- |
| `debiased_lasso` | Gaussian de-biased lasso; also routes binomial outcomes to the logistic analogue |
| `debiased_logistic_lasso` | explicit binomial/logistic de-biased lasso path |
| `debiased_cox_lasso` | Cox survival B-path approximation after sparse selection |
| `refit` | active-set refit on selected mediators |
| `bootstrap` | bootstrap B-path p-values and intervals |

Examples:

```r
b_inference = "debiased_lasso"
b_inference = "debiased_logistic_lasso"
b_inference = "debiased_cox_lasso"
b_inference = "refit"
b_inference = "bootstrap"
```

For survival outcomes, the practical options are:

```r
b_inference = "refit"
b_inference = "debiased_cox_lasso"
b_inference = "bootstrap"
```

`debiased_cox_lasso` currently reports a Cox active-set Wald approximation. Use
`bootstrap` when resampling-based B-path uncertainty is desired and runtime is
acceptable.

## Primary Evidence and FDR

The default primary mediator evidence combines A-path and B-path evidence. The
final mediator table includes:

- `p_primary`: primary mediator p-value
- `q_primary_global`: FDR across all mediators from all views
- `q_primary_within_view`: FDR separately within each view
- `q_primary`: selected q-value family based on `fdr_scope`

Choose the correction method:

```r
fdr_method = "BH"  # Benjamini-Hochberg
fdr_method = "BY"  # Benjamini-Yekutieli
```

Choose the final q-value family:

```r
fdr_scope = "global"
fdr_scope = "within_view"
```

Global FDR treats the full multiview mediator search as one family. Within-view
FDR treats each assay view as its own family.

## Bootstrap

Bootstrap can be used for B-path inference and uncertainty summaries:

```r
fit_boot <- fit_multiview_parallel_zentangler(
  mae,
  x_var = "X",
  y_var = "Y",
  b_inference = "bootstrap",
  bootstrap_repeats = 100
)
```

Bootstrap columns include:

- `p_b_bootstrap`
- `b_boot_mean`, `b_boot_sd`, `b_boot_low`, `b_boot_high`
- `score_boot_mean`, `score_boot_sd`, `score_boot_low`, `score_boot_high`
- `score_boot_selection_freq`

For clustered or repeated-measures data:

```r
fit_cluster_boot <- fit_multiview_parallel_zentangler(
  mae,
  x_var = "X",
  y_var = "Y",
  bootstrap_repeats = 100,
  bootstrap_id = "SubjectID"
)
```

## Method Presets

Presets provide coherent starting configurations while keeping every low-level
option available.

```r
zentangler_method_presets()
```

Current presets:

| Preset | Description |
| --- | --- |
| `custom` | preserve user-supplied options |
| `fast_lasso` | LM A-stage, SIS screening, early-fusion lasso |
| `elastic_net` | LM A-stage, SIS screening, early-fusion elastic net |
| `longitudinal_maaslin2` | MaAsLin2 A-stage, SIS screening, early-fusion lasso |
| `full_exploratory` | LM A-stage, no hard screening, early-fusion elastic net |

## Output Object

The fit object is a list. The most commonly used fields are:

| Field | Contents |
| --- | --- |
| `settings` | model settings and selected options |
| `diagnostics` | sample counts, feature counts, runtime, screening counts |
| `views` | per-view mediator results |
| `combined_mediators` | stacked mediator table across all views |
| `mediators_active` | active mediator table using default q <= 0.25 |
| `mediators_top` | top mediators by absolute score |
| `view_summary` | tested, screened, and active mediator counts by view |
| `model_summary` | one-row model/run summary |
| `effect_decomposition` | direct, indirect, total-effect-style summaries |
| `bootstrap` | bootstrap matrices and failures, if bootstrap was enabled |
| `fits` | fitted model objects when `return_fits = TRUE` |

Inspect the full mediator table:

```r
all_mediators <- zentangler_all_mediators(fit)
head(all_mediators)
```

Common mediator-table columns:

| Column | Meaning |
| --- | --- |
| `omics` | assay view |
| `mediator` | feature name |
| `a` | exposure-to-mediator estimate |
| `p_a`, `q_a` | A-path p-value and q-value |
| `b`, `p_b` | B-path estimate and p-value |
| `score`, `abs_score` | mediation score and absolute score |
| `p_primary` | primary mediator evidence p-value |
| `q_primary_global` | global FDR q-value |
| `q_primary_within_view` | within-view FDR q-value |
| `q_primary` | selected q-value based on `fdr_scope` |

Convenience helpers:

```r
zentangler_model_summary(fit)
zentangler_view_summary(fit)
zentangler_top_mediators(fit, n = 20)
zentangler_active_mediators(fit, q_threshold = 0.25)

summary <- summarize_zentangler(
  fit,
  q_threshold = seq(0.05, 0.25, by = 0.05)
)
summary$threshold_summary
```

## Simulation Examples

Zentangler includes SIMMBA/InterSIM-style simulation helpers for checking
mediator recovery when the truth is known.

### Continuous Outcome

```r
sim <- gen_simmba(
  nsample = 100,
  nrep = 1,
  outcome.type = "continuous",
  ygen.mode = "LM",
  seed = 1
)

fit_continuous <- fit_multiview_parallel_zentangler(
  sim$trainMae[[1]],
  x_var = "A",
  y_var = "Y",
  y_family = "gaussian",
  fusion_mode = "early",
  b_inference = "debiased_lasso",
  sis_n = 50,
  screen_method = "sis",
  seed = 1
)

head(zentangler_all_mediators(fit_continuous))
```

### Binary Outcome

```r
sim_binary <- gen_simmba(
  nsample = 100,
  nrep = 1,
  outcome.type = "binary",
  ygen.mode = "LM",
  seed = 1
)

fit_binary <- fit_multiview_parallel_zentangler(
  sim_binary$trainMae[[1]],
  x_var = "A",
  y_var = "Y",
  y_family = "binomial",
  fusion_mode = "early",
  b_inference = "debiased_logistic_lasso",
  sis_n = 50,
  screen_method = "sis",
  seed = 1
)

head(zentangler_all_mediators(fit_binary))
```

### Survival Outcome

```r
sim_survival <- gen_simmba(
  nsample = 100,
  nrep = 1,
  outcome.type = "survival",
  ygen.mode = "LM",
  seed = 1
)

fit_survival_refit <- fit_multiview_parallel_zentangler(
  sim_survival$trainMae[[1]],
  x_var = "A",
  y_family = "survival",
  survival_time_var = "time",
  survival_event_var = "status",
  fusion_mode = "early",
  b_inference = "refit",
  sis_n = 50,
  screen_method = "sis",
  seed = 1
)

fit_survival_cox <- fit_multiview_parallel_zentangler(
  sim_survival$trainMae[[1]],
  x_var = "A",
  y_family = "survival",
  survival_time_var = "time",
  survival_event_var = "status",
  fusion_mode = "early",
  b_inference = "debiased_cox_lasso",
  sis_n = 50,
  screen_method = "sis",
  seed = 1
)

head(zentangler_all_mediators(fit_survival_refit))
head(zentangler_all_mediators(fit_survival_cox))
```

### Benchmark Runner

`run_intersim_zentangler()` evaluates recovery across simulation replicates,
fusion modes, q-value thresholds, and FDR scopes.

```r
res <- run_intersim_zentangler(
  nrep = 10,
  nsample = 1000,
  fusion_modes = c("early", "intermediate", "late"),
  outcome.type = "continuous",
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

When `q_threshold` is a vector, each simulated dataset is fitted once and then
evaluated at each threshold. This is useful for HPC runs because thresholds do
not require separate model fits.

Common recovery metrics include:

| Metric | Meaning |
| --- | --- |
| `n_active` | number of selected active mediators |
| `true_active` | selected mediators that are truly active |
| `false_active` | selected mediators that are false positives |
| `precision` | true_active / n_active |
| `recall` | true_active / total true mediators |
| `fdr` | false_active / n_active |
| `top50_true` | number of true mediators among top 50 ranked mediators |
| `top50_precision` | top50_true / 50 |

## Package Functions

- `fit_multiview_parallel_zentangler()`: fit the multiview mediation model
- `zentangler_method_presets()`: inspect named presets
- `zentangler_all_mediators()`: return the full mediator table
- `zentangler_top_mediators()`: inspect top-ranked mediators
- `zentangler_active_mediators()`: inspect active mediators under a q threshold
- `zentangler_view_summary()`: summarize mediator counts by view
- `zentangler_model_summary()`: summarize settings and diagnostics
- `summarize_zentangler()`: return model, view, threshold, top, and active summaries
- `gen_simmba()`: generate MAE simulation data with mediation truth
- `run_intersim_zentangler()`: run simulation benchmarks across settings
- `zentangler_truth_key()`: standardize simulation truth tables
- `zentangler_score_truth_recovery()`: score mediator recovery against truth
- `trigger_InterSIM()`: generate the null InterSIM data object used by the simulator
