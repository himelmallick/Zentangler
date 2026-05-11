# Zentangler parallel multiview mediation package
#
# The original monolithic implementation has been split across focused files in
# this R/ directory. The public API remains:
#
#   fit_multiview_parallel_zentangler(mae, x_var, y_var, ...)
#
# File map:
#   00-utils.R                  validation, preprocessing, FDR helpers
#   01-mae-input.R              MultiAssayExperiment extraction/alignment
#   02-presets-output-helpers.R presets and output inspection helpers
#   03-a-stage-screening.R      A-path estimation and SIS/MaAsLin2 screening
#   04-y-stage-fusion.R         early/intermediate/late Y-stage fusion models
#   05-b-path-inference.R       B-path refit and de-biased lasso inference
#   06-causal-effects.R         causal direct/indirect/total effect summaries
#   07-bootstrap.R              subject/sample bootstrap utilities
#   08-fit-blocks.R             internal matrix/block fitting engine
#   09-fit-mae.R                public MAE-first fitting interface
