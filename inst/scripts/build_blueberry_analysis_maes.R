#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(SummarizedExperiment)
})

subset_blueberry_mae_views <- function(mae, view_names) {
  all_views <- names(MultiAssayExperiment::experiments(mae))
  missing_views <- setdiff(view_names, all_views)
  if (length(missing_views) > 0) {
    stop("Views not found in MAE: ", paste(missing_views, collapse = ", "), call. = FALSE)
  }

  exps <- MultiAssayExperiment::experiments(mae)[view_names]
  sm <- as.data.frame(MultiAssayExperiment::sampleMap(mae), stringsAsFactors = FALSE)
  sm <- sm[sm$assay %in% view_names, , drop = FALSE]
  sm$assay <- factor(sm$assay, levels = view_names)

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(exps),
    colData = S4Vectors::DataFrame(MultiAssayExperiment::colData(mae)),
    sampleMap = S4Vectors::DataFrame(
      assay = sm$assay,
      primary = sm$primary,
      colname = sm$colname
    )
  )
}

save_mae_pair <- function(mae, base_path_no_ext, object_name = "mae") {
  saveRDS(mae, paste0(base_path_no_ext, ".rds"))
  assign(object_name, mae, envir = environment())
  save(list = object_name, file = paste0(base_path_no_ext, ".RData"), envir = environment())
  invisible(TRUE)
}

build_blueberry_analysis_maes <- function(
  input_mae_path = file.path("data", "Blueberry_MAE.rds")
) {
  if (!file.exists(input_mae_path)) {
    stop("Input MAE not found: ", input_mae_path, call. = FALSE)
  }

  mae <- readRDS(input_mae_path)

  maes <- list(
    pathways = subset_blueberry_mae_views(
      mae,
      c("species", "dna_pathways", "rna_pathways")
    ),
    ecs = subset_blueberry_mae_views(
      mae,
      c("species", "dna_ecs", "rna_ecs")
    ),
    all_omics = subset_blueberry_mae_views(
      mae,
      c("species", "dna_pathways", "rna_pathways", "dna_ecs", "rna_ecs")
    )
  )

  out_bases <- c(
    pathways = file.path("data", "Blueberry_MAE_pathways"),
    ecs = file.path("data", "Blueberry_MAE_ecs"),
    all_omics = file.path("data", "Blueberry_MAE_all_omics")
  )

  for (nm in names(maes)) {
    save_mae_pair(maes[[nm]], out_bases[[nm]], object_name = "mae")
  }

  maes
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  input_mae_path <- if (length(args) >= 1L) args[[1L]] else file.path("data", "Blueberry_MAE.rds")
  maes <- build_blueberry_analysis_maes(input_mae_path = input_mae_path)

  for (nm in names(maes)) {
    exp_names <- names(MultiAssayExperiment::experiments(maes[[nm]]))
    message(
      nm,
      ": ",
      paste(exp_names, collapse = ", "),
      " | samples in colData = ",
      nrow(as.data.frame(MultiAssayExperiment::colData(maes[[nm]])))
    )
  }
}
