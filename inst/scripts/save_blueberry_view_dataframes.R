#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})

extract_view_df <- function(mae, view_name) {
  se <- MultiAssayExperiment::experiments(mae)[[view_name]]
  mat <- SummarizedExperiment::assay(se, 1)
  out <- as.data.frame(t(mat), check.names = FALSE, stringsAsFactors = FALSE)
  out$sample <- rownames(out)
  out <- out[, c("sample", setdiff(colnames(out), "sample")), drop = FALSE]
  rownames(out) <- out$sample
  out
}

save_blueberry_view_dataframes <- function(
  mae_path = file.path("data", "Blueberry_MAE_all_omics.rds")
) {
  mae <- readRDS(mae_path)

  species_df <- extract_view_df(mae, "species")
  dna_pathways_df <- extract_view_df(mae, "dna_pathways")
  rna_pathways_df <- extract_view_df(mae, "rna_pathways")
  dna_ecs_df <- extract_view_df(mae, "dna_ecs")
  rna_ecs_df <- extract_view_df(mae, "rna_ecs")

  pathways_dfs <- list(
    species = species_df,
    dna_pathways = dna_pathways_df,
    rna_pathways = rna_pathways_df
  )

  ecs_dfs <- list(
    species = species_df,
    dna_ecs = dna_ecs_df,
    rna_ecs = rna_ecs_df
  )

  all_omics_dfs <- list(
    species = species_df,
    dna_pathways = dna_pathways_df,
    rna_pathways = rna_pathways_df,
    dna_ecs = dna_ecs_df,
    rna_ecs = rna_ecs_df
  )

  save(pathways_dfs, file = file.path("data", "Blueberry_pathways_dfs.RData"))
  save(ecs_dfs, file = file.path("data", "Blueberry_ecs_dfs.RData"))
  save(all_omics_dfs, file = file.path("data", "Blueberry_all_omics_dfs.RData"))

  invisible(list(
    pathways_dfs = pathways_dfs,
    ecs_dfs = ecs_dfs,
    all_omics_dfs = all_omics_dfs
  ))
}

if (sys.nframe() == 0L) {
  save_blueberry_view_dataframes()
}
