#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(SummarizedExperiment)
  library(readxl)
})

normalize_blueberry_metadata <- function(metadata_df, sample_col = "sample") {
  stopifnot(sample_col %in% colnames(metadata_df))
  metadata_df <- as.data.frame(metadata_df, check.names = FALSE, stringsAsFactors = FALSE)
  metadata_df[[sample_col]] <- trimws(as.character(metadata_df[[sample_col]]))
  metadata_df <- metadata_df[nzchar(metadata_df[[sample_col]]), , drop = FALSE]
  metadata_df$sample <- metadata_df[[sample_col]]

  group_map <- c(
    X = "13g_blueberry_plus_13g_placebo",
    Y = "26g_blueberry",
    Z = "26g_placebo"
  )

  blueberry_dose_map <- c(X = 13, Y = 26, Z = 0)
  placebo_dose_map <- c(X = 13, Y = 0, Z = 26)

  metadata_df$group <- trimws(as.character(metadata_df$group))
  metadata_df$group_label <- unname(group_map[metadata_df$group])
  metadata_df$blueberry_dose_g <- unname(blueberry_dose_map[metadata_df$group])
  metadata_df$placebo_dose_g <- unname(placebo_dose_map[metadata_df$group])
  metadata_df$treatment_any <- ifelse(metadata_df$group %in% c("X", "Y"), 1L, 0L)
  metadata_df$treatment_full <- ifelse(
    metadata_df$group %in% c("Y", "Z"),
    ifelse(metadata_df$group == "Y", 1L, 0L),
    NA_integer_
  )
  metadata_df$month6 <- if ("timepoint" %in% colnames(metadata_df)) {
    suppressWarnings(as.numeric(metadata_df$timepoint))
  } else {
    NA_real_
  }

  rownames(metadata_df) <- metadata_df$sample
  metadata_df
}

build_month6_outcome_metadata <- function(longitudinal_df) {
  needed <- c("Sample", "subject", "timepoint", "group", "cGMPpmolmL")
  missing_needed <- setdiff(needed, colnames(longitudinal_df))
  if (length(missing_needed) > 0) {
    stop(
      "Longitudinal metadata is missing required columns: ",
      paste(missing_needed, collapse = ", "),
      call. = FALSE
    )
  }

  df <- normalize_blueberry_metadata(longitudinal_df, sample_col = "Sample")
  df$timepoint_num <- suppressWarnings(as.numeric(as.character(df$timepoint)))
  df$cGMPpmolmL <- suppressWarnings(as.numeric(as.character(df$cGMPpmolmL)))

  baseline <- df[df$timepoint_num == 0, c("subject", "cGMPpmolmL"), drop = FALSE]
  month6 <- df[df$timepoint_num == 6, , drop = FALSE]

  if (nrow(month6) == 0) {
    stop("No month-6 rows were found in the longitudinal metadata.", call. = FALSE)
  }

  baseline_idx <- match(month6$subject, baseline$subject)
  month6$cGMP_baseline <- baseline$cGMPpmolmL[baseline_idx]
  month6$cGMP_month6 <- month6$cGMPpmolmL
  month6$cGMP_change <- month6$cGMP_month6 - month6$cGMP_baseline

  # Restore sample IDs as the primary row key.
  rownames(month6) <- month6$sample
  month6
}

numeric_matrix_from_df <- function(df, row_ids, matrix_name) {
  out <- as.matrix(data.frame(
    lapply(df, function(x) suppressWarnings(as.numeric(as.character(x)))),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ))
  rownames(out) <- row_ids
  storage.mode(out) <- "numeric"
  if (!all(is.finite(out[!is.na(out)]))) {
    warning(matrix_name, " contains non-finite values; coercion may have introduced NAs.")
  }
  out
}

build_species_experiment <- function(xlsx_path, metadata_samples) {
  df <- readxl::read_excel(
    xlsx_path,
    sheet = "Table 4 - species abundance",
    .name_repair = "minimal"
  )
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)

  colnames(df)[1] <- "sample"
  sample_col <- "sample"
  sample_ids <- trimws(as.character(df[[sample_col]]))
  keep <- sample_ids %in% metadata_samples
  df <- df[keep, , drop = FALSE]
  sample_ids <- sample_ids[keep]

  assay_df <- df[, -1, drop = FALSE]
  mat_sample_feature <- numeric_matrix_from_df(
    assay_df,
    row_ids = sample_ids,
    matrix_name = "species abundance"
  )
  mat_feature_sample <- t(mat_sample_feature)

  rownames(mat_feature_sample) <- make.unique(colnames(assay_df))
  colnames(mat_feature_sample) <- sample_ids

  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = mat_feature_sample),
    rowData = S4Vectors::DataFrame(
      feature_id = rownames(mat_feature_sample),
      source_sheet = "Table 4 - species abundance",
      row.names = rownames(mat_feature_sample)
    )
  )
}

build_feature_by_sample_experiment <- function(
  xlsx_path,
  sheet_name,
  metadata_samples,
  keep_unstratified_only = FALSE
) {
  df <- readxl::read_excel(
    xlsx_path,
    sheet = sheet_name,
    .name_repair = "minimal"
  )
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)

  colnames(df)[1] <- "feature_id"
  feature_col <- "feature_id"
  feature_ids <- trimws(as.character(df[[feature_col]]))
  sample_cols <- colnames(df)[-1]

  if (keep_unstratified_only) {
    keep_feature <- !grepl("\\|", feature_ids, fixed = FALSE)
    df <- df[keep_feature, , drop = FALSE]
    feature_ids <- feature_ids[keep_feature]
  }

  keep_samples <- sample_cols %in% metadata_samples
  assay_df <- df[, c(FALSE, keep_samples), drop = FALSE]
  kept_sample_cols <- sample_cols[keep_samples]

  mat <- numeric_matrix_from_df(
    assay_df,
    row_ids = make.unique(feature_ids),
    matrix_name = sheet_name
  )
  colnames(mat) <- kept_sample_cols

  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = mat),
    rowData = S4Vectors::DataFrame(
      feature_id = rownames(mat),
      source_sheet = sheet_name,
      row.names = rownames(mat)
    )
  )
}

build_blueberry_mae <- function(
  xlsx_path = file.path("data", "Blueberry_Multiomics.xlsx"),
  metadata_path = file.path("data", "metadata_emma_longitudinal_only.csv"),
  existing_mae_path = file.path("data", "Blueberry_MAE.rds")
) {
  if (file.exists(metadata_path)) {
    metadata_df <- utils::read.csv(
      metadata_path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    metadata_df <- build_month6_outcome_metadata(metadata_df)
  } else {
    metadata_df <- readxl::read_excel(
      xlsx_path,
      sheet = "Table 3 - working metadata"
    )
    metadata_df <- normalize_blueberry_metadata(metadata_df)
  }
  metadata_samples <- rownames(metadata_df)

  if (file.exists(xlsx_path)) {
    experiments <- list(
      species = build_species_experiment(xlsx_path, metadata_samples),
      dna_pathways = build_feature_by_sample_experiment(
        xlsx_path = xlsx_path,
        sheet_name = "Table 5 - pathways DNA",
        metadata_samples = metadata_samples,
        keep_unstratified_only = TRUE
      ),
      rna_pathways = build_feature_by_sample_experiment(
        xlsx_path = xlsx_path,
        sheet_name = "Table 6 - pathways RNA",
        metadata_samples = metadata_samples,
        keep_unstratified_only = TRUE
      ),
      dna_ecs = build_feature_by_sample_experiment(
        xlsx_path = xlsx_path,
        sheet_name = "Table 7 - ECs DNA",
        metadata_samples = metadata_samples,
        keep_unstratified_only = TRUE
      ),
      rna_ecs = build_feature_by_sample_experiment(
        xlsx_path = xlsx_path,
        sheet_name = "Table 8 - ECs RNA",
        metadata_samples = metadata_samples,
        keep_unstratified_only = TRUE
      )
    )

    sample_map_rows <- lapply(names(experiments), function(view) {
      sample_ids <- colnames(experiments[[view]])
      data.frame(
        assay = view,
        primary = sample_ids,
        colname = sample_ids,
        stringsAsFactors = FALSE
      )
    })
    sample_map <- do.call(rbind, sample_map_rows)
    sample_map <- S4Vectors::DataFrame(
      assay = factor(sample_map$assay, levels = names(experiments)),
      primary = sample_map$primary,
      colname = sample_map$colname
    )
  } else if (file.exists(existing_mae_path)) {
    existing_mae <- readRDS(existing_mae_path)
    experiments <- MultiAssayExperiment::experiments(existing_mae)
    sample_map <- MultiAssayExperiment::sampleMap(existing_mae)
  } else {
    stop(
      "Neither workbook nor existing MAE was found. Looked for: ",
      xlsx_path,
      " and ",
      existing_mae_path,
      call. = FALSE
    )
  }

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(experiments),
    colData = S4Vectors::DataFrame(
      metadata_df,
      row.names = rownames(metadata_df)
    ),
    sampleMap = sample_map
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  xlsx_path <- if (length(args) >= 1L) args[[1L]] else file.path("data", "Blueberry_Multiomics.xlsx")
  out_path <- if (length(args) >= 2L) args[[2L]] else file.path("data", "Blueberry_MAE.rds")
  metadata_path <- if (length(args) >= 3L) args[[3L]] else file.path("data", "metadata_emma_longitudinal_only.csv")
  existing_mae_path <- if (length(args) >= 4L) args[[4L]] else file.path("data", "Blueberry_MAE.rds")

  mae <- build_blueberry_mae(
    xlsx_path = xlsx_path,
    metadata_path = metadata_path,
    existing_mae_path = existing_mae_path
  )
  saveRDS(mae, out_path)

  message("Saved MultiAssayExperiment to: ", normalizePath(out_path, mustWork = FALSE))
  message("Views: ", paste(names(MultiAssayExperiment::experiments(mae)), collapse = ", "))
  message("Samples in colData: ", nrow(SummarizedExperiment::colData(mae)))
}
