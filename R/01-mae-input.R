# -----------------------------------------------------------------------------
# MultiAssayExperiment helpers
# -----------------------------------------------------------------------------

zentangler_require_mae <- function() {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("Package 'MultiAssayExperiment' is required for Zentangler MAE input.", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required for Zentangler MAE input.", call. = FALSE)
  }
  invisible(TRUE)
}

mae_fill_nonfinite_zero <- function(df_like) {
  M <- as.matrix(df_like)
  M[!is.finite(M)] <- 0
  out <- as.data.frame(M, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- rownames(df_like)
  out
}

mae_rows_mean_by_id <- function(df_like) {
  if (nrow(df_like) == 0) return(df_like)
  M <- as.matrix(df_like)
  grp <- rownames(df_like)
  rs <- rowsum(M, group = grp, reorder = FALSE)
  n_grp <- as.numeric(table(grp)[rownames(rs)])
  out <- rs / n_grp
  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

zentangler_mae_coldata <- function(mae) {
  zentangler_require_mae()
  cd <- as.data.frame(SummarizedExperiment::colData(mae), stringsAsFactors = FALSE)
  if (is.null(rownames(cd)) || any(!nzchar(rownames(cd)))) {
    stop("MAE colData must have rownames as primary sample IDs.", call. = FALSE)
  }
  cd
}

zentangler_mae_extract_block <- function(
  mae,
  experiment_name,
  assay_name = NULL,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
) {
  # Return one MAE experiment as samples x features.
  #
  # SummarizedExperiment assays are Bioconductor-standard features x samples, so
  # we transpose them. Plain matrices/data.frames are accepted too; if their
  # columns match MAE samples/sampleMap, they are treated as features x samples,
  # otherwise rows are treated as samples.
  duplicate_primary <- match.arg(duplicate_primary)
  zentangler_require_mae()

  exps <- MultiAssayExperiment::experiments(mae)
  if (!(experiment_name %in% names(exps))) {
    stop("experiment_name not found in MAE experiments: ", experiment_name, call. = FALSE)
  }

  exp_obj <- exps[[experiment_name]]
  if (inherits(exp_obj, "SummarizedExperiment")) {
    a_names <- SummarizedExperiment::assayNames(exp_obj)
    if (length(a_names) == 0) stop("No assays found in experiment: ", experiment_name, call. = FALSE)
    if (is.null(assay_name)) assay_name <- a_names[1]
    if (!(assay_name %in% a_names)) {
      stop("assay_name '", assay_name, "' not found in experiment '", experiment_name, "'.", call. = FALSE)
    }
    mat <- as.matrix(SummarizedExperiment::assay(exp_obj, assay_name))
    if (is.null(colnames(mat))) {
      stop("Experiment '", experiment_name, "' assay has no sample colnames.", call. = FALSE)
    }
    if (is.null(rownames(mat))) rownames(mat) <- paste0(experiment_name, "_feature_", seq_len(nrow(mat)))
    block <- as.data.frame(t(mat), check.names = FALSE, stringsAsFactors = FALSE)
    rownames(block) <- colnames(mat)
  } else if (is.matrix(exp_obj) || is.data.frame(exp_obj)) {
    mat <- as.matrix(exp_obj)
    if (is.null(rownames(mat))) rownames(mat) <- paste0(experiment_name, "_row_", seq_len(nrow(mat)))

    primary_ids <- rownames(zentangler_mae_coldata(mae))
    sm <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
    mapped_cols <- character(0)
    if (all(c("assay", "colname") %in% colnames(sm))) {
      mapped_cols <- as.character(sm$colname[sm$assay == experiment_name])
    }
    looks_feature_by_sample <- !is.null(colnames(mat)) && (
      any(colnames(mat) %in% primary_ids) || any(colnames(mat) %in% mapped_cols)
    )

    if (looks_feature_by_sample) {
      block <- as.data.frame(t(mat), check.names = FALSE, stringsAsFactors = FALSE)
      rownames(block) <- colnames(mat)
    } else {
      if (is.null(rownames(mat))) {
        stop("Plain matrix experiment '", experiment_name, "' needs rownames as samples or colnames as samples.", call. = FALSE)
      }
      block <- as.data.frame(mat, check.names = FALSE, stringsAsFactors = FALSE)
    }
  } else {
    stop(
      "Unsupported experiment class for '", experiment_name, "': ",
      paste(class(exp_obj), collapse = ", "),
      ". Use SummarizedExperiment-like or matrix/data.frame.",
      call. = FALSE
    )
  }

  sm <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
  if (all(c("assay", "colname", "primary") %in% colnames(sm))) {
    smi <- sm[sm$assay == experiment_name, c("colname", "primary"), drop = FALSE]
    if (nrow(smi) > 0) {
      idx <- match(rownames(block), as.character(smi$colname))
      pri <- as.character(smi$primary[idx])
      use_map <- !is.na(pri) & nzchar(trimws(pri))
      rownames(block)[use_map] <- pri[use_map]
    }
  }

  keep_rows <- !is.na(rownames(block)) & nzchar(trimws(rownames(block)))
  block <- block[keep_rows, , drop = FALSE]
  if (nrow(block) == 0) stop("No valid sample rows after extracting experiment: ", experiment_name, call. = FALSE)

  if (anyDuplicated(rownames(block)) > 0) {
    if (identical(duplicate_primary, "first")) {
      block <- block[!duplicated(rownames(block)), , drop = FALSE]
    } else {
      block <- mae_rows_mean_by_id(block)
    }
  }

  raw <- as.matrix(block)
  num <- suppressWarnings(matrix(
    as.numeric(raw),
    nrow = nrow(raw),
    ncol = ncol(raw),
    dimnames = dimnames(raw)
  ))
  block_num <- as.data.frame(num, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(block_num) <- rownames(block)
  if (fill_nonfinite_zero) block_num <- mae_fill_nonfinite_zero(block_num)
  block_num
}

zentangler_mae_to_blocks <- function(
  mae,
  view_names = NULL,
  assay_names = NULL,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
) {
  duplicate_primary <- match.arg(duplicate_primary)
  zentangler_require_mae()

  exps <- MultiAssayExperiment::experiments(mae)
  if (is.null(view_names)) view_names <- names(exps)
  if (length(view_names) == 0) stop("MAE contains no experiments/views.", call. = FALSE)
  missing_views <- setdiff(view_names, names(exps))
  if (length(missing_views) > 0) {
    stop("Requested views not found in MAE: ", paste(missing_views, collapse = ", "), call. = FALSE)
  }

  if (is.null(assay_names)) {
    assay_names <- setNames(rep(NA_character_, length(view_names)), view_names)
  } else {
    if (is.null(names(assay_names))) {
      if (length(assay_names) != length(view_names)) {
        stop("Unnamed assay_names must have the same length as view_names.", call. = FALSE)
      }
      names(assay_names) <- view_names
    }
    assay_names <- assay_names[view_names]
  }

  safe_names <- sanitize_view_names(view_names)
  blocks <- vector("list", length(view_names))
  names(blocks) <- safe_names
  view_map <- data.frame(
    experiment_name = view_names,
    view = safe_names,
    assay_name = unname(assay_names),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(view_names)) {
    an <- assay_names[[view_names[i]]]
    if (is.na(an) || !nzchar(an)) an <- NULL
    blocks[[safe_names[i]]] <- zentangler_mae_extract_block(
      mae = mae,
      experiment_name = view_names[i],
      assay_name = an,
      duplicate_primary = duplicate_primary,
      fill_nonfinite_zero = fill_nonfinite_zero
    )
  }

  list(blocks = blocks, pheno_df = zentangler_mae_coldata(mae), view_map = view_map)
}

