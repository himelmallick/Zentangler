# -----------------------------------------------------------------------------
# A-stage screening
# -----------------------------------------------------------------------------

hima_sis_default_n <- function(n, p) {
  max(1L, min(p, floor(n / max(log(n), 1.0))))
}

hima_screen_mediators <- function(
  M,
  X,
  C = NULL,
  sis_n = NULL,
  rank_by = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  fdr_method = c("BH", "BY")
) {
  rank_by <- match.arg(rank_by)
  screen_method <- match.arg(screen_method)
  fdr_method <- validate_fdr_method(fdr_method)

  n <- nrow(M)
  p <- ncol(M)

  tab <- data.frame(
    mediator = colnames(M),
    a_hat = NA_real_,
    p_a = NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(p)) {
    dat <- data.frame(y = M[, j], X = X)
    if (!is.null(C) && ncol(C) > 0) dat <- cbind(dat, as.data.frame(C))

    fit <- try(lm(y ~ ., data = dat), silent = TRUE)
    if (inherits(fit, "try-error")) next

    co <- coef(summary(fit))
    if ("X" %in% rownames(co)) {
      tab$a_hat[j] <- co["X", "Estimate"]
      tab$p_a[j] <- co["X", "Pr(>|t|)"]
    }
  }

  tab$q_a <- p.adjust(tab$p_a, method = fdr_method)
  ord <- if (rank_by == "abs_a") {
    order(abs(tab$a_hat), decreasing = TRUE, na.last = TRUE)
  } else {
    order(tab$p_a, na.last = TRUE)
  }

  if (screen_method == "none") {
    # Keep every mediator that survived sample alignment and numeric filtering.
    # A-path estimates are still computed because they are part of the mediation score.
    sis_n_used <- p
    selected <- tab$mediator
  } else {
    if (is.null(sis_n)) sis_n <- hima_sis_default_n(n = n, p = p)
    sis_n_used <- max(1L, min(p, as.integer(sis_n)))
    selected <- tab$mediator[ord][seq_len(sis_n_used)]
  }

  tab$selected <- tab$mediator %in% selected
  tab <- tab[order(tab$p_a, na.last = TRUE), ]
  rownames(tab) <- NULL

  list(table = tab, selected = selected, sis_n = sis_n_used, screen_method = screen_method)
}

pick_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

hima_screen_mediators_maaslin2 <- function(
  M,
  pheno_df,
  x_var,
  covariates = NULL,
  random_effects = NULL,
  sis_n = NULL,
  rank_by = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  standardize = FALSE,
  output_dir = NULL,
  fdr_method = c("BH", "BY")
) {
  # MaAsLin2 A-stage:
  #   mediator_j ~ X + fixed covariates + optional random effects
  #
  # This is mainly for longitudinal/repeated-measures settings where the
  # exposure-to-mediator leg should account for subject-level clustering.
  rank_by <- match.arg(rank_by)
  screen_method <- match.arg(screen_method)
  fdr_method <- validate_fdr_method(fdr_method)

  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    stop("Package 'Maaslin2' is required when a_stage_model = 'maaslin2'.")
  }

  if (is.null(rownames(M)) || is.null(rownames(pheno_df))) {
    stop("M and pheno_df must have sample IDs in rownames for MaAsLin2 screening.")
  }

  ids <- rownames(M)
  if (!all(ids %in% rownames(pheno_df))) {
    stop("Some mediator sample IDs are missing in pheno_df for MaAsLin2 screening.")
  }

  n <- nrow(M)
  p <- ncol(M)

  covariates <- unique(covariates)
  random_effects <- unique(random_effects)
  needed <- unique(c(x_var, covariates, random_effects))
  missing_vars <- setdiff(needed, colnames(pheno_df))
  if (length(missing_vars) > 0) {
    stop("Missing MaAsLin2 metadata columns in pheno_df: ", paste(missing_vars, collapse = ", "))
  }

  fixed_effects <- unique(c(x_var, covariates))
  if (length(random_effects) > 0) fixed_effects <- setdiff(fixed_effects, random_effects)
  if (!(x_var %in% fixed_effects)) fixed_effects <- c(x_var, fixed_effects)

  meta <- pheno_df[ids, unique(c(fixed_effects, random_effects)), drop = FALSE]

  out_dir <- output_dir
  if (is.null(out_dir) || !nzchar(out_dir)) {
    out_dir <- tempfile(pattern = "zentangler_maaslin2_a_")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fit <- Maaslin2::Maaslin2(
    input_data = as.data.frame(M, check.names = FALSE, stringsAsFactors = FALSE),
    input_metadata = as.data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE),
    output = out_dir,
    fixed_effects = fixed_effects,
    random_effects = if (length(random_effects) > 0) random_effects else NULL,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    standardize = standardize,
    max_significance = 1,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    save_models = FALSE
  )

  res <- fit$results
  if (is.null(res) || nrow(res) == 0) {
    all_res_path <- file.path(out_dir, "all_results.tsv")
    if (file.exists(all_res_path)) {
      res <- utils::read.delim(all_res_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  if (is.null(res) || nrow(res) == 0) stop("MaAsLin2 returned no association results.")

  feature_col <- pick_first_col(res, c("feature", "name"))
  meta_col <- pick_first_col(res, c("metadata", "variable", "covariate"))
  coef_col <- pick_first_col(res, c("coef", "estimate", "value", "beta"))
  p_col <- pick_first_col(res, c("pval", "pvalue", "p"))
  q_col <- pick_first_col(res, c("qval", "qvalue", "q"))

  if (is.null(feature_col) || is.null(meta_col) || is.null(coef_col) || is.null(p_col)) {
    stop("Unexpected MaAsLin2 result schema; could not find feature, metadata, coefficient, and p-value columns.")
  }

  res_x <- res[as.character(res[[meta_col]]) == x_var, , drop = FALSE]

  tab <- data.frame(
    mediator = colnames(M),
    a_hat = NA_real_,
    p_a = NA_real_,
    q_a = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(res_x) > 0) {
    feat <- as.character(res_x[[feature_col]])
    idx <- match(feat, tab$mediator)
    keep <- !is.na(idx)
    if (any(keep)) {
      tab$a_hat[idx[keep]] <- as.numeric(res_x[[coef_col]][keep])
      tab$p_a[idx[keep]] <- as.numeric(res_x[[p_col]][keep])
      if (!is.null(q_col) && q_col %in% colnames(res_x)) {
        tab$q_a[idx[keep]] <- as.numeric(res_x[[q_col]][keep])
      }
    }
  }

  # Recompute A-stage q-values with the requested FDR procedure so LM and
  # MaAsLin2 A-stage outputs are comparable inside Zentangler.
  tab$q_a <- p.adjust(tab$p_a, method = fdr_method)

  ord <- if (rank_by == "abs_a") {
    order(abs(tab$a_hat), decreasing = TRUE, na.last = TRUE)
  } else {
    order(tab$p_a, na.last = TRUE)
  }

  if (screen_method == "none") {
    sis_n_used <- p
    selected <- tab$mediator
  } else {
    if (is.null(sis_n)) sis_n <- hima_sis_default_n(n = n, p = p)
    sis_n_used <- max(1L, min(p, as.integer(sis_n)))
    selected <- tab$mediator[ord][seq_len(sis_n_used)]
  }

  tab$selected <- tab$mediator %in% selected
  tab <- tab[order(tab$p_a, na.last = TRUE), ]
  rownames(tab) <- NULL

  list(
    table = tab,
    selected = selected,
    sis_n = sis_n_used,
    screen_method = screen_method,
    a_stage_model = "maaslin2",
    output_dir = out_dir,
    fixed_effects = fixed_effects,
    random_effects = random_effects
  )
}

