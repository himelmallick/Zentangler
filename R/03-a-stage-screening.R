# -----------------------------------------------------------------------------
# A-stage screening
# -----------------------------------------------------------------------------

hima_sis_default_n <- function(n, p) {
  max(1L, min(p, floor(n / max(log(n), 1.0))))
}

zentangler_a_stage_lm <- function(M, X, C = NULL) {
  # Vectorized A-path regressions M_j ~ X + C for all mediators at once. This is
  # numerically identical to one lm() per mediator but uses a single QR.
  M <- as.matrix(M)
  W <- cbind(`(Intercept)` = 1, X = as.numeric(X))
  if (!is.null(C) && ncol(C) > 0) W <- cbind(W, as.matrix(C))
  n <- nrow(W)
  qr_w <- qr(W)
  rank <- qr_w$rank
  df <- n - rank
  out <- data.frame(
    mediator = colnames(M),
    a_hat = NA_real_,
    se_a = NA_real_,
    p_a = NA_real_,
    stringsAsFactors = FALSE
  )
  if (rank < ncol(W) || df <= 0) {
    warning("A-stage design (intercept, X, covariates) is rank deficient or has no residual df; A-path estimates are NA.")
    return(out)
  }
  coef <- qr.coef(qr_w, M)
  coef <- matrix(coef, nrow = ncol(W))
  resid <- M - W %*% coef
  sigma2 <- colSums(resid^2) / df
  xtx_inv <- chol2inv(qr.R(qr_w))
  x_pos <- 2L
  se <- sqrt(xtx_inv[x_pos, x_pos] * sigma2)
  a_hat <- coef[x_pos, ]
  tstat <- a_hat / se
  out$a_hat <- as.numeric(a_hat)
  out$se_a <- as.numeric(se)
  out$p_a <- as.numeric(2 * stats::pt(abs(tstat), df = df, lower.tail = FALSE))
  out$p_a[!is.finite(out$p_a)] <- NA_real_
  out
}

zentangler_screen_order <- function(tab, M, X, rank_by = c("abs_a", "pvalue"), log_link = FALSE) {
  # SIS ranking. "abs_a" ranks by the standardized A-path coefficient
  # |a_j| * sd(X) / sd(M_j), so features measured on large scales do not crowd
  # out the others. For log-link models (Tweedieverse) a_j is already a log
  # fold change, so the standardized effect is |a_j| * sd(X). "pvalue" ranks by
  # the A-path p-value.
  rank_by <- match.arg(rank_by)
  if (identical(rank_by, "pvalue")) return(order(tab$p_a, na.last = TRUE))
  sd_x <- stats::sd(as.numeric(X))
  if (isTRUE(log_link)) return(order(abs(tab$a_hat) * sd_x, decreasing = TRUE, na.last = TRUE))
  sd_m <- apply(as.matrix(M), 2, stats::sd)
  std_a <- abs(tab$a_hat) * sd_x / sd_m[match(tab$mediator, colnames(M))]
  order(std_a, decreasing = TRUE, na.last = TRUE)
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

  tab <- zentangler_a_stage_lm(M = M, X = X, C = C)
  tab$se_a <- NULL

  tab$q_a <- p.adjust(tab$p_a, method = fdr_method)
  ord <- zentangler_screen_order(tab, M = M, X = X, rank_by = rank_by)

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

zentangler_maaslin2_defaults <- function() {
  # Zentangler's defaults for Maaslin2::Maaslin2(); any of them, and any other
  # Maaslin2 argument, can be overridden through maaslin2_args.
  list(
    normalization = "NONE",
    transform = "NONE",
    analysis_method = "LM",
    standardize = FALSE,
    max_significance = 1,
    min_prevalence = 0,
    min_abundance = -Inf,
    min_variance = 0,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    save_models = FALSE,
    cores = 1
  )
}

hima_screen_mediators_maaslin2 <- function(
  M,
  pheno_df,
  x_var,
  covariates = NULL,
  sis_n = NULL,
  rank_by = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  fdr_method = c("BH", "BY"),
  maaslin2_args = list()
) {
  # MaAsLin2 A-stage: mediator_j ~ X + fixed covariates (+ random effects).
  # maaslin2_args is passed to Maaslin2::Maaslin2() as is, on top of
  # zentangler_maaslin2_defaults(). Zentangler supplies input_data,
  # input_metadata and output; the exposure is always added to fixed_effects.
  rank_by <- match.arg(rank_by)
  screen_method <- match.arg(screen_method)
  fdr_method <- validate_fdr_method(fdr_method)

  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    stop("Package 'Maaslin2' is required when a_stage_model = 'maaslin2'.", call. = FALSE)
  }
  ids <- rownames(M)
  if (is.null(ids) || !all(ids %in% rownames(pheno_df))) {
    stop("Mediator sample IDs must match pheno_df rownames for MaAsLin2 screening.", call. = FALSE)
  }
  n <- nrow(M)
  p <- ncol(M)

  user <- as.list(maaslin2_args %||% list())
  allowed <- setdiff(names(formals(Maaslin2::Maaslin2)), c("input_data", "input_metadata"))
  unknown <- setdiff(names(user), c(allowed, "output_dir"))
  if (length(unknown) > 0L) {
    stop("Unknown Maaslin2 argument(s) in maaslin2_args: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (!is.null(user$output_dir) && is.null(user$output)) user$output <- user$output_dir
  user$output_dir <- NULL
  random_effects <- unique(as.character(user$random_effects %||% character(0)))
  fixed_effects <- unique(c(x_var, as.character(user$fixed_effects %||% covariates)))
  fixed_effects <- c(x_var, setdiff(fixed_effects, c(x_var, random_effects)))
  meta_cols <- unique(c(fixed_effects, random_effects, as.character(user$reference %||% character(0))))
  meta_cols <- intersect(meta_cols, colnames(pheno_df))
  missing_vars <- setdiff(c(fixed_effects, random_effects), colnames(pheno_df))
  if (length(missing_vars) > 0L) {
    stop("Missing MaAsLin2 metadata columns: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  out_dir <- user$output
  if (is.null(out_dir) || !nzchar(out_dir)) out_dir <- tempfile(pattern = "zentangler_maaslin2_a_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  call_args <- utils::modifyList(zentangler_maaslin2_defaults(), user)
  call_args$input_data <- as.data.frame(M, check.names = FALSE, stringsAsFactors = FALSE)
  call_args$input_metadata <- as.data.frame(pheno_df[ids, meta_cols, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  call_args$output <- out_dir
  call_args$fixed_effects <- fixed_effects
  call_args$random_effects <- if (length(random_effects) > 0L) random_effects else NULL
  # MaAsLin2 logs every step to the console; keep it quiet (its log file is
  # still written to the output folder).
  fit <- NULL
  utils::capture.output(
    fit <- suppressMessages(suppressWarnings(do.call(Maaslin2::Maaslin2, call_args))),
    type = "output"
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
    # MaAsLin2 may rewrite feature names with make.names(); map those back.
    miss <- is.na(idx)
    if (any(miss)) idx[miss] <- match(feat[miss], make.names(tab$mediator))
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
  n_missing_a <- sum(!is.finite(tab$a_hat))
  if (n_missing_a > 0) {
    warning(
      "MaAsLin2 A-stage returned no estimate for ", n_missing_a, " of ", p,
      " features; they are treated as untested (p_a = NA) and cannot be screened in."
    )
  }
  tab$q_a <- p.adjust(tab$p_a, method = fdr_method)

  ord <- zentangler_screen_order(tab, M = M, X = pheno_df[ids, x_var], rank_by = rank_by)

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
    random_effects = random_effects,
    maaslin2_args = call_args[setdiff(names(call_args), c("input_data", "input_metadata"))]
  )
}


# -----------------------------------------------------------------------------
# Tweedieverse A-stage
#
# Per-feature Tweedie GLMs, mediator_j ~ X + covariates (+ random effects), on
# the untransformed values of the view. Tweedieverse is technology-agnostic:
# counts, relative abundances and continuous non-negative data (with or without
# zeros) use the compound Poisson range 1 < p < 2 with a log link; data with
# negative values resolve to p = 0 (Gaussian, identity link). Only screening
# and p_a come from this model; the B-path and the causal-effect stage use the
# (transformed) feature scale, so effects are always re-estimated with linear
# mediator models.
# -----------------------------------------------------------------------------

zentangler_tweedieverse_defaults <- function() {
  # Zentangler's defaults for Tweedieverse::Tweedieverse(); any of them, and any
  # other Tweedieverse argument, can be overridden through tweedieverse_args.
  list(
    output = NULL,
    abd_threshold = -Inf,
    prev_threshold = 0,
    var_threshold = 0,
    entropy_threshold = 0,
    max_significance = 1,
    standardize = FALSE,
    median_comparison = FALSE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
}

hima_screen_mediators_tweedieverse <- function(
  M,
  pheno_df,
  x_var,
  covariates = NULL,
  sis_n = NULL,
  rank_by = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  fdr_method = c("BH", "BY"),
  tweedieverse_args = list(),
  view = NULL
) {
  rank_by <- match.arg(rank_by)
  screen_method <- match.arg(screen_method)
  fdr_method <- validate_fdr_method(fdr_method)
  if (!requireNamespace("Tweedieverse", quietly = TRUE)) {
    stop("Package 'Tweedieverse' is required when a_stage_model = 'tweedieverse'. ",
         "Install it with remotes::install_github('himelmallick/Tweedieverse').", call. = FALSE)
  }
  M <- as.matrix(M)
  ids <- rownames(M)
  if (is.null(ids) || !all(ids %in% rownames(pheno_df))) {
    stop("Mediator sample IDs must match pheno_df rownames for Tweedieverse screening.", call. = FALSE)
  }
  lab <- if (is.null(view)) "" else paste0(" (view '", view, "')")
  n <- nrow(M)
  p <- ncol(M)

  user <- as.list(tweedieverse_args %||% list())
  allowed <- setdiff(names(formals(Tweedieverse::Tweedieverse)), c("input_features", "input_metadata", "assay_name"))
  unknown <- setdiff(names(user), allowed)
  if (length(unknown) > 0L) {
    stop("Unknown Tweedieverse argument(s) in tweedieverse_args: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  random_effects <- unique(as.character(user$random_effects %||% character(0)))
  fixed_effects <- unique(c(x_var, as.character(user$fixed_effects %||% covariates)))
  fixed_effects <- c(x_var, setdiff(fixed_effects, c(x_var, random_effects)))
  extra_cols <- as.character(c(user$scale_factor %||% character(0)))
  missing_vars <- setdiff(c(fixed_effects, random_effects, extra_cols), colnames(pheno_df))
  if (length(missing_vars) > 0L) {
    stop("Missing Tweedieverse metadata columns: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  meta_cols <- unique(c(fixed_effects, random_effects, extra_cols))
  if (!is.null(user$output) && nzchar(user$output)) dir.create(user$output, recursive = TRUE, showWarnings = FALSE)

  call_args <- utils::modifyList(zentangler_tweedieverse_defaults(), user, keep.null = TRUE)
  meta <- as.data.frame(pheno_df[ids, meta_cols, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  call_args$input_features <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(M)),
    colData = meta
  )
  call_args$input_metadata <- NULL
  call_args$assay_name <- "counts"
  call_args$fixed_effects <- fixed_effects
  call_args$random_effects <- if (length(random_effects) > 0L) random_effects else NULL
  call_args$max_significance <- 1  # every feature is needed for screening

  # Tweedieverse logs every feature fit to the console and prints the errors of
  # failed fits through try(); keep it quiet.
  res <- NULL
  null_con <- file(nullfile(), open = "w")
  old_try <- options(try.outFile = null_con)
  on.exit({ options(old_try); close(null_con) }, add = TRUE)
  utils::capture.output(
    res <- suppressMessages(suppressWarnings(do.call(Tweedieverse::Tweedieverse, call_args))),
    type = "output"
  )
  if (is.null(res) || !is.data.frame(res) || nrow(res) == 0L) stop("Tweedieverse returned no association results", lab, ".", call. = FALSE)

  res_x <- res[as.character(res$metadata) == x_var, , drop = FALSE]
  tab <- data.frame(mediator = colnames(M), a_hat = NA_real_, p_a = NA_real_, q_a = NA_real_, stringsAsFactors = FALSE)
  tweedie_index <- setNames(rep(NA_real_, p), colnames(M))
  if (nrow(res_x) > 0L) {
    feat <- as.character(res_x$feature)
    idx <- match(feat, tab$mediator)
    miss <- is.na(idx)
    if (any(miss)) idx[miss] <- match(feat[miss], make.names(tab$mediator))
    keep <- !is.na(idx)
    tab$a_hat[idx[keep]] <- as.numeric(res_x$coef[keep])
    tab$p_a[idx[keep]] <- as.numeric(res_x$pval[keep])
    if ("tweedie.index" %in% colnames(res_x)) tweedie_index[idx[keep]] <- as.numeric(res_x$tweedie.index[keep])
  }
  bad <- !is.finite(tab$a_hat) | !is.finite(tab$p_a)
  tab$a_hat[bad] <- NA_real_
  tab$p_a[bad] <- NA_real_
  if (any(bad)) {
    warning("Tweedieverse A-stage", lab, " returned no estimate for ", sum(bad), " of ", p,
            " features (fit failures or filtering); they are treated as untested (p_a = NA) and cannot be screened in.", call. = FALSE)
  }
  tab$q_a <- p.adjust(tab$p_a, method = fdr_method)

  # With p = 0 (data with negative values, or tweedie_p = 0) Tweedieverse uses
  # the identity link, so a_j is then on the feature's own scale.
  p_used <- tweedie_index[is.finite(tweedie_index)]
  gaussian_p0 <- length(p_used) > 0L && all(p_used == 0)
  log_link <- !gaussian_p0 && identical(as.character(call_args$link %||% "log"), "log")
  ord <- zentangler_screen_order(tab, M = M, X = pheno_df[ids, x_var], rank_by = rank_by, log_link = log_link)
  if (screen_method == "none") {
    sis_n_used <- p
    selected <- tab$mediator[is.finite(tab$p_a)]
  } else {
    if (is.null(sis_n)) sis_n <- hima_sis_default_n(n = n, p = p)
    sis_n_used <- max(0L, min(sum(is.finite(tab$p_a)), as.integer(sis_n)))
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
    a_stage_model = "tweedieverse",
    a_scale = if (log_link) "log fold change per unit X" else if (gaussian_p0) "identity (Gaussian, p = 0)" else paste0(call_args$link, " link"),
    tweedie_index = tweedie_index,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    tweedieverse_args = call_args[setdiff(names(call_args), c("input_features", "input_metadata"))]
  )
}

# Raw (untransformed) values travel inside the blocks list under a reserved
# name so that bootstrap resampling keeps them aligned with the transformed
# block; the engines split them off before sample alignment.
zentangler_raw_key <- function(view) if (length(view)) paste0(".raw::", view) else character(0)

zentangler_attach_raw_blocks <- function(blocks, raw_blocks, a_stage) {
  views <- names(a_stage$model)[a_stage$model == "tweedieverse"]
  for (v in intersect(views, names(raw_blocks))) blocks[[zentangler_raw_key(v)]] <- raw_blocks[[v]]
  blocks
}

zentangler_split_raw_blocks <- function(blocks) {
  is_raw <- startsWith(names(blocks), ".raw::")
  raw <- blocks[is_raw]
  names(raw) <- sub("^\\.raw::", "", names(raw))
  list(blocks = blocks[!is_raw], raw = raw)
}

zentangler_align_raw_blocks <- function(raw, B) {
  # Subset raw values to the aligned samples and retained features of B.
  out <- list()
  for (v in intersect(names(raw), names(B))) {
    R <- as.matrix(raw[[v]])
    ids <- rownames(B[[v]])
    cols <- colnames(B[[v]])
    if (!all(ids %in% rownames(R)) || !all(cols %in% colnames(R))) {
      stop("Raw values for view '", v, "' do not match the transformed block.", call. = FALSE)
    }
    R <- R[ids, cols, drop = FALSE]
    R[!is.finite(R)] <- NA_real_
    out[[v]] <- R
  }
  out
}

# -----------------------------------------------------------------------------
# Per-view A-stage models
#
# a_stage_model is either one model for every view ("lm", "maaslin2" or
# "tweedieverse") or a vector named by view, e.g.
# c(species = "tweedieverse", metabolites = "lm"); views not named fall back to
# "lm". maaslin2_args / tweedieverse_args are each either one list of
# arguments shared by every view using that model, or a list named by view
# whose elements are such lists. Names may be MAE experiment names or the
# sanitized view names Zentangler uses internally. This only affects the
# A-path (exposure -> mediator) stage.
# -----------------------------------------------------------------------------

zentangler_view_name_map <- function(views, view_map = NULL) {
  # Map experiment names and sanitized names to sanitized view names.
  map <- setNames(views, views)
  if (!is.null(view_map) && nrow(view_map) > 0L) {
    map[as.character(view_map$experiment_name)] <- as.character(view_map$view)
  }
  map
}

zentangler_a_stage_models <- function() c("lm", "maaslin2", "tweedieverse")

zentangler_resolve_method_args <- function(args, views, name_map, to_view, label) {
  # One list shared by every view, or a list named by view of such lists.
  args <- as.list(args %||% list())
  per_view <- length(args) > 0L && !is.null(names(args)) &&
    all(names(args) %in% names(name_map)) && all(vapply(args, is.list, logical(1)))
  out <- setNames(rep(list(list()), length(views)), views)
  if (per_view) {
    for (nm in names(args)) out[[to_view(nm, label)]] <- args[[nm]]
  } else {
    for (v in views) out[[v]] <- args
  }
  out
}

zentangler_resolve_a_stage <- function(a_stage_model, maaslin2_args, views, view_map = NULL, tweedieverse_args = list()) {
  name_map <- zentangler_view_name_map(views, view_map)
  to_view <- function(nm, label) {
    bad <- setdiff(nm, names(name_map))
    if (length(bad) > 0L) stop(label, " names views that are not analysed: ", paste(bad, collapse = ", "), call. = FALSE)
    unname(name_map[nm])
  }
  a_stage_model <- a_stage_model %||% "lm"
  models <- setNames(rep("lm", length(views)), views)
  if (is.null(names(a_stage_model))) {
    if (length(a_stage_model) != 1L) {
      stop("a_stage_model must be a single model or a vector named by view, e.g. c(species = 'tweedieverse', metabolites = 'lm').", call. = FALSE)
    }
    models[] <- as.character(a_stage_model)
  } else {
    models[to_view(names(a_stage_model), "a_stage_model")] <- as.character(a_stage_model)
  }
  bad <- setdiff(unique(models), zentangler_a_stage_models())
  if (length(bad) > 0L) {
    stop("Unknown a_stage_model: ", paste(bad, collapse = ", "), ". Use ",
         paste0("'", zentangler_a_stage_models(), "'", collapse = ", "), ".", call. = FALSE)
  }

  maaslin2_by_view <- zentangler_resolve_method_args(maaslin2_args, views, name_map, to_view, "maaslin2_args")
  tweedie_by_view <- zentangler_resolve_method_args(tweedieverse_args, views, name_map, to_view, "tweedieverse_args")
  re_of <- function(v) {
    switch(models[[v]],
      maaslin2 = as.character(maaslin2_by_view[[v]]$random_effects %||% character(0)),
      tweedieverse = as.character(c(tweedie_by_view[[v]]$random_effects, tweedie_by_view[[v]]$scale_factor) %||% character(0)),
      character(0)
    )
  }
  random_effects <- unique(unlist(lapply(views, re_of)))
  list(model = models, maaslin2_args = maaslin2_by_view, tweedieverse_args = tweedie_by_view, random_effects = random_effects)
}

# Point every MaAsLin2 / Tweedieverse output folder of a bootstrap replicate
# to its own sub-folder.
zentangler_a_stage_bootstrap_dirs <- function(a_stage, b) {
  if (is.null(a_stage)) return(a_stage)
  for (slot in c("maaslin2_args", "tweedieverse_args")) {
    for (v in names(a_stage[[slot]])) {
      ma <- a_stage[[slot]][[v]]
      base_out <- ma$output %||% ma$output_dir
      if (!is.null(base_out) && nzchar(base_out)) {
        ma$output_dir <- NULL
        ma$output <- file.path(base_out, paste0("bootstrap_", b))
        a_stage[[slot]][[v]] <- ma
      }
    }
  }
  a_stage
}

zentangler_a_stage_screen <- function(view, M, X, C, pheno_df, x_var, covariates, model, maaslin2_args,
                                      sis_n, rank_by, screen_method, fdr_method,
                                      tweedieverse_args = list(), M_raw = NULL) {
  if (identical(model, "tweedieverse")) {
    ta <- tweedieverse_args %||% list()
    if (!is.null(ta$output) && nzchar(ta$output)) ta$output <- file.path(ta$output, paste0("a_stage_", view))
    out <- hima_screen_mediators_tweedieverse(
      M = M_raw %||% M, pheno_df = pheno_df, x_var = x_var, covariates = covariates,
      sis_n = sis_n, rank_by = rank_by, screen_method = screen_method,
      fdr_method = fdr_method, tweedieverse_args = ta, view = view
    )
    # Keep the table aligned with the (transformed) block used downstream.
    out$table <- out$table[out$table$mediator %in% colnames(M), , drop = FALSE]
    out$selected <- intersect(out$selected, colnames(M))
    return(out)
  }
  if (identical(model, "maaslin2")) {
    ma <- maaslin2_args
    base_out <- ma$output %||% ma$output_dir
    ma$output_dir <- NULL
    ma$output <- if (!is.null(base_out) && nzchar(base_out)) file.path(base_out, paste0("a_stage_", view)) else NULL
    return(hima_screen_mediators_maaslin2(
      M = M, pheno_df = pheno_df, x_var = x_var, covariates = covariates,
      sis_n = sis_n, rank_by = rank_by, screen_method = screen_method,
      fdr_method = fdr_method, maaslin2_args = ma
    ))
  }
  out <- hima_screen_mediators(M, X, C = C, sis_n = sis_n, rank_by = rank_by, screen_method = screen_method, fdr_method = fdr_method)
  out$a_stage_model <- "lm"
  out
}
