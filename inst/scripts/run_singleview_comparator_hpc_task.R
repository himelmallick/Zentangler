#!/usr/bin/env Rscript
# Run one row of the Zentangler simulation grid through single-view comparator methods.
#
# Comparator logic:
# - Generate the same SIMMBA/InterSIM MAE data as the Zentangler grid row.
# - For each replicate and each omics view, run the single-view method separately.
# - Stack mediator-level results across views.
# - Apply one final FDR correction across stacked mediation-effect p-values.
# - Score recovery against the known SIMMBA mediator truth.
#
# Supported methods:
#   hima        : HIMA package, single-view high-dimensional mediation
#   hima2       : HIMA2-style call from the HIMA package when available
#   multimedia  : multimedia package pathwise indirect effects; bootstrap p-values optional
#
# Example:
#   Rscript inst/scripts/run_singleview_comparator_hpc_task.R \
#     --grid intersim_zentangler_grid.csv \
#     --task-id 1 \
#     --out-dir singleview_comparator_results \
#     --package-dir /home/fs01/naa4050/Zentangler \
#     --methods hima,hima2,multimedia

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!grepl("^--", key)) {
      i <- i + 1L
      next
    }
    nm <- sub("^--", "", key)
    val <- if (i + 1L <= length(args) && !grepl("^--", args[[i + 1L]])) args[[i + 1L]] else TRUE
    out[[nm]] <- val
    i <- i + if (isTRUE(val)) 1L else 2L
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(default)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

as_chr_or_null <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(NULL)
  as.character(x)
}

as_int_or_null <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(NULL)
  as.integer(x)
}

as_num_or_null <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(NULL)
  as.numeric(x)
}

as_chr_vec <- function(x, default = character(0)) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(default)
  vals <- unlist(strsplit(as.character(x), "[,;|[:space:]]+"))
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0) default else vals
}

as_num_vec <- function(x, default = 0.25) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(default)
  vals <- unlist(strsplit(as.character(x), "[,;|[:space:]]+"))
  vals <- vals[nzchar(vals)]
  out <- suppressWarnings(as.numeric(vals))
  out <- out[is.finite(out)]
  if (length(out) == 0) default else out
}

row_value <- function(row, name, default = NULL) {
  if (!(name %in% colnames(row))) return(default)
  val <- row[[name]][1]
  if (length(val) == 0 || is.na(val) || !nzchar(as.character(val))) return(default)
  val
}

safe_tag <- function(x) {
  x <- ifelse(is.na(x) | !nzchar(as.character(x)), "NA", as.character(x))
  gsub("[^A-Za-z0-9]+", "", x)
}

load_zentangler_sources <- function(package_dir = NULL) {
  if (!is.null(package_dir) && nzchar(package_dir) && file.exists(file.path(package_dir, "DESCRIPTION"))) {
    r_files <- sort(list.files(file.path(package_dir, "R"), pattern = "[.]R$", full.names = TRUE))
    for (f in r_files) source(f)
    return(invisible(TRUE))
  }
  suppressPackageStartupMessages(library(Zentangler))
  invisible(TRUE)
}

zt_get <- function(name) {
  if (exists(name, mode = "function", inherits = TRUE)) {
    return(get(name, mode = "function", inherits = TRUE))
  }
  if (requireNamespace("Zentangler", quietly = TRUE) &&
      exists(name, envir = asNamespace("Zentangler"), inherits = FALSE)) {
    return(get(name, envir = asNamespace("Zentangler")))
  }
  stop("Could not find Zentangler helper function: ", name, call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

pick_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) hit[1] else NA_character_
}

standardize_method_table <- function(x, all_mediators, view, method, rep_id, p_default = 1) {
  all_mediators <- as.character(all_mediators)
  out <- data.frame(
    rep = rep_id,
    method = method,
    omics = view,
    mediator = all_mediators,
    a = NA_real_,
    b = 0,
    score = 0,
    abs_score = 0,
    p_primary = p_default,
    p_source = NA_character_,
    stringsAsFactors = FALSE
  )

  if (is.null(x)) return(out)
  x <- try(as.data.frame(x, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(x, "try-error") || nrow(x) == 0) return(out)

  id_col <- pick_col(x, c("mediator", "Mediator", "ID", "id", "feature", "Feature", "name", "Name"))
  med <- if (!is.na(id_col)) as.character(x[[id_col]]) else rownames(x)
  med[is.na(med) | !nzchar(med)] <- rownames(x)[is.na(med) | !nzchar(med)]

  a_col <- pick_col(x, c("a", "alpha", "alpha_hat", "alpha.hat", "ahat", "a_hat", "effect.a", "beta.xm"))
  b_col <- pick_col(x, c("b", "beta", "beta_hat", "beta.hat", "bhat", "b_hat", "effect.b", "beta.my"))
  score_col <- pick_col(x, c("score", "indirect", "indirect_effect", "indirect.effect", "alpha_beta", "alpha*beta", "ab", "IDE", "NIE"))
  p_col <- pick_col(x, c("p_primary", "p.value", "p_value", "pvalue", "p", "pval", "p.val", "p.joint", "p_joint", "pmax", "p.max", "raw.p", "p_raw"))
  q_col <- pick_col(x, c("q", "q.value", "q_value", "qvalue", "FDR", "BH.FDR", "adj.p", "adjusted.p", "padj", "qval"))

  tmp <- data.frame(mediator = med, stringsAsFactors = FALSE)
  tmp$a <- if (!is.na(a_col)) suppressWarnings(as.numeric(x[[a_col]])) else NA_real_
  tmp$b <- if (!is.na(b_col)) suppressWarnings(as.numeric(x[[b_col]])) else NA_real_
  tmp$score <- if (!is.na(score_col)) suppressWarnings(as.numeric(x[[score_col]])) else tmp$a * tmp$b
  tmp$p_primary <- if (!is.na(p_col)) suppressWarnings(as.numeric(x[[p_col]])) else NA_real_
  if (all(!is.finite(tmp$p_primary)) && !is.na(q_col)) {
    tmp$p_primary <- suppressWarnings(as.numeric(x[[q_col]]))
    p_src <- q_col
  } else {
    p_src <- p_col
  }

  tmp <- tmp[!duplicated(tmp$mediator), , drop = FALSE]
  idx <- match(out$mediator, tmp$mediator)
  ok <- !is.na(idx)

  out$a[ok] <- tmp$a[idx[ok]]
  out$b[ok] <- tmp$b[idx[ok]]
  out$score[ok] <- tmp$score[idx[ok]]
  out$p_primary[ok] <- tmp$p_primary[idx[ok]]
  out$p_source[ok] <- p_src

  out$b[!is.finite(out$b)] <- 0
  out$score[!is.finite(out$score)] <- 0
  out$abs_score <- abs(out$score)
  out$p_primary[!is.finite(out$p_primary)] <- p_default
  out$p_primary <- pmin(1, pmax(0, out$p_primary))
  out
}

run_hima_singleview <- function(M, pheno, x_var = "A", y_var = "Y", y_family = "gaussian", penalty = "lasso") {
  if (!requireNamespace("HIMA", quietly = TRUE)) {
    stop("Package 'HIMA' is required for method = 'hima'.")
  }
  data_pheno <- as.data.frame(pheno, stringsAsFactors = FALSE)
  data_M <- as.data.frame(M, check.names = FALSE, stringsAsFactors = FALSE)
  form <- stats::as.formula(paste(y_var, "~", x_var))

  f <- get("hima", envir = asNamespace("HIMA"))
  argn <- names(formals(f))

  call_args <- list()
  if ("formula" %in% argn) call_args$formula <- form
  if ("data.pheno" %in% argn) call_args$data.pheno <- data_pheno
  if ("data.M" %in% argn) call_args$data.M <- data_M
  if ("X" %in% argn) call_args$X <- data_pheno[[x_var]]
  if ("Y" %in% argn) call_args$Y <- data_pheno[[y_var]]
  if ("M" %in% argn) call_args$M <- data_M
  if ("COV.XM" %in% argn) call_args$COV.XM <- NULL
  if ("COV.MY" %in% argn) call_args$COV.MY <- NULL
  if ("mediator.type" %in% argn) call_args$mediator.type <- "gaussian"
  if ("M.family" %in% argn) call_args$M.family <- "gaussian"
  if ("outcome.family" %in% argn) call_args$outcome.family <- if (identical(y_family, "binomial")) "binomial" else "gaussian"
  if ("Y.family" %in% argn) call_args$Y.family <- if (identical(y_family, "binomial")) "binomial" else "gaussian"
  if ("penalty" %in% argn) call_args$penalty <- penalty
  if ("topN" %in% argn) call_args$topN <- ncol(data_M)
  if ("sigcut" %in% argn) call_args$sigcut <- 1
  if ("Sigcut" %in% argn) call_args$Sigcut <- 1
  if ("scale" %in% argn) call_args$scale <- FALSE
  if ("verbose" %in% argn) call_args$verbose <- FALSE
  if ("parallel" %in% argn) call_args$parallel <- FALSE

  do.call(f, call_args)
}

run_hima2_singleview <- function(M, pheno, x_var = "A", y_var = "Y", y_family = "gaussian", penalty = "DBlasso") {
  if (!requireNamespace("HIMA", quietly = TRUE)) {
    stop("Package 'HIMA' is required for method = 'hima2'.")
  }
  ns <- asNamespace("HIMA")
  if (!exists("hima2", envir = ns, inherits = FALSE)) {
    stop("Function HIMA::hima2 is not available in the installed HIMA package.")
  }

  data_pheno <- as.data.frame(pheno, stringsAsFactors = FALSE)
  data_M <- as.data.frame(M, check.names = FALSE, stringsAsFactors = FALSE)
  form <- stats::as.formula(paste(y_var, "~", x_var))

  f <- get("hima2", envir = ns)
  argn <- names(formals(f))

  call_args <- list()
  if ("formula" %in% argn) call_args$formula <- form
  if ("data.pheno" %in% argn) call_args$data.pheno <- data_pheno
  if ("data.M" %in% argn) call_args$data.M <- data_M
  if ("X" %in% argn) call_args$X <- data_pheno[[x_var]]
  if ("Y" %in% argn) call_args$Y <- data_pheno[[y_var]]
  if ("M" %in% argn) call_args$M <- data_M
  if ("mediator.family" %in% argn) call_args$mediator.family <- "gaussian"
  if ("M.family" %in% argn) call_args$M.family <- "gaussian"
  if ("outcome.family" %in% argn) call_args$outcome.family <- if (identical(y_family, "binomial")) "binomial" else "gaussian"
  if ("Y.family" %in% argn) call_args$Y.family <- if (identical(y_family, "binomial")) "binomial" else "gaussian"
  if ("penalty" %in% argn) call_args$penalty <- penalty
  if ("topN" %in% argn) call_args$topN <- ncol(data_M)
  if ("sigcut" %in% argn) call_args$sigcut <- 1
  if ("Sigcut" %in% argn) call_args$Sigcut <- 1
  if ("scale" %in% argn) call_args$scale <- FALSE
  if ("verbose" %in% argn) call_args$verbose <- FALSE
  if ("parallel" %in% argn) call_args$parallel <- FALSE

  do.call(f, call_args)
}

multimedia_indirect_from_object <- function(x) {
  xdf <- try(as.data.frame(x, stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(xdf, "try-error") && nrow(xdf) > 0) return(xdf)
  if (is.list(x)) {
    for (nm in names(x)) {
      xdf <- try(as.data.frame(x[[nm]], stringsAsFactors = FALSE), silent = TRUE)
      if (!inherits(xdf, "try-error") && nrow(xdf) > 0) return(xdf)
    }
  }
  data.frame()
}

multimedia_indirect_pathwise_safe <- function(fit, exper) {
  f <- get("indirect_pathwise", envir = asNamespace("multimedia"))
  argn <- names(formals(f))
  call_args <- list(fit, exper)
  if ("progress" %in% argn) call_args$progress <- FALSE
  do.call(f, call_args)
}

run_multimedia_singleview <- function(M,
                                      pheno,
                                      x_var = "A",
                                      y_var = "Y",
                                      y_family = "gaussian",
                                      model = c("lm", "rf"),
                                      bootstrap_repeats = 0L,
                                      rf_trees = 500L,
                                      seed = 1L) {
  if (!requireNamespace("multimedia", quietly = TRUE)) {
    stop("Package 'multimedia' is required for method = 'multimedia'.")
  }
  model <- match.arg(model)
  if (identical(model, "rf") && !requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required for multimedia random forest models.")
  }
  if (!requireNamespace("tidyselect", quietly = TRUE)) {
    stop("Package 'tidyselect' is required by this multimedia wrapper.")
  }

  med_names <- colnames(M)
  df <- data.frame(Y = as.numeric(pheno[[y_var]]), A = as.numeric(pheno[[x_var]]), M, check.names = FALSE)

  exper <- multimedia::mediation_data(df, "Y", "A", tidyselect::all_of(med_names))
  base_model <- if (identical(model, "rf")) {
    f <- get("rf_model", envir = asNamespace("multimedia"))
    argn <- names(formals(f))
    if ("num.trees" %in% argn) do.call(f, list(num.trees = as.integer(rf_trees))) else f()
  } else {
    multimedia::lm_model()
  }

  fit <- multimedia::multimedia(exper, base_model)
  if (exists("estimate", envir = asNamespace("multimedia"), inherits = FALSE)) {
    fit <- multimedia::estimate(fit, exper)
  }
  ind <- multimedia_indirect_pathwise_safe(fit, exper)
  tab <- multimedia_indirect_from_object(ind)
  if (nrow(tab) == 0) return(data.frame())

  med_col <- pick_col(tab, c("mediator", "Mediator", "feature", "Feature", "term", "variable", "name"))
  score_col <- pick_col(tab, c("indirect", "indirect_effect", "indirect.effect", "effect", "estimate", "value"))
  if (is.na(med_col) || is.na(score_col)) {
    stop("Could not identify mediator/effect columns in multimedia::indirect_pathwise output.")
  }

  out <- data.frame(
    mediator = as.character(tab[[med_col]]),
    score = suppressWarnings(as.numeric(tab[[score_col]])),
    stringsAsFactors = FALSE
  )
  out <- stats::aggregate(score ~ mediator, data = out, FUN = mean, na.rm = TRUE)
  out$a <- NA_real_
  out$b <- out$score
  out$p_primary <- NA_real_

  bootstrap_repeats <- max(0L, as.integer(bootstrap_repeats))
  if (bootstrap_repeats > 0L) {
    # multimedia exposes bootstrap inference through repeated model refits on
    # resampled mediation_data objects. We implement that explicitly here so
    # that the output is a per-mediator p-value vector compatible with the
    # cross-view FDR comparison used for HIMA/HIMA2/Zentangler.
    set.seed(seed)
    boot_scores <- matrix(
      NA_real_,
      nrow = bootstrap_repeats,
      ncol = length(out$mediator),
      dimnames = list(NULL, out$mediator)
    )

    for (bb in seq_len(bootstrap_repeats)) {
      idx <- sample(seq_len(nrow(df)), size = nrow(df), replace = TRUE)
      df_b <- df[idx, , drop = FALSE]

      boot_tab <- try({
        exper_b <- multimedia::mediation_data(df_b, "Y", "A", tidyselect::all_of(med_names))
        base_model_b <- if (identical(model, "rf")) {
          f <- get("rf_model", envir = asNamespace("multimedia"))
          argn <- names(formals(f))
          if ("num.trees" %in% argn) do.call(f, list(num.trees = as.integer(rf_trees))) else f()
        } else {
          multimedia::lm_model()
        }
        fit_b <- multimedia::multimedia(exper_b, base_model_b)
        if (exists("estimate", envir = asNamespace("multimedia"), inherits = FALSE)) {
          fit_b <- multimedia::estimate(fit_b, exper_b)
        }
        ind_b <- multimedia_indirect_pathwise_safe(fit_b, exper_b)
        multimedia_indirect_from_object(ind_b)
      }, silent = TRUE)

      if (inherits(boot_tab, "try-error") || nrow(boot_tab) == 0) next
      b_med_col <- pick_col(boot_tab, c("mediator", "Mediator", "feature", "Feature", "term", "variable", "name"))
      b_score_col <- pick_col(boot_tab, c("indirect", "indirect_effect", "indirect.effect", "effect", "estimate", "value"))
      if (is.na(b_med_col) || is.na(b_score_col)) next

      tmp <- data.frame(
        mediator = as.character(boot_tab[[b_med_col]]),
        score = suppressWarnings(as.numeric(boot_tab[[b_score_col]])),
        stringsAsFactors = FALSE
      )
      tmp <- stats::aggregate(score ~ mediator, data = tmp, FUN = mean, na.rm = TRUE)
      hit <- match(tmp$mediator, colnames(boot_scores))
      ok <- !is.na(hit)
      boot_scores[bb, hit[ok]] <- tmp$score[ok]
    }

    pvals <- vapply(out$mediator, function(med) {
      obs <- out$score[out$mediator == med][1]
      vals <- boot_scores[, med]
      vals <- vals[is.finite(vals)]
      if (!is.finite(obs) || length(vals) < 2L) return(NA_real_)
      centered <- vals - obs
      (1 + sum(abs(centered) >= abs(obs), na.rm = TRUE)) / (length(centered) + 1)
    }, numeric(1))
    if (all(!is.finite(pvals))) {
      stop("multimedia bootstrap requested, but no per-mediator bootstrap p-values could be computed.")
    }
    out$p_primary <- pvals
  }

  out
}

apply_comparator_fdr <- function(tab, fdr_method = c("BH", "BY")) {
  fdr_method <- match.arg(fdr_method)
  if (is.null(tab) || nrow(tab) == 0) return(tab)
  tab$q_primary_global <- NA_real_
  ok <- is.finite(tab$p_primary)
  tab$q_primary_global[ok] <- p.adjust(tab$p_primary[ok], method = fdr_method)
  tab$q_primary_within_view <- NA_real_
  for (idx in split(seq_len(nrow(tab)), paste(tab$method, tab$omics, sep = "::"))) {
    ok_i <- is.finite(tab$p_primary[idx])
    if (any(ok_i)) tab$q_primary_within_view[idx[ok_i]] <- p.adjust(tab$p_primary[idx[ok_i]], method = fdr_method)
  }
  tab$q_primary <- tab$q_primary_global
  tab
}

score_truth_recovery <- function(tab, truth_key, rep, method, q_threshold, q_col, fdr_scope, top_n) {
  if (is.null(tab) || nrow(tab) == 0) {
    return(data.frame(
      rep = rep,
      method = method,
      fdr_scope = fdr_scope,
      q_col = q_col,
      q_threshold = q_threshold,
      n_active = 0,
      true_active = 0,
      false_active = 0,
      n_true = sum(truth_key$is_true, na.rm = TRUE),
      precision = NA_real_,
      recall = 0,
      fdr = NA_real_,
      top50_true = 0,
      top50_precision = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  key <- paste(as.character(tab$omics), as.character(tab$mediator), sep = "::")
  truth_map <- stats::setNames(truth_key$is_true, truth_key$key)
  is_true <- as.logical(truth_map[key])
  is_true[is.na(is_true)] <- FALSE

  active <- is.finite(tab[[q_col]]) & tab[[q_col]] <= q_threshold & is.finite(tab$score) & tab$score != 0
  n_active <- sum(active, na.rm = TRUE)
  true_active <- sum(active & is_true, na.rm = TRUE)
  false_active <- sum(active & !is_true, na.rm = TRUE)
  n_true <- sum(truth_key$is_true, na.rm = TRUE)

  ord <- order(tab$abs_score, decreasing = TRUE, na.last = NA)
  top_idx <- head(ord, min(as.integer(top_n), length(ord)))
  top_true <- sum(is_true[top_idx], na.rm = TRUE)
  top_precision <- if (length(top_idx) > 0) top_true / length(top_idx) else NA_real_

  data.frame(
    rep = rep,
    method = method,
    fdr_scope = fdr_scope,
    q_col = q_col,
    q_threshold = q_threshold,
    n_active = n_active,
    true_active = true_active,
    false_active = false_active,
    n_true = n_true,
    precision = if (n_active > 0) true_active / n_active else NA_real_,
    recall = if (n_true > 0) true_active / n_true else NA_real_,
    fdr = if (n_active > 0) false_active / n_active else NA_real_,
    top50_true = top_true,
    top50_precision = top_precision,
    stringsAsFactors = FALSE
  )
}

summarize_recovery <- function(detail) {
  if (is.null(detail) || nrow(detail) == 0) return(data.frame())
  group_cols <- intersect(c("method", "fdr_scope", "q_col", "q_threshold"), colnames(detail))
  key <- interaction(detail[, group_cols, drop = FALSE], drop = TRUE, sep = "\r")
  rows <- lapply(split(seq_len(nrow(detail)), key), function(idx) {
    d <- detail[idx, , drop = FALSE]
    g <- d[1L, group_cols, drop = FALSE]
    rownames(g) <- NULL
    data.frame(
      g,
      n_active = mean(d$n_active, na.rm = TRUE),
      true_active = mean(d$true_active, na.rm = TRUE),
      false_active = mean(d$false_active, na.rm = TRUE),
      precision = mean(d$precision, na.rm = TRUE),
      recall = mean(d$recall, na.rm = TRUE),
      fdr = mean(d$fdr, na.rm = TRUE),
      top50_true = mean(d$top50_true, na.rm = TRUE),
      top50_precision = mean(d$top50_precision, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$method, out$fdr_scope, out$q_threshold), , drop = FALSE]
}

write_if_nonempty <- function(x, file) {
  if (is.null(x)) return(invisible(FALSE))
  if (is.data.frame(x) && nrow(x) == 0L) return(invisible(FALSE))
  utils::write.csv(x, file, row.names = FALSE)
  invisible(TRUE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$grid)) stop("Missing --grid")
if (is.null(args$`task-id`)) stop("Missing --task-id")

out_dir <- args$`out-dir` %||% args$out_dir %||% "singleview_comparator_results"
package_dir <- args$`package-dir` %||% args$package_dir %||% NULL
sim_dir <- args$`sim-dir` %||% args$sim_dir %||% NULL
methods <- tolower(as_chr_vec(args$methods %||% "hima,hima2,multimedia", default = c("hima", "hima2", "multimedia")))
methods <- intersect(methods, c("hima", "hima2", "multimedia"))
if (length(methods) == 0) stop("No supported methods requested.")

grid <- utils::read.csv(args$grid, stringsAsFactors = FALSE, check.names = FALSE)
task_offset <- as.integer(args$`task-offset` %||% args$task_offset %||% 0L)
array_task_id <- as.integer(args$`task-id`)
task_id <- array_task_id + task_offset
if (!is.finite(task_id) || task_id < 1L || task_id > nrow(grid)) {
  stop("effective task-id must be between 1 and ", nrow(grid), ". Got ", task_id)
}
row <- grid[task_id, , drop = FALSE]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

q_thresholds <- as_num_vec(row_value(row, "q_thresholds", row_value(row, "q_threshold", 0.25)), default = 0.25)
fdr_method <- as.character(row_value(row, "fdr_method", "BH"))
fdr_method <- match.arg(fdr_method, choices = c("BH", "BY"))
fdr_scope <- as.character(row_value(row, "fdr_scope", "both"))
if (!fdr_scope %in% c("global", "within_view", "both")) fdr_scope <- "both"
# Comparator runs should always report both correction families. This lets
# HIMA/HIMA2/multimedia be compared against Zentangler's global and
# within-view FDR summaries without rerunning the same model fit.
eval_scopes <- c("global", "within_view")
TIE <- as_num_or_null(row_value(row, "TIE", 1))
if (is.null(TIE)) TIE <- 1
outcome_type <- as.character(row_value(row, "outcome_type", "continuous"))
y_family <- if (identical(outcome_type, "binary")) "binomial" else "gaussian"

tag <- sprintf(
  "job%04d_singleview_%s_n%s_q%s_scope%s_snr%s_path%s_tie%s",
  as.integer(row_value(row, "job_id", task_id)),
  paste(safe_tag(methods), collapse = ""),
  safe_tag(row_value(row, "nsample", NA)),
  paste0(safe_tag(fdr_method), "grid", gsub("\\.", "p", as.character(min(q_thresholds))), "to", gsub("\\.", "p", as.character(max(q_thresholds)))),
  safe_tag(paste(eval_scopes, collapse = "-")),
  safe_tag(row_value(row, "snr", NA)),
  safe_tag(row_value(row, "n_pathways", NA)),
  gsub("\\.", "p", as.character(row_value(row, "TIE", 1)))
)

out_rds <- file.path(out_dir, paste0(tag, ".rds"))
out_summary <- file.path(out_dir, paste0(tag, "_summary.csv"))
out_detail <- file.path(out_dir, paste0(tag, "_detail.csv"))
out_mediators <- file.path(out_dir, paste0(tag, "_mediators.csv"))
out_errors <- file.path(out_dir, paste0(tag, "_errors.txt"))
out_settings <- file.path(out_dir, paste0(tag, "_settings.csv"))
out_grid_row <- file.path(out_dir, paste0(tag, "_grid_row.csv"))
out_session <- file.path(out_dir, paste0(tag, "_sessionInfo.txt"))
out_log <- file.path(out_dir, paste0(tag, "_log.txt"))
out_sim_source <- file.path(out_dir, paste0(tag, "_simulation_source.csv"))

sim_cache_file <- function(sim_dir, row, task_id) {
  if (is.null(sim_dir) || !nzchar(sim_dir)) return(NULL)
  sim_id <- as.integer(row_value(row, "sim_id", row_value(row, "job_id", task_id)))
  file.path(sim_dir, sprintf("sim%04d_sim.rds", sim_id))
}

log_msg <- function(..., .sep = "") {
  msg <- paste(..., sep = .sep)
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = out_log, append = TRUE)
  invisible(line)
}

log_block <- function(title, x) {
  lines <- capture.output(print(x))
  log_msg(title)
  for (ln in lines) log_msg("  ", ln)
  invisible(lines)
}

log_msg("Running single-view comparator HPC task")
log_msg("Array task ID: ", array_task_id)
log_msg("Task offset: ", task_offset)
log_msg("Effective grid row: ", task_id, " of ", nrow(grid))
log_msg("Methods: ", paste(methods, collapse = ", "))
log_msg("Grid file: ", normalizePath(args$grid, mustWork = FALSE))
log_msg("Output directory: ", normalizePath(out_dir, mustWork = FALSE))
log_msg("Package directory: ", ifelse(is.null(package_dir), "<installed package>", normalizePath(package_dir, mustWork = FALSE)))
log_msg("Simulation cache directory: ", ifelse(is.null(sim_dir), "<none; generate internally>", normalizePath(sim_dir, mustWork = FALSE)))
log_msg("Log file: ", normalizePath(out_log, mustWork = FALSE))
log_block("Grid row:", row)

log_msg("Loading Zentangler sources/package")
load_zentangler_sources(package_dir)
mae_to_blocks <- zt_get("zentangler_mae_to_blocks")
truth_key_fun <- zt_get("zentangler_truth_key")
utils::write.csv(row, out_grid_row, row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_session)
log_msg("Wrote grid row: ", out_grid_row)
log_msg("Wrote session info: ", out_session)

sim_path <- sim_cache_file(sim_dir, row, task_id)
if (!is.null(sim_path)) {
  if (!file.exists(sim_path)) stop("Requested simulation cache does not exist: ", sim_path)
  log_msg("Loading paired simulation cache: ", sim_path)
  sim <- readRDS(sim_path)
  log_msg("Finished loading paired simulation cache")
} else {
  log_msg("No simulation cache supplied; generating SIMMBA data internally")
  sim <- gen_simmba(
    nsample = as.integer(row_value(row, "nsample", 600L)),
    nrep = as.integer(row_value(row, "nrep", 10L)),
    seed = as.integer(row_value(row, "sim_seed", row_value(row, "seed", 1L))),
    p.train = 1,
    ygen.mode = as.character(row_value(row, "ygen_mode", "LM")),
    outcome.type = outcome_type,
    snr = as.numeric(row_value(row, "snr", 1)),
    n.pathways = as.integer(row_value(row, "n_pathways", 10L)),
    TIE = rep(TIE, 3)
  )
  log_msg("Finished SIMMBA generation")
}
utils::write.csv(
  data.frame(
    task_id = task_id,
    job_id = as.integer(row_value(row, "job_id", task_id)),
    sim_id = as.integer(row_value(row, "sim_id", row_value(row, "job_id", task_id))),
    sim_path = ifelse(is.null(sim_path), NA_character_, sim_path),
    used_cache = !is.null(sim_path),
    stringsAsFactors = FALSE
  ),
  out_sim_source,
  row.names = FALSE
)
log_msg("Wrote simulation source: ", out_sim_source)

nrep <- as.integer(row_value(row, "nrep", 10L))
top_n <- as.integer(row_value(row, "top_n", 50L))
hima_penalty <- as.character(args$`hima-penalty` %||% args$hima_penalty %||% "DBlasso")
hima2_penalty <- as.character(args$`hima2-penalty` %||% args$hima2_penalty %||% "DBlasso")
multimedia_model <- as.character(args$`multimedia-model` %||% args$multimedia_model %||% "lm")
multimedia_bootstrap <- as.integer(args$`multimedia-bootstrap` %||% args$multimedia_bootstrap %||% 100L)
multimedia_rf_trees <- as.integer(args$`multimedia-rf-trees` %||% args$multimedia_rf_trees %||% 500L)

all_mediator_rows <- list()
detail_rows <- list()
errors <- character(0)

for (rep_i in seq_len(nrep)) {
  log_msg("Starting replicate ", rep_i, " of ", nrep)
  mae_i <- sim$trainMae[[rep_i]]
  truth_i <- sim$truthDat[[rep_i]]
  if (is.null(mae_i)) {
    err <- paste0("rep_", rep_i, ": missing trainMae")
    errors <- c(errors, err)
    log_msg("ERROR: ", err)
    next
  }

  inputs <- mae_to_blocks(mae_i, duplicate_primary = "mean", fill_nonfinite_zero = FALSE)
  blocks <- inputs$blocks
  pheno <- inputs$pheno_df
  truth_key <- truth_key_fun(truth_i)
  log_msg(
    "Replicate ", rep_i, " contains ", length(blocks), " view(s): ",
    paste(paste0(names(blocks), "=", vapply(blocks, ncol, integer(1)), " features"), collapse = "; ")
  )

  for (method in methods) {
    method_start <- Sys.time()
    log_msg("Replicate ", rep_i, ": starting method ", method)
    method_tabs <- list()
    for (view in names(blocks)) {
      M <- blocks[[view]]
      view_start <- Sys.time()
      log_msg(
        "Replicate ", rep_i, ", method ", method, ", view ", view,
        ": fitting ", nrow(M), " samples x ", ncol(M), " mediators"
      )
      res <- try(
        switch(
          method,
          hima = run_hima_singleview(M, pheno, y_family = y_family, penalty = hima_penalty),
          hima2 = run_hima2_singleview(M, pheno, y_family = y_family, penalty = hima2_penalty),
          multimedia = run_multimedia_singleview(
            M,
            pheno,
            y_family = y_family,
            model = multimedia_model,
            bootstrap_repeats = multimedia_bootstrap,
            rf_trees = multimedia_rf_trees,
            seed = as.integer(row_value(row, "seed", row_value(row, "sim_seed", 1L))) + rep_i
          )
        ),
        silent = TRUE
      )

      key <- paste0("rep_", rep_i, "__", method, "__", view)
      if (inherits(res, "try-error")) {
        err <- paste0(key, ": ", conditionMessage(attr(res, "condition")))
        errors <- c(errors, err)
        log_msg("ERROR: ", err)
        tab <- standardize_method_table(NULL, colnames(M), view, method, rep_i)
      } else {
        tab <- standardize_method_table(res, colnames(M), view, method, rep_i)
        log_msg(
          "Replicate ", rep_i, ", method ", method, ", view ", view,
          ": finished with ", sum(is.finite(tab$p_primary) & tab$p_primary < 1, na.rm = TRUE),
          " mediator(s) carrying method p-values; runtime ",
          round(as.numeric(difftime(Sys.time(), view_start, units = "mins")), 3), " min"
        )
      }
      method_tabs[[view]] <- tab
    }

    combined <- do.call(rbind, method_tabs)
    rownames(combined) <- NULL
    combined <- apply_comparator_fdr(combined, fdr_method = fdr_method)
    all_mediator_rows[[length(all_mediator_rows) + 1L]] <- combined

    for (scope_i in eval_scopes) {
      q_col_i <- if (identical(scope_i, "within_view")) "q_primary_within_view" else "q_primary_global"
      for (threshold_i in q_thresholds) {
        detail_rows[[length(detail_rows) + 1L]] <- score_truth_recovery(
          tab = combined,
          truth_key = truth_key,
          rep = rep_i,
          method = method,
          q_threshold = threshold_i,
          q_col = q_col_i,
          fdr_scope = scope_i,
          top_n = top_n
        )
      }
    }
    log_msg(
      "Replicate ", rep_i, ": finished method ", method,
      "; runtime ", round(as.numeric(difftime(Sys.time(), method_start, units = "mins")), 3), " min"
    )
  }
  log_msg("Finished replicate ", rep_i, " of ", nrep)
}

log_msg("Combining mediator/detail rows")
mediators <- if (length(all_mediator_rows) > 0) do.call(rbind, all_mediator_rows) else data.frame()
detail <- if (length(detail_rows) > 0) do.call(rbind, detail_rows) else data.frame()
if (nrow(detail) > 0) rownames(detail) <- NULL
summary <- summarize_recovery(detail)

settings <- list(
  methods = methods,
  hima_penalty = hima_penalty,
  hima2_penalty = hima2_penalty,
  multimedia_model = multimedia_model,
  multimedia_bootstrap = multimedia_bootstrap,
  multimedia_rf_trees = multimedia_rf_trees,
  nrep = nrep,
  top_n = top_n,
  q_threshold = q_thresholds,
  fdr_method = fdr_method,
  fdr_scope_requested = fdr_scope,
  fdr_scopes_evaluated = eval_scopes,
  grid_row = row
)

res <- list(
  summary = summary,
  detail = detail,
  mediators = mediators,
  truth = sim$truthDat,
  errors = errors,
  settings = settings
)

log_msg("Writing output files")
saveRDS(res, out_rds)
utils::write.csv(summary, out_summary, row.names = FALSE)
utils::write.csv(detail, out_detail, row.names = FALSE)
write_if_nonempty(mediators, out_mediators)
settings_df <- data.frame(
  setting = names(settings)[names(settings) != "grid_row"],
  value = vapply(settings[names(settings) != "grid_row"], function(x) paste(as.character(x), collapse = ","), character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(settings_df, out_settings, row.names = FALSE)
if (length(errors) > 0) writeLines(errors, out_errors) else writeLines(character(0), out_errors)

log_msg("Finished single-view comparator task ", task_id)
log_msg("Saved RDS: ", out_rds)
log_msg("Saved summary: ", out_summary)
log_msg("Saved detail: ", out_detail)
log_msg("Saved mediators: ", out_mediators)
log_msg("Saved errors: ", out_errors)
log_msg("Number of recorded errors: ", length(errors))
log_block("Summary:", summary)
