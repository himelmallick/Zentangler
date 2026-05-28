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

empty_standardized_method_table <- function() {
  data.frame(
    rep = integer(0),
    method = character(0),
    omics = character(0),
    mediator = character(0),
    a = numeric(0),
    b = numeric(0),
    score = numeric(0),
    abs_score = numeric(0),
    p_primary = numeric(0),
    p_source = character(0),
    stringsAsFactors = FALSE
  )
}

standardize_method_table <- function(x, all_mediators, view, method, rep_id, p_default = 1) {
  all_mediators <- as.character(all_mediators)

  # HIMA/HIMA2 usually return a selected mediator result table, not a full
  # feature-level table. For a fair comparison with Zentangler, build a full
  # feature table first and merge the returned HIMA/HIMA2 p-values onto it.
  # Features not returned by HIMA/HIMA2 are retained as non-selected features
  # with p_primary = 1, so the downstream FDR correction is applied over all
  # tested features rather than only the selected mediator subset.
  if (method %in% c("hima", "hima2")) {
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
      p_source = "not_returned_by_method",
      stringsAsFactors = FALSE
    )

    if (is.null(x)) return(out)

    # HIMA/HIMA2 objects can have class "hima" while internally behaving
    # like a named list of same-length vectors (ID, alpha, beta, alpha*beta,
    # p-value, etc.). as.data.frame() is not reliable for all HIMA versions,
    # so coerce that list structure explicitly before parsing columns.
    if (is.list(x) && !is.data.frame(x)) {
      keep <- vapply(x, function(z) is.atomic(z) && length(z) > 0, logical(1))
      lens <- vapply(x[keep], length, integer(1))
      if (length(lens) > 0) {
        target_len <- max(lens)
        keep_names <- names(lens)[lens == target_len]
        x <- try(as.data.frame(x[keep_names], check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
      } else {
        x <- data.frame()
      }
    } else {
      x <- try(as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
    }

    if (inherits(x, "try-error") || is.null(x) || nrow(x) == 0) return(out)

    raw_names <- colnames(x)
    safe_names <- make.names(raw_names, unique = TRUE)
    colnames(x) <- safe_names

    pick_hima_col <- function(candidates) {
      candidate_safe <- make.names(candidates, unique = FALSE)
      hit <- safe_names[safe_names %in% candidate_safe]
      if (length(hit) > 0) return(hit[1])
      hit_raw <- raw_names[raw_names %in% candidates]
      if (length(hit_raw) > 0) return(make.names(hit_raw[1], unique = FALSE))
      NA_character_
    }

    num_hima_col <- function(candidates) {
      nm <- pick_hima_col(candidates)
      if (is.na(nm) || !(nm %in% colnames(x))) return(rep(NA_real_, nrow(x)))
      suppressWarnings(as.numeric(x[[nm]]))
    }

    id_col <- pick_hima_col(c("mediator", "Mediator", "ID", "Index", "id", "feature", "Feature", "name", "Name"))
    med <- if (!is.na(id_col) && id_col %in% colnames(x)) as.character(x[[id_col]]) else rownames(x)
    med[is.na(med) | !nzchar(med)] <- rownames(x)[is.na(med) | !nzchar(med)]

    a <- num_hima_col(c("a", "alpha", "alpha_hat", "alpha.hat", "ahat", "a_hat", "effect.a", "beta.xm"))
    b <- num_hima_col(c("b", "beta", "beta_hat", "beta.hat", "bhat", "b_hat", "effect.b", "beta.my"))
    score_from_method <- num_hima_col(c("score", "indirect", "indirect_effect", "indirect.effect", "alpha_beta", "alpha*beta", "ab", "IDE", "NIE"))
    score <- ifelse(is.finite(score_from_method), score_from_method, a * b)

    p_primary <- num_hima_col(c("p_primary", "p-value", "p.value", "p_value", "pvalue", "p", "pval", "p.val", "p.joint", "p_joint", "pmax", "p.max", "raw.p", "p_raw"))
    p_col <- pick_hima_col(c("p_primary", "p-value", "p.value", "p_value", "pvalue", "p", "pval", "p.val", "p.joint", "p_joint", "pmax", "p.max", "raw.p", "p_raw"))
    q_col <- pick_hima_col(c("q", "q.value", "q_value", "qvalue", "FDR", "BH.FDR", "adj.p", "adjusted.p", "padj", "qval"))
    if (all(!is.finite(p_primary)) && !is.na(q_col)) {
      p_primary <- suppressWarnings(as.numeric(x[[q_col]]))
      p_src <- q_col
    } else {
      p_src <- p_col
    }

    tmp <- data.frame(
      omics = view,
      mediator = med,
      a = a,
      b = b,
      score = score,
      abs_score = abs(score),
      p_primary = p_primary,
      p_source = p_src,
      stringsAsFactors = FALSE
    )
    tmp <- tmp[!duplicated(tmp$mediator), , drop = FALSE]
    idx <- match(out$mediator, tmp$mediator)
    ok <- !is.na(idx)

    out$a[ok] <- tmp$a[idx[ok]]
    out$b[ok] <- tmp$b[idx[ok]]
    out$score[ok] <- tmp$score[idx[ok]]
    out$p_primary[ok] <- tmp$p_primary[idx[ok]]
    out$p_source[ok] <- tmp$p_source[idx[ok]]

    out$b[!is.finite(out$b)] <- 0
    out$score[!is.finite(out$score)] <- 0
    out$abs_score <- abs(out$score)
    out$p_primary[!is.finite(out$p_primary)] <- p_default
    out$p_primary <- pmin(1, pmax(0, out$p_primary))
    return(out)
  }

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
    q_primary_global = NA_real_,
    q_primary_within_view = NA_real_,
    q_primary = NA_real_,
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
  q_col <- pick_col(x, c("q_primary", "q_primary_global", "q_primary_within_view", "multimedia_fdr_hat", "fdr_hat", "q", "q.value", "q_value", "qvalue", "FDR", "BH.FDR", "adj.p", "adjusted.p", "padj", "qval"))
  q_global_col <- pick_col(x, c("q_primary_global", "multimedia_fdr_hat", "fdr_hat"))
  q_within_col <- pick_col(x, c("q_primary_within_view", "multimedia_fdr_hat", "fdr_hat"))
  q_primary_col <- pick_col(x, c("q_primary", "multimedia_fdr_hat", "fdr_hat"))

  tmp <- data.frame(mediator = med, stringsAsFactors = FALSE)
  tmp$a <- if (!is.na(a_col)) suppressWarnings(as.numeric(x[[a_col]])) else NA_real_
  tmp$b <- if (!is.na(b_col)) suppressWarnings(as.numeric(x[[b_col]])) else NA_real_
  tmp$score <- if (!is.na(score_col)) suppressWarnings(as.numeric(x[[score_col]])) else tmp$a * tmp$b
  tmp$p_primary <- if (!is.na(p_col)) suppressWarnings(as.numeric(x[[p_col]])) else NA_real_
  tmp$q_primary_global <- if (!is.na(q_global_col)) suppressWarnings(as.numeric(x[[q_global_col]])) else NA_real_
  tmp$q_primary_within_view <- if (!is.na(q_within_col)) suppressWarnings(as.numeric(x[[q_within_col]])) else NA_real_
  tmp$q_primary <- if (!is.na(q_primary_col)) suppressWarnings(as.numeric(x[[q_primary_col]])) else NA_real_
  if (all(!is.finite(tmp$p_primary)) && !is.na(q_col)) {
    if (method %in% c("hima", "hima2")) {
      tmp$p_primary <- suppressWarnings(as.numeric(x[[q_col]]))
    }
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
  out$q_primary_global[ok] <- tmp$q_primary_global[idx[ok]]
  out$q_primary_within_view[ok] <- tmp$q_primary_within_view[idx[ok]]
  out$q_primary[ok] <- tmp$q_primary[idx[ok]]
  out$p_source[ok] <- p_src

  out$b[!is.finite(out$b)] <- 0
  out$score[!is.finite(out$score)] <- 0
  out$abs_score <- abs(out$score)
  out$p_primary[!is.finite(out$p_primary)] <- if (method %in% c("hima", "hima2")) p_default else NA_real_
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
  if ("FDRcut" %in% argn) call_args$FDRcut <- 1
  if ("Bonfcut" %in% argn) call_args$Bonfcut <- 1
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
  if ("FDRcut" %in% argn) call_args$FDRcut <- 1
  if ("Bonfcut" %in% argn) call_args$Bonfcut <- 1
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

multimedia_null_contrast_safe <- function(fit, exper) {
  f <- get("null_contrast", envir = asNamespace("multimedia"))
  argn <- names(formals(f))
  call_args <- list(model = fit, exper = exper)
  if ("nullification" %in% argn) call_args$nullification <- "M->Y"
  if ("f" %in% argn) call_args$f <- get("indirect_pathwise", envir = asNamespace("multimedia"))
  do.call(f, call_args)
}

multimedia_fdr_summary_safe <- function(contrast, q_value = 1) {
  f <- get("fdr_summary", envir = asNamespace("multimedia"))
  argn <- names(formals(f))
  call_args <- list(contrast)
  if ("effect" %in% argn) call_args$effect <- "indirect_pathwise"
  if ("q_value" %in% argn) call_args$q_value <- q_value
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
                                      cores = 1L,
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
    multimedia::lm_model(progress = FALSE)
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
  out$q_primary_global <- NA_real_
  out$q_primary_within_view <- NA_real_
  out$q_primary <- NA_real_
  out$multimedia_fdr_hat <- NA_real_
  out$multimedia_keep_q1 <- NA

  native_fdr <- try({
    contrast <- multimedia_null_contrast_safe(fit, exper)
    multimedia_fdr_summary_safe(contrast, q_value = 1)
  }, silent = TRUE)

  if (!inherits(native_fdr, "try-error") && nrow(as.data.frame(native_fdr)) > 0) {
    fdr_tab <- as.data.frame(native_fdr, stringsAsFactors = FALSE)
    med_fdr_col <- pick_col(fdr_tab, c("mediator", "Mediator", "feature", "Feature", "term", "variable", "name"))
    fdr_col <- pick_col(fdr_tab, c("fdr_hat", "fdr", "q_value", "q", "FDR"))
    keep_col <- pick_col(fdr_tab, c("keep", "selected", "active"))
    source_col <- pick_col(fdr_tab, c("source", "Source"))
    if (!is.na(source_col)) {
      fdr_tab <- fdr_tab[as.character(fdr_tab[[source_col]]) == "real", , drop = FALSE]
    }
    if (!is.na(med_fdr_col) && !is.na(fdr_col) && nrow(fdr_tab) > 0) {
      fdr_tab <- stats::aggregate(
        fdr_tab[[fdr_col]],
        by = list(mediator = as.character(fdr_tab[[med_fdr_col]])),
        FUN = min,
        na.rm = TRUE
      )
      colnames(fdr_tab)[2] <- "fdr_hat"
      hit <- match(out$mediator, fdr_tab$mediator)
      ok <- !is.na(hit)
      out$multimedia_fdr_hat[ok] <- suppressWarnings(as.numeric(fdr_tab$fdr_hat[hit[ok]]))
      out$q_primary_global[ok] <- out$multimedia_fdr_hat[ok]
      out$q_primary_within_view[ok] <- out$multimedia_fdr_hat[ok]
      out$q_primary[ok] <- out$multimedia_fdr_hat[ok]
    }
    if (!is.na(keep_col) && !is.na(med_fdr_col) && nrow(fdr_tab) > 0) {
      keep_map <- stats::aggregate(
        as.logical(native_fdr[[keep_col]]),
        by = list(mediator = as.character(native_fdr[[med_fdr_col]])),
        FUN = any,
        na.rm = TRUE
      )
      hit <- match(out$mediator, keep_map$mediator)
      ok <- !is.na(hit)
      out$multimedia_keep_q1[ok] <- keep_map$x[hit[ok]]
    }
  }

  bootstrap_repeats <- max(0L, as.integer(bootstrap_repeats))
  if (all(!is.finite(out$q_primary)) && bootstrap_repeats > 0L) {
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

    fit_bootstrap_one <- function(bb) {
      idx <- sample(seq_len(nrow(df)), size = nrow(df), replace = TRUE)
      df_b <- df[idx, , drop = FALSE]

      boot_tab <- try({
        exper_b <- multimedia::mediation_data(df_b, "Y", "A", tidyselect::all_of(med_names))
        base_model_b <- if (identical(model, "rf")) {
          f <- get("rf_model", envir = asNamespace("multimedia"))
          argn <- names(formals(f))
          if ("num.trees" %in% argn) do.call(f, list(num.trees = as.integer(rf_trees))) else f()
        } else {
          multimedia::lm_model(progress = FALSE)
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
      out_b <- setNames(rep(NA_real_, length(out$mediator)), out$mediator)
      hit <- match(tmp$mediator, names(out_b))
      ok <- !is.na(hit)
      out_b[hit[ok]] <- tmp$score[ok]
      out_b
    }

    cores <- max(1L, as.integer(cores))
    if (.Platform$OS.type != "windows" && cores > 1L && bootstrap_repeats > 1L) {
      boot_list <- parallel::mclapply(
        seq_len(bootstrap_repeats),
        fit_bootstrap_one,
        mc.cores = min(cores, bootstrap_repeats),
        mc.set.seed = TRUE
      )
      for (bb in seq_along(boot_list)) {
        if (is.numeric(boot_list[[bb]]) && length(boot_list[[bb]]) == ncol(boot_scores)) {
          boot_scores[bb, ] <- boot_list[[bb]]
        }
      }
    } else {
      for (bb in seq_len(bootstrap_repeats)) {
        boot_scores[bb, ] <- fit_bootstrap_one(bb)
      }
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

apply_comparator_fdr <- function(tab, fdr_method = c("BH", "BY"), tested_by_view = NULL) {
  fdr_method <- match.arg(fdr_method)
  if (is.null(tab) || nrow(tab) == 0) return(tab)
  q_global_existing <- if ("q_primary_global" %in% colnames(tab)) tab$q_primary_global else rep(NA_real_, nrow(tab))
  q_within_existing <- if ("q_primary_within_view" %in% colnames(tab)) tab$q_primary_within_view else rep(NA_real_, nrow(tab))
  q_existing <- if ("q_primary" %in% colnames(tab)) tab$q_primary else rep(NA_real_, nrow(tab))

  # HIMA/HIMA2 return only a selected mediator table rather than one row for
  # every tested mediator. For fair FDR correction, adjust over the full
  # single-view testing universe using p.adjust(..., n = m), while keeping the
  # output table clean and limited to returned mediators.
  if (is.null(tested_by_view)) {
    tested_by_view <- stats::setNames(
      as.integer(table(tab$omics)),
      names(table(tab$omics))
    )
  } else {
    tested_names <- names(tested_by_view)
    tested_by_view <- as.integer(tested_by_view)
    names(tested_by_view) <- tested_names
  }

  tested_n_for_view <- function(view) {
    view <- as.character(view)[1]
    if (!view %in% names(tested_by_view)) {
      return(sum(tab$omics == view, na.rm = TRUE))
    }
    n <- tested_by_view[[view]]
    if (is.null(n) || !is.finite(n) || n < 1L) {
      sum(tab$omics == view, na.rm = TRUE)
    } else {
      as.integer(n)
    }
  }

  tab$q_primary_global <- q_global_existing
  for (method_i in unique(tab$method)) {
    idx_m <- which(tab$method == method_i)
    ok_m <- is.finite(tab$p_primary[idx_m]) & !is.finite(tab$q_primary_global[idx_m])
    if (any(ok_m)) {
      n_global <- sum(vapply(unique(tab$omics[idx_m]), tested_n_for_view, integer(1)), na.rm = TRUE)
      n_global <- max(n_global, sum(ok_m))
      tab$q_primary_global[idx_m[ok_m]] <- p.adjust(tab$p_primary[idx_m[ok_m]], method = fdr_method, n = n_global)
    }
  }

  tab$q_primary_within_view <- q_within_existing
  for (idx in split(seq_len(nrow(tab)), paste(tab$method, tab$omics, sep = "::"))) {
    ok_i <- is.finite(tab$p_primary[idx]) & !is.finite(tab$q_primary_within_view[idx])
    if (any(ok_i)) {
      n_view <- tested_n_for_view(tab$omics[idx][1])
      n_view <- max(n_view, sum(ok_i))
      tab$q_primary_within_view[idx[ok_i]] <- p.adjust(tab$p_primary[idx[ok_i]], method = fdr_method, n = n_view)
    }
  }
  tab$q_primary <- ifelse(is.finite(q_existing), q_existing, tab$q_primary_global)
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
      precision = 0,
      recall = 0,
      fdr = 0,
      top50_true = 0,
      top50_precision = 0,
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
  top_precision <- if (length(top_idx) > 0) top_true / length(top_idx) else 0

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
    precision = if (n_active > 0) true_active / n_active else 0,
    recall = if (n_true > 0) true_active / n_true else 0,
    fdr = if (n_active > 0) false_active / n_active else 0,
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
    mean_metric <- function(x) {
      x <- ifelse(is.finite(x), x, 0)
      mean(x)
    }
    data.frame(
      g,
      n_active = mean_metric(d$n_active),
      true_active = mean_metric(d$true_active),
      false_active = mean_metric(d$false_active),
      n_true = mean_metric(d$n_true),
      precision = mean_metric(d$precision),
      precision_discovery_only = if (sum(d$n_active, na.rm = TRUE) > 0) {
        sum(d$true_active, na.rm = TRUE) / sum(d$n_active, na.rm = TRUE)
      } else {
        NA_real_
      },
      recall = mean_metric(d$recall),
      fdr = mean_metric(d$fdr),
      fdr_discovery_only = if (sum(d$n_active, na.rm = TRUE) > 0) {
        sum(d$false_active, na.rm = TRUE) / sum(d$n_active, na.rm = TRUE)
      } else {
        NA_real_
      },
      n_reps = length(idx),
      n_reps_with_discovery = sum(d$n_active > 0, na.rm = TRUE),
      top50_true = mean_metric(d$top50_true),
      top50_precision = mean_metric(d$top50_precision),
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
multimedia_cores <- as.integer(args$`multimedia-cores` %||% args$multimedia_cores %||% Sys.getenv("MULTIMEDIA_CORES", unset = "1"))
if (!is.finite(multimedia_cores) || multimedia_cores < 1L) multimedia_cores <- 1L

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
            cores = multimedia_cores,
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
    combined <- apply_comparator_fdr(
      combined,
      fdr_method = fdr_method,
      tested_by_view = vapply(blocks, ncol, integer(1))
    )
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
  multimedia_cores = multimedia_cores,
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
