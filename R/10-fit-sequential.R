# -----------------------------------------------------------------------------
# K-layer sequential mediation
# -----------------------------------------------------------------------------

zentangler_validate_stage_views <- function(stage_views, available_views) {
  if (is.null(stage_views)) {
    stage_views <- as.list(available_views)
    names(stage_views) <- paste0("stage", seq_along(stage_views))
  }
  if (!is.list(stage_views) || length(stage_views) < 2L) {
    stop("stage_views must be a list with at least two mediator layers.", call. = FALSE)
  }

  stage_views <- lapply(stage_views, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) stop("Each stage_views entry must contain at least one view.", call. = FALSE)
    x
  })
  if (is.null(names(stage_views)) || any(!nzchar(names(stage_views)))) {
    names(stage_views) <- paste0("stage", seq_along(stage_views))
  }

  missing <- setdiff(unique(unlist(stage_views, use.names = FALSE)), available_views)
  if (length(missing) > 0L) {
    stop("stage_views contains views not present in the input: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  duplicated_views <- unique(unlist(stage_views, use.names = FALSE)[
    duplicated(unlist(stage_views, use.names = FALSE))
  ])
  if (length(duplicated_views) > 0L) {
    stop("Views can appear in only one sequential layer: ", paste(duplicated_views, collapse = ", "), call. = FALSE)
  }

  stage_views
}

zentangler_map_stage_views <- function(stage_views, view_map) {
  if (is.null(stage_views) || is.null(view_map) || nrow(view_map) == 0L) return(stage_views)
  original_to_safe <- setNames(as.character(view_map$view), as.character(view_map$experiment_name))
  safe_to_safe <- setNames(as.character(view_map$view), as.character(view_map$view))
  lapply(stage_views, function(x) {
    x <- as.character(x)
    mapped <- original_to_safe[x]
    missing <- is.na(mapped)
    mapped[missing] <- safe_to_safe[x[missing]]
    mapped
  })
}

zentangler_feature_key <- function(view, mediator) {
  paste(view, mediator, sep = "::")
}

zentangler_residualize_vector <- function(y, Z = NULL) {
  y <- as.numeric(y)
  if (is.null(Z) || ncol(Z) == 0L) return(y)
  fit <- lm.fit(x = cbind(`(Intercept)` = 1, as.matrix(Z)), y = y)
  as.numeric(fit$residuals)
}

zentangler_cor_pvalue <- function(r, n, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  if (!is.finite(r) || n < 4L || abs(r) >= 1) {
    if (is.finite(r) && abs(r) >= 1) return(0)
    return(NA_real_)
  }
  if (identical(method, "pearson")) {
    stat <- r * sqrt((n - 2) / max(1e-12, 1 - r^2))
    2 * stats::pt(abs(stat), df = n - 2, lower.tail = FALSE)
  } else {
    # Large-sample approximation used for speed in high-dimensional screens.
    stat <- r * sqrt((n - 2) / max(1e-12, 1 - r^2))
    2 * stats::pt(abs(stat), df = n - 2, lower.tail = FALSE)
  }
}

zentangler_standardize_for_cor <- function(M, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  M <- as.matrix(M)
  dn <- dimnames(M)
  if (identical(method, "spearman")) {
    M <- apply(M, 2, rank, ties.method = "average", na.last = "keep")
    M <- as.matrix(M)
    dimnames(M) <- dn
  }
  center <- colMeans(M, na.rm = TRUE)
  scale <- apply(M, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale <= 0] <- NA_real_
  out <- sweep(sweep(M, 2, center, "-"), 2, scale, "/")
  dimnames(out) <- dn
  out
}

zentangler_screen_cor_links <- function(
  from_df,
  to_df,
  X,
  C = NULL,
  cor_method = c("spearman", "pearson"),
  cor_fdr_method = c("BH", "BY"),
  min_abs_cor = 0.3,
  cor_q_threshold = 0.25,
  residualize_links = TRUE,
  max_links = 5000L
) {
  cor_method <- match.arg(cor_method)
  cor_fdr_method <- validate_fdr_method(cor_fdr_method)
  max_links <- max(1L, as.integer(max_links))

  if (nrow(from_df) != nrow(to_df)) stop("from_df and to_df must have the same samples.", call. = FALSE)
  if (ncol(from_df) == 0L || ncol(to_df) == 0L) return(data.frame())

  Z <- NULL
  if (isTRUE(residualize_links)) {
    Z <- cbind(X = as.numeric(X))
    if (!is.null(C) && ncol(C) > 0L) Z <- cbind(Z, as.matrix(C))
  }

  A <- as.matrix(from_df)
  B <- as.matrix(to_df)
  if (!is.null(Z)) {
    a_dn <- dimnames(A)
    b_dn <- dimnames(B)
    A <- apply(A, 2, zentangler_residualize_vector, Z = Z)
    B <- apply(B, 2, zentangler_residualize_vector, Z = Z)
    A <- as.matrix(A)
    B <- as.matrix(B)
    dimnames(A) <- a_dn
    dimnames(B) <- b_dn
  }

  A <- zentangler_standardize_for_cor(A, method = cor_method)
  B <- zentangler_standardize_for_cor(B, method = cor_method)
  keep_a <- colSums(!is.finite(A)) == 0L
  keep_b <- colSums(!is.finite(B)) == 0L
  A <- A[, keep_a, drop = FALSE]
  B <- B[, keep_b, drop = FALSE]
  if (ncol(A) == 0L || ncol(B) == 0L) return(data.frame())

  R <- crossprod(A, B) / max(1L, nrow(A) - 1L)
  idx <- which(abs(R) >= min_abs_cor, arr.ind = TRUE)
  if (nrow(idx) == 0L) return(data.frame())

  cor_val <- R[idx]
  p_val <- vapply(cor_val, zentangler_cor_pvalue, numeric(1), n = nrow(A), method = cor_method)
  q_val <- stats::p.adjust(p_val, method = cor_fdr_method)
  out <- data.frame(
    from_key = colnames(A)[idx[, 1]],
    to_key = colnames(B)[idx[, 2]],
    cor = as.numeric(cor_val),
    abs_cor = abs(as.numeric(cor_val)),
    p_cor = as.numeric(p_val),
    q_cor = as.numeric(q_val),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$q_cor) & out$q_cor <= cor_q_threshold, , drop = FALSE]
  if (nrow(out) == 0L) return(out)
  out <- out[order(out$q_cor, -out$abs_cor), , drop = FALSE]
  out <- out[seq_len(min(max_links, nrow(out))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

zentangler_estimate_link_effects <- function(links, data_by_key, X, C = NULL) {
  if (is.null(links) || nrow(links) == 0L) return(links)
  d_hat <- p_d <- rep(NA_real_, nrow(links))
  for (i in seq_len(nrow(links))) {
    from <- data_by_key[[links$from_key[i]]]
    to <- data_by_key[[links$to_key[i]]]
    dat <- data.frame(to = as.numeric(to), X = as.numeric(X), from = as.numeric(from), check.names = FALSE)
    if (!is.null(C) && ncol(C) > 0L) dat <- cbind(dat, as.data.frame(C, check.names = FALSE))
    fit <- try(stats::lm(to ~ ., data = dat), silent = TRUE)
    if (inherits(fit, "try-error")) next
    co <- summary(fit)$coefficients
    if ("from" %in% rownames(co)) {
      d_hat[i] <- co["from", "Estimate"]
      p_d[i] <- co["from", "Pr(>|t|)"]
    }
  }
  links$d_hat <- d_hat
  links$p_d <- p_d
  links
}

zentangler_fit_terminal_effects <- function(
  Y,
  X,
  terminal_df,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  outcome_model = c("joint", "marginal")
) {
  y_family <- match.arg(y_family)
  outcome_model <- match.arg(outcome_model)
  terminal_keys <- colnames(terminal_df)
  out <- data.frame(
    terminal_key = terminal_keys,
    b = NA_real_,
    p_b = NA_real_,
    outcome_model = outcome_model,
    stringsAsFactors = FALSE
  )
  if (length(terminal_keys) == 0L) return(out)

  use_joint <- identical(outcome_model, "joint")
  n_extra <- 1L + if (!is.null(C)) ncol(C) else 0L
  if (use_joint && ncol(terminal_df) >= (nrow(terminal_df) - n_extra - 2L)) {
    warning("Too many terminal mediators for a stable joint outcome refit; using marginal terminal effects.")
    use_joint <- FALSE
    out$outcome_model <- "marginal"
  }

  fit_one <- function(key) {
    dat <- data.frame(X = as.numeric(X), M = as.numeric(terminal_df[, key]), check.names = FALSE)
    if (!is.null(C) && ncol(C) > 0L) dat <- cbind(dat, as.data.frame(C, check.names = FALSE))
    if (identical(y_family, "survival")) {
      fit <- try(suppressWarnings(survival::coxph(Y ~ ., data = dat)), silent = TRUE)
    } else if (identical(y_family, "binomial")) {
      fit <- try(suppressWarnings(stats::glm(as.numeric(Y) ~ ., data = dat, family = stats::binomial())), silent = TRUE)
    } else {
      fit <- try(stats::lm(as.numeric(Y) ~ ., data = dat), silent = TRUE)
    }
    if (inherits(fit, "try-error")) return(c(b = NA_real_, p_b = NA_real_))
    sm <- summary(fit)$coefficients
    if (!("M" %in% rownames(sm))) return(c(b = NA_real_, p_b = NA_real_))
    p_col <- grep("^Pr\\(", colnames(sm), value = TRUE)
    if (length(p_col) == 0L) p_col <- colnames(sm)[ncol(sm)]
    c(b = unname(sm["M", "Estimate"]), p_b = unname(sm["M", p_col[[1L]]]))
  }

  if (!use_joint) {
    vals <- t(vapply(terminal_keys, fit_one, numeric(2)))
    out$b <- vals[, "b"]
    out$p_b <- vals[, "p_b"]
    return(out)
  }

  raw_names <- c("X", if (!is.null(C) && ncol(C) > 0L) colnames(C) else character(0), terminal_keys)
  safe_names <- make.names(raw_names, unique = TRUE)
  dat <- data.frame(X = as.numeric(X), check.names = FALSE)
  if (!is.null(C) && ncol(C) > 0L) dat <- cbind(dat, as.data.frame(C, check.names = FALSE))
  dat <- cbind(dat, as.data.frame(terminal_df, check.names = FALSE))
  colnames(dat) <- safe_names

  if (identical(y_family, "survival")) {
    fit <- try(suppressWarnings(survival::coxph(Y ~ ., data = dat)), silent = TRUE)
  } else if (identical(y_family, "binomial")) {
    fit <- try(suppressWarnings(stats::glm(as.numeric(Y) ~ ., data = dat, family = stats::binomial())), silent = TRUE)
  } else {
    fit <- try(stats::lm(as.numeric(Y) ~ ., data = dat), silent = TRUE)
  }
  if (inherits(fit, "try-error")) {
    warning("Joint terminal outcome refit failed; using marginal terminal effects.")
    return(zentangler_fit_terminal_effects(
      Y = Y, X = X, terminal_df = terminal_df, C = C, y_family = y_family, outcome_model = "marginal"
    ))
  }

  sm <- summary(fit)$coefficients
  p_col <- grep("^Pr\\(", colnames(sm), value = TRUE)
  if (length(p_col) == 0L) p_col <- colnames(sm)[ncol(sm)]
  name_map <- setNames(raw_names, safe_names)
  for (rn in rownames(sm)) {
    key <- name_map[[rn]]
    if (is.null(key) || !(key %in% terminal_keys)) next
    ii <- match(key, out$terminal_key)
    out$b[ii] <- sm[rn, "Estimate"]
    out$p_b[ii] <- sm[rn, p_col[[1L]]
    ]
  }
  out
}

zentangler_extend_paths <- function(paths, links, max_paths = 10000L) {
  if (length(paths) == 0L || is.null(links) || nrow(links) == 0L) return(list())
  by_from <- split(seq_len(nrow(links)), links$from_key)
  out <- list()
  ii <- 0L
  for (path in paths) {
    last_key <- tail(path$keys, 1L)
    hit <- by_from[[last_key]]
    if (is.null(hit)) next
    for (j in hit) {
      ii <- ii + 1L
      out[[ii]] <- list(
        keys = c(path$keys, links$to_key[j]),
        link_rows = c(path$link_rows, j)
      )
      if (ii >= max_paths) return(out)
    }
  }
  out
}

#' Fit k-layer sequential Zentangler mediation
#'
#' Fits a sequential mediation screen for ordered mediator layers such as
#' `X -> M1 -> M2 -> ... -> Mk -> Y`. The first layer is screened by association
#' with the exposure. Adjacent mediator layers are screened by mediator-to-
#' mediator correlation, then retained links are refit by regression to estimate
#' directed link coefficients.
#'
#' @param mae A `MultiAssayExperiment`.
#' @param x_var Exposure column in `colData(mae)`.
#' @param y_var Outcome column in `colData(mae)` for Gaussian/binomial outcomes.
#' @param stage_views List of ordered mediator layers. Each element is a
#'   character vector of experiment/view names. If `NULL`, `view_names` order is
#'   used as one view per layer.
#' @param view_names Optional subset/order of MAE experiments.
#' @param assay_names Optional assay names passed to MAE extraction.
#' @param covariates Optional phenotype covariates.
#' @param y_family Outcome family: `"gaussian"`, `"binomial"`, or `"survival"`.
#' @param survival_time_var,survival_event_var Survival outcome columns.
#' @param sis_n Number of first-layer mediators to retain per first-layer view.
#' @param sis_rank Ranking for first-layer screen.
#' @param screen_method First-layer screen method.
#' @param cor_method Correlation method for adjacent mediator-layer screens.
#' @param min_abs_cor Minimum absolute correlation for adjacent-layer links.
#' @param cor_q_threshold FDR cutoff for adjacent-layer correlation links.
#' @param cor_fdr_method FDR method for adjacent-layer correlation links.
#' @param residualize_links Whether to residualize mediator-to-mediator
#'   correlations on exposure and covariates before screening.
#' @param max_links_per_transition Maximum retained links per adjacent transition.
#' @param max_paths Maximum full paths to keep while expanding across k layers.
#' @param outcome_model `"joint"` fits terminal mediators jointly when feasible;
#'   `"marginal"` fits one terminal mediator at a time.
#' @param fdr_method FDR method for final path-level p-values.
#' @param fdr_scope Currently `"global"` only for path-level q-values.
#' @param seed Random seed.
#' @param duplicate_primary,fill_nonfinite_zero MAE extraction options.
#'
#' @return A list with settings, transition links, terminal effects, and
#'   `sequential_paths`.
#' @export
fit_sequential_zentangler <- function(
  mae,
  x_var,
  y_var = NULL,
  stage_views = NULL,
  view_names = NULL,
  assay_names = NULL,
  covariates = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  survival_time_var = NULL,
  survival_event_var = NULL,
  sis_n = 50L,
  sis_rank = c("abs_a", "pvalue"),
  screen_method = c("sis", "none"),
  cor_method = c("spearman", "pearson"),
  min_abs_cor = 0.3,
  cor_q_threshold = 0.25,
  cor_fdr_method = c("BH", "BY"),
  residualize_links = TRUE,
  max_links_per_transition = 5000L,
  max_paths = 10000L,
  outcome_model = c("joint", "marginal"),
  fdr_method = c("BH", "BY"),
  fdr_scope = c("global"),
  seed = 1,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE
) {
  start_time <- Sys.time()
  duplicate_primary <- match.arg(duplicate_primary)
  y_family <- match.arg(y_family)
  sis_rank <- match.arg(sis_rank)
  screen_method <- match.arg(screen_method)
  cor_method <- match.arg(cor_method)
  cor_fdr_method <- validate_fdr_method(cor_fdr_method)
  outcome_model <- match.arg(outcome_model)
  fdr_method <- validate_fdr_method(fdr_method)
  fdr_scope <- match.arg(fdr_scope)
  set.seed(seed)

  if (is.null(x_var) || length(x_var) != 1L || !nzchar(as.character(x_var))) {
    stop("x_var is required.", call. = FALSE)
  }
  if (!identical(y_family, "survival") && (is.null(y_var) || length(y_var) != 1L || !nzchar(as.character(y_var)))) {
    stop("y_var is required unless y_family = 'survival'.", call. = FALSE)
  }

  inputs <- zentangler_mae_to_blocks(
    mae = mae,
    view_names = view_names,
    assay_names = assay_names,
    duplicate_primary = duplicate_primary,
    fill_nonfinite_zero = fill_nonfinite_zero
  )
  stage_views <- zentangler_map_stage_views(stage_views, inputs$view_map)
  stage_views <- zentangler_validate_stage_views(stage_views, names(inputs$blocks))
  needed <- unique(c(x_var, y_var, survival_time_var, survival_event_var, covariates))
  al <- align_samples_multiview_blocks(inputs$blocks, inputs$pheno_df, needed_pheno_cols = needed)
  blocks <- al$blocks
  pheno_df <- al$pheno
  X <- pheno_df[[x_var]]
  if (!is.numeric(X)) stop("x_var must be numeric.", call. = FALSE)
  if (identical(y_family, "survival")) {
    if (is.null(survival_time_var) || is.null(survival_event_var)) {
      stop("survival_time_var and survival_event_var are required for survival outcomes.", call. = FALSE)
    }
    time_raw <- suppressWarnings(as.numeric(pheno_df[[survival_time_var]]))
    event_raw <- suppressWarnings(as.numeric(pheno_df[[survival_event_var]]))
    Y <- survival::Surv(time_raw, event_raw)
  } else {
    Y <- pheno_df[[y_var]]
    if (!is.numeric(Y)) stop("y_var must be numeric.", call. = FALSE)
  }
  C <- make_covariate_matrix(pheno_df, covariates)

  data_by_key <- list()
  feature_meta <- data.frame(
    key = character(), stage = integer(), stage_name = character(),
    view = character(), mediator = character(), stringsAsFactors = FALSE
  )
  for (stage_i in seq_along(stage_views)) {
    for (view in stage_views[[stage_i]]) {
      M <- blocks[[view]]
      keys <- zentangler_feature_key(view, colnames(M))
      for (j in seq_along(keys)) data_by_key[[keys[j]]] <- M[, j]
      feature_meta <- rbind(
        feature_meta,
        data.frame(
          key = keys,
          stage = stage_i,
          stage_name = names(stage_views)[stage_i],
          view = view,
          mediator = colnames(M),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  layer1_views <- stage_views[[1L]]
  first_screens <- lapply(layer1_views, function(view) {
    hima_screen_mediators(
      M = blocks[[view]],
      X = X,
      C = C,
      sis_n = sis_n,
      rank_by = sis_rank,
      screen_method = screen_method,
      fdr_method = fdr_method
    )
  })
  names(first_screens) <- layer1_views

  first_selected <- unlist(lapply(names(first_screens), function(view) {
    zentangler_feature_key(view, first_screens[[view]]$selected)
  }), use.names = FALSE)
  if (length(first_selected) == 0L) stop("No first-layer mediators were selected.", call. = FALSE)

  a_lookup <- do.call(rbind, lapply(names(first_screens), function(view) {
    tab <- first_screens[[view]]$table
    data.frame(
      key = zentangler_feature_key(view, tab$mediator),
      a = tab$a_hat,
      p_a = tab$p_a,
      q_a = tab$q_a,
      selected = tab$selected,
      stringsAsFactors = FALSE
    )
  }))

  transition_links <- list()
  current_keys <- first_selected
  paths <- lapply(first_selected, function(key) list(keys = key, link_rows = integer(0)))

  for (stage_i in seq_len(length(stage_views) - 1L)) {
    to_views <- stage_views[[stage_i + 1L]]
    to_keys <- feature_meta$key[feature_meta$stage == stage_i + 1L]
    from_df <- as.data.frame(data_by_key[current_keys], check.names = FALSE)
    to_df <- as.data.frame(data_by_key[to_keys], check.names = FALSE)
    rownames(from_df) <- rownames(pheno_df)
    rownames(to_df) <- rownames(pheno_df)

    links <- zentangler_screen_cor_links(
      from_df = from_df,
      to_df = to_df,
      X = X,
      C = C,
      cor_method = cor_method,
      cor_fdr_method = cor_fdr_method,
      min_abs_cor = min_abs_cor,
      cor_q_threshold = cor_q_threshold,
      residualize_links = residualize_links,
      max_links = max_links_per_transition
    )
    if (nrow(links) == 0L) {
      warning("No correlation-screened links retained for transition ", stage_i, " -> ", stage_i + 1L, ".")
      transition_links[[paste0("stage", stage_i, "_to_stage", stage_i + 1L)]] <- links
      paths <- list()
      break
    }

    from_meta <- feature_meta[match(links$from_key, feature_meta$key), , drop = FALSE]
    to_meta <- feature_meta[match(links$to_key, feature_meta$key), , drop = FALSE]
    links$transition <- stage_i
    links$from_stage <- stage_i
    links$to_stage <- stage_i + 1L
    links$from_view <- from_meta$view
    links$from_mediator <- from_meta$mediator
    links$to_view <- to_meta$view
    links$to_mediator <- to_meta$mediator
    links <- zentangler_estimate_link_effects(links, data_by_key = data_by_key, X = X, C = C)
    transition_links[[paste0("stage", stage_i, "_to_stage", stage_i + 1L)]] <- links

    paths <- zentangler_extend_paths(paths, links, max_paths = max_paths)
    current_keys <- unique(vapply(paths, function(path) tail(path$keys, 1L), character(1)))
    if (length(paths) == 0L || length(current_keys) == 0L) break
  }

  if (length(paths) == 0L) {
    sequential_paths <- data.frame()
    terminal_effects <- data.frame()
  } else {
    terminal_keys <- unique(vapply(paths, function(path) tail(path$keys, 1L), character(1)))
    terminal_df <- as.data.frame(data_by_key[terminal_keys], check.names = FALSE)
    rownames(terminal_df) <- rownames(pheno_df)
    terminal_effects <- zentangler_fit_terminal_effects(
      Y = Y,
      X = X,
      terminal_df = terminal_df,
      C = C,
      y_family = y_family,
      outcome_model = outcome_model
    )

    link_all <- do.call(rbind, transition_links)
    path_rows <- vector("list", length(paths))
    for (i in seq_along(paths)) {
      keys <- paths[[i]]$keys
      terminal_key <- tail(keys, 1L)
      link_hits <- lapply(seq_len(length(keys) - 1L), function(j) {
        hit <- which(link_all$from_key == keys[j] & link_all$to_key == keys[j + 1L])
        if (length(hit) == 0L) NA_integer_ else hit[[1L]]
      })
      link_hits <- unlist(link_hits)
      link_dat <- link_all[link_hits, , drop = FALSE]
      a_row <- a_lookup[match(keys[[1L]], a_lookup$key), , drop = FALSE]
      b_row <- terminal_effects[match(terminal_key, terminal_effects$terminal_key), , drop = FALSE]
      path_meta <- feature_meta[match(keys, feature_meta$key), , drop = FALSE]
      d_vals <- link_dat$d_hat
      p_vals <- c(a_row$p_a, link_dat$p_d, b_row$p_b)
      score <- as.numeric(a_row$a) * prod(d_vals) * as.numeric(b_row$b)
      path_rows[[i]] <- data.frame(
        path_id = i,
        n_layers = length(keys),
        mediator_path = paste(paste(path_meta$view, path_meta$mediator, sep = ":"), collapse = " -> "),
        key_path = paste(keys, collapse = " -> "),
        first_view = path_meta$view[[1L]],
        first_mediator = path_meta$mediator[[1L]],
        terminal_view = tail(path_meta$view, 1L),
        terminal_mediator = tail(path_meta$mediator, 1L),
        a = as.numeric(a_row$a),
        d_product = prod(d_vals),
        b = as.numeric(b_row$b),
        sequential_score = score,
        abs_sequential_score = abs(score),
        p_a = as.numeric(a_row$p_a),
        max_p_d = if (any(is.finite(link_dat$p_d))) max(link_dat$p_d, na.rm = TRUE) else NA_real_,
        p_b = as.numeric(b_row$p_b),
        p_primary = if (any(is.finite(p_vals))) max(p_vals, na.rm = TRUE) else NA_real_,
        min_abs_cor = min(link_dat$abs_cor, na.rm = TRUE),
        max_q_cor = max(link_dat$q_cor, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    sequential_paths <- do.call(rbind, path_rows)
    sequential_paths$q_primary_global <- stats::p.adjust(sequential_paths$p_primary, method = fdr_method)
    sequential_paths$q_primary <- sequential_paths$q_primary_global
    sequential_paths <- sequential_paths[order(sequential_paths$abs_sequential_score, decreasing = TRUE), , drop = FALSE]
    rownames(sequential_paths) <- NULL
  }

  end_time <- Sys.time()
  out <- list(
    settings = list(
      input_container = "MultiAssayExperiment",
      model = "sequential",
      x_var = x_var,
      y_var = y_var,
      y_family = y_family,
      survival_time_var = survival_time_var,
      survival_event_var = survival_event_var,
      stage_views = stage_views,
      covariates = covariates,
      sis_n = sis_n,
      sis_rank = sis_rank,
      screen_method = screen_method,
      cor_method = cor_method,
      min_abs_cor = min_abs_cor,
      cor_q_threshold = cor_q_threshold,
      cor_fdr_method = cor_fdr_method,
      residualize_links = residualize_links,
      outcome_model = outcome_model,
      fdr_method = fdr_method,
      fdr_scope = fdr_scope
    ),
    diagnostics = list(
      started = start_time,
      finished = end_time,
      runtime_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
      n_samples = nrow(pheno_df),
      n_layers = length(stage_views),
      n_paths = nrow(sequential_paths)
    ),
    sample_ids = rownames(pheno_df),
    feature_metadata = feature_meta,
    first_layer_screens = first_screens,
    transition_links = transition_links,
    terminal_effects = terminal_effects,
    sequential_paths = sequential_paths,
    paths = sequential_paths
  )
  class(out) <- c("zentangler_sequential_fit", "list")
  out
}

#' Extract sequential Zentangler paths
#'
#' @param fit Output from `fit_sequential_zentangler()`.
#' @return Data frame of sequential paths.
#' @export
zentangler_sequential_paths <- function(fit) {
  out <- fit$sequential_paths %||% fit$paths
  if (is.null(out)) return(data.frame())
  rownames(out) <- NULL
  out
}
