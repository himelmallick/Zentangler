# -----------------------------------------------------------------------------
# K-layer sequential mediation
#
# Route X -> M1 -> M2 -> ... -> Mk -> Y over k ordered views.
#
# 1. Layer-1 A-stage: M1_j ~ X + C for every layer-1 feature; SIS keeps the
#    top sis_n.
# 2. Transition screens: residual correlations (on X and C) between retained
#    layer-s features and every layer-(s+1) feature. BH/BY is applied over all
#    p_s x p_(s+1) tested pairs before the |r| and q filters are applied.
# 3. Link coefficients d: each retained layer-(s+1) feature is regressed on X,
#    C, all of its retained parents and all retained upstream ancestors
#    (sequential-ignorability mediator model). OLS is used when the model is
#    low-dimensional; otherwise a de-biased lasso.
# 4. Terminal B-stage: Y ~ X + C + terminal features + all retained upstream
#    features. The terminal b is the conditional effect of the terminal
#    mediator given X and the upstream route, de-biased (b_inference =
#    "debiased") or bootstrapped (b_inference = "bootstrap").
# 5. Path score a * prod(d) * b and joint-significance p = max(p_a, p_d, p_b),
#    with BH/BY over all enumerated paths.
# -----------------------------------------------------------------------------

zentangler_feature_key <- function(view, mediator) {
  paste(view, mediator, sep = "::")
}

zentangler_resolve_route <- function(route, mae) {
  # Map route entries given either as MAE experiment names or as sanitized
  # view names back to MAE experiment names.
  exp_names <- names(MultiAssayExperiment::experiments(mae))
  safe <- sanitize_view_names(exp_names)
  route <- as.character(route)
  route <- route[!is.na(route) & nzchar(route)]
  out <- vapply(route, function(v) {
    if (v %in% exp_names) return(v)
    hit <- match(v, safe)
    if (!is.na(hit)) return(exp_names[hit])
    NA_character_
  }, character(1))
  if (anyNA(out)) {
    stop("route contains views not present in the MAE: ", paste(route[is.na(out)], collapse = ", "), call. = FALSE)
  }
  unname(out)
}

zentangler_check_route <- function(route) {
  if (is.list(route)) {
    stop("route is now a character vector, e.g. route = c('species', 'metabolites') (path_templates was removed).", call. = FALSE)
  }
  route <- as.character(route)
  route <- route[!is.na(route) & nzchar(route)]
  if (length(route) < 2L) {
    stop(
      "fit_sequential_zentangler() requires at least two modalities in its route. ",
      "Use fit_multiview_parallel_zentangler() for pure parallel mediation.",
      call. = FALSE
    )
  }
  if (anyDuplicated(route) > 0L) {
    stop("A sequential route cannot visit the same view twice: ", paste(route, collapse = " -> "), call. = FALSE)
  }
  invisible(route)
}

zentangler_residualize_columns <- function(M, Z = NULL) {
  M <- as.matrix(M)
  if (is.null(Z) || ncol(Z) == 0L) return(sweep(M, 2, colMeans(M), "-"))
  W <- cbind(`(Intercept)` = 1, as.matrix(Z))
  M - W %*% qr.coef(qr(W), M)
}

zentangler_standardize_for_cor <- function(M, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  M <- as.matrix(M)
  dn <- dimnames(M)
  n0 <- nrow(M)
  if (identical(method, "spearman")) {
    M <- matrix(apply(M, 2, rank, ties.method = "average", na.last = "keep"), nrow = n0)
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
  from_mat,
  to_mat,
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
  empty <- data.frame(from_key = character(0), to_key = character(0), cor = numeric(0), abs_cor = numeric(0), p_cor = numeric(0), q_cor = numeric(0), stringsAsFactors = FALSE)
  if (ncol(from_mat) == 0L || ncol(to_mat) == 0L) return(list(links = empty, n_tests = 0L))

  n <- nrow(from_mat)
  Z <- NULL
  if (isTRUE(residualize_links)) {
    Z <- cbind(X = as.numeric(X))
    if (!is.null(C) && ncol(C) > 0L) Z <- cbind(Z, as.matrix(C))
  }
  # Rank first (Spearman), then residualize the ranks on X and C.
  A <- if (identical(cor_method, "spearman")) zentangler_standardize_for_cor(from_mat, "spearman") else as.matrix(from_mat)
  B <- if (identical(cor_method, "spearman")) zentangler_standardize_for_cor(to_mat, "spearman") else as.matrix(to_mat)
  A <- zentangler_standardize_for_cor(zentangler_residualize_columns(A, Z), "pearson")
  B <- zentangler_standardize_for_cor(zentangler_residualize_columns(B, Z), "pearson")
  keep_a <- colSums(!is.finite(A)) == 0L
  keep_b <- colSums(!is.finite(B)) == 0L
  A <- A[, keep_a, drop = FALSE]
  B <- B[, keep_b, drop = FALSE]
  if (ncol(A) == 0L || ncol(B) == 0L) return(list(links = empty, n_tests = 0L))

  R <- crossprod(A, B) / (n - 1)
  R[R > 1] <- 1
  R[R < -1] <- -1
  df <- n - 2L - (if (is.null(Z)) 0L else ncol(Z))
  if (df < 2L) stop("Too few samples for transition correlation tests.", call. = FALSE)
  tstat <- R * sqrt(df / pmax(1e-12, 1 - R^2))
  P <- 2 * stats::pt(abs(tstat), df = df, lower.tail = FALSE)
  # FDR over every tested pair, not only those passing the |r| filter.
  Q <- matrix(stats::p.adjust(as.numeric(P), method = cor_fdr_method), nrow = nrow(P), dimnames = dimnames(P))
  idx <- which(abs(R) >= min_abs_cor & Q <= cor_q_threshold, arr.ind = TRUE)
  n_tests <- length(P)
  if (nrow(idx) == 0L) return(list(links = empty, n_tests = n_tests))
  out <- data.frame(
    from_key = colnames(A)[idx[, 1]],
    to_key = colnames(B)[idx[, 2]],
    cor = R[idx],
    abs_cor = abs(R[idx]),
    p_cor = P[idx],
    q_cor = Q[idx],
    stringsAsFactors = FALSE
  )
  out <- out[order(out$q_cor, -out$abs_cor), , drop = FALSE]
  if (nrow(out) > max_links) {
    warning(
      nrow(out), " links passed the transition screen; keeping the ", max_links,
      " strongest (max_links_per_transition).",
      call. = FALSE
    )
    out <- out[seq_len(max_links), , drop = FALSE]
  }
  rownames(out) <- NULL
  list(links = out, n_tests = n_tests)
}

zentangler_fit_node <- function(y, W, targets, adjust = NULL, mode = c("debiased", "none"), lambda_choice = "lambda.min") {
  # Mediator model y ~ W (unpenalized: X, C) + targets + adjust. Returns
  # estimates and p-values for the target columns. OLS when low-dimensional,
  # otherwise de-biased lasso (mode = "debiased") or plain lasso estimates
  # without p-values (mode = "none", used inside bootstrap replicates).
  mode <- match.arg(mode)
  P <- cbind(targets, adjust)
  n <- length(y)
  k <- ncol(W) + ncol(P) + 1L
  tn <- colnames(targets)
  if (k <= n - 10L && k <= n / 2) {
    D <- cbind(`(Intercept)` = 1, W, P)
    fit <- stats::lm.fit(D, y)
    co <- fit$coefficients
    df <- n - fit$rank
    if (fit$rank < ncol(D) || df <= 0) {
      return(data.frame(name = tn, estimate = NA_real_, p = NA_real_, stringsAsFactors = FALSE))
    }
    sigma2 <- sum(fit$residuals^2) / df
    xtx_inv <- chol2inv(qr.R(fit$qr))
    piv <- fit$qr$pivot
    se <- rep(NA_real_, ncol(D))
    se[piv] <- sqrt(diag(xtx_inv) * sigma2)
    pos <- match(tn, colnames(D))
    est <- co[pos]
    p <- 2 * stats::pt(abs(est / se[pos]), df = df, lower.tail = FALSE)
    return(data.frame(name = tn, estimate = as.numeric(est), p = as.numeric(p), stringsAsFactors = FALSE))
  }
  if (identical(mode, "debiased")) {
    db <- zentangler_debiased_lasso(Y = y, W = W, P = P, targets = tn, y_family = "gaussian", lambda_choice = lambda_choice, max_targets = length(tn))
    out <- data.frame(name = tn, estimate = NA_real_, p = NA_real_, stringsAsFactors = FALSE)
    hit <- match(db$table$name, tn)
    out$estimate[hit] <- db$table$estimate
    out$p[hit] <- db$table$p
    return(out)
  }
  Zp <- cbind(W, P)
  pf <- c(rep(0, ncol(W)), rep(1, ncol(P)))
  cv <- glmnet::cv.glmnet(Zp, y, alpha = 1, penalty.factor = pf, foldid = zentangler_make_foldid(n))
  co <- as.matrix(stats::coef(cv, s = lambda_choice))[, 1]
  data.frame(name = tn, estimate = as.numeric(co[tn]), p = NA_real_, stringsAsFactors = FALSE)
}

zentangler_extend_paths <- function(paths, links, max_paths = 10000L) {
  if (length(paths) == 0L || is.null(links) || nrow(links) == 0L) return(list(paths = list(), truncated = FALSE))
  by_from <- split(seq_len(nrow(links)), links$from_key)
  out <- list()
  ii <- 0L
  for (path in paths) {
    last_key <- utils::tail(path, 1L)
    hit <- by_from[[last_key]]
    if (is.null(hit)) next
    for (j in hit) {
      ii <- ii + 1L
      out[[ii]] <- c(path, links$to_key[j])
      if (ii >= max_paths) return(list(paths = out, truncated = TRUE))
    }
  }
  list(paths = out, truncated = FALSE)
}

zentangler_sequential_core <- function(
  data_by_view,
  route,
  X,
  Y,
  C = NULL,
  y_family = "gaussian",
  sis_n = 50L,
  sis_rank = "abs_a",
  screen_method = "sis",
  cor_method = "spearman",
  cor_fdr_method = "BH",
  min_abs_cor = 0.3,
  cor_q_threshold = 0.25,
  residualize_links = TRUE,
  max_links_per_transition = 5000L,
  max_paths = 10000L,
  lambda_choice = "lambda.min",
  glmnet_alpha = 1,
  engine = "multiview",
  rho = NULL,
  rho_grid = zentangler_default_rho_grid(),
  b_inference = "debiased",
  debias_max_targets = 200L,
  fdr_method = "BH",
  pheno_df = NULL,
  x_var = NULL,
  covariates = NULL,
  a_stage = NULL,
  verbose = FALSE
) {
  # Route discovery + estimation on one dataset. b_inference "none" skips all
  # p-value machinery (bootstrap replicates).
  mode <- if (identical(b_inference, "debiased")) "debiased" else "none"
  n <- length(X)
  warn <- character(0)
  key_mat <- function(view, cols = NULL) {
    M <- data_by_view[[view]]
    if (!is.null(cols)) M <- M[, cols, drop = FALSE]
    colnames(M) <- zentangler_feature_key(view, colnames(M))
    M
  }
  all_mat <- do.call(cbind, lapply(route, key_mat))
  W <- zentangler_y_design(X, C)

  # 1. Layer-1 screen.
  v1 <- route[[1L]]
  # Only layer 1 has an A-path; its model may be "lm" or "maaslin2".
  screen1 <- zentangler_a_stage_screen(
    view = v1, M = data_by_view[[v1]], X = X, C = C, pheno_df = pheno_df, x_var = x_var, covariates = covariates,
    model = a_stage$model[[v1]] %||% "lm", maaslin2_args = a_stage$maaslin2_args[[v1]] %||% list(),
    sis_n = sis_n, rank_by = sis_rank, screen_method = screen_method, fdr_method = fdr_method
  )
  a_lookup <- data.frame(
    key = zentangler_feature_key(v1, screen1$table$mediator),
    a = screen1$table$a_hat,
    p_a = screen1$table$p_a,
    q_a = screen1$table$q_a,
    selected = screen1$table$selected,
    stringsAsFactors = FALSE
  )
  zentangler_log(verbose, sprintf("Layer 1 (%s) [%s]: %d features, %d with p_a < 0.05, %d screened in",
                                  v1, a_stage$model[[v1]] %||% "lm", nrow(screen1$table),
                                  sum(screen1$table$p_a < 0.05, na.rm = TRUE), length(screen1$selected)))
  current <- zentangler_feature_key(v1, screen1$selected)
  ancestors <- setNames(lapply(current, function(k) character(0)), current)
  paths <- as.list(current)
  link_tables <- list()
  n_link_tests <- 0L
  paths_truncated <- FALSE

  # 2-3. Transitions.
  for (s in seq_len(length(route) - 1L)) {
    to_view <- route[[s + 1L]]
    scr <- zentangler_screen_cor_links(
      from_mat = all_mat[, current, drop = FALSE],
      to_mat = key_mat(to_view),
      X = X, C = C,
      cor_method = cor_method, cor_fdr_method = cor_fdr_method,
      min_abs_cor = min_abs_cor, cor_q_threshold = cor_q_threshold,
      residualize_links = residualize_links, max_links = max_links_per_transition
    )
    links <- scr$links
    n_link_tests <- n_link_tests + scr$n_tests
    zentangler_log(verbose, sprintf("Transition %d (%s -> %s): %d pairs tested, %d links kept (|r| >= %g, q <= %g)",
                                    s, route[[s]], to_view, scr$n_tests, nrow(links), min_abs_cor, cor_q_threshold))
    tr_name <- paste0("stage", s, "_to_stage", s + 1L)
    if (nrow(links) == 0L) {
      warn <- c(warn, paste0("No links passed the transition screen ", s, " -> ", s + 1L, "."))
      link_tables[[tr_name]] <- links
      paths <- list()
      break
    }
    links$transition <- s
    links$from_stage <- s
    links$to_stage <- s + 1L
    links$d_hat <- NA_real_
    links$p_d <- NA_real_
    new_anc <- list()
    for (to_key in unique(links$to_key)) {
      parents <- unique(links$from_key[links$to_key == to_key])
      anc <- unique(unlist(ancestors[parents], use.names = FALSE))
      anc <- setdiff(anc, parents)
      fit_d <- try(zentangler_fit_node(
        y = as.numeric(all_mat[, to_key]),
        W = W,
        targets = all_mat[, parents, drop = FALSE],
        adjust = if (length(anc) > 0L) all_mat[, anc, drop = FALSE] else NULL,
        mode = mode,
        lambda_choice = lambda_choice
      ), silent = TRUE)
      if (!inherits(fit_d, "try-error")) {
        rows <- which(links$to_key == to_key)
        hit <- match(links$from_key[rows], fit_d$name)
        links$d_hat[rows] <- fit_d$estimate[hit]
        links$p_d[rows] <- fit_d$p[hit]
      }
      new_anc[[to_key]] <- unique(c(parents, anc))
    }
    link_tables[[tr_name]] <- links
    ext <- zentangler_extend_paths(paths, links, max_paths = max_paths)
    paths <- ext$paths
    if (isTRUE(ext$truncated)) {
      paths_truncated <- TRUE
      warn <- c(warn, paste0("Path enumeration hit max_paths = ", max_paths, " at transition ", s, "; remaining paths were not enumerated."))
    }
    current <- unique(vapply(paths, function(p) utils::tail(p, 1L), character(1)))
    ancestors <- new_anc[current]
    zentangler_log(verbose, sprintf("%d paths so far, %d distinct %s mediators", length(paths), length(current), to_view), indent = 1L)
    if (length(paths) == 0L) break
  }

  empty_out <- list(
    first_layer_screen = screen1, a_lookup = a_lookup, transition_links = link_tables,
    terminal_effects = data.frame(), sequential_paths = data.frame(),
    n_link_tests = n_link_tests, paths_truncated = paths_truncated, warnings = warn,
    n_first_layer_selected = length(screen1$selected)
  )
  if (length(paths) == 0L) return(empty_out)

  # 4. Terminal model conditional on X, C and all retained upstream features.
  # The retained mediators of each route layer form one view, so the terminal
  # model is the same cooperative-learning (multiview) fit as the parallel API.
  terminal_keys <- unique(vapply(paths, function(p) utils::tail(p, 1L), character(1)))
  upstream_keys <- setdiff(unique(unlist(lapply(paths, function(p) utils::head(p, -1L)), use.names = FALSE)), terminal_keys)
  T_mat <- all_mat[, terminal_keys, drop = FALSE]
  U_mat <- if (length(upstream_keys) > 0L) all_mat[, upstream_keys, drop = FALSE] else NULL
  key_view <- function(k) vapply(strsplit(k, "::", fixed = TRUE), `[`, character(1), 1L)
  model_keys <- c(terminal_keys, upstream_keys)
  term_blocks <- lapply(route, function(v) {
    ks <- model_keys[key_view(model_keys) == v]
    if (length(ks) == 0L) return(NULL)
    all_mat[, ks, drop = FALSE]
  })
  names(term_blocks) <- route
  term_blocks <- term_blocks[!vapply(term_blocks, is.null, logical(1))]
  zentangler_log(verbose, sprintf("Terminal outcome model: %d terminal + %d upstream mediators, %s ...",
                                  length(terminal_keys), length(upstream_keys),
                                  if (identical(engine, "multiview")) "multiview (tuning rho by CV)" else paste0("IntegratedLearner-style ", sub("il_", "", engine), " fusion")))
  yfit <- fit_y_multiview_stage(
    Y = Y, X = X, blocks = term_blocks, C = C, engine = engine, y_family = y_family,
    lambda_choice = lambda_choice, glmnet_alpha = glmnet_alpha,
    rho = rho, rho_grid = rho_grid
  )
  if (!is.null(yfit$rho)) zentangler_log(verbose, sprintf("rho chosen = %g", yfit$rho), indent = 1L)
  coef_all <- unlist(unname(yfit$coef_by_view))
  b_fusion <- setNames(as.numeric(coef_all[terminal_keys]), terminal_keys)
  b_fusion[!is.finite(b_fusion)] <- 0
  terminal_effects <- data.frame(
    terminal_key = terminal_keys,
    b = b_fusion,
    se_b = NA_real_,
    p_b = NA_real_,
    b_fusion = b_fusion,
    n_upstream_adjusted = length(upstream_keys),
    stringsAsFactors = FALSE
  )
  b_method <- "none"
  if (identical(mode, "debiased")) {
    init <- c(yfit$w_coef, coef_all)
    db <- try(zentangler_debiased_lasso(
      Y = Y, W = W, P = cbind(T_mat, U_mat), targets = terminal_keys, y_family = y_family,
      lambda_choice = lambda_choice, max_targets = debias_max_targets, beta_init = init
    ), silent = TRUE)
    terminal_effects$b <- NA_real_
    if (inherits(db, "try-error")) {
      warn <- c(warn, paste0("Terminal de-biased inference failed: ", conditionMessage(attr(db, "condition"))))
      b_method <- "debiased_failed"
    } else {
      hit <- match(db$table$name, terminal_keys)
      terminal_effects$b[hit] <- db$table$estimate
      terminal_effects$se_b[hit] <- db$table$se
      terminal_effects$p_b[hit] <- db$table$p
      b_method <- db$method
      if (length(db$dropped) > 0L) warn <- c(warn, paste0(length(db$dropped), " terminal mediators exceeded debias_max_targets; they get p = 1."))
    }
    terminal_effects$b[!is.finite(terminal_effects$b)] <- 0
  }
  terminal_effects$b_inference_method <- b_method

  # 5. Path table.
  link_all <- do.call(rbind, link_tables)
  link_id <- paste(link_all$from_key, link_all$to_key, sep = " => ")
  path_rows <- lapply(seq_along(paths), function(i) {
    keys <- paths[[i]]
    lid <- paste(utils::head(keys, -1L), keys[-1L], sep = " => ")
    ld <- link_all[match(lid, link_id), , drop = FALSE]
    a_row <- a_lookup[match(keys[[1L]], a_lookup$key), , drop = FALSE]
    t_row <- terminal_effects[match(utils::tail(keys, 1L), terminal_effects$terminal_key), , drop = FALSE]
    parts <- strsplit(keys, "::", fixed = TRUE)
    views <- vapply(parts, `[`, character(1), 1L)
    meds <- vapply(parts, function(x) paste(x[-1L], collapse = "::"), character(1))
    d_prod <- prod(ld$d_hat)
    score <- as.numeric(a_row$a) * d_prod * as.numeric(t_row$b)
    pv <- c(a_row$p_a, ld$p_d, t_row$p_b)
    data.frame(
      path_id = i,
      n_layers = length(keys),
      route_type = "sequential",
      mediator_path = paste(paste(views, meds, sep = ":"), collapse = " -> "),
      key_path = paste(keys, collapse = " -> "),
      first_view = views[[1L]],
      first_mediator = meds[[1L]],
      terminal_view = utils::tail(views, 1L),
      terminal_mediator = utils::tail(meds, 1L),
      terminal_key = utils::tail(keys, 1L),
      a = as.numeric(a_row$a),
      d_product = d_prod,
      b = as.numeric(t_row$b),
      b_fusion = as.numeric(t_row$b_fusion),
      sequential_score = score,
      abs_sequential_score = abs(score),
      p_a = as.numeric(a_row$p_a),
      max_p_d = if (any(is.finite(ld$p_d))) max(ld$p_d, na.rm = TRUE) else NA_real_,
      p_b = as.numeric(t_row$p_b),
      p_primary = if (all(is.finite(pv))) max(pv) else 1,
      min_abs_cor = min(ld$abs_cor),
      max_q_cor = max(ld$q_cor),
      stringsAsFactors = FALSE
    )
  })
  sequential_paths <- do.call(rbind, path_rows)
  sequential_paths$p_primary_model_based <- sequential_paths$p_primary
  sequential_paths$q_primary_global <- stats::p.adjust(sequential_paths$p_primary, method = fdr_method)
  sequential_paths$q_primary <- sequential_paths$q_primary_global

  out <- empty_out
  out$terminal_effects <- terminal_effects
  out$sequential_paths <- sequential_paths
  out$upstream_keys <- upstream_keys
  out$terminal_rho <- yfit$rho %||% NA_real_
  out$terminal_model_method <- yfit$method
  out
}

#' Fit one k-layer sequential Zentangler mediation route
#'
#' Fits `X -> M1 -> M2 -> ... -> Mk -> Y` over an ordered route of views. The
#' first layer is screened on its A-path; adjacent layers are linked by
#' residual-correlation screens (FDR over every tested pair); link
#' coefficients come from mediator models that condition on the exposure,
#' covariates, all retained parents and upstream ancestors; and the terminal
#' outcome model conditions on the exposure, covariates and every retained
#' upstream mediator, fitted with multiview cooperative learning (each route
#' layer is one view; `rho` in \[0, 1\] tuned by cross-validation).
#'
#' @param mae A `MultiAssayExperiment`.
#' @param x_var Numeric exposure column in `colData(mae)`.
#' @param y_var Outcome column (Gaussian / 0-1 binomial outcomes).
#' @param route Character vector of at least two distinct views, in order.
#' @param covariates Optional covariate columns.
#' @param y_family `"gaussian"`, `"binomial"` or `"survival"`.
#' @param survival_time_var,survival_event_var Survival columns.
#' @param assay_names Optional assay names per view.
#' @param transform,pseudocount Per-view `"none"` / `"log"` / `"clr"`
#'   transformation and pseudocount (see [fit_multiview_parallel_zentangler()]).
#' @param duplicate_primary,fill_nonfinite_zero MAE extraction options.
#' @param sis_n,sis_rank,screen_method First-layer screen.
#' @param a_stage_model,maaslin2_args A-path model for the first layer (the
#'   only layer with an exposure -> mediator path): `"lm"` or `"maaslin2"`,
#'   given once or named by view as in [fit_multiview_parallel_zentangler()];
#'   entries for other route views are accepted and ignored.
#' @param cor_method,min_abs_cor,cor_q_threshold,cor_fdr_method,residualize_links
#'   Transition-link screen.
#' @param max_links_per_transition,max_paths Caps (warned when hit).
#' @param fusion Terminal outcome model: `"multiview"` (default; each route
#'   layer is one view) or the supplementary `"integrated_learner"`.
#' @param rho,rho_grid Fixed multiview `rho` in \[0, 1\], or the CV grid.
#' @param integrated_learner_fusion `"late"` or `"early"` for
#'   `fusion = "integrated_learner"`.
#' @param lambda_choice,glmnet_alpha Penalty choices.
#' @param b_inference `"debiased"` or `"bootstrap"`.
#' @param debias_max_targets Maximum terminal mediators to de-bias.
#' @param path_inference `"model_based"` (max of path-leg p-values) or
#'   `"bootstrap_score"` (pipeline bootstrap of the path score).
#' @param fdr_method `"BH"` or `"BY"` for path q-values.
#' @param causal_inference `"none"` or `"bootstrap"`.
#' @param effect_q_threshold Paths at or below this q-value form the effect set.
#' @param effect_x0,effect_x1 Exposure contrast for causal effects.
#' @param bootstrap_repeats,bootstrap_ci_level,bootstrap_seed,bootstrap_id Bootstrap settings.
#' @param seed Random seed (the caller's RNG state is restored on exit).
#' @param verbose Print timestamped progress messages.
#'
#' @return An object of class `zentangler_sequential_fit`: `settings`,
#'   `diagnostics`, `paths`, `edges` (retained transition links), `terminals`
#'   (terminal B-stage), `first_layer` (layer-1 screen), `effects` and
#'   `bootstrap`. Use [zentangler_paths()], [zentangler_effects()] and `summary()`.
#' @export
fit_sequential_zentangler <- function(
  mae,
  x_var,
  y_var = NULL,
  route,
  covariates = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  survival_time_var = NULL,
  survival_event_var = NULL,
  assay_names = NULL,
  transform = "none",
  pseudocount = NULL,
  duplicate_primary = c("mean", "first"),
  fill_nonfinite_zero = FALSE,
  a_stage_model = "lm",
  maaslin2_args = list(),
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
  fusion = c("multiview", "integrated_learner"),
  rho = NULL,
  rho_grid = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  integrated_learner_fusion = c("late", "early"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  b_inference = c("debiased", "bootstrap"),
  debias_max_targets = 200L,
  path_inference = c("model_based", "bootstrap_score"),
  fdr_method = c("BH", "BY"),
  causal_inference = c("none", "bootstrap"),
  effect_q_threshold = 0.25,
  effect_x0 = NULL,
  effect_x1 = NULL,
  bootstrap_repeats = 0L,
  bootstrap_ci_level = 0.95,
  bootstrap_seed = NULL,
  bootstrap_id = NULL,
  seed = 1,
  verbose = TRUE
) {
  start_time <- Sys.time()
  restore_rng <- zentangler_save_rng()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)
  y_family <- match.arg(y_family)
  fusion <- match.arg(fusion)
  integrated_learner_fusion <- match.arg(integrated_learner_fusion)
  engine <- zentangler_resolve_fusion(fusion, integrated_learner_fusion)
  b_inference <- zentangler_validate_b_inference(b_inference)
  if (identical(b_inference, "none")) stop("Use b_inference = 'debiased' or 'bootstrap'.", call. = FALSE)
  path_inference <- match.arg(path_inference)
  causal_inference <- match.arg(causal_inference)
  fdr_method <- validate_fdr_method(fdr_method)
  bootstrap_repeats <- max(0L, as.integer(bootstrap_repeats))
  bootstrap_seed <- bootstrap_seed %||% (as.integer(seed) + 10000L)
  needs_boot <- c(b_inference = identical(b_inference, "bootstrap"), path_inference = identical(path_inference, "bootstrap_score"), causal_inference = identical(causal_inference, "bootstrap"))
  if (any(needs_boot) && bootstrap_repeats < 1L) {
    stop(paste(names(needs_boot)[needs_boot], collapse = ", "), " = bootstrap requires bootstrap_repeats > 0.", call. = FALSE)
  }
  x_var <- zentangler_require_string(x_var, "x_var")
  o <- list(y_family = y_family, y_var = y_var, survival_time_var = survival_time_var, survival_event_var = survival_event_var)
  if (identical(y_family, "survival")) {
    o$survival_time_var <- zentangler_require_string(survival_time_var, "survival_time_var")
    o$survival_event_var <- zentangler_require_string(survival_event_var, "survival_event_var")
    o$y_var <- NULL
  } else {
    o$y_var <- zentangler_require_string(y_var, "y_var")
    o$survival_time_var <- o$survival_event_var <- NULL
  }

  zentangler_check_route(route)
  zentangler_require_mae()
  route_exp <- zentangler_resolve_route(route, mae)
  route_assays <- if (!is.null(assay_names) && !is.null(names(assay_names))) assay_names[intersect(names(assay_names), route_exp)] else NULL
  inputs <- zentangler_mae_to_blocks(
    mae = mae, view_names = route_exp, assay_names = route_assays,
    duplicate_primary = match.arg(duplicate_primary), fill_nonfinite_zero = fill_nonfinite_zero,
    transform = transform, pseudocount = pseudocount
  )
  route <- inputs$view_map$view[match(route_exp, inputs$view_map$experiment_name)]
  a_stage <- zentangler_resolve_a_stage(a_stage_model, maaslin2_args, route, inputs$view_map)
  if (!identical(a_stage$model[[route[[1L]]]], "maaslin2")) a_stage$random_effects <- character(0)

  needed <- unique(c(x_var, o$y_var, o$survival_time_var, o$survival_event_var, covariates, bootstrap_id, a_stage$random_effects))
  al <- align_samples_multiview_blocks(inputs$blocks, inputs$pheno_df, needed_pheno_cols = needed)
  blocks <- al$blocks
  pheno_df <- al$pheno
  X <- pheno_df[[x_var]]
  if (!is.numeric(X)) stop("x_var must be numeric.", call. = FALSE)
  Y <- zentangler_make_outcome(pheno_df, o)
  C <- make_covariate_matrix(pheno_df, covariates)
  zentangler_log(verbose, sprintf("Sequential mediation: X -> %s -> Y", paste(route, collapse = " -> ")))
  tr <- inputs$view_map$transform
  if (any(tr != "none")) zentangler_log(verbose, "Transforms: ", paste(paste0(inputs$view_map$view, " = ", tr)[tr != "none"], collapse = ", "), indent = 1L)
  zentangler_log(verbose, sprintf("Data: %d of %d samples complete across route views; outcome %s%s.",
                                  nrow(pheno_df), nrow(inputs$pheno_df), y_family,
                                  if (length(covariates)) paste0("; covariates: ", paste(covariates, collapse = ", ")) else ""))
  contrast <- infer_effect_contrast(X, x0 = effect_x0, x1 = effect_x1)
  exposure_outcome <- zentangler_exposure_outcome_test(Y, X, C, y_family, contrast[["x0"]], contrast[["x1"]])
  zentangler_log_exposure_outcome(verbose, exposure_outcome)

  core_args <- list(
    route = route, y_family = y_family, sis_n = sis_n, sis_rank = match.arg(sis_rank),
    screen_method = match.arg(screen_method), cor_method = match.arg(cor_method),
    cor_fdr_method = validate_fdr_method(cor_fdr_method), min_abs_cor = min_abs_cor,
    cor_q_threshold = cor_q_threshold, residualize_links = residualize_links,
    max_links_per_transition = max_links_per_transition, max_paths = max_paths,
    lambda_choice = match.arg(lambda_choice), glmnet_alpha = validate_glmnet_alpha(glmnet_alpha),
    engine = engine, rho = rho, rho_grid = rho_grid,
    debias_max_targets = max(1L, as.integer(debias_max_targets)), fdr_method = fdr_method,
    x_var = x_var, covariates = covariates, a_stage = a_stage
  )
  core <- zentangler_muffle_cox_ties(do.call(zentangler_sequential_core, c(list(data_by_view = blocks, X = X, Y = Y, C = C, pheno_df = pheno_df, b_inference = b_inference, verbose = verbose), core_args)))
  for (w in core$warnings) warning(w, call. = FALSE)
  paths <- core$sequential_paths
  terminals <- core$terminal_effects
  zentangler_log(verbose, sprintf("%d complete paths; terminal B-stage: %s", nrow(paths),
                                  if (nrow(terminals)) terminals$b_inference_method[[1L]] else "none"))

  # Pipeline bootstrap: the whole route discovery is re-run per replicate.
  bootstrap <- list()
  if (bootstrap_repeats > 0L && nrow(paths) > 0L && (needs_boot[["b_inference"]] || needs_boot[["path_inference"]])) {
    if (is.finite(core$terminal_rho %||% NA_real_)) core_args$rho <- core$terminal_rho  # rho fixed in replicates
    zentangler_log(verbose, sprintf("Pipeline bootstrap: %d replicates ...", bootstrap_repeats))
    boot <- zentangler_bootstrap_sequential_pipeline(
      verbose = verbose,
      blocks = blocks, pheno_df = pheno_df, X = X, Y = Y, C = C, core_args = core_args,
      path_keys = paths$key_path, terminal_keys = terminals$terminal_key,
      repeats = bootstrap_repeats, seed = bootstrap_seed, bootstrap_id = bootstrap_id
    )
    paths <- zentangler_boot_summary_columns(paths, boot$score_matrix, "sequential_score", bootstrap_ci_level)
    paths$path_boot_selection_freq <- colMeans(boot$score_matrix != 0, na.rm = TRUE)
    paths$p_sequential_score_boot <- bootstrap_sign_pvalues(boot$score_matrix)
    terminals <- zentangler_boot_summary_columns(terminals, boot$b_matrix, "b", bootstrap_ci_level)
    terminals$p_b_bootstrap <- bootstrap_sign_pvalues(boot$b_matrix)
    if (identical(b_inference, "bootstrap")) {
      terminals$p_b <- terminals$p_b_bootstrap
      paths$p_b <- terminals$p_b[match(paths$terminal_key, terminals$terminal_key)]
      pv <- cbind(paths$p_a, paths$max_p_d, paths$p_b)
      paths$p_primary <- ifelse(rowSums(!is.finite(pv)) == 0, apply(pv, 1, max), 1)
    }
    if (identical(path_inference, "bootstrap_score")) {
      paths$p_primary <- paths$p_sequential_score_boot
      paths$p_primary[!is.finite(paths$p_primary)] <- 1
    }
    paths$q_primary <- stats::p.adjust(paths$p_primary, method = fdr_method)
    bootstrap$pipeline <- boot
  }
  paths$q_primary_global <- NULL
  paths$p_primary_model_based <- NULL
  if (verbose && nrow(paths) > 0L) {
    act <- vapply(c(0.05, 0.1, 0.25), function(q) sum(zentangler_is_active(paths, q)), integer(1))
    zentangler_log(TRUE, sprintf("Active paths: %d at q <= 0.05, %d at q <= 0.10, %d at q <= 0.25", act[1], act[2], act[3]), indent = 1L)
  }

  # Post-selection causal effects for the fixed set of significant paths.
  data_by_key <- list()
  for (view in route) {
    M <- blocks[[view]]
    for (j in seq_len(ncol(M))) data_by_key[[zentangler_feature_key(view, colnames(M)[j])]] <- M[, j]
  }
  effects <- data.frame()
  if (nrow(paths) > 0L) {
    effect_paths <- zentangler_sequential_effect_paths(paths, effect_q_threshold)
    effects <- compute_effect_decomposition_sequential(
      sequential_paths = effect_paths, Y = Y, X = X, data_by_key = data_by_key, C = C,
      y_family = y_family, x0 = contrast[["x0"]], x1 = contrast[["x1"]], n_layers = length(route)
    )
    zentangler_log(verbose, sprintf(
      "Causal effects (%d paths with q <= %g; %s, X %.3g -> %.3g): NDE %.4g, NIE %.4g, TE %.4g, proportion mediated %.3g",
      nrow(effect_paths), effect_q_threshold, gsub("_", " ", effects$effect_scale[[1L]]), contrast[["x0"]], contrast[["x1"]],
      effects$nde[[1L]], effects$nie_total[[1L]], effects$te[[1L]], effects$prop_mediated[[1L]]
    ))
    if (identical(causal_inference, "bootstrap")) {
      zentangler_log(verbose, sprintf("Causal-effect bootstrap: %d replicates (fixed path set) ...", bootstrap_repeats))
      eb <- zentangler_bootstrap_sequential_effects(
        verbose = verbose,
        effect_paths = effect_paths, data_by_key = data_by_key, pheno_df = pheno_df, X = X, Y = Y, C = C,
        y_family = y_family, x0 = contrast[["x0"]], x1 = contrast[["x1"]], n_layers = length(route),
        repeats = bootstrap_repeats, seed = bootstrap_seed, bootstrap_id = bootstrap_id
      )
      effects <- add_effect_bootstrap_columns(effects, eb$effect_decomposition, bootstrap_ci_level)
      bootstrap$effects <- eb
    }
    paths$in_effect_set <- paths$key_path %in% effect_paths$key_path
    paths <- paths[order(paths$abs_sequential_score, decreasing = TRUE), , drop = FALSE]
    rownames(paths) <- NULL
  }
  edges <- do.call(rbind, core$transition_links)
  if (is.null(edges)) edges <- data.frame()
  rownames(edges) <- NULL

  end_time <- Sys.time()
  out <- list(
    settings = c(core_args, list(
      x_var = x_var, y_var = o$y_var, survival_time_var = o$survival_time_var,
      survival_event_var = o$survival_event_var, covariates = covariates,
      a_stage_model_layer1 = a_stage$model[[route[[1L]]]],
      b_inference = b_inference,
      b_inference_method = if (nrow(terminals) > 0L) terminals$b_inference_method[[1L]] else NA_character_,
      fusion = fusion, integrated_learner_fusion = integrated_learner_fusion,
      rho_used = core$terminal_rho %||% NA_real_,
      terminal_model_method = core$terminal_model_method %||% NA_character_,
      path_inference = path_inference, causal_inference = causal_inference,
      effect_q_threshold = effect_q_threshold, effect_x0 = contrast[["x0"]], effect_x1 = contrast[["x1"]],
      bootstrap_repeats = bootstrap_repeats, bootstrap_ci_level = bootstrap_ci_level,
      bootstrap_seed = bootstrap_seed, bootstrap_id = bootstrap_id, seed = seed
    )),
    diagnostics = list(
      runtime_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
      n_samples = nrow(pheno_df),
      sample_ids = rownames(pheno_df),
      n_first_layer_selected = core$n_first_layer_selected,
      n_transition_link_tests = core$n_link_tests,
      n_transition_links = nrow(edges),
      paths_truncated = core$paths_truncated,
      n_upstream_adjusted_in_terminal = length(core$upstream_keys %||% character(0)),
      exposure_outcome = exposure_outcome,
      bootstrap_successful = c(pipeline = bootstrap$pipeline$repeats_successful %||% NA_integer_, effects = bootstrap$effects$repeats_successful %||% NA_integer_),
      bootstrap_failures = c(bootstrap$pipeline$failures, bootstrap$effects$failures)
    ),
    paths = paths,
    edges = edges,
    terminals = terminals,
    first_layer = core$first_layer_screen$table,
    effects = effects,
    bootstrap = if (length(bootstrap) > 0L) bootstrap else NULL
  )
  class(out) <- c("zentangler_sequential_fit", "list")
  zentangler_log(verbose, sprintf("Done in %.1f s.", out$diagnostics$runtime_seconds))
  out
}

zentangler_bootstrap_sequential_pipeline <- function(blocks, pheno_df, X, Y, C, core_args, path_keys, terminal_keys, repeats, seed, bootstrap_id = NULL, verbose = FALSE) {
  if (is.null(seed) || !is.finite(as.numeric(seed))) seed <- 10001L
  seed <- as.integer(seed)
  score_mat <- matrix(NA_real_, nrow = repeats, ncol = length(path_keys), dimnames = list(NULL, path_keys))
  b_mat <- matrix(NA_real_, nrow = repeats, ncol = length(terminal_keys), dimnames = list(NULL, terminal_keys))
  failures <- character(0)
  n_ok <- 0L
  for (b in seq_len(repeats)) {
    bt <- try(make_bootstrap_sample(blocks, pheno_df, bootstrap_id = bootstrap_id, seed = seed + b), silent = TRUE)
    if (inherits(bt, "try-error")) {
      failures <- c(failures, paste0("bootstrap_", b, ": ", conditionMessage(attr(bt, "condition"))))
      next
    }
    idx <- bt$idx
    args_b <- core_args
    for (v in names(args_b$a_stage$maaslin2_args)) {
      base_out <- args_b$a_stage$maaslin2_args[[v]]$output
      if (!is.null(base_out) && nzchar(base_out)) args_b$a_stage$maaslin2_args[[v]]$output <- file.path(base_out, paste0("bootstrap_", b))
    }
    core_b <- try(suppressWarnings(zentangler_with_fold_groups(bt$fold_group, do.call(
      zentangler_sequential_core,
      c(list(
        data_by_view = bt$blocks, pheno_df = bt$pheno_df, X = X[idx], Y = zentangler_subset_y(Y, idx),
        C = if (is.null(C)) NULL else C[idx, , drop = FALSE], b_inference = "none"
      ), args_b)
    ))), silent = TRUE)
    if (inherits(core_b, "try-error")) {
      failures <- c(failures, paste0("bootstrap_", b, ": ", conditionMessage(attr(core_b, "condition"))))
      next
    }
    s_row <- setNames(rep(0, length(path_keys)), path_keys)
    b_row <- setNames(rep(0, length(terminal_keys)), terminal_keys)
    sp <- core_b$sequential_paths
    if (nrow(sp) > 0L) {
      hit <- match(path_keys, sp$key_path)
      ok <- !is.na(hit)
      s_row[ok] <- sp$sequential_score[hit[ok]]
      te <- core_b$terminal_effects
      hit_t <- match(terminal_keys, te$terminal_key)
      ok_t <- !is.na(hit_t)
      b_row[ok_t] <- te$b_fusion[hit_t[ok_t]]
    }
    s_row[!is.finite(s_row)] <- 0
    score_mat[b, ] <- s_row
    b_mat[b, ] <- b_row
    n_ok <- n_ok + 1L
    if (b %% zentangler_progress_every(repeats) == 0L || b == repeats) {
      zentangler_log(verbose, sprintf("bootstrap %d / %d (%d succeeded)", b, repeats, n_ok), indent = 1L)
    }
  }
  zentangler_check_bootstrap_success(n_ok, repeats, "Sequential pipeline bootstrap")
  list(score_matrix = score_mat, b_matrix = b_mat, repeats_requested = repeats, repeats_successful = n_ok, failures = failures)
}

zentangler_bootstrap_sequential_effects <- function(effect_paths, data_by_key, pheno_df, X, Y, C, y_family, x0, x1, n_layers, repeats, seed, bootstrap_id = NULL, verbose = FALSE) {
  if (is.null(seed) || !is.finite(as.numeric(seed))) seed <- 10001L
  seed <- as.integer(seed) + 500000L
  dummy <- list(.d = matrix(0, nrow = nrow(pheno_df), ncol = 1L, dimnames = list(rownames(pheno_df), ".d")))
  rows <- list()
  failures <- character(0)
  for (b in seq_len(repeats)) {
    bt <- try(make_bootstrap_sample(dummy, pheno_df, bootstrap_id = bootstrap_id, seed = seed + b), silent = TRUE)
    if (inherits(bt, "try-error")) {
      failures <- c(failures, paste0("effect_bootstrap_", b, ": ", conditionMessage(attr(bt, "condition"))))
      next
    }
    idx <- bt$idx
    eff_b <- try(suppressWarnings(compute_effect_decomposition_sequential(
      sequential_paths = effect_paths, Y = zentangler_subset_y(Y, idx), X = X[idx],
      data_by_key = lapply(data_by_key, function(v) v[idx]),
      C = if (is.null(C)) NULL else C[idx, , drop = FALSE],
      y_family = y_family, x0 = x0, x1 = x1, n_layers = n_layers
    )), silent = TRUE)
    if (inherits(eff_b, "try-error") || !any(is.finite(eff_b$nie_total))) {
      failures <- c(failures, paste0("effect_bootstrap_", b, ": effect model failed"))
      next
    }
    rows[[length(rows) + 1L]] <- eff_b
    if (b %% zentangler_progress_every(repeats) == 0L || b == repeats) {
      zentangler_log(verbose, sprintf("effect bootstrap %d / %d (%d succeeded)", b, repeats, length(rows)), indent = 1L)
    }
  }
  zentangler_check_bootstrap_success(length(rows), repeats, "Sequential causal-effect bootstrap")
  list(effect_decomposition = zentangler_bind_rows_fill(rows), repeats_requested = repeats, repeats_successful = length(rows), failures = failures)
}

