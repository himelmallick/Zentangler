infer_effect_contrast <- function(X, x0 = NULL, x1 = NULL) {
  x_fin <- as.numeric(X[is.finite(X)])
  if (length(x_fin) == 0) stop("Cannot infer effect contrast from non-finite X.")

  if (is.null(x0) || !is.finite(as.numeric(x0))) {
    ux <- sort(unique(x_fin))
    x0 <- if (length(ux) == 2L) ux[1] else as.numeric(stats::quantile(x_fin, probs = 0.25, names = FALSE))
  }
  if (is.null(x1) || !is.finite(as.numeric(x1))) {
    ux <- sort(unique(x_fin))
    x1 <- if (length(ux) == 2L) ux[2] else as.numeric(stats::quantile(x_fin, probs = 0.75, names = FALSE))
  }
  x0 <- as.numeric(x0)
  x1 <- as.numeric(x1)
  if (!is.finite(x0) || !is.finite(x1)) stop("effect_x0 and effect_x1 must be finite.")
  if (isTRUE(all.equal(x0, x1))) stop("effect_x0 and effect_x1 must be different.")
  c(x0 = x0, x1 = x1, delta = x1 - x0)
}

estimate_reduced_x_effect <- function(Y, X, C = NULL, y_family = c("gaussian", "binomial", "survival"), x0 = NULL, x1 = NULL) {
  y_family <- match.arg(y_family)
  contrast <- infer_effect_contrast(X, x0 = x0, x1 = x1)

  W <- cbind(`(Intercept)` = 1, X = as.numeric(X))
  if (!is.null(C) && ncol(C) > 0) W <- cbind(W, as.matrix(C))
  W <- as.matrix(W)

  fit <- try(
    if (identical(y_family, "survival")) {
      dat <- as.data.frame(W[, -1, drop = FALSE], check.names = TRUE)
      survival::coxph(Y ~ ., data = dat)
    } else if (identical(y_family, "gaussian")) {
      lm.fit(x = W, y = as.numeric(Y))
    } else {
      suppressWarnings(glm.fit(x = W, y = as.numeric(Y), family = binomial(), intercept = FALSE))
    },
    silent = TRUE
  )
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NA_real_)

  out <- as.numeric(fit$coefficients["X"])
  if (!is.finite(out)) return(NA_real_)
  if (identical(y_family, "survival")) return(out * contrast[["delta"]])
  if (identical(y_family, "gaussian")) return(out * contrast[["delta"]])

  co <- fit$coefficients
  co[!is.finite(co)] <- 0
  W0 <- W
  W1 <- W
  W0[, "X"] <- contrast[["x0"]]
  W1[, "X"] <- contrast[["x1"]]
  mean(stats::plogis(as.numeric(W1 %*% co)) - stats::plogis(as.numeric(W0 %*% co)), na.rm = TRUE)
}

make_active_mediator_matrix <- function(combined_mediators, blocks) {
  if (is.null(combined_mediators) || nrow(combined_mediators) == 0) {
    return(list(M = matrix(nrow = nrow(blocks[[1]]), ncol = 0), map = data.frame()))
  }
  active <- is.finite(combined_mediators$a) &
    is.finite(combined_mediators$b) &
    combined_mediators$b != 0
  active_tab <- combined_mediators[active, , drop = FALSE]
  if (nrow(active_tab) == 0) {
    return(list(M = matrix(nrow = nrow(blocks[[1]]), ncol = 0), map = data.frame()))
  }

  mats <- list()
  map <- data.frame(
    omics = character(),
    mediator = character(),
    design_name = character(),
    a = numeric(),
    b_original = numeric(),
    stringsAsFactors = FALSE
  )
  for (ii in seq_len(nrow(active_tab))) {
    view <- as.character(active_tab$omics[ii])
    med <- as.character(active_tab$mediator[ii])
    if (!(view %in% names(blocks)) || !(med %in% colnames(blocks[[view]]))) next
    design_name <- paste0(view, "::", med)
    vals <- as.matrix(blocks[[view]][, med, drop = FALSE])
    colnames(vals) <- design_name
    mats[[length(mats) + 1L]] <- vals
    map <- rbind(
      map,
      data.frame(
        omics = view,
        mediator = med,
        design_name = design_name,
        a = as.numeric(active_tab$a[ii]),
        b_original = as.numeric(active_tab$b[ii]),
        stringsAsFactors = FALSE
      )
    )
  }
  if (length(mats) == 0) {
    return(list(M = matrix(nrow = nrow(blocks[[1]]), ncol = 0), map = data.frame()))
  }
  list(M = do.call(cbind, mats), map = map)
}

fit_active_outcome_effect_model <- function(Y, X, M_active, C = NULL, y_family = c("gaussian", "binomial", "survival")) {
  y_family <- match.arg(y_family)
  Z <- cbind(`(Intercept)` = 1, X = as.numeric(X))
  if (!is.null(C) && ncol(C) > 0) {
    C2 <- as.matrix(C)
    colnames(C2) <- paste0("C::", colnames(C2))
    Z <- cbind(Z, C2)
  }
  if (!is.null(M_active) && ncol(M_active) > 0) Z <- cbind(Z, as.matrix(M_active))
  Z <- as.matrix(Z)

  fit <- try(
    if (identical(y_family, "survival")) {
      dat <- as.data.frame(Z[, -1, drop = FALSE], check.names = TRUE)
      survival::coxph(Y ~ ., data = dat)
    } else if (identical(y_family, "gaussian")) {
      stats::lm.fit(x = Z, y = as.numeric(Y))
    } else {
      suppressWarnings(stats::glm.fit(x = Z, y = as.numeric(Y), family = stats::binomial(), intercept = FALSE))
    },
    silent = TRUE
  )
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) {
    co <- setNames(rep(0, ncol(Z)), colnames(Z))
  } else {
    co <- fit$coefficients
    names(co) <- colnames(Z)
    co[!is.finite(co)] <- 0
  }
  list(coefficients = co, design = Z)
}

predict_effect_mean <- function(co, Z_template, X_value, M_value = NULL, y_family = c("gaussian", "binomial", "survival")) {
  y_family <- match.arg(y_family)
  Z <- Z_template
  Z[, "X"] <- as.numeric(X_value)
  if (!is.null(M_value) && ncol(M_value) > 0) {
    hit <- intersect(colnames(M_value), colnames(Z))
    Z[, hit] <- M_value[, hit, drop = FALSE]
  }
  eta <- as.numeric(Z %*% co[colnames(Z)])
  if (identical(y_family, "binomial")) stats::plogis(eta) else eta
}

compute_effect_decomposition_multiview <- function(
  combined_mediators,
  direct_effect,
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  x0 = NULL,
  x1 = NULL
) {
  y_family <- match.arg(y_family)
  contrast <- infer_effect_contrast(X, x0 = x0, x1 = x1)
  x0 <- contrast[["x0"]]
  x1 <- contrast[["x1"]]
  delta_x <- contrast[["delta"]]

  active_info <- make_active_mediator_matrix(combined_mediators, blocks)
  M_obs <- active_info$M
  active_map <- active_info$map

  if (identical(y_family, "survival")) {
    direct_linear_product <- as.numeric(direct_effect) * delta_x
    indirect_linear_product <- if (nrow(active_map) > 0) sum(active_map$a * active_map$b_original * delta_x, na.rm = TRUE) else 0
    reduced_te <- estimate_reduced_x_effect(Y = Y, X = X, C = C, y_family = y_family, x0 = x0, x1 = x1)
    return(data.frame(
      effect_method = "cox_log_hazard_product_summary",
      effect_scale = "log_hazard_ratio",
      x0 = x0,
      x1 = x1,
      delta_x = delta_x,
      nde = direct_linear_product,
      nie_total = indirect_linear_product,
      te = direct_linear_product + indirect_linear_product,
      prop_mediated = if (is.finite(reduced_te) && abs(reduced_te) > 1e-08) indirect_linear_product / reduced_te else NA_real_,
      reduced_total_effect = reduced_te,
      decomposition_gap = NA_real_,
      direct_effect = direct_linear_product,
      indirect_total = indirect_linear_product,
      total_decomp = direct_linear_product + indirect_linear_product,
      direct_coefficient_active_model = as.numeric(direct_effect),
      direct_effect_y_model_coefficient = as.numeric(direct_effect),
      direct_effect_y_model_product = direct_linear_product,
      indirect_product_sum = indirect_linear_product,
      n_active_total = nrow(active_map),
      stringsAsFactors = FALSE
    ))
  }

  effect_fit <- fit_active_outcome_effect_model(
    Y = Y,
    X = X,
    M_active = M_obs,
    C = C,
    y_family = y_family
  )
  co <- effect_fit$coefficients
  Z_template <- effect_fit$design

  M_x0 <- M_obs
  M_x1 <- M_obs
  if (ncol(M_obs) > 0 && nrow(active_map) > 0) {
    for (jj in seq_len(nrow(active_map))) {
      dn <- active_map$design_name[jj]
      a_j <- active_map$a[jj]
      if (!is.finite(a_j) || !(dn %in% colnames(M_obs))) next
      M_x0[, dn] <- M_obs[, dn] + a_j * (x0 - as.numeric(X))
      M_x1[, dn] <- M_obs[, dn] + a_j * (x1 - as.numeric(X))
    }
  }

  mu_00 <- predict_effect_mean(co, Z_template, X_value = x0, M_value = M_x0, y_family = y_family)
  mu_10 <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_x0, y_family = y_family)
  mu_11 <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_x1, y_family = y_family)

  nde <- mean(mu_10 - mu_00, na.rm = TRUE)
  nie <- mean(mu_11 - mu_10, na.rm = TRUE)
  te <- mean(mu_11 - mu_00, na.rm = TRUE)
  reduced_te <- estimate_reduced_x_effect(Y = Y, X = X, C = C, y_family = y_family, x0 = x0, x1 = x1)

  direct_coef <- if ("X" %in% names(co)) as.numeric(co[["X"]]) else NA_real_
  direct_linear_product <- as.numeric(direct_effect) * delta_x
  indirect_linear_product <- if (nrow(active_map) > 0) sum(active_map$a * active_map$b_original * delta_x, na.rm = TRUE) else 0

  out <- data.frame(
    effect_method = if (identical(y_family, "binomial")) "parametric_g_computation_logistic" else "parametric_g_computation_linear",
    effect_scale = if (identical(y_family, "binomial")) "risk_difference" else "mean_difference",
    x0 = x0,
    x1 = x1,
    delta_x = delta_x,
    nde = nde,
    nie_total = nie,
    te = te,
    prop_mediated = if (is.finite(te) && abs(te) > 1e-08) nie / te else NA_real_,
    reduced_total_effect = reduced_te,
    decomposition_gap = te - (nde + nie),
    direct_effect = nde,
    indirect_total = nie,
    total_decomp = te,
    direct_coefficient_active_model = direct_coef,
    direct_effect_y_model_coefficient = as.numeric(direct_effect),
    direct_effect_y_model_product = direct_linear_product,
    indirect_product_sum = indirect_linear_product,
    n_active_total = nrow(active_map),
    stringsAsFactors = FALSE
  )

  views <- if (nrow(active_map) > 0) unique(as.character(active_map$omics)) else character(0)
  mu_base <- mu_10
  for (view in views) {
    safe <- make.names(view)
    M_view <- M_x0
    idx <- active_map$omics == view
    view_cols <- active_map$design_name[idx]
    M_view[, view_cols] <- M_x1[, view_cols, drop = FALSE]
    mu_view <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_view, y_family = y_family)
    out[[paste0("nie_", safe)]] <- mean(mu_view - mu_base, na.rm = TRUE)
    out[[paste0("n_active_", safe)]] <- sum(idx)
    out[[paste0("indirect_", safe)]] <- out[[paste0("nie_", safe)]]
  }

  out
}

build_causal_stage_summary <- function(effect_decomposition, combined_mediators, causal_inference = "none", bootstrap = NULL) {
  if (is.null(effect_decomposition) || nrow(effect_decomposition) == 0) return(list())
  active_count <- if (!is.null(combined_mediators) && nrow(combined_mediators) > 0) {
    sum(is.finite(combined_mediators$b) & combined_mediators$b != 0, na.rm = TRUE)
  } else {
    0L
  }
  eff <- effect_decomposition[1, , drop = FALSE]
  list(
    stage = "post_selection_causal_effects",
    stage_method = as.character(eff$effect_method[[1]]),
    effect_scale = as.character(eff$effect_scale[[1]]),
    causal_inference = as.character(causal_inference)[1],
    bootstrap_used = isTRUE(identical(causal_inference, "bootstrap")) && !is.null(bootstrap),
    x0 = as.numeric(eff$x0[[1]]),
    x1 = as.numeric(eff$x1[[1]]),
    delta_x = as.numeric(eff$delta_x[[1]]),
    n_active_mediators = as.integer(active_count),
    nde = as.numeric(eff$nde[[1]]),
    nie_total = as.numeric(eff$nie_total[[1]]),
    te = as.numeric(eff$te[[1]]),
    prop_mediated = as.numeric(eff$prop_mediated[[1]])
  )
}

zentangler_sequential_active_keys <- function(paths) {
  if (is.null(paths) || nrow(paths) == 0L || !("key_path" %in% colnames(paths))) return(character(0))
  keep <- is.finite(paths$b) & paths$b != 0
  if (any("q_primary" %in% colnames(paths))) {
    keep <- keep & is.finite(paths$q_primary)
  }
  paths_use <- paths[keep, , drop = FALSE]
  if (nrow(paths_use) == 0L) return(character(0))
  unique(unlist(strsplit(as.character(paths_use$key_path), " -> ", fixed = TRUE), use.names = FALSE))
}

zentangler_sequential_edge_table <- function(paths) {
  if (is.null(paths) || nrow(paths) == 0L || !("key_path" %in% colnames(paths))) return(data.frame())
  keep <- is.finite(paths$b) & paths$b != 0
  paths_use <- paths[keep, , drop = FALSE]
  if (nrow(paths_use) == 0L) return(data.frame())
  rows <- list()
  ii <- 0L
  for (path in as.character(paths_use$key_path)) {
    keys <- strsplit(path, " -> ", fixed = TRUE)[[1L]]
    if (length(keys) < 2L) next
    for (j in seq_len(length(keys) - 1L)) {
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        from_key = keys[j],
        to_key = keys[j + 1L],
        stage = j + 1L,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) return(data.frame())
  out <- unique(do.call(rbind, rows))
  rownames(out) <- NULL
  out
}

zentangler_fit_sequential_mediator_models <- function(paths, data_by_key, X, C = NULL) {
  active_keys <- zentangler_sequential_active_keys(paths)
  if (length(active_keys) == 0L) return(list(models = list(), active_keys = character(0), edge_table = data.frame()))
  edge_table <- zentangler_sequential_edge_table(paths)
  keys_by_stage <- split(c(edge_table$from_key[edge_table$stage == 2L], edge_table$to_key), c(rep(1L, sum(edge_table$stage == 2L)), edge_table$stage))
  if (length(edge_table) == 0L) {
    keys_by_stage <- list(active_keys)
    names(keys_by_stage) <- "1"
  } else {
    first_keys <- unique(edge_table$from_key[edge_table$stage == 2L])
    keys_by_stage <- split(edge_table$to_key, edge_table$stage)
    keys_by_stage <- c(list(`1` = first_keys), keys_by_stage)
    keys_by_stage <- lapply(keys_by_stage, unique)
  }
  keys_by_stage <- keys_by_stage[order(as.integer(names(keys_by_stage)))]
  models <- list()
  for (stage_name in names(keys_by_stage)) {
    stage <- as.integer(stage_name)
    for (key in keys_by_stage[[stage_name]]) {
      y <- as.numeric(data_by_key[[key]])
      if (stage <= 1L) {
        dat <- data.frame(y = y, X = as.numeric(X), check.names = FALSE)
        if (!is.null(C) && ncol(C) > 0L) dat <- cbind(dat, as.data.frame(C, check.names = FALSE))
      } else {
        parent_keys <- unique(edge_table$from_key[edge_table$to_key == key])
        parent_df <- as.data.frame(lapply(parent_keys, function(pk) as.numeric(data_by_key[[pk]])), check.names = FALSE)
        colnames(parent_df) <- parent_keys
        dat <- data.frame(y = y, X = as.numeric(X), check.names = FALSE)
        if (ncol(parent_df) > 0L) dat <- cbind(dat, parent_df)
        if (!is.null(C) && ncol(C) > 0L) dat <- cbind(dat, as.data.frame(C, check.names = FALSE))
      }
      fit <- try(stats::lm(y ~ ., data = dat), silent = TRUE)
      if (inherits(fit, "try-error")) next
      co <- stats::coef(fit)
      co[!is.finite(co)] <- 0
      models[[key]] <- list(stage = stage, parents = if (stage <= 1L) character(0) else parent_keys, coefficients = co)
    }
  }
  list(models = models, active_keys = active_keys, edge_table = edge_table)
}

zentangler_predict_sequential_mediators <- function(model_info, X_value, X_obs, data_by_key, C = NULL) {
  models <- model_info$models %||% list()
  if (length(models) == 0L) return(matrix(nrow = length(X_obs), ncol = 0))
  order_keys <- names(sort(vapply(models, function(x) x$stage %||% 1L, integer(1))))
  out <- matrix(NA_real_, nrow = length(X_obs), ncol = length(order_keys))
  colnames(out) <- order_keys
  for (key in order_keys) {
    info <- models[[key]]
    co <- info$coefficients
    val <- rep(co[["(Intercept)"]] %||% 0, length(X_obs))
    if ("X" %in% names(co)) val <- val + co[["X"]] * as.numeric(X_value)
    if (length(info$parents) > 0L) {
      for (pk in info$parents) {
        if (pk %in% colnames(out)) {
          val <- val + (co[[pk]] %||% 0) * out[, pk]
        } else {
          val <- val + (co[[pk]] %||% 0) * as.numeric(data_by_key[[pk]])
        }
      }
    }
    if (!is.null(C) && ncol(C) > 0L) {
      for (cn in colnames(C)) {
        if (cn %in% names(co)) val <- val + co[[cn]] * as.numeric(C[, cn])
      }
    }
    out[, key] <- val
  }
  out
}

compute_effect_decomposition_sequential <- function(
  sequential_paths,
  Y,
  X,
  data_by_key,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  x0 = NULL,
  x1 = NULL
) {
  y_family <- match.arg(y_family)
  contrast <- infer_effect_contrast(X, x0 = x0, x1 = x1)
  x0 <- contrast[["x0"]]
  x1 <- contrast[["x1"]]
  delta_x <- contrast[["delta"]]
  active_keys <- zentangler_sequential_active_keys(sequential_paths)
  if (length(active_keys) == 0L) {
    reduced_te <- estimate_reduced_x_effect(Y = Y, X = X, C = C, y_family = y_family, x0 = x0, x1 = x1)
    return(data.frame(
      effect_method = if (identical(y_family, "binomial")) "parametric_g_computation_logistic_sequential" else if (identical(y_family, "survival")) "cox_log_hazard_product_summary_sequential" else "parametric_g_computation_linear_sequential",
      effect_scale = if (identical(y_family, "binomial")) "risk_difference" else if (identical(y_family, "survival")) "log_hazard_ratio" else "mean_difference",
      x0 = x0, x1 = x1, delta_x = delta_x,
      nde = reduced_te, nie_total = 0, te = reduced_te,
      prop_mediated = 0, reduced_total_effect = reduced_te,
      decomposition_gap = 0, direct_effect = reduced_te, indirect_total = 0,
      total_decomp = reduced_te, n_active_paths = 0L, n_active_mediators = 0L,
      n_active_edges = 0L, indirect_product_sum = 0,
      stringsAsFactors = FALSE
    ))
  }
  model_info <- zentangler_fit_sequential_mediator_models(sequential_paths, data_by_key = data_by_key, X = X, C = C)
  M_obs <- as.matrix(as.data.frame(lapply(active_keys, function(k) as.numeric(data_by_key[[k]])), check.names = FALSE))
  colnames(M_obs) <- active_keys
  effect_fit <- fit_active_outcome_effect_model(Y = Y, X = X, M_active = M_obs, C = C, y_family = y_family)
  co <- effect_fit$coefficients
  Z_template <- effect_fit$design
  if (identical(y_family, "survival")) {
    reduced_te <- estimate_reduced_x_effect(Y = Y, X = X, C = C, y_family = y_family, x0 = x0, x1 = x1)
    path_scores <- sequential_paths$sequential_score[is.finite(sequential_paths$sequential_score) & is.finite(sequential_paths$b) & sequential_paths$b != 0]
    indirect_linear_product <- sum(path_scores, na.rm = TRUE) * delta_x
    direct_linear_product <- if ("X" %in% names(co)) as.numeric(co[["X"]]) * delta_x else reduced_te - indirect_linear_product
    return(data.frame(
      effect_method = "cox_log_hazard_product_summary_sequential",
      effect_scale = "log_hazard_ratio",
      x0 = x0, x1 = x1, delta_x = delta_x,
      nde = direct_linear_product, nie_total = indirect_linear_product, te = direct_linear_product + indirect_linear_product,
      prop_mediated = if (is.finite(reduced_te) && abs(reduced_te) > 1e-08) indirect_linear_product / reduced_te else NA_real_,
      reduced_total_effect = reduced_te, decomposition_gap = NA_real_,
      direct_effect = direct_linear_product, indirect_total = indirect_linear_product, total_decomp = direct_linear_product + indirect_linear_product,
      n_active_paths = nrow(sequential_paths), n_active_mediators = length(active_keys), n_active_edges = nrow(model_info$edge_table),
      indirect_product_sum = indirect_linear_product,
      stringsAsFactors = FALSE
    ))
  }
  M_x0 <- zentangler_predict_sequential_mediators(model_info, X_value = x0, X_obs = X, data_by_key = data_by_key, C = C)
  M_x1 <- zentangler_predict_sequential_mediators(model_info, X_value = x1, X_obs = X, data_by_key = data_by_key, C = C)
  mu_00 <- predict_effect_mean(co, Z_template, X_value = x0, M_value = M_x0, y_family = y_family)
  mu_10 <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_x0, y_family = y_family)
  mu_11 <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_x1, y_family = y_family)
  nde <- mean(mu_10 - mu_00, na.rm = TRUE)
  nie <- mean(mu_11 - mu_10, na.rm = TRUE)
  te <- mean(mu_11 - mu_00, na.rm = TRUE)
  reduced_te <- estimate_reduced_x_effect(Y = Y, X = X, C = C, y_family = y_family, x0 = x0, x1 = x1)
  stage_counts <- table(vapply(model_info$models, function(x) x$stage %||% 1L, integer(1)))
  out <- data.frame(
    effect_method = if (identical(y_family, "binomial")) "parametric_g_computation_logistic_sequential" else "parametric_g_computation_linear_sequential",
    effect_scale = if (identical(y_family, "binomial")) "risk_difference" else "mean_difference",
    x0 = x0, x1 = x1, delta_x = delta_x,
    nde = nde, nie_total = nie, te = te,
    prop_mediated = if (is.finite(te) && abs(te) > 1e-08) nie / te else NA_real_,
    reduced_total_effect = reduced_te, decomposition_gap = te - (nde + nie),
    direct_effect = nde, indirect_total = nie, total_decomp = te,
    n_active_paths = nrow(sequential_paths), n_active_mediators = length(active_keys), n_active_edges = nrow(model_info$edge_table),
    indirect_product_sum = sum(sequential_paths$sequential_score[is.finite(sequential_paths$sequential_score) & is.finite(sequential_paths$b) & sequential_paths$b != 0], na.rm = TRUE) * delta_x,
    stringsAsFactors = FALSE
  )
  for (st in names(stage_counts)) {
    keys_st <- names(which(vapply(model_info$models, function(x) identical(as.integer(x$stage %||% 1L), as.integer(st)), logical(1))))
    M_stage <- M_x0
    hit <- intersect(keys_st, colnames(M_stage))
    if (length(hit) > 0L) M_stage[, hit] <- M_x1[, hit, drop = FALSE]
    mu_stage <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_stage, y_family = y_family)
    safe <- paste0("stage", st)
    out[[paste0("nie_", safe)]] <- mean(mu_stage - mu_10, na.rm = TRUE)
    out[[paste0("n_active_", safe)]] <- as.integer(stage_counts[[st]])
  }
  out
}
