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

estimate_reduced_x_effect <- function(Y, X, C = NULL, y_family = c("gaussian", "binomial"), x0 = NULL, x1 = NULL) {
  y_family <- match.arg(y_family)
  contrast <- infer_effect_contrast(X, x0 = x0, x1 = x1)

  W <- cbind(`(Intercept)` = 1, X = as.numeric(X))
  if (!is.null(C) && ncol(C) > 0) W <- cbind(W, as.matrix(C))
  W <- as.matrix(W)

  fit <- try(
    if (identical(y_family, "gaussian")) {
      lm.fit(x = W, y = as.numeric(Y))
    } else {
      suppressWarnings(glm.fit(x = W, y = as.numeric(Y), family = binomial(), intercept = FALSE))
    },
    silent = TRUE
  )
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) return(NA_real_)

  out <- as.numeric(fit$coefficients["X"])
  if (!is.finite(out)) return(NA_real_)
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

fit_active_outcome_effect_model <- function(Y, X, M_active, C = NULL, y_family = c("gaussian", "binomial")) {
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
    if (identical(y_family, "gaussian")) {
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

predict_effect_mean <- function(co, Z_template, X_value, M_value = NULL, y_family = c("gaussian", "binomial")) {
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
  y_family = c("gaussian", "binomial"),
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

