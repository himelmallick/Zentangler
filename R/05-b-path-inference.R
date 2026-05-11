# -----------------------------------------------------------------------------
# B-stage p-values and effect summaries
# -----------------------------------------------------------------------------

infer_p_b_multiview_refit <- function(Y, X, blocks, C = NULL, coef_by_view, y_family = c("gaussian", "binomial")) {
  y_family <- match.arg(y_family)
  out <- lapply(blocks, function(M) setNames(rep(NA_real_, ncol(M)), colnames(M)))

  active_design <- character(0)
  active_map <- data.frame(view = character(), mediator = character(), design_name = character())

  for (view in names(blocks)) {
    b <- coef_by_view[[view]]
    active <- names(b)[is.finite(b) & b != 0]
    if (length(active) == 0) next

    dn <- paste0(view, "::", active)
    active_design <- c(active_design, dn)
    active_map <- rbind(
      active_map,
      data.frame(view = view, mediator = active, design_name = dn, stringsAsFactors = FALSE)
    )
  }
  if (length(active_design) == 0) return(out)

  active_views <- unique(active_map$view)
  active_blocks <- lapply(active_views, function(view) {
    blocks[[view]][, active_map$mediator[active_map$view == view], drop = FALSE]
  })
  names(active_blocks) <- active_views

  pref <- prefix_multiview_blocks(active_blocks)
  M_active <- pref$X[, intersect(colnames(pref$X), active_design), drop = FALSE]

  Z <- cbind(`(Intercept)` = 1, X = X)
  if (!is.null(C) && ncol(C) > 0) Z <- cbind(Z, as.matrix(C))
  Z <- cbind(Z, M_active)

  if (identical(y_family, "gaussian")) {
    fit <- lm.fit(x = as.matrix(Z), y = as.numeric(Y))
    rank <- fit$rank
    df <- nrow(Z) - rank
    if (df <= 0) return(out)

    rss <- sum(fit$residuals^2)
    sigma2 <- rss / df
    R <- try(qr.R(fit$qr)[seq_len(rank), seq_len(rank), drop = FALSE], silent = TRUE)
    if (inherits(R, "try-error")) return(out)
    xtx_inv <- try(chol2inv(R), silent = TRUE)
    if (inherits(xtx_inv, "try-error")) return(out)

    beta <- fit$coefficients
    se <- rep(NA_real_, length(beta))
    se[fit$qr$pivot[seq_len(rank)]] <- sqrt(diag(xtx_inv) * sigma2)
    names(se) <- names(beta)
    pvals <- 2 * stats::pt(abs(beta / se), df = df, lower.tail = FALSE)
  } else {
    dat <- as.data.frame(Z, check.names = TRUE)
    safe_names <- colnames(dat)
    fit <- try(suppressWarnings(glm(Y ~ . - 1, data = cbind(Y = Y, dat), family = binomial())), silent = TRUE)
    if (inherits(fit, "try-error")) return(out)
    sm <- summary(fit)$coefficients
    pvals <- sm[, ncol(sm)]
    names(pvals) <- colnames(Z)[match(rownames(sm), safe_names)]
  }

  for (i in seq_len(nrow(active_map))) {
    view <- active_map$view[i]
    med <- active_map$mediator[i]
    dn <- active_map$design_name[i]
    if (dn %in% names(pvals)) out[[view]][med] <- as.numeric(pvals[dn])
  }

  out
}

pick_lambda_from_cv <- function(cvfit, lambda_choice = c("lambda.1se", "lambda.min")) {
  lambda_choice <- match.arg(lambda_choice)
  lam <- suppressWarnings(as.numeric(cvfit[[lambda_choice]]))
  if (!is.finite(lam) || lam <= 0) {
    lam_grid <- suppressWarnings(as.numeric(cvfit$lambda))
    lam <- if (length(lam_grid) > 0) max(min(lam_grid), 1e-06) else 1e-04
  }
  lam
}

debiased_lasso_pvals_multiview_gaussian <- function(
  Y,
  X,
  blocks,
  C = NULL,
  coef_by_view,
  lambda_choice = c("lambda.1se", "lambda.min"),
  max_targets = 200L
) {
  # HIMA-style de-biased lasso for the B-stage.
  #
  # This uses the whole SIS-selected multiview design, then de-biases the
  # nonzero lasso targets. It is more defensible than active-set refit because
  # it corrects lasso shrinkage, but it is still not a complete selective
  # inference theorem for the full screening + fusion pipeline.
  lambda_choice <- match.arg(lambda_choice)

  out <- lapply(blocks, function(M) setNames(rep(NA_real_, ncol(M)), colnames(M)))

  active_map <- data.frame(view = character(), mediator = character(), design_name = character())
  for (view in names(blocks)) {
    b <- coef_by_view[[view]]
    active <- names(b)[is.finite(b) & b != 0]
    if (length(active) == 0) next
    active_map <- rbind(
      active_map,
      data.frame(
        view = view,
        mediator = active,
        design_name = paste0(view, "::", active),
        stringsAsFactors = FALSE
      )
    )
  }
  if (nrow(active_map) == 0) {
    return(list(p_b = out, method = "debiased_lasso_gaussian_no_active_targets"))
  }

  pref <- prefix_multiview_blocks(blocks)
  dat <- data.frame(y = as.numeric(Y), X = as.numeric(X), check.names = FALSE)
  x_cols <- "X"
  if (!is.null(C) && ncol(C) > 0) {
    C_df <- as.data.frame(C, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(C_df) <- paste0("C::", colnames(C_df))
    dat <- cbind(dat, C_df)
    x_cols <- c(x_cols, colnames(C_df))
  }
  dat <- cbind(dat, as.data.frame(pref$X, check.names = FALSE, stringsAsFactors = FALSE))

  cc <- complete.cases(dat)
  dat <- dat[cc, , drop = FALSE]
  if (nrow(dat) < 20) {
    return(list(p_b = out, method = "debiased_lasso_gaussian_too_few_samples"))
  }

  y <- as.numeric(dat$y)
  Z_raw <- as.matrix(dat[, setdiff(colnames(dat), "y"), drop = FALSE])
  n <- nrow(Z_raw)
  p <- ncol(Z_raw)
  if (p < 2) {
    return(list(p_b = out, method = "debiased_lasso_gaussian_design_too_small"))
  }

  z_center <- colMeans(Z_raw)
  z_scale <- apply(Z_raw, 2, stats::sd)
  z_scale[!is.finite(z_scale) | z_scale < 1e-08] <- 1
  Z <- sweep(sweep(Z_raw, 2, z_center, "-"), 2, z_scale, "/")
  y0 <- y - mean(y)

  pf <- rep(1, p)
  names(pf) <- colnames(Z)
  pf[colnames(Z) %in% x_cols] <- 0

  cv_main <- glmnet::cv.glmnet(
    x = Z,
    y = y0,
    family = "gaussian",
    alpha = 1,
    penalty.factor = pf,
    standardize = FALSE,
    intercept = FALSE
  )
  lam_main <- pick_lambda_from_cv(cv_main, lambda_choice = lambda_choice)
  fit_main <- glmnet::glmnet(
    x = Z,
    y = y0,
    family = "gaussian",
    alpha = 1,
    lambda = lam_main,
    penalty.factor = pf,
    standardize = FALSE,
    intercept = FALSE
  )
  beta_hat <- as.numeric(glmnet::coef.glmnet(fit_main, s = lam_main))[-1]
  names(beta_hat) <- colnames(Z)

  resid <- y0 - as.numeric(Z %*% beta_hat)
  df_eff <- sum(abs(beta_hat) > 0) + 1L
  sigma2_hat <- sum(resid^2) / max(1, n - df_eff)
  sigma2_hat <- max(sigma2_hat, 1e-10)

  target_names <- active_map$design_name[active_map$design_name %in% colnames(Z)]
  if (length(target_names) == 0) {
    return(list(p_b = out, method = "debiased_lasso_gaussian_no_valid_targets"))
  }

  max_targets <- max(1L, as.integer(max_targets))
  if (length(target_names) > max_targets) {
    ord <- order(abs(beta_hat[target_names]), decreasing = TRUE, na.last = TRUE)
    target_names <- target_names[ord[seq_len(max_targets)]]
  }

  p_map <- setNames(rep(NA_real_, length(target_names)), target_names)
  for (nm in target_names) {
    j <- match(nm, colnames(Z))
    if (is.na(j)) next

    z_j <- Z[, j]
    z_m <- Z[, -j, drop = FALSE]
    if (ncol(z_m) == 0) next

    cv_j <- glmnet::cv.glmnet(
      x = z_m,
      y = z_j,
      family = "gaussian",
      alpha = 1,
      standardize = FALSE,
      intercept = FALSE
    )
    lam_j <- pick_lambda_from_cv(cv_j, lambda_choice = lambda_choice)
    fit_j <- glmnet::glmnet(
      x = z_m,
      y = z_j,
      family = "gaussian",
      alpha = 1,
      lambda = lam_j,
      standardize = FALSE,
      intercept = FALSE
    )
    gamma_j <- as.numeric(glmnet::coef.glmnet(fit_j, s = lam_j))[-1]
    r_j <- z_j - as.numeric(z_m %*% gamma_j)
    tau_j <- mean(z_j * r_j)
    if (!is.finite(tau_j) || abs(tau_j) < 1e-08) next

    w_j <- r_j / tau_j
    beta_db <- beta_hat[[nm]] + mean(w_j * resid)
    se_db <- sqrt(sigma2_hat * mean(w_j^2) / n)
    if (!is.finite(se_db) || se_db <= 0) next

    p_map[[nm]] <- 2 * stats::pnorm(-abs(beta_db / se_db))
  }

  for (i in seq_len(nrow(active_map))) {
    view <- active_map$view[i]
    med <- active_map$mediator[i]
    dn <- active_map$design_name[i]
    if (dn %in% names(p_map)) out[[view]][med] <- as.numeric(p_map[[dn]])
  }

  list(p_b = out, method = "debiased_lasso_gaussian")
}

debiased_lasso_pvals_multiview_logistic <- function(
  Y,
  X,
  blocks,
  C = NULL,
  coef_by_view,
  lambda_choice = c("lambda.1se", "lambda.min"),
  max_targets = 200L
) {
  # Approximate de-biased logistic lasso for the B-stage.
  #
  # This is the binomial analogue of the Gaussian de-biasing routine above. It
  # uses a penalized logistic outcome model followed by weighted nodewise lasso
  # corrections for active mediator coefficients. It provides an approximate
  # Wald-style p-value for binary outcomes; it is not a full selective-inference
  # theorem for screening + fusion.
  lambda_choice <- match.arg(lambda_choice)

  out <- lapply(blocks, function(M) setNames(rep(NA_real_, ncol(M)), colnames(M)))

  y_raw <- as.numeric(Y)
  y_unique <- sort(unique(y_raw[is.finite(y_raw)]))
  if (length(y_unique) != 2L) {
    return(list(p_b = out, method = "debiased_lasso_logistic_invalid_binary_y"))
  }
  y01 <- as.numeric(y_raw == max(y_unique))

  active_map <- data.frame(view = character(), mediator = character(), design_name = character())
  for (view in names(blocks)) {
    b <- coef_by_view[[view]]
    active <- names(b)[is.finite(b) & b != 0]
    if (length(active) == 0) next
    active_map <- rbind(
      active_map,
      data.frame(
        view = view,
        mediator = active,
        design_name = paste0(view, "::", active),
        stringsAsFactors = FALSE
      )
    )
  }
  if (nrow(active_map) == 0) {
    return(list(p_b = out, method = "debiased_lasso_logistic_no_active_targets"))
  }

  pref <- prefix_multiview_blocks(blocks)
  dat <- data.frame(y = y01, X = as.numeric(X), check.names = FALSE)
  x_cols <- "X"
  if (!is.null(C) && ncol(C) > 0) {
    C_df <- as.data.frame(C, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(C_df) <- paste0("C::", colnames(C_df))
    dat <- cbind(dat, C_df)
    x_cols <- c(x_cols, colnames(C_df))
  }
  dat <- cbind(dat, as.data.frame(pref$X, check.names = FALSE, stringsAsFactors = FALSE))

  cc <- complete.cases(dat)
  dat <- dat[cc, , drop = FALSE]
  if (nrow(dat) < 30 || length(unique(dat$y)) < 2L) {
    return(list(p_b = out, method = "debiased_lasso_logistic_too_few_or_one_class"))
  }

  y <- as.numeric(dat$y)
  Z_raw <- as.matrix(dat[, setdiff(colnames(dat), "y"), drop = FALSE])
  n <- nrow(Z_raw)
  p <- ncol(Z_raw)
  if (p < 2) {
    return(list(p_b = out, method = "debiased_lasso_logistic_design_too_small"))
  }

  z_center <- colMeans(Z_raw)
  z_scale <- apply(Z_raw, 2, stats::sd)
  z_scale[!is.finite(z_scale) | z_scale < 1e-08] <- 1
  Z <- sweep(sweep(Z_raw, 2, z_center, "-"), 2, z_scale, "/")

  pf <- rep(1, p)
  names(pf) <- colnames(Z)
  pf[colnames(Z) %in% x_cols] <- 0

  cv_main <- glmnet::cv.glmnet(
    x = Z,
    y = y,
    family = "binomial",
    alpha = 1,
    penalty.factor = pf,
    standardize = FALSE,
    intercept = TRUE
  )
  lam_main <- pick_lambda_from_cv(cv_main, lambda_choice = lambda_choice)
  fit_main <- glmnet::glmnet(
    x = Z,
    y = y,
    family = "binomial",
    alpha = 1,
    lambda = lam_main,
    penalty.factor = pf,
    standardize = FALSE,
    intercept = TRUE
  )
  coef_main <- as.numeric(glmnet::coef.glmnet(fit_main, s = lam_main))
  beta0_hat <- coef_main[1]
  beta_hat <- coef_main[-1]
  names(beta_hat) <- colnames(Z)

  eta <- as.numeric(beta0_hat + Z %*% beta_hat)
  mu <- stats::plogis(eta)
  w <- pmax(mu * (1 - mu), 1e-05)
  score_resid <- y - mu

  target_names <- active_map$design_name[active_map$design_name %in% colnames(Z)]
  if (length(target_names) == 0) {
    return(list(p_b = out, method = "debiased_lasso_logistic_no_valid_targets"))
  }

  max_targets <- max(1L, as.integer(max_targets))
  if (length(target_names) > max_targets) {
    ord <- order(abs(beta_hat[target_names]), decreasing = TRUE, na.last = TRUE)
    target_names <- target_names[ord[seq_len(max_targets)]]
  }

  p_map <- setNames(rep(NA_real_, length(target_names)), target_names)
  sqrt_w <- sqrt(w)
  for (nm in target_names) {
    j <- match(nm, colnames(Z))
    if (is.na(j)) next

    z_j <- Z[, j]
    z_m <- Z[, -j, drop = FALSE]
    if (ncol(z_m) == 0) next

    z_j_w <- z_j * sqrt_w
    z_m_w <- sweep(z_m, 1, sqrt_w, "*")

    cv_j <- glmnet::cv.glmnet(
      x = z_m_w,
      y = z_j_w,
      family = "gaussian",
      alpha = 1,
      standardize = FALSE,
      intercept = FALSE
    )
    lam_j <- pick_lambda_from_cv(cv_j, lambda_choice = lambda_choice)
    fit_j <- glmnet::glmnet(
      x = z_m_w,
      y = z_j_w,
      family = "gaussian",
      alpha = 1,
      lambda = lam_j,
      standardize = FALSE,
      intercept = FALSE
    )
    gamma_j <- as.numeric(glmnet::coef.glmnet(fit_j, s = lam_j))[-1]
    r_j <- z_j - as.numeric(z_m %*% gamma_j)
    tau_j <- mean(w * z_j * r_j)
    if (!is.finite(tau_j) || abs(tau_j) < 1e-08) next

    beta_db <- beta_hat[[nm]] + mean(r_j * score_resid) / tau_j
    se_db <- sqrt(mean(w * r_j^2) / (n * tau_j^2))
    if (!is.finite(se_db) || se_db <= 0) next

    p_map[[nm]] <- 2 * stats::pnorm(-abs(beta_db / se_db))
  }

  for (i in seq_len(nrow(active_map))) {
    view <- active_map$view[i]
    med <- active_map$mediator[i]
    dn <- active_map$design_name[i]
    if (dn %in% names(p_map)) out[[view]][med] <- as.numeric(p_map[[dn]])
  }

  list(p_b = out, method = "debiased_lasso_logistic")
}

infer_p_b_multiview <- function(
  Y,
  X,
  blocks,
  C = NULL,
  coef_by_view,
  y_family = c("gaussian", "binomial"),
  method = c("debiased_lasso", "debiased_logistic_lasso", "refit"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  max_debias_targets = 200L
) {
  y_family <- match.arg(y_family)
  method <- match.arg(method)
  lambda_choice <- match.arg(lambda_choice)

  if (identical(method, "debiased_lasso") && identical(y_family, "gaussian")) {
    db_try <- try(
      debiased_lasso_pvals_multiview_gaussian(
        Y = Y,
        X = X,
        blocks = blocks,
        C = C,
        coef_by_view = coef_by_view,
        lambda_choice = lambda_choice,
        max_targets = max_debias_targets
      ),
      silent = TRUE
    )
    if (!inherits(db_try, "try-error")) return(db_try)

    msg <- conditionMessage(attr(db_try, "condition"))
    warning("Debiased-lasso B-path inference failed; falling back to refit p-values. Reason: ", msg)
  } else if (
    (identical(method, "debiased_lasso") && identical(y_family, "binomial")) ||
      identical(method, "debiased_logistic_lasso")
  ) {
    if (!identical(y_family, "binomial")) {
      warning("debiased_logistic_lasso is intended for binomial Y; falling back to refit p-values.")
    } else {
      db_try <- try(
        debiased_lasso_pvals_multiview_logistic(
          Y = Y,
          X = X,
          blocks = blocks,
          C = C,
          coef_by_view = coef_by_view,
          lambda_choice = lambda_choice,
          max_targets = max_debias_targets
        ),
        silent = TRUE
      )
      if (!inherits(db_try, "try-error")) return(db_try)

      msg <- conditionMessage(attr(db_try, "condition"))
      warning("Debiased-logistic-lasso B-path inference failed; falling back to refit p-values. Reason: ", msg)
    }
  }

  list(
    p_b = infer_p_b_multiview_refit(
      Y = Y,
      X = X,
      blocks = blocks,
      C = C,
      coef_by_view = coef_by_view,
      y_family = y_family
    ),
    method = "active_set_refit"
  )
}

assemble_mediator_table <- function(screen_tab, b_vec, p_b_vec, fdr_method = c("BH", "BY")) {
  fdr_method <- validate_fdr_method(fdr_method)
  a_map <- setNames(screen_tab$a_hat, screen_tab$mediator)
  p_a_map <- setNames(screen_tab$p_a, screen_tab$mediator)
  q_a_map <- setNames(screen_tab$q_a, screen_tab$mediator)
  sel_map <- setNames(screen_tab$selected, screen_tab$mediator)
  med <- names(a_map)
  selected <- as.logical(sel_map[med])

  out <- data.frame(
    mediator = med,
    a = a_map[med],
    p_a = p_a_map[med],
    q_a = q_a_map[med],
    b = b_vec[med],
    p_b = p_b_vec[med],
    selected_by_screen = selected,
    selected_by_sis = selected,
    stringsAsFactors = FALSE
  )

  out$score <- out$a * out$b
  out$abs_score <- abs(out$score)
  both_legs <- is.finite(out$p_a) & is.finite(out$p_b)
  out$joint_p_ab <- NA_real_
  out$joint_p_ab[both_legs] <- pmax(out$p_a[both_legs], out$p_b[both_legs])
  out$p_primary <- out$joint_p_ab
  out$q_primary_within_view <- p.adjust(out$p_primary, method = fdr_method)
  out$q_primary <- out$q_primary_within_view

  out <- out[order(out$abs_score, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}

