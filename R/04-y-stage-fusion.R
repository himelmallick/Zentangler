# -----------------------------------------------------------------------------
# Fusion / Y-stage helpers
# -----------------------------------------------------------------------------

validate_glmnet_alpha <- function(glmnet_alpha) {
  if (length(glmnet_alpha) != 1L || !is.numeric(glmnet_alpha) || !is.finite(glmnet_alpha)) {
    stop("glmnet_alpha must be one finite numeric value between 0 and 1.")
  }
  if (glmnet_alpha < 0 || glmnet_alpha > 1) {
    stop("glmnet_alpha must be between 0 and 1. Use 1 for lasso, 0.5 for elastic net, and 0 for ridge.")
  }
  as.numeric(glmnet_alpha)
}

glmnet_alpha_label <- function(glmnet_alpha) {
  if (isTRUE(all.equal(glmnet_alpha, 1))) return("lasso")
  if (isTRUE(all.equal(glmnet_alpha, 0))) return("ridge")
  paste0("elastic net alpha=", format(glmnet_alpha, trim = TRUE, scientific = FALSE))
}

choose_block_lambda <- function(Y, W, M, y_family, lambda_choice) {
  Z <- cbind(W, M)
  pf <- c(rep(0, ncol(W)), rep(1, ncol(M)))

  cvfit <- glmnet::cv.glmnet(
    x = as.matrix(Z),
    y = Y,
    family = y_family,
    alpha = 1,
    penalty.factor = pf,
    standardize = TRUE,
    intercept = FALSE
  )

  lam <- as.numeric(cvfit[[lambda_choice]])
  if (!is.finite(lam)) lam <- as.numeric(cvfit$lambda[min(10, length(cvfit$lambda))])
  lam <- max(lam, 1e-8)

  list(lambda = lam, cvfit = cvfit)
}

update_gamma_gaussian <- function(Y, W, offset) {
  fit <- lm.fit(x = W, y = Y - offset)
  gamma <- as.numeric(fit$coefficients)
  gamma[!is.finite(gamma)] <- 0
  names(gamma) <- colnames(W)
  gamma
}

fit_block_lasso_for_late <- function(
  Y,
  X,
  M,
  C = NULL,
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  prefix = "B::"
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)

  M2 <- as.matrix(M)
  colnames(M2) <- paste0(prefix, colnames(M))

  Z <- cbind(X = X)
  if (!is.null(C) && ncol(C) > 0) Z <- cbind(Z, C)
  Z <- cbind(Z, M2)

  pf <- rep(1, ncol(Z))
  pf[colnames(Z) == "X"] <- 0
  if (!is.null(C) && ncol(C) > 0) pf[colnames(Z) %in% colnames(C)] <- 0

  cvfit <- glmnet::cv.glmnet(
    x = as.matrix(Z),
    y = Y,
    family = y_family,
    alpha = glmnet_alpha,
    penalty.factor = pf,
    standardize = TRUE
  )

  beta <- as.matrix(coef(cvfit, s = lambda_choice))[, 1]
  names(beta) <- gsub("`", "", names(beta), fixed = TRUE)

  pred <- as.numeric(predict(cvfit, newx = as.matrix(Z), s = lambda_choice, type = "response"))

  x_coef <- beta["X"]
  if (!is.finite(x_coef)) x_coef <- 0

  med_coef <- setNames(rep(0, ncol(M)), colnames(M))
  med_tag <- paste0(prefix, colnames(M))
  hit <- intersect(med_tag, names(beta))
  if (length(hit) > 0) med_coef[sub(paste0("^", prefix), "", hit)] <- beta[hit]

  list(cvfit = cvfit, pred = pred, x_coef = as.numeric(x_coef), med_coef = med_coef)
}

fit_y_multiview_glmnet_lasso <- function(
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)

  pref <- prefix_multiview_blocks(blocks)
  Z <- cbind(X = X)
  if (!is.null(C) && ncol(C) > 0) Z <- cbind(Z, C)
  Z <- cbind(Z, pref$X)

  pf <- rep(1, ncol(Z))
  pf[colnames(Z) == "X"] <- 0
  if (!is.null(C) && ncol(C) > 0) pf[colnames(Z) %in% colnames(C)] <- 0

  cvfit <- glmnet::cv.glmnet(
    x = as.matrix(Z),
    y = Y,
    family = y_family,
    alpha = glmnet_alpha,
    penalty.factor = pf,
    standardize = TRUE
  )

  beta <- as.matrix(coef(cvfit, s = lambda_choice))[, 1]
  names(beta) <- gsub("`", "", names(beta), fixed = TRUE)

  coef_by_view <- lapply(names(blocks), function(view) {
    dn <- paste0(view, "::", colnames(blocks[[view]]))
    out <- setNames(rep(0, length(dn)), colnames(blocks[[view]]))
    hit <- intersect(dn, names(beta))
    if (length(hit) > 0) out[sub(paste0("^", view, "::"), "", hit)] <- beta[hit]
    out
  })
  names(coef_by_view) <- names(blocks)

  x_coef <- beta["X"]
  if (!is.finite(x_coef)) x_coef <- 0

  list(
    method = paste0("K-view early fusion glmnet (", glmnet_alpha_label(glmnet_alpha), ")"),
    coef_by_view = coef_by_view,
    x_coef = as.numeric(x_coef),
    fit = cvfit,
    feature_map = pref$map
  )
}

fit_y_multiview_cooperative_intermediate <- function(
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  rho = 0.2,
  maxit = 100,
  tol = 1e-04
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  if (!identical(y_family, "gaussian")) {
    stop("K-view cooperative intermediate fusion currently supports y_family='gaussian'.")
  }
  if (rho < 0) stop("rho must be non-negative.")

  W <- cbind(`(Intercept)` = 1, X = X)
  if (!is.null(C) && ncol(C) > 0) {
    C2 <- as.matrix(C)
    colnames(C2) <- paste0("C::", colnames(C2))
    W <- cbind(W, C2)
  }
  W <- as.matrix(W)
  blocks <- lapply(blocks, as.matrix)

  n <- length(Y)
  if (any(vapply(blocks, nrow, integer(1)) != n) || nrow(W) != n) {
    stop("Dimension mismatch in K-view cooperative fit.")
  }

  lambda_info <- lapply(blocks, function(M) {
    choose_block_lambda(Y = Y, W = W, M = M, y_family = y_family, lambda_choice = lambda_choice)
  })
  lambdas <- vapply(lambda_info, function(x) x$lambda, numeric(1))

  betas <- lapply(blocks, function(M) setNames(rep(0, ncol(M)), colnames(M)))
  gamma <- update_gamma_gaussian(Y = Y, W = W, offset = rep(0, n))
  converged <- FALSE
  iter <- 0L
  sqrt_rho <- sqrt(rho)

  for (it in seq_len(maxit)) {
    iter <- it
    betas_old <- betas
    gamma_old <- gamma

    eta_list <- Map(function(M, b) as.numeric(M %*% b), blocks, betas)
    offset_old <- Reduce(`+`, eta_list)
    gamma <- update_gamma_gaussian(Y = Y, W = W, offset = offset_old)
    r <- as.numeric(Y - W %*% gamma)

    for (view in names(blocks)) {
      M <- blocks[[view]]
      other_views <- setdiff(names(blocks), view)
      other_eta <- if (length(other_views) > 0) {
        Reduce(`+`, Map(function(M2, b2) as.numeric(M2 %*% b2), blocks[other_views], betas[other_views]))
      } else {
        rep(0, n)
      }
      target_main <- r - other_eta
      target_agree <- if (length(other_views) > 0) other_eta / length(other_views) else rep(0, n)

      if (rho > 0) {
        x_aug <- rbind(M, sqrt_rho * M)
        y_aug <- c(target_main, sqrt_rho * target_agree)
      } else {
        x_aug <- M
        y_aug <- target_main
      }

      fit_k <- glmnet::glmnet(
        x = as.matrix(x_aug),
        y = y_aug,
        family = "gaussian",
        alpha = 1,
        lambda = lambdas[view],
        intercept = FALSE,
        standardize = TRUE
      )
      b_new <- as.matrix(coef(fit_k, s = lambdas[view]))[-1, 1]
      names(b_new) <- colnames(M)
      betas[[view]] <- b_new
    }

    offset_new <- Reduce(`+`, Map(function(M, b) as.numeric(M %*% b), blocks, betas))
    gamma <- update_gamma_gaussian(Y = Y, W = W, offset = offset_new)

    max_db <- max(unlist(Map(function(a, b) max(abs(a - b)), betas, betas_old)))
    max_dg <- if (length(gamma) > 0) max(abs(gamma - gamma_old)) else 0
    delta <- max(max_db, max_dg)
    if (is.finite(delta) && delta < tol) {
      converged <- TRUE
      break
    }
  }

  x_coef <- if ("X" %in% names(gamma)) as.numeric(gamma["X"]) else 0
  if (!is.finite(x_coef)) x_coef <- 0

  list(
    method = "K-view intermediate fusion cooperative lasso",
    coef_by_view = betas,
    x_coef = x_coef,
    fit = list(
      gamma = gamma,
      lambda = lambdas,
      rho = rho,
      iterations = iter,
      converged = converged,
      lambda_cv = lapply(lambda_info, `[[`, "cvfit")
    )
  )
}

fit_y_multiview_late_fusion_lasso <- function(
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)

  block_fits <- lapply(names(blocks), function(view) {
    fit_block_lasso_for_late(
      Y = Y,
      X = X,
      M = blocks[[view]],
      C = C,
      y_family = y_family,
      lambda_choice = lambda_choice,
      glmnet_alpha = glmnet_alpha,
      prefix = paste0(view, "::")
    )
  })
  names(block_fits) <- names(blocks)

  meta_dat <- as.data.frame(lapply(block_fits, `[[`, "pred"), check.names = FALSE)
  colnames(meta_dat) <- paste0("pred_", names(blocks))
  meta_dat$Y <- Y
  meta_formula <- as.formula(paste("Y ~", paste(setdiff(colnames(meta_dat), "Y"), collapse = " + ")))

  meta_fit <- if (identical(y_family, "gaussian")) {
    lm(meta_formula, data = meta_dat)
  } else {
    glm(meta_formula, data = meta_dat, family = binomial())
  }

  meta_co <- coef(meta_fit)
  weights <- setNames(rep(0, length(blocks)), names(blocks))
  for (view in names(blocks)) {
    nm <- paste0("pred_", view)
    if (nm %in% names(meta_co) && is.finite(meta_co[nm])) weights[view] <- as.numeric(meta_co[nm])
  }

  coef_by_view <- lapply(names(blocks), function(view) block_fits[[view]]$med_coef * weights[view])
  names(coef_by_view) <- names(blocks)
  x_coef <- sum(vapply(names(blocks), function(view) block_fits[[view]]$x_coef * weights[view], numeric(1)))

  list(
    method = paste0("K-view late fusion two-stage glmnet (", glmnet_alpha_label(glmnet_alpha), ")"),
    coef_by_view = coef_by_view,
    x_coef = as.numeric(x_coef),
    fit = list(blocks = lapply(block_fits, `[[`, "cvfit"), meta = meta_fit),
    meta_weights = weights
  )
}

fit_y_multiview_stage <- function(
  Y,
  X,
  blocks,
  C = NULL,
  fusion_mode = c("early", "intermediate", "late"),
  y_family = c("gaussian", "binomial"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  coop_rho = 0.2,
  coop_maxit = 100,
  coop_tol = 1e-04
) {
  fusion_mode <- match.arg(fusion_mode)
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)

  if (identical(fusion_mode, "intermediate")) {
    if (!isTRUE(all.equal(glmnet_alpha, 1))) {
      warning(
        "glmnet_alpha is currently applied only to early and late fusion. ",
        "Intermediate fusion uses the cooperative lasso-style update with alpha=1."
      )
    }
    return(fit_y_multiview_cooperative_intermediate(
      Y = Y,
      X = X,
      blocks = blocks,
      C = C,
      y_family = y_family,
      lambda_choice = lambda_choice,
      rho = coop_rho,
      maxit = coop_maxit,
      tol = coop_tol
    ))
  }

  if (identical(fusion_mode, "late")) {
    return(fit_y_multiview_late_fusion_lasso(
      Y = Y,
      X = X,
      blocks = blocks,
      C = C,
      y_family = y_family,
      lambda_choice = lambda_choice,
      glmnet_alpha = glmnet_alpha
    ))
  }

  fit_y_multiview_glmnet_lasso(
    Y = Y,
    X = X,
    blocks = blocks,
    C = C,
    y_family = y_family,
    lambda_choice = lambda_choice,
    glmnet_alpha = glmnet_alpha
  )
}

