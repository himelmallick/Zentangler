#' Generate simulated multi-omics data for benchmarking
#'
#' Generates simulated multi-omics datasets with specified parameters,
#' including sample size, signal-to-noise ratio, DE probabilities,
#' and response variable generation mode. It also splits the data into
#' training and testing sets.
#'
#' @param nsample Sample size
#' @param snr Signal to noise ratio (for continuous outcomes)
#' @param p.train Train-test split ratio
#' @param de.prob DE probability across all modalities (vector of length 3)
#' @param de.downProb Down-regulation probability (vector of length 3)
#' @param de.facLoc DE factor location (vector of length 3)
#' @param de.facScale DE factor scale (vector of length 3)
#' @param mediation.spike Logical. If TRUE, use mediation-style spike-in
#'   effects instead of Splatter DE factors.
#' @param TIE Total indirect effect per modality (vector of length 3)
#' @param n.pathways Number of true mediating features per modality. If NULL,
#'   uses round(de.prob * number of features) per modality.
#' @param dirPathways_a Direction of exposure-to-feature paths (vector of length 3)
#' @param dirPathways_b Direction of feature-to-outcome paths (vector of length 3)
#' @param effect_type Mediation truth type per modality. One of "spiked",
#'   "null_a", "null_b", or "null_ab".
#' @param ygen.mode Y generation mode: "LM", "Friedman", or "Friedman2"
#' @param outcome.type Outcome type: "continuous", "binary", or "survival"
#' @param surv.hscale Multiplicative scale on hazard for survival outcomes
#'   (hazard_i = surv.hscale * exp(eta_i))
#' @param cens.lower Lower bound for uniform censoring time
#' @param cens.upper Upper bound for uniform censoring time
#' @param nrep Number of repetitions
#' @param seed Random seed
#' @param return.pcl Logical. If TRUE, also return legacy pcl-like InterSIM
#'   objects for debugging. The default FALSE keeps the simulation output
#'   MAE-first for Bioconductor-style workflows.
#'
#' @return List containing training and testing datasets and simulation parameters:
#'   \itemize{
#'     \item trainMae: list of length nrep, each a MultiAssayExperiment with training samples
#'     \item testMae: list of length nrep, each a MultiAssayExperiment with test samples
#'     \item truthDat: per-feature simulation truth tables
#'     \item snr, p.train, de.prob, de.downProb, de.facLoc, de.facScale, nrep, seed
#'   }
#' @export
#' @examples
#' simulated_data <- gen_simmba(nsample = 100)
gen_simmba <- function(nsample,                      # Sample size
                       snr = 1,                      # Signal to noise ratio (continuous)
                       p.train = 0.7,                # Train-test split ratio
                       de.prob = rep(0.1, 3),        # DE probability across modalities
                       de.downProb = rep(0.5, 3),    # Down-regulation probability
                       de.facLoc = rep(1, 3),        # DE factor location
                       de.facScale = rep(0.4, 3),    # DE factor scale
                       mediation.spike = TRUE,        # Use mediation-style spike-in
                       TIE = rep(1 / 3, 3),           # Total indirect effect per modality
                       n.pathways = NULL,             # True mediators per modality
                       dirPathways_a = rep(1, 3),     # Direction of A -> M
                       dirPathways_b = rep(1, 3),     # Direction of M -> Y
                       effect_type = rep("spiked", 3),
                       ygen.mode = c("LM", "Friedman", "Friedman2"),
                       outcome.type = c("continuous", "binary", "survival"),
                       surv.hscale = 1,              # hazard scale for survival
                       cens.lower = 1,               # censoring lower bound
                       cens.upper = 3,               # censoring upper bound
                       nrep = 100,                   # Number of repetitions
                       seed = 1234,                  # Random seed
                       return.pcl = FALSE) {         # Return legacy raw objects?

  set.seed(seed)

  ygen.mode    <- match.arg(ygen.mode)
  outcome.type <- match.arg(outcome.type)
  effect_type <- match.arg(effect_type,
                           choices = c("spiked", "null_a", "null_b", "null_ab"),
                           several.ok = TRUE)
  effect_type <- rep(effect_type, length.out = 3)
  TIE <- rep(TIE, length.out = 3)
  de.prob <- rep(de.prob, length.out = 3)
  de.downProb <- rep(de.downProb, length.out = 3)
  de.facLoc <- rep(de.facLoc, length.out = 3)
  de.facScale <- rep(de.facScale, length.out = 3)
  dirPathways_a <- rep(dirPathways_a, length.out = 3)
  dirPathways_b <- rep(dirPathways_b, length.out = 3)

  # Initialize lists to store training and testing datasets.
  # The MAE lists are the Bioconductor-facing simulation output.
  trainDat <- testDat <- vector("list", nrep)
  names(trainDat) <- names(testDat) <- paste("Rep", 1:nrep, sep = "_")
  trainMae <- testMae <- vector("list", nrep)
  names(trainMae) <- names(testMae) <- paste("Rep", 1:nrep, sep = "_")
  truthDat <- vector("list", nrep)
  names(truthDat) <- paste("Rep", 1:nrep, sep = "_")

  # Loop over the number of repetitions
  for (k in seq_len(nrep)) {

    ## 1. Generate feature data using InterSIM
    pcl <- trigger_InterSIM(n = nsample)  # user-defined wrapper
    feature_types <- unique(as.character(pcl$feature_metadata$featureType))
    nfeature <- table(pcl$feature_metadata$featureType)[feature_types]

    ## Use InterSIM's sample clustering as the binary exposure. This keeps the
    ## InterSIM data structure while replacing its response with a mediation
    ## outcome generated below.
    exposure <- as.numeric(factor(pcl$sample_metadata$Y)) - 1
    pcl$sample_metadata$A <- exposure

    ## 2. Construct mediator truth and optionally spike features by exposure.
    if (isTRUE(mediation.spike)) {
      spike <- simmba_apply_mediation_spikes(
        pcl = pcl,
        exposure = exposure,
        feature_types = feature_types,
        TIE = TIE,
        n.pathways = n.pathways,
        de.prob = de.prob,
        dirPathways_a = dirPathways_a,
        dirPathways_b = dirPathways_b,
        effect_type = effect_type
      )
      pcl <- spike$pcl
      beta0 <- spike$beta0
      truth <- spike$truth
    } else {
      ## Original SIMMBA-style beta generation using Splatter factors.
      de.facs <- vector("list", 3)
      for (i in 1:3) {
        de.facs[[i]] <- getFromNamespace("getLNormFactors", "splatter")(
          n.facs   = nfeature[i],
          sel.prob = de.prob[i],
          neg.prob = de.downProb[i],
          fac.loc  = de.facLoc[i],
          fac.scale = de.facScale[i]
        )
      }

      # Convert to log fold changes (LFCs): this is your beta
      beta0 <- log2(unlist(de.facs))
      names(beta0) <- rownames(pcl$feature_table)
      truth <- data.frame(
        featureID = rownames(pcl$feature_table),
        featureType = pcl$feature_metadata$featureType,
        is_mediator = beta0 != 0,
        a_strength = NA_real_,
        b_strength = beta0,
        TIE_true = NA_real_,
        effect_type = "splatter",
        stringsAsFactors = FALSE
      )
    }

    X <- t(as.matrix(pcl$feature_table))
    beta0 <- beta0[colnames(X)]

    ## 3. Base linear predictor from beta (this is what you want to preserve)
    eta_lin <- as.numeric(X %*% beta0)

    ## Optional non-linear Friedman component (beta0 itself is unchanged)
    if (ygen.mode %in% c("Friedman", "Friedman2")) {

      # Choose 5 non-zero beta features for the Friedman function
      nonzero_index <- which(beta0 != 0)
      if (length(nonzero_index) < 5) {
        stop("Not enough non-zero coefficients to select 5 Friedman features.")
      }
      friedman_index <- sample(nonzero_index, 5)
      X.friedman <- X[, friedman_index, drop = FALSE]

      # Friedman function (embedded here for self-containment)
      f <- function(x) {
        10 * sin(pi * x[, 1] * x[, 2]) +
          20 * (x[, 3] - 0.5)^2 +
          10 * x[, 4] +
          5  * x[, 5]
      }

      Xbeta.friedman <- f(X.friedman)
    }

    ## 4. Final linear predictor eta, used for *all* outcome types
    if (ygen.mode == "LM") {
      eta <- eta_lin
    } else if (ygen.mode == "Friedman") {
      eta <- Xbeta.friedman
    } else if (ygen.mode == "Friedman2") {
      eta <- eta_lin + Xbeta.friedman
    }

    # Store Xbeta in metadata for reference
    pcl$sample_metadata$Xbeta <- eta
    pcl$sample_metadata$TIE_total <- sum(unique(truth$TIE_true[truth$is_mediator]), na.rm = TRUE)
    truth_match <- truth[
      match(pcl$feature_metadata$featureID, truth$featureID),
      c("is_mediator", "a_strength", "b_strength", "TIE_true", "effect_type")
    ]
    pcl$feature_metadata <- cbind(pcl$feature_metadata, truth_match)
    rownames(pcl$feature_metadata) <- pcl$feature_metadata$featureID

    ###################################
    # 5. Generate outcome from eta    #
    ###################################

    if (outcome.type == "continuous") {

      # Y = eta + noise, SNR controlled by snr
      sigma2 <- as.vector(var(eta) / snr)
      Y <- eta + rnorm(nsample) * sqrt(sigma2)
      pcl$sample_metadata$Y <- as.vector(Y)

    } else if (outcome.type == "binary") {

      # Bernoulli with probability p = logit^{-1}(eta)
      p <- plogis(eta)
      Ybin <- rbinom(nsample, size = 1, prob = p)
      pcl$sample_metadata$Y <- Ybin

    } else if (outcome.type == "survival") {

      # Simple survival model with exponential event times:
      # hazard_i = surv.hscale * exp(eta_i)
      h  <- as.vector(surv.hscale * exp(eta))   # rate parameter for rexp
      X0 <- rexp(nsample, rate = h)             # true event times

      # Censoring times
      C <- runif(nsample, cens.lower, cens.upper)

      # Event indicator and observed time
      delta <- ifelse(C >= X0, 1L, 0L)
      Tobs  <- ifelse(C >= X0, X0, C)           # NOTE: uses X0, not X

      pcl$sample_metadata$time   <- Tobs
      pcl$sample_metadata$status <- delta
    }

    #################################
    # 6. Train / test split         #
    #################################

    train <- test <- pcl
    tr.row <- sample.int(nsample, size = round(nsample * p.train), replace = FALSE)

    train$sample_metadata <- pcl$sample_metadata[tr.row, , drop = FALSE]
    test$sample_metadata  <- pcl$sample_metadata[-tr.row, , drop = FALSE]

    train$feature_table <- pcl$feature_table[, tr.row, drop = FALSE]
    test$feature_table  <- pcl$feature_table[, -tr.row, drop = FALSE]

    # Store in lists
    trainDat[[k]] <- train
    testDat[[k]]  <- test
    trainMae[[k]] <- simmba_pcl_to_mae(train)
    testMae[[k]]  <- simmba_pcl_to_mae(test)
    truthDat[[k]] <- truth
  }

  # Return synthetic data and simulation parameters. By default, keep this
  # object MAE-specific; raw pcl-like InterSIM objects are optional debugging
  # payloads, not the public-facing simulation contract.
  out <- list(
    trainMae    = trainMae,
    testMae     = testMae,
    snr         = snr,
    p.train     = p.train,
    de.prob     = de.prob,
    de.downProb = de.downProb,
    de.facLoc   = de.facLoc,
    de.facScale = de.facScale,
    nrep        = nrep,
    seed        = seed,
    mediation.spike = mediation.spike,
    TIE         = TIE,
    n.pathways  = n.pathways,
    dirPathways_a = dirPathways_a,
    dirPathways_b = dirPathways_b,
    effect_type = effect_type,
    truthDat    = truthDat,
    ygen.mode   = ygen.mode,
    outcome.type = outcome.type,
    surv.hscale = surv.hscale,
    cens.lower  = cens.lower,
    cens.upper  = cens.upper
  )
  if (isTRUE(return.pcl)) {
    out$trainDat <- trainDat
    out$testDat <- testDat
  }
  return(out)
}

simmba_sanitize_view_names <- function(x) {
  out <- make.names(as.character(x), unique = TRUE)
  out <- gsub("\\.+", "_", out)
  out <- gsub("^_+|_+$", "", out)
  out[!nzchar(out)] <- paste0("view", seq_len(sum(!nzchar(out))))
  make.unique(out, sep = "_")
}

simmba_pcl_to_mae <- function(pcl, assay_name = "abundance") {
  # Convert the InterSIM/SIMMBA pcl-like object into a MultiAssayExperiment.
  #
  # MAE is the package-facing representation we want for Zentangler:
  # - one experiment per modality/view
  # - features in rows and samples in columns
  # - sample-level X/Y/covariates in colData
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("Package 'MultiAssayExperiment' is required to create SIMMBA MAE output.", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required to create SIMMBA MAE output.", call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required to create SIMMBA MAE output.", call. = FALSE)
  }

  feature_table <- as.matrix(pcl$feature_table)
  feature_metadata <- as.data.frame(pcl$feature_metadata, stringsAsFactors = FALSE)
  sample_metadata <- as.data.frame(pcl$sample_metadata, stringsAsFactors = FALSE)

  if (is.null(rownames(feature_table)) || is.null(colnames(feature_table))) {
    stop("pcl$feature_table must have feature IDs as rownames and sample IDs as colnames.", call. = FALSE)
  }
  if (!all(c("featureID", "featureType") %in% colnames(feature_metadata))) {
    stop("pcl$feature_metadata must contain featureID and featureType columns.", call. = FALSE)
  }

  sample_ids <- colnames(feature_table)
  if (is.null(rownames(sample_metadata)) || !all(sample_ids %in% rownames(sample_metadata))) {
    rownames(sample_metadata) <- sample_ids
  }
  sample_metadata <- sample_metadata[sample_ids, , drop = FALSE]

  feature_metadata <- feature_metadata[
    match(rownames(feature_table), as.character(feature_metadata$featureID)),
    , drop = FALSE
  ]
  rownames(feature_metadata) <- rownames(feature_table)

  raw_views <- unique(as.character(feature_metadata$featureType))
  raw_views <- raw_views[!is.na(raw_views) & nzchar(raw_views)]
  safe_views <- simmba_sanitize_view_names(raw_views)

  experiments <- vector("list", length(raw_views))
  names(experiments) <- safe_views
  for (i in seq_along(raw_views)) {
    feature_ids <- rownames(feature_table)[as.character(feature_metadata$featureType) == raw_views[i]]
    if (length(feature_ids) == 0) next

    assay_mat <- feature_table[feature_ids, sample_ids, drop = FALSE]
    row_dat <- S4Vectors::DataFrame(feature_metadata[feature_ids, , drop = FALSE])

    experiments[[safe_views[i]]] <- SummarizedExperiment::SummarizedExperiment(
      assays = stats::setNames(list(assay_mat), assay_name),
      rowData = row_dat
    )
  }

  col_dat <- S4Vectors::DataFrame(sample_metadata)
  rownames(col_dat) <- sample_ids

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(experiments),
    colData = col_dat
  )
}

simmba_apply_mediation_spikes <- function(pcl,
                                          exposure,
                                          feature_types,
                                          TIE,
                                          n.pathways,
                                          de.prob,
                                          dirPathways_a,
                                          dirPathways_b,
                                          effect_type) {
  feature_table <- as.matrix(pcl$feature_table)
  beta0 <- rep(0, nrow(feature_table))
  names(beta0) <- rownames(feature_table)

  truth_list <- vector("list", length(feature_types))

  for (i in seq_along(feature_types)) {
    type_i <- feature_types[i]
    feature_ids <- rownames(pcl$feature_metadata)[
      as.character(pcl$feature_metadata$featureType) == type_i
    ]
    k <- simmba_n_pathways(
      n_features = length(feature_ids),
      n.pathways = n.pathways,
      de.prob = de.prob[i],
      modality_index = i
    )

    selected <- if (k > 0) sample(feature_ids, size = k, replace = FALSE) else character(0)
    base_strength <- if (k > 0 && TIE[i] != 0) sqrt(abs(TIE[i]) / k) else 0
    a_strength <- base_strength * dirPathways_a[i]
    b_strength <- base_strength * dirPathways_b[i]

    has_a <- effect_type[i] %in% c("spiked", "null_b")
    has_b <- effect_type[i] %in% c("spiked", "null_a")
    TIE_true <- if (has_a && has_b) TIE[i] * dirPathways_a[i] * dirPathways_b[i] else 0

    if (has_a && length(selected) > 0) {
      feature_table[selected, ] <- simmba_additive_spike(
        feature_table[selected, , drop = FALSE],
        exposure = exposure,
        a_strength = a_strength
      )
    }

    if (has_b && length(selected) > 0) {
      beta0[selected] <- b_strength
    }

    truth_list[[i]] <- data.frame(
      featureID = feature_ids,
      featureType = type_i,
      is_mediator = feature_ids %in% selected,
      a_strength = ifelse(feature_ids %in% selected, a_strength, 0),
      b_strength = ifelse(feature_ids %in% selected, ifelse(has_b, b_strength, 0), 0),
      TIE_requested = TIE[i],
      TIE_true = ifelse(feature_ids %in% selected, TIE_true, 0),
      effect_type = effect_type[i],
      stringsAsFactors = FALSE
    )
  }

  pcl$feature_table <- as.data.frame(feature_table)
  list(
    pcl = pcl,
    beta0 = beta0,
    truth = do.call(rbind, truth_list)
  )
}

simmba_n_pathways <- function(n_features, n.pathways, de.prob, modality_index) {
  if (is.null(n.pathways)) {
    return(max(1, round(n_features * de.prob)))
  }
  n.pathways <- rep(n.pathways, length.out = 3)
  k <- n.pathways[modality_index]
  if (k > n_features) {
    stop("n.pathways cannot exceed the number of features in modality ", modality_index, ".")
  }
  k
}

simmba_additive_spike <- function(feature_block, exposure, a_strength) {
  exposure_scaled <- as.numeric(scale(exposure))
  exposure_scaled[is.na(exposure_scaled)] <- 0

  out <- feature_block
  for (j in seq_len(nrow(out))) {
    sd_j <- stats::sd(out[j, ], na.rm = TRUE)
    if (is.na(sd_j) || sd_j == 0) {
      sd_j <- 1
    }
    out[j, ] <- out[j, ] + exposure_scaled * a_strength * sd_j
  }
  out
}
