suppressPackageStartupMessages({
  get_script_path <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg) > 0) return(sub("^--file=", "", file_arg[1]))
    if (!is.null(sys.frames()[[1]]$ofile)) return(sys.frames()[[1]]$ofile)
    NULL
  }

  script_path <- get_script_path()
  pkg_dir <- if (!is.null(script_path)) {
    normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }

  if (
    requireNamespace("devtools", quietly = TRUE) &&
      file.exists(file.path(pkg_dir, "DESCRIPTION")) &&
      file.exists(file.path(pkg_dir, "R", "zentangler.R"))
  ) {
    devtools::load_all(pkg_dir, quiet = TRUE)
  } else {
    library(Zentangler)
  }

  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
})

pcl_to_mae_local <- function(pcl) {
  feature_table <- as.matrix(pcl$feature_table)
  feature_metadata <- as.data.frame(pcl$feature_metadata)
  sample_metadata <- as.data.frame(pcl$sample_metadata)
  
  sample_ids <- colnames(feature_table)
  
  if (is.null(rownames(sample_metadata)) ||
      !all(sample_ids %in% rownames(sample_metadata))) {
    rownames(sample_metadata) <- sample_ids
  }
  
  sample_metadata <- sample_metadata[sample_ids, , drop = FALSE]
  
  feature_types <- unique(as.character(feature_metadata$featureType))
  experiments <- list()
  
  for (ft in feature_types) {
    feature_ids <- rownames(feature_metadata)[feature_metadata$featureType == ft]
    
    assay_mat <- feature_table[feature_ids, sample_ids, drop = FALSE]
    
    row_data <- DataFrame(
      featureID = feature_ids,
      featureType = ft
    )
    
    experiments[[make.names(ft)]] <- SummarizedExperiment(
      assays = list(abundance = assay_mat),
      rowData = row_data
    )
  }
  
  col_data <- DataFrame(sample_metadata)
  rownames(col_data) <- sample_ids
  
  MultiAssayExperiment(
    experiments = ExperimentList(experiments),
    colData = col_data
  )
}

simulate_intersim_mediation_mae <- function(
    nsample = 150,
    p.train = 0.7,
    n_mediators_per_view = 10,
    a_strength = 0.5,
    b_strength = 0.5,
    direct_effect = 0.2,
    snr = 1,
    seed = 1
) {
  set.seed(seed)
  
  pcl <- trigger_InterSIM(nsample)
  
  feature_table <- as.matrix(pcl$feature_table)
  feature_metadata <- as.data.frame(pcl$feature_metadata)
  sample_metadata <- as.data.frame(pcl$sample_metadata)
  
  sample_ids <- colnames(feature_table)
  
  A <- as.numeric(factor(sample_metadata$Y)) - 1
  sample_metadata$A <- A
  
  beta <- rep(0, nrow(feature_table))
  names(beta) <- rownames(feature_table)
  
  truth_list <- list()
  feature_types <- unique(as.character(feature_metadata$featureType))
  
  for (ft in feature_types) {
    feature_ids <- rownames(feature_metadata)[feature_metadata$featureType == ft]
    
    n_med <- min(n_mediators_per_view, length(feature_ids))
    true_features <- sample(feature_ids, n_med)
    
    a_vec <- rep(0, length(feature_ids))
    b_vec <- rep(0, length(feature_ids))
    names(a_vec) <- feature_ids
    names(b_vec) <- feature_ids
    
    signs_a <- sample(c(-1, 1), n_med, replace = TRUE)
    signs_b <- sample(c(-1, 1), n_med, replace = TRUE)
    
    a_vec[true_features] <- a_strength * signs_a
    b_vec[true_features] <- b_strength * signs_b
    
    for (feature in true_features) {
      sd_j <- sd(feature_table[feature, ], na.rm = TRUE)
      if (!is.finite(sd_j) || sd_j == 0) sd_j <- 1
      
      feature_table[feature, ] <- feature_table[feature, ] +
        scale(A)[, 1] * a_vec[feature] * sd_j
      
      beta[feature] <- b_vec[feature]
    }
    
    truth_list[[ft]] <- data.frame(
      featureType = ft,
      featureID = feature_ids,
      is_mediator = feature_ids %in% true_features,
      a_true = a_vec,
      b_true = b_vec,
      indirect_true = a_vec * b_vec,
      stringsAsFactors = FALSE
    )
  }
  
  eta <- as.numeric(t(feature_table) %*% beta) + direct_effect * A
  sigma <- sqrt(var(eta) / snr)
  Y <- eta + rnorm(nsample, sd = sigma)
  
  sample_metadata$Y <- Y
  rownames(sample_metadata) <- sample_ids
  
  pcl$feature_table <- feature_table
  pcl$sample_metadata <- sample_metadata
  pcl$feature_metadata <- feature_metadata
  
  train_idx <- sample(seq_len(nsample), size = round(nsample * p.train))
  
  pcl_train <- pcl
  pcl_train$feature_table <- feature_table[, train_idx, drop = FALSE]
  pcl_train$sample_metadata <- sample_metadata[train_idx, , drop = FALSE]
  
  mae_train <- pcl_to_mae_local(pcl_train)
  
  list(
    mae = mae_train,
    truth = do.call(rbind, truth_list)
  )
}

run_intersim_zentangler <- function(
    nrep = 20,
    nsample = 150,
    fusion_modes = c("early", "intermediate", "late"),
    sis_n = 50,
    sis_rank = c("abs_a", "pvalue"),
    screen_method = c("sis", "none"),
    q_threshold = 0.25,
    out_prefix = "zentangler_intersim"
) {
  sis_rank <- match.arg(sis_rank)
  screen_method <- match.arg(screen_method)

  results <- list()
  
  for (r in seq_len(nrep)) {
    message("InterSIM rep ", r, " / ", nrep)
    
    sim <- simulate_intersim_mediation_mae(
      nsample = nsample,
      seed = 10000 + r
    )
    
    mae <- sim$mae
    truth <- sim$truth
    truth$key <- paste(make.names(truth$featureType), truth$featureID, sep = "::")
    
    for (fm in fusion_modes) {
      message("  fusion = ", fm)
      
      fit <- try(
        fit_multiview_parallel_zentangler(
          mae = mae,
          x_var = "A",
          y_var = "Y",
          sis_n = sis_n,
          sis_rank = sis_rank,
          screen_method = screen_method,
          fusion_mode = fm,
          y_family = "gaussian",
          b_inference = "debiased_lasso",
          bootstrap_repeats = 0,
          seed = 20000 + r
        ),
        silent = TRUE
      )
      
      if (inherits(fit, "try-error")) {
        results[[length(results) + 1]] <- data.frame(
          rep = r,
          fusion_mode = fm,
          error = as.character(fit),
          stringsAsFactors = FALSE
        )
        next
      }
      
      tab <- fit$combined_mediators
      tab$key <- paste(tab$omics, tab$mediator, sep = "::")
      tab$is_true <- tab$key %in% truth$key[truth$is_mediator]
      
      active <- is.finite(tab$q_primary) &
        tab$q_primary <= q_threshold &
        is.finite(tab$b) &
        tab$b != 0
      
      top50 <- head(tab[order(tab$abs_score, decreasing = TRUE), ], 50)
      
      results[[length(results) + 1]] <- data.frame(
        rep = r,
        fusion_mode = fm,
        n_tested = nrow(tab),
        n_true = sum(tab$is_true),
        n_active = sum(active),
        true_active = sum(active & tab$is_true),
        false_active = sum(active & !tab$is_true),
        precision = ifelse(sum(active) > 0, sum(active & tab$is_true) / sum(active), NA_real_),
        recall = ifelse(sum(tab$is_true) > 0, sum(active & tab$is_true) / sum(tab$is_true), NA_real_),
        fdr = ifelse(sum(active) > 0, sum(active & !tab$is_true) / sum(active), NA_real_),
        top50_true = sum(top50$is_true),
        top50_precision = sum(top50$is_true) / nrow(top50),
        error = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  
  raw <- do.call(rbind, results)
  
  ok <- is.na(raw$error)

  metric_cols <- c(
    "n_active",
    "true_active",
    "false_active",
    "precision",
    "recall",
    "fdr",
    "top50_true",
    "top50_precision"
  )

  if (any(ok)) {
    summary_list <- lapply(split(raw[ok, , drop = FALSE], raw$fusion_mode[ok]), function(d) {
      vals <- vapply(metric_cols, function(nm) {
        mean(as.numeric(d[[nm]]), na.rm = TRUE)
      }, numeric(1))
      data.frame(
        fusion_mode = unique(d$fusion_mode)[1],
        as.data.frame(as.list(vals), check.names = FALSE),
        stringsAsFactors = FALSE
      )
    })
    summary <- do.call(rbind, summary_list)
    rownames(summary) <- NULL
  } else {
    summary <- data.frame()
  }
  
  out_prefix <- paste0(out_prefix, "_", screen_method)
  write.csv(raw, paste0(out_prefix, "_raw_results.csv"), row.names = FALSE)
  write.csv(summary, paste0(out_prefix, "_summary_results.csv"), row.names = FALSE)
  
  list(
    raw = raw,
    summary = summary,
    settings = list(
      nrep = nrep,
      nsample = nsample,
      fusion_modes = fusion_modes,
      sis_n = sis_n,
      sis_rank = sis_rank,
      screen_method = screen_method,
      q_threshold = q_threshold,
      out_prefix = out_prefix
    )
  )
}

res_sis <- run_intersim_zentangler(
  nrep = 10,
  nsample = 600,
  fusion_modes = c("early", "intermediate", "late"),
  sis_n = 200,
  screen_method = "sis",
  q_threshold = 0.25
)

res_none <- run_intersim_zentangler(
  nrep = 10,
  nsample = 600,
  fusion_modes = c("early", "intermediate", "late"),
  screen_method = "none",
  q_threshold = 0.25
)

res_sis$summary
res_none$summary
