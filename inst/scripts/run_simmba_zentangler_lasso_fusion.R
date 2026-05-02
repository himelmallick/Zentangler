# SIMMBA / InterSIM runner for the Zentangler package
#
# This script is installed with the package under inst/scripts. It generates a
# MultiAssayExperiment simulation and runs early, intermediate, and late fusion.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Zentangler)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})

NSAMPLE <- 180L
NREP <- 1L
REP_INDEX <- 1L
P_TRAIN <- 0.70
SEED <- 1234L
YGEN_MODE <- "LM"
OUTCOME_TYPE <- "continuous"
FUSION_MODES <- c("early", "intermediate", "late")
SIS_N_PER_VIEW <- 50L
SIS_RANK <- "abs_a"
SCREEN_METHOD <- "sis" # "sis" or "none"; "none" passes all filtered mediators to the B-stage
LAMBDA_CHOICE <- "lambda.1se"
COOP_RHO <- 0.20
B_INFERENCE <- "debiased_lasso"
DEBIAS_MAX_TARGETS <- 200L
BOOTSTRAP_REPEATS <- 10L
BOOTSTRAP_CI_LEVEL <- 0.95
BOOTSTRAP_SEED <- 9001L
Q_THRESHOLD <- 0.25
TOP_N <- 50L
OUT_DIR <- file.path(getwd(), "zentangler_simmba_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

slugify <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", tolower(as.character(x)))
  gsub("^_+|_+$", "", out)
}

mae_view_map <- function(mae) {
  exps <- MultiAssayExperiment::experiments(mae)
  out <- lapply(names(exps), function(view) {
    feature_type <- view
    exp_obj <- exps[[view]]
    if (inherits(exp_obj, "SummarizedExperiment")) {
      rd <- as.data.frame(SummarizedExperiment::rowData(exp_obj), stringsAsFactors = FALSE)
      if ("featureType" %in% colnames(rd)) {
        ft <- unique(as.character(rd$featureType))
        ft <- ft[!is.na(ft) & nzchar(ft)]
        if (length(ft) > 0) feature_type <- ft[1]
      }
    }
    data.frame(featureType = feature_type, view = view, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

feature_counts <- function(mae) {
  exps <- MultiAssayExperiment::experiments(mae)
  out <- vapply(exps, function(exp_obj) {
    if (inherits(exp_obj, "SummarizedExperiment")) {
      nrow(SummarizedExperiment::assay(exp_obj, 1))
    } else {
      ncol(as.matrix(exp_obj))
    }
  }, integer(1))
  out
}

add_primary_inference <- function(tab) {
  if (is.null(tab) || nrow(tab) == 0) return(tab)
  if (!("p_primary" %in% colnames(tab))) {
    tab$p_primary <- if ("joint_p_ab" %in% colnames(tab)) tab$joint_p_ab else NA_real_
  }
  if (!("q_primary" %in% colnames(tab))) {
    tab$q_primary <- stats::p.adjust(tab$p_primary, method = "BH")
  }
  tab
}

is_active_mediator <- function(tab, q_threshold = Q_THRESHOLD) {
  if (is.null(tab) || nrow(tab) == 0) return(logical(0))
  is.finite(tab$b) & tab$b != 0 & is.finite(tab$q_primary) & tab$q_primary <= q_threshold
}

annotate_truth <- function(tab, truth, view_map) {
  tab <- add_primary_inference(tab)
  truth2 <- merge(as.data.frame(truth), view_map, by = "featureType", all.x = TRUE, sort = FALSE)
  truth_keep <- truth2[, intersect(
    c("view", "featureID", "is_mediator", "a_strength", "b_strength", "TIE_true", "effect_type"),
    colnames(truth2)
  ), drop = FALSE]
  out <- merge(
    tab,
    truth_keep,
    by.x = c("omics", "mediator"),
    by.y = c("view", "featureID"),
    all.x = TRUE,
    sort = FALSE
  )
  out$is_mediator[is.na(out$is_mediator)] <- FALSE
  out
}

summarize_truth_recovery <- function(tab, fusion_mode, q_threshold = Q_THRESHOLD, top_n = TOP_N) {
  tab <- add_primary_inference(tab)
  active <- is_active_mediator(tab, q_threshold = q_threshold)
  truth_flag <- if ("is_mediator" %in% colnames(tab)) as.logical(tab$is_mediator) else rep(FALSE, nrow(tab))
  top_tab <- tab[order(tab$abs_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  top_tab <- top_tab[seq_len(min(top_n, nrow(top_tab))), , drop = FALSE]
  data.frame(
    fusion_mode = fusion_mode,
    n_tested = nrow(tab),
    n_active = sum(active),
    n_truth = sum(truth_flag, na.rm = TRUE),
    true_active = sum(active & truth_flag, na.rm = TRUE),
    false_active = sum(active & !truth_flag, na.rm = TRUE),
    top_n = top_n,
    top_n_true = sum(as.logical(top_tab$is_mediator), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

cat("Generating SIMMBA MAE data...\n")
sim <- gen_simmba(
  nsample = NSAMPLE,
  p.train = P_TRAIN,
  ygen.mode = YGEN_MODE,
  outcome.type = OUTCOME_TYPE,
  nrep = NREP,
  seed = SEED
)

mae_train <- sim$trainMae[[REP_INDEX]]
truth <- sim$truthDat[[REP_INDEX]]
view_names <- names(MultiAssayExperiment::experiments(mae_train))
view_map <- mae_view_map(mae_train)

write_csv(view_map, file.path(OUT_DIR, "simmba_view_map.csv"))
write_csv(truth, file.path(OUT_DIR, "simmba_truth.csv"))
saveRDS(mae_train, file.path(OUT_DIR, "simmba_train_mae.rds"))

cat("Views detected: ", paste(view_names, collapse = ", "), "\n", sep = "")
cat("Features per view:\n")
print(feature_counts(mae_train))

summary_list <- list()
for (fusion_mode in FUSION_MODES) {
  cat("\nRunning Zentangler fusion mode: ", fusion_mode, "\n", sep = "")
  fit <- fit_multiview_parallel_zentangler(
    mae = mae_train,
    x_var = "A",
    y_var = "Y",
    view_names = view_names,
    sis_n = SIS_N_PER_VIEW,
    sis_rank = SIS_RANK,
    screen_method = SCREEN_METHOD,
    fusion_mode = fusion_mode,
    y_family = if (OUTCOME_TYPE == "binary") "binomial" else "gaussian",
    lambda_choice = LAMBDA_CHOICE,
    b_inference = B_INFERENCE,
    debias_max_targets = DEBIAS_MAX_TARGETS,
    coop_rho = COOP_RHO,
    seed = SEED,
    return_fits = TRUE,
    bootstrap_repeats = BOOTSTRAP_REPEATS,
    bootstrap_ci_level = BOOTSTRAP_CI_LEVEL,
    bootstrap_seed = BOOTSTRAP_SEED
  )

  all_tab <- annotate_truth(fit$combined_mediators, truth = truth, view_map = view_map)
  all_tab <- all_tab[order(all_tab$q_primary, -all_tab$abs_score, na.last = TRUE), , drop = FALSE]
  active_tab <- all_tab[is_active_mediator(all_tab), , drop = FALSE]
  top_tab <- all_tab[order(all_tab$abs_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  top_tab <- top_tab[seq_len(min(TOP_N, nrow(top_tab))), , drop = FALSE]

  prefix <- paste0("simmba_", slugify(fusion_mode))
  saveRDS(fit, file.path(OUT_DIR, paste0(prefix, "_fit.rds")))
  write_csv(all_tab, file.path(OUT_DIR, paste0(prefix, "_all_mediators.csv")))
  write_csv(active_tab, file.path(OUT_DIR, paste0(prefix, "_active_q", Q_THRESHOLD, ".csv")))
  write_csv(top_tab, file.path(OUT_DIR, paste0(prefix, "_top", TOP_N, ".csv")))
  write_csv(fit$effect_decomposition, file.path(OUT_DIR, paste0(prefix, "_effect_decomposition.csv")))

  summary_list[[fusion_mode]] <- summarize_truth_recovery(all_tab, fusion_mode = fusion_mode)
}

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL
write_csv(summary_df, file.path(OUT_DIR, "simmba_fusion_truth_recovery_summary.csv"))
cat("\nFinished. Results folder: ", OUT_DIR, "\n", sep = "")
print(summary_df)
