#!/usr/bin/env Rscript
# Run ASD sequential Zentangler negative-control permutations.
#
# This script runs observed fits plus optional negative controls:
#   - permute_y: breaks terminal mediator -> outcome alignment
#   - permute_x: breaks exposure -> first mediator alignment
#   - permute_first: permutes the first route assay across samples
#   - permute_terminal: permutes the terminal route assay across samples
#   - permute_all_assays: permutes every assay independently across samples
#
# Example:
# Rscript inst/scripts/run_asd_sequential_permutation_controls.R \
#   --out-dir sequential_asd_permutation_controls \
#   --contrasts T1 \
#   --symptoms ABC \
#   --routes kos_to_fecal,species_to_fecal \
#   --null-types permute_y,permute_terminal \
#   --n-perm 100 \
#   --sis-n 100 \
#   --min-abs-cor 0.05 \
#   --cor-q-threshold 0.75 \
#   --fusion-mode early \
#   --glmnet-alpha 1 \
#   --lambda-choice lambda.min

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

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

split_arg <- function(x, default = character()) {
  if (is.null(x) || length(x) == 0L || !nzchar(as.character(x))) return(default)
  z <- unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE)
  trimws(z[nzchar(trimws(z))])
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(as.character(x))) return(default)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

as_num_vec <- function(x, default) {
  z <- split_arg(x)
  if (length(z) == 0L) return(default)
  as.numeric(z)
}

as_int_vec <- function(x, default) {
  as.integer(as_num_vec(x, default = default))
}

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

safe_num <- function(x) gsub("[.]", "p", as.character(x))

bind_fill <- function(dfs) {
  dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, dfs)
  if (length(dfs) == 0L) return(data.frame())
  all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(x) {
    miss <- setdiff(all_cols, names(x))
    for (m in miss) x[[m]] <- NA
    x[, all_cols, drop = FALSE]
  })
  do.call(rbind, dfs)
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

source_asd_helpers <- function(package_dir = NULL) {
  candidates <- c(
    if (!is.null(package_dir) && nzchar(package_dir)) file.path(package_dir, "vignettes", "asd-data-helpers.R") else character(),
    file.path(getwd(), "vignettes", "asd-data-helpers.R"),
    file.path(dirname(getwd()), "vignettes", "asd-data-helpers.R")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("Could not find vignettes/asd-data-helpers.R. Run from the Zentangler repo or pass --package-dir /path/to/Zentangler.")
  }
  source(hit[[1L]])
  invisible(hit[[1L]])
}

make_default_routes <- function() {
  list(
    species_to_kos = c("species", "kos"),
    species_to_fecal = c("species", "fecal_metabolites"),
    species_to_plasma = c("species", "plasma_metabolites"),
    kos_to_fecal = c("kos", "fecal_metabolites"),
    kos_to_plasma = c("kos", "plasma_metabolites"),
    fecal_to_plasma = c("fecal_metabolites", "plasma_metabolites"),
    species_to_kos_to_fecal = c("species", "kos", "fecal_metabolites"),
    species_to_kos_to_plasma = c("species", "kos", "plasma_metabolites"),
    species_to_fecal_to_plasma = c("species", "fecal_metabolites", "plasma_metabolites"),
    kos_to_fecal_to_plasma = c("kos", "fecal_metabolites", "plasma_metabolites")
  )
}

get_q_values <- function(paths) {
  if ("q_primary_global" %in% names(paths)) return(suppressWarnings(as.numeric(paths$q_primary_global)))
  if ("q_primary" %in% names(paths)) return(suppressWarnings(as.numeric(paths$q_primary)))
  rep(NA_real_, nrow(paths))
}

summarize_paths <- function(paths, q_threshold = 0.25) {
  if (!is.data.frame(paths) || nrow(paths) == 0L) {
    return(data.frame(
      n_paths = 0L,
      n_significant_paths = 0L,
      best_q = NA_real_,
      best_p = NA_real_,
      max_abs_sequential_score = NA_real_
    ))
  }
  q <- get_q_values(paths)
  p <- if ("p_primary" %in% names(paths)) suppressWarnings(as.numeric(paths$p_primary)) else rep(NA_real_, nrow(paths))
  score <- if ("abs_sequential_score" %in% names(paths)) {
    suppressWarnings(as.numeric(paths$abs_sequential_score))
  } else {
    rep(NA_real_, nrow(paths))
  }
  data.frame(
    n_paths = nrow(paths),
    n_significant_paths = sum(is.finite(q) & q <= q_threshold),
    best_q = {
      v <- suppressWarnings(min(q, na.rm = TRUE))
      if (is.infinite(v)) NA_real_ else v
    },
    best_p = {
      v <- suppressWarnings(min(p, na.rm = TRUE))
      if (is.infinite(v)) NA_real_ else v
    },
    max_abs_sequential_score = {
      v <- suppressWarnings(max(score, na.rm = TRUE))
      if (is.infinite(v)) NA_real_ else v
    }
  )
}

permute_mae_coldata <- function(mae, var, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mae_null <- mae
  cd <- MultiAssayExperiment::colData(mae_null)
  cd[[var]] <- sample(as.vector(cd[[var]]))
  MultiAssayExperiment::colData(mae_null) <- cd
  mae_null
}

permute_one_assay_samples <- function(mae, assay_name, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mae_null <- mae
  se <- mae_null[[assay_name]]
  assay_names <- SummarizedExperiment::assayNames(se)
  first_assay <- assay_names[[1L]]
  mat <- SummarizedExperiment::assay(se, first_assay)
  old_colnames <- colnames(mat)
  mat <- mat[, sample(seq_len(ncol(mat))), drop = FALSE]
  colnames(mat) <- old_colnames
  SummarizedExperiment::assay(se, first_assay) <- mat
  mae_null[[assay_name]] <- se
  mae_null
}

permute_all_assays_samples <- function(mae, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mae_null <- mae
  for (assay_name in names(MultiAssayExperiment::experiments(mae_null))) {
    mae_null <- permute_one_assay_samples(mae_null, assay_name = assay_name)
  }
  mae_null
}

make_null_mae <- function(mae, null_type, route_views, x_var = "X", y_var = "Y", seed = NULL) {
  if (identical(null_type, "observed")) return(mae)
  if (identical(null_type, "permute_y")) return(permute_mae_coldata(mae, y_var, seed = seed))
  if (identical(null_type, "permute_x")) return(permute_mae_coldata(mae, x_var, seed = seed))
  if (identical(null_type, "permute_first")) return(permute_one_assay_samples(mae, route_views[[1L]], seed = seed))
  if (identical(null_type, "permute_terminal")) return(permute_one_assay_samples(mae, tail(route_views, 1L), seed = seed))
  if (identical(null_type, "permute_all_assays")) return(permute_all_assays_samples(mae, seed = seed))
  stop("Unknown null_type: ", null_type)
}

make_iteration_dir <- function(base_dir, contrast, symptom, route_name, sis_n,
                               min_abs_cor, cor_q_threshold, fusion_mode,
                               glmnet_alpha, lambda_choice, null_type, perm_id) {
  iter_name <- paste(
    safe_name(contrast),
    safe_name(symptom),
    safe_name(route_name),
    paste0("sis", sis_n),
    paste0("cor", safe_num(min_abs_cor)),
    paste0("corq", safe_num(cor_q_threshold)),
    safe_name(fusion_mode),
    paste0("alpha", safe_num(glmnet_alpha)),
    safe_name(lambda_choice),
    safe_name(null_type),
    paste0("perm", sprintf("%04d", as.integer(perm_id))),
    sep = "__"
  )
  iter_dir <- file.path(base_dir, iter_name)
  dir.create(iter_dir, recursive = TRUE, showWarnings = FALSE)
  iter_dir
}

run_one <- function(mae_obs, route_views, meta, args) {
  iter_dir <- make_iteration_dir(
    base_dir = args$out_dir,
    contrast = meta$contrast,
    symptom = meta$symptom,
    route_name = meta$route_name,
    sis_n = meta$sis_n,
    min_abs_cor = meta$min_abs_cor,
    cor_q_threshold = meta$cor_q_threshold,
    fusion_mode = meta$fusion_mode,
    glmnet_alpha = meta$glmnet_alpha,
    lambda_choice = meta$lambda_choice,
    null_type = meta$null_type,
    perm_id = meta$perm_id
  )
  meta$iteration_dir <- iter_dir
  utils::write.csv(meta, file.path(iter_dir, "grid_settings.csv"), row.names = FALSE)

  cat(
    "Running:", meta$contrast, "|", meta$symptom, "|", meta$route_name,
    "|", meta$null_type, "perm", meta$perm_id,
    "| sis=", meta$sis_n,
    "| cor>=", meta$min_abs_cor,
    "| corq<=", meta$cor_q_threshold,
    "|", meta$fusion_mode,
    "| alpha=", meta$glmnet_alpha, "\n"
  )

  mae_run <- make_null_mae(
    mae_obs,
    null_type = meta$null_type,
    route_views = route_views,
    x_var = args$x_var,
    y_var = args$y_var,
    seed = meta$seed
  )

  fit <- fit_sequential_zentangler(
    mae_run,
    x_var = args$x_var,
    y_var = args$y_var,
    path_templates = stats::setNames(list(route_views), meta$route_name),
    sis_n = meta$sis_n,
    cor_method = args$cor_method,
    min_abs_cor = meta$min_abs_cor,
    cor_q_threshold = meta$cor_q_threshold,
    fusion_mode = meta$fusion_mode,
    lambda_choice = meta$lambda_choice,
    glmnet_alpha = meta$glmnet_alpha,
    y_family = args$y_family,
    seed = meta$seed
  )

  paths <- zentangler_sequential_paths(fit)
  terminals <- zentangler_sequential_terminals(fit)
  model_summary <- summarize_sequential_zentangler(fit)$model_summary
  path_summary <- summarize_paths(paths, q_threshold = args$q_threshold)

  if (nrow(paths) > 0L) paths <- cbind(meta[rep(1L, nrow(paths)), , drop = FALSE], paths)
  if (nrow(terminals) > 0L) terminals <- cbind(meta[rep(1L, nrow(terminals)), , drop = FALSE], terminals)
  if (nrow(model_summary) > 0L) {
    model_summary <- cbind(meta[rep(1L, nrow(model_summary)), , drop = FALSE], model_summary)
  }
  path_summary <- cbind(meta, path_summary)

  utils::write.csv(paths, file.path(iter_dir, "sequential_paths.csv"), row.names = FALSE)
  utils::write.csv(terminals, file.path(iter_dir, "terminal_effects.csv"), row.names = FALSE)
  utils::write.csv(model_summary, file.path(iter_dir, "model_summary.csv"), row.names = FALSE)
  utils::write.csv(path_summary, file.path(iter_dir, "path_summary.csv"), row.names = FALSE)

  if (isTRUE(args$save_fit)) {
    saveRDS(
      list(fit = fit, paths = paths, terminals = terminals, model_summary = model_summary, path_summary = path_summary),
      file.path(iter_dir, "fit_result.rds")
    )
  }

  list(paths = paths, terminals = terminals, model_summary = model_summary, summary = path_summary, error = data.frame())
}

args_in <- parse_args(commandArgs(trailingOnly = TRUE))

package_dir <- args_in$`package-dir` %||% args_in$package_dir %||% getwd()
load_zentangler_sources(package_dir = package_dir)
source_asd_helpers(package_dir = package_dir)

data_dir <- args_in$`data-dir` %||% args_in$data_dir %||% Sys.getenv("ZENTANGLER_DATA_DIR", unset = "")
asd <- if (nzchar(data_dir)) load_asd_example_data(data_dir = data_dir) else load_asd_example_data()
omics <- asd$omics
symptoms <- asd$symptoms

all_routes <- make_default_routes()
route_names <- split_arg(args_in$routes, default = names(all_routes))
missing_routes <- setdiff(route_names, names(all_routes))
if (length(missing_routes) > 0L) stop("Unknown route(s): ", paste(missing_routes, collapse = ", "))
routes <- all_routes[route_names]

symptom_names <- split_arg(args_in$symptoms, default = sort(unique(as.character(symptoms$Symptom))))
contrasts <- split_arg(args_in$contrasts, default = c("T1"))
null_types <- split_arg(args_in$`null-types` %||% args_in$null_types, default = c("permute_y", "permute_terminal"))
allowed_null_types <- c("permute_y", "permute_x", "permute_first", "permute_terminal", "permute_all_assays")
unknown_nulls <- setdiff(null_types, allowed_null_types)
if (length(unknown_nulls) > 0L) stop("Unknown null type(s): ", paste(unknown_nulls, collapse = ", "))

script_args <- list(
  out_dir = args_in$`out-dir` %||% args_in$out_dir %||% "sequential_asd_permutation_controls",
  x_var = args_in$`x-var` %||% args_in$x_var %||% "X",
  y_var = args_in$`y-var` %||% args_in$y_var %||% "Y",
  y_family = args_in$`y-family` %||% args_in$y_family %||% "gaussian",
  cor_method = args_in$`cor-method` %||% args_in$cor_method %||% "spearman",
  q_threshold = as.numeric(args_in$`q-threshold` %||% args_in$q_threshold %||% 0.25),
  save_fit = as_bool(args_in$`save-fit` %||% args_in$save_fit, default = FALSE)
)
dir.create(script_args$out_dir, recursive = TRUE, showWarnings = FALSE)

n_perm <- as.integer(args_in$`n-perm` %||% args_in$n_perm %||% 100L)
base_seed <- as.integer(args_in$seed %||% 1L)
sis_n_values <- as_int_vec(args_in$`sis-n` %||% args_in$sis_n, default = 100L)
min_abs_cor_values <- as_num_vec(args_in$`min-abs-cor` %||% args_in$min_abs_cor, default = 0.05)
cor_q_values <- as_num_vec(args_in$`cor-q-threshold` %||% args_in$cor_q_threshold, default = 0.75)
fusion_modes <- split_arg(args_in$`fusion-mode` %||% args_in$fusion_mode, default = "early")
glmnet_alpha_values <- as_num_vec(args_in$`glmnet-alpha` %||% args_in$glmnet_alpha, default = 1)
lambda_choices <- split_arg(args_in$`lambda-choice` %||% args_in$lambda_choice, default = "lambda.min")

grid <- expand.grid(
  contrast = contrasts,
  symptom = symptom_names,
  route_name = names(routes),
  sis_n = sis_n_values,
  min_abs_cor = min_abs_cor_values,
  cor_q_threshold = cor_q_values,
  fusion_mode = fusion_modes,
  glmnet_alpha = glmnet_alpha_values,
  lambda_choice = lambda_choices,
  stringsAsFactors = FALSE
)

utils::write.csv(grid, file.path(script_args$out_dir, "run_grid_base.csv"), row.names = FALSE)

all_results <- list()
rr <- 1L
for (ii in seq_len(nrow(grid))) {
  row <- grid[ii, , drop = FALSE]
  route_views <- routes[[row$route_name]]
  mae_obs <- tryCatch(
    build_asd_mae(omics, symptoms, contrast = row$contrast, symptom = row$symptom),
    error = function(e) e
  )
  if (inherits(mae_obs, "error")) {
    err <- data.frame(row, null_type = "build_mae", perm_id = NA_integer_, seed = NA_integer_, error = conditionMessage(mae_obs))
    all_results[[rr]] <- list(paths = data.frame(), terminals = data.frame(), model_summary = data.frame(), summary = data.frame(), error = err)
    rr <- rr + 1L
    next
  }

  run_plan <- rbind(
    data.frame(row, null_type = "observed", perm_id = 0L, seed = base_seed, stringsAsFactors = FALSE),
    do.call(rbind, lapply(null_types, function(nt) {
      data.frame(row, null_type = nt, perm_id = seq_len(n_perm), seed = base_seed + seq_len(n_perm), stringsAsFactors = FALSE)
    }))
  )

  for (jj in seq_len(nrow(run_plan))) {
    meta <- run_plan[jj, , drop = FALSE]
    meta$view_route <- paste(route_views, collapse = " -> ")
    one <- tryCatch(
      run_one(mae_obs = mae_obs, route_views = route_views, meta = meta, args = script_args),
      error = function(e) {
        iter_dir <- make_iteration_dir(
          base_dir = script_args$out_dir,
          contrast = meta$contrast,
          symptom = meta$symptom,
          route_name = meta$route_name,
          sis_n = meta$sis_n,
          min_abs_cor = meta$min_abs_cor,
          cor_q_threshold = meta$cor_q_threshold,
          fusion_mode = meta$fusion_mode,
          glmnet_alpha = meta$glmnet_alpha,
          lambda_choice = meta$lambda_choice,
          null_type = meta$null_type,
          perm_id = meta$perm_id
        )
        meta$iteration_dir <- iter_dir
        err <- cbind(meta, data.frame(error = conditionMessage(e), stringsAsFactors = FALSE))
        utils::write.csv(err, file.path(iter_dir, "error.csv"), row.names = FALSE)
        list(paths = data.frame(), terminals = data.frame(), model_summary = data.frame(), summary = data.frame(), error = err)
      }
    )
    all_results[[rr]] <- one
    rr <- rr + 1L
  }
}

all_paths <- bind_fill(lapply(all_results, `[[`, "paths"))
all_terminals <- bind_fill(lapply(all_results, `[[`, "terminals"))
all_model_summaries <- bind_fill(lapply(all_results, `[[`, "model_summary"))
all_summaries <- bind_fill(lapply(all_results, `[[`, "summary"))
all_errors <- bind_fill(lapply(all_results, `[[`, "error"))

utils::write.csv(all_paths, file.path(script_args$out_dir, "ALL_sequential_paths.csv"), row.names = FALSE)
utils::write.csv(all_terminals, file.path(script_args$out_dir, "ALL_terminal_effects.csv"), row.names = FALSE)
utils::write.csv(all_model_summaries, file.path(script_args$out_dir, "ALL_model_summaries.csv"), row.names = FALSE)
utils::write.csv(all_summaries, file.path(script_args$out_dir, "ALL_path_summaries.csv"), row.names = FALSE)
utils::write.csv(all_errors, file.path(script_args$out_dir, "ALL_errors.csv"), row.names = FALSE)

if (nrow(all_summaries) > 0L) {
  observed <- all_summaries[all_summaries$null_type == "observed", , drop = FALSE]
  nulls <- all_summaries[all_summaries$null_type != "observed", , drop = FALSE]
  empirical <- list()
  ee <- 1L
  key_cols <- c(
    "contrast", "symptom", "route_name", "sis_n", "min_abs_cor",
    "cor_q_threshold", "fusion_mode", "glmnet_alpha", "lambda_choice"
  )
  for (i in seq_len(nrow(observed))) {
    obs <- observed[i, , drop = FALSE]
    idx_base <- rep(TRUE, nrow(nulls))
    for (cc in key_cols) idx_base <- idx_base & as.character(nulls[[cc]]) == as.character(obs[[cc]])
    for (nt in unique(nulls$null_type[idx_base])) {
      nn <- nulls[idx_base & nulls$null_type == nt, , drop = FALSE]
      empirical[[ee]] <- cbind(
        obs[, key_cols, drop = FALSE],
        data.frame(
          null_type = nt,
          obs_n_significant_paths = obs$n_significant_paths,
          obs_best_q = obs$best_q,
          n_null = nrow(nn),
          null_mean_significant_paths = mean(nn$n_significant_paths, na.rm = TRUE),
          null_max_significant_paths = max(nn$n_significant_paths, na.rm = TRUE),
          empirical_p_count = mean(nn$n_significant_paths >= obs$n_significant_paths, na.rm = TRUE),
          empirical_p_best_q = mean(nn$best_q <= obs$best_q, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      )
      ee <- ee + 1L
    }
  }
  empirical_tab <- bind_fill(empirical)
  utils::write.csv(empirical_tab, file.path(script_args$out_dir, "ALL_empirical_negative_control_summary.csv"), row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(script_args$out_dir, "sessionInfo.txt"))
cat("Done. Output directory:", normalizePath(script_args$out_dir, mustWork = FALSE), "\n")
