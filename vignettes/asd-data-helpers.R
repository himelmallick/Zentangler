find_asd_data_dir <- function() {
  env_dir <- Sys.getenv("ZENTANGLER_DATA_DIR", unset = "")
  candidates <- c(
    env_dir,
    "data",
    file.path("..", "data"),
    file.path(getwd(), "data"),
    file.path(dirname(getwd()), "data")
  )
  candidates <- candidates[nzchar(candidates)]
  hits <- candidates[dir.exists(candidates)]
  if (length(hits) == 0) {
    stop(
      "Could not find the ASD data directory. Set ZENTANGLER_DATA_DIR to ",
      "the folder containing ASD_Species.RData, ASD_KOs.RData, ",
      "ASD_fecalMetabolites.RData, ASD_plasmaMetabolites.RData, and ",
      "ASD_Symptoms.RData.",
      call. = FALSE
    )
  }
  normalizePath(hits[[1]], mustWork = TRUE)
}

read_asd_pcl <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("pcl", envir = e, inherits = FALSE)) {
    stop("Expected object `pcl` in ", path, call. = FALSE)
  }
  get("pcl", envir = e)
}

load_asd_example_data <- function(data_dir = find_asd_data_dir()) {
  omics_files <- c(
    species = "ASD_Species.RData",
    kos = "ASD_KOs.RData",
    fecal_metabolites = "ASD_fecalMetabolites.RData",
    plasma_metabolites = "ASD_plasmaMetabolites.RData"
  )

  omics <- lapply(file.path(data_dir, omics_files), read_asd_pcl)
  names(omics) <- names(omics_files)
  symptoms <- read_asd_pcl(file.path(data_dir, "ASD_Symptoms.RData"))$features
  list(data_dir = data_dir, omics = omics, symptoms = symptoms)
}

normalize_asd_omics_time <- function(x) {
  z <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(z))
  out[z %in% c("1", "0", "T0", "BASELINE")] <- "T0"
  out[z %in% c("3", "10", "T1", "10WK", "WEEK10", "MTT_10WK")] <- "T1"
  out[z %in% c("2YR", "2YR_FOLLOWUP", "T2", "YEAR2", "2Y")] <- "T2"
  out
}

normalize_asd_symptom_time <- function(x) {
  z <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(z))
  out[z %in% c("T0", "BASELINE", "0")] <- "T0"
  out[z %in% c("T1", "10", "10WK", "WEEK10", "1")] <- "T1"
  out[z %in% c("T2", "2YR", "2YR_FOLLOWUP", "2")] <- "T2"
  out
}

make_asd_primary_id <- function(subject_id, timepoint) {
  paste(as.character(subject_id), as.character(timepoint), sep = "__")
}

asd_numeric_feature_matrix <- function(x, fill_zero = TRUE, scale_features = TRUE) {
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(data.frame(lapply(x, function(z) {
    suppressWarnings(as.numeric(as.character(z)))
  }), check.names = FALSE))
  colnames(mat) <- colnames(x)
  rownames(mat) <- rownames(x)
  storage.mode(mat) <- "numeric"
  if (fill_zero) mat[!is.finite(mat)] <- 0
  keep <- apply(mat, 2, function(z) stats::var(z, na.rm = TRUE) > 0)
  mat <- mat[, keep, drop = FALSE]
  if (scale_features && ncol(mat) > 0) {
    mat <- scale(mat)
    mat[!is.finite(mat)] <- 0
  }
  mat
}

build_asd_mae <- function(
  omics,
  symptoms,
  contrast = c("T1", "T2"),
  symptom = "ABC",
  views = names(omics),
  fill_zero = TRUE,
  scale_features = TRUE
) {
  contrast <- match.arg(contrast)
  keep_times <- c("T0", contrast)

  sym <- as.data.frame(symptoms, check.names = FALSE, stringsAsFactors = FALSE)
  sym$time_norm <- normalize_asd_symptom_time(sym$Timepoint)
  sym <- sym[sym$Symptom == symptom & sym$time_norm %in% keep_times, , drop = FALSE]
  sym$primary <- make_asd_primary_id(sym$ID, sym$time_norm)
  sym$X <- as.integer(sym$time_norm == contrast)
  sym$Y <- suppressWarnings(as.numeric(as.character(sym$Score)))
  sym <- sym[is.finite(sym$Y), c("primary", "ID", "time_norm", "X", "Y"), drop = FALSE]
  rownames(sym) <- sym$primary

  experiments <- list()
  sample_map_rows <- list()

  for (view in views) {
    pcl <- omics[[view]]
    meta <- as.data.frame(pcl$metadata, check.names = FALSE, stringsAsFactors = FALSE)
    feat <- as.data.frame(pcl$features, check.names = FALSE, stringsAsFactors = FALSE)

    meta$time_norm <- normalize_asd_omics_time(meta$Time_point)
    meta$primary <- make_asd_primary_id(meta[["Subject ID"]], meta$time_norm)
    rownames(feat) <- meta$primary

    idx <- meta$time_norm %in% keep_times & meta$primary %in% sym$primary
    if (!any(idx)) next

    mat <- asd_numeric_feature_matrix(
      feat[idx, , drop = FALSE],
      fill_zero = fill_zero,
      scale_features = scale_features
    )
    if (nrow(mat) < 4 || ncol(mat) < 2) next

    experiments[[view]] <- SummarizedExperiment::SummarizedExperiment(
      assays = list(abundance = t(mat))
    )
    sample_map_rows[[view]] <- data.frame(
      assay = view,
      primary = rownames(mat),
      colname = rownames(mat),
      stringsAsFactors = FALSE
    )
  }

  if (length(experiments) == 0) {
    stop("No usable omics views were available for contrast=", contrast, " and symptom=", symptom, call. = FALSE)
  }

  primary_ids <- Reduce(union, lapply(sample_map_rows, `[[`, "primary"))
  pheno <- sym[intersect(rownames(sym), primary_ids), , drop = FALSE]
  sample_map <- do.call(rbind, sample_map_rows)
  sample_map <- sample_map[sample_map$primary %in% rownames(pheno), , drop = FALSE]
  sample_map <- S4Vectors::DataFrame(
    assay = sample_map$assay,
    primary = sample_map$primary,
    colname = sample_map$colname
  )

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = experiments,
    colData = S4Vectors::DataFrame(
      X = pheno$X,
      Y = pheno$Y,
      subject_id = pheno$ID,
      timepoint = pheno$time_norm,
      row.names = rownames(pheno)
    ),
    sampleMap = sample_map
  )
}
