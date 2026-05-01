#' Trigger InterSIM to Generate Null Multi-Omics Dataset
#'
#' This function triggers the `InterSIM` library to generate a null multi-omics dataset with specified parameters.
#' It processes the simulated data and returns it as a list containing feature data, sample metadata, and feature metadata.
#'
#' @param n Integer. The number of samples to generate.
#' @return List. A list containing:
#' \item{feature_table}{Data frame. Combined feature data (methylation, expression, protein) with formatted column names.}
#' \item{sample_metadata}{Data frame. Sample metadata including subject IDs and clustering assignments.}
#' \item{feature_metadata}{Data frame. Feature metadata including feature IDs and types.}
#' @export
#' @examples
#' simulated_data <- trigger_InterSIM(100)
#'
trigger_InterSIM <- function(n) {

  # Generate null multi-omics dataset with specified parameters
  sim.data <- InterSIM::InterSIM(n.sample = n,
                                 cluster.sample.prop = c(0.51, 0.49), # Sample proportion for clustering
                                 delta.methyl = 0, # Set delta to zero for null data simulation
                                 delta.expr = 0, # Set delta to zero for null data simulation
                                 delta.protein = 0) # Set delta to zero for null data simulation

  # Process features into a list
  names <- c('methyl', 'expr', 'protein')  # Define feature types
  list_features <- vector("list", length(names))  # Create an empty list for features
  names(list_features) <- names  # Name the list elements

  # Transpose and assign simulated data to the list
  list_features[[1]] <- t(sim.data$dat.methyl)
  list_features[[2]] <- t(sim.data$dat.expr)
  list_features[[3]] <- t(sim.data$dat.protein)

  # Combine the features into a single data frame
  feature_table <- as.data.frame(Reduce(rbind, list_features))
  colnames(feature_table) <- stringr::str_to_title(colnames(feature_table))  # Format column names

  # Extract and format sample metadata
  sample_metadata <- sim.data$clustering.assignment
  colnames(sample_metadata) <- c('subjectID', 'Y')
  sample_metadata$subjectID <- stringr::str_to_title(sample_metadata$subjectID)
  rownames(sample_metadata) <- sample_metadata$subjectID

  # Create feature metadata
  rowID <- rep(names, sapply(list_features, nrow))
  feature_metadata <- cbind.data.frame(featureID = rownames(feature_table), featureType = rowID)
  rownames(feature_metadata) <- feature_metadata$featureID

  # Save datasets as a list and return it
  pcl <- list(feature_table = feature_table,
              sample_metadata = sample_metadata,
              feature_metadata = feature_metadata)

  return(pcl)
}
