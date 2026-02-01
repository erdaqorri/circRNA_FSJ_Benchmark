#' Analyze a batch of null simulation replicates
#'
#' @param data_dir Directory containing all simulation files
#' @param dataset_name Name of dataset (e.g., "BC", "LC") for output naming
#' @param file_pattern Pattern to match count files (e.g., "BC_autofilter_counts_DE0_rep")
#' @param pvalue_thresholds pvalue thresholds to evaluate (default: c(0.01, 0.05, 0.10))
#' @param output_dir Directory to save results (default: same as data_dir)
#'
#' @return Dataframe with all replicates' results
#'
analyze_null_simulation_batch <- function(data_dir,
                                          dataset_name,
                                          file_pattern,
                                          pvalue_thresholds = c(0.01, 0.05, 0.10),
                                          output_dir = NULL) {
  
  cat("\n========================================\n")
  cat(sprintf("Starting batch analysis for dataset: %s\n", dataset_name))
  cat("========================================\n\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- data_dir
  }
  
  # Ensure data_dir ends with a slash
  if (!endsWith(data_dir, "/")) {
    data_dir <- paste0(data_dir, "/")
  }
  
  #### Step 1: Discover all count files ####
  cat("Discovering simulation files...\n")
  
  # Find all count files matching the pattern
  all_files <- list.files(data_dir, pattern = paste0(file_pattern, ".*_counts\\.csv$"), full.names = FALSE)
  
  if (length(all_files) == 0) {
    stop(sprintf("No files found matching pattern: %s*_counts.csv", file_pattern))
  }
  
  cat(sprintf("Found %d count files\n\n", length(all_files)))
  
  #### Step 2: Extract rep numbers and pair with design files ####
  # Extract rep numbers from filenames
  rep_numbers <- gsub(paste0(".*", file_pattern, "(\\d+)_counts\\.csv"), "\\1", all_files)
  rep_numbers <- as.numeric(rep_numbers)
  
  # Sort by rep number
  file_order <- order(rep_numbers)
  all_files <- all_files[file_order]
  rep_numbers <- rep_numbers[file_order]
  
  # Create corresponding design file names
  design_files <- gsub("_counts\\.csv", "_design.csv", all_files)
  
  # Verify design files exist
  design_exist <- file.exists(paste0(data_dir, design_files))
  if (!all(design_exist)) {
    missing <- design_files[!design_exist]
    stop(sprintf("Missing design files: %s", paste(missing, collapse = ", ")))
  }
  
  #### Step 3: Initialize results storage ####
  results_list <- vector("list", length = length(all_files))
  
  #### Step 4: Process each replicate ####
  n_total <- length(all_files)
  n_success <- 0
  n_failed <- 0
  
  for (i in seq_along(all_files)) {
    cat(sprintf("Processing replicate %d/%d (rep_%d)...\n", i, n_total, rep_numbers[i]))
    
    # Build full paths
    counts_path <- paste0(data_dir, all_files[i])
    design_path <- paste0(data_dir, design_files[i])
    
    # Run analysis
    result <- run_DESeq2_DE0(
      bsj_counts_path = counts_path,
      metadata_path = design_path,
      pvalue_thresholds = pvalue_thresholds,
      rep_id = rep_numbers[i]
    )
    
    # Store result
    results_list[[i]] <- result
    
    # Track success/failure
    if (!is.null(result) && !is.na(result$n_genes_tested[1])) {
      n_success <- n_success + 1
    } else {
      n_failed <- n_failed + 1
    }
  }
  
  #### Step 5: Combine results ####
  cat("\nCombining results...\n")
  results_df <- do.call(rbind, results_list)
  
  #### Step 6: Add dataset metadata ####
  results_df <- results_df %>%
    mutate(
      dataset = dataset_name,
      .before = 1,
      norm_method = "Median of Ratios"
    )
  
  #### Step 7: Save results ####
  output_filename <- sprintf("%s_DE0_DESeq2_summary.csv", dataset_name)
  output_path <- file.path(output_dir, output_filename)
  
  write.csv(results_df, output_path, row.names = FALSE)
  
  #### Step 8: Print summary ####
  cat("\n========================================\n")
  cat("BATCH ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(sprintf("Dataset: %s\n", dataset_name))
  cat(sprintf("Total replicates: %d\n", n_total))
  cat(sprintf("Successful: %d\n", n_success))
  cat(sprintf("Failed: %d\n", n_failed))
  cat(sprintf("\nResults saved to:\n%s\n", output_path))
  
  # Print some summary statistics if we have successful runs
  if (n_success > 0) {
    valid_results <- results_df %>% filter(!is.na(n_genes_tested))
    
    cat("\n--- Summary Statistics ---\n")
    cat(sprintf("Mean genes after filtering: %.0f (±%.0f)\n", 
                mean(valid_results$n_genes_tested, na.rm = TRUE),
                sd(valid_results$n_genes_tested, na.rm = TRUE)))
    
    for (pvalue in pvalue_thresholds) {
      pvalue_label <- gsub("\\.", "", sprintf("%.3f", pvalue))
      col_name <- paste0("n_sig_", pvalue_label)
      
      if (col_name %in% colnames(valid_results)) {
        cat(sprintf("Mean significant genes (pvalue<%.2f): %.1f (±%.1f)\n",
                    pvalue,
                    mean(valid_results[[col_name]], na.rm = TRUE),
                    sd(valid_results[[col_name]], na.rm = TRUE)))
      }
    }
    
    cat(sprintf("Mean p-value: %.3f (±%.3f) [Expected ~0.5]\n",
                mean(valid_results$mean_pval, na.rm = TRUE),
                sd(valid_results$mean_pval, na.rm = TRUE)))
    cat(sprintf("Mean runtime: %.2f seconds\n", 
                mean(valid_results$runtime_sec, na.rm = TRUE)))
  }
  
  cat("========================================\n\n")
  
  return(results_df)
}


# ==============================================================================
# USAGE
# ==============================================================================

# Run batch analysis on BC dataset
bc_results <- analyze_null_simulation_batch(
  data_dir = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/",
  dataset_name = "BC",
  file_pattern = "BC_autofilter_counts_DE0_rep",
  output_dir = "/Users/alexakouri/Documents/DE_benchmarking_paper/"
)

head(bc_results, 11)
