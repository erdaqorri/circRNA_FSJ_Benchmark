#### edgeR Pipeline for Null Simulations (DE=0) ####
#### Purpose: Flexible edgeR analysis for multiple datasets with null hypothesis

library(dplyr)
library(edgeR)

#' Run edgeR pipeline on a single simulation replicate
#'
#' @param bsj_counts_path Path to counts CSV file
#' @param metadata_path Path to design/metadata CSV file
#' @param edgeR_norm Normalization method (default: "TMM", options: "TMM", "TMMwsp)
#' @param pvalue_thresholds Vector of pvalue thresholds to evaluate, the default is set to 3 but you can change it to as many as you want to test (default: c(0.01, 0.05, 0.10))
#' @param rep_id Replicate identifier (for tracking the simulated datasets in the results table)
#' 
#' @return Single-row dataframe with summary metrics
#' 
#' 

###### Automated function ###### 
run_edgeR_DE0 <- function(bsj_counts_path,
                          metadata_path,
                          edgeR_norm = "TMM",
                          pvalue_thresholds = c(0.01, 0.05, 0.10),
                          rep_id = NA) {
  # Start timer
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata #
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    rownames(bsj_fc) = bsj_fc$X
    bsj_fc$X = NULL
    #head(bsj_fc)

    # print the initial gene number (should be the same across all the simulated datasets)
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    bsj_fc_md$sample_id = colnames(bsj_fc)
    #head(bsj_fc_md)
    
    # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
    if(!identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      # Try to reorder the metadata to match count columns
      bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), ]
      
      # Verify again
      if(!identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
        stop("Sample IDs in the BSJ counts matrix and the metadata files do not match!")
      }
    
    }
    
    # Step 2: Create DGEList and filter #
    cat(" Starting Step 2: Creating the DGEList and filtering... \n")
    
    # first define groups for the dgelist
    group = as.factor(bsj_fc_md$Group)
    #print(group)
    
    # Check that we have exactly 2 conditions
    if (length(levels(group)) != 2) {
      stop(paste("Expected 2 groups, found:", length(levels(group))))
    }
    
    # Create DGEList
    dgelist = DGEList(counts = bsj_fc, group = group)
    #head(dgelist)
    
    # Apply automated filtering
    keep = filterByExpr(dgelist, group = group)
    dgelist_filtered = dgelist[keep, , keep.lib.sizes = FALSE]
    #head(dgelist_filtered)
    
    # these are our truth set that are being tested in the DE for which we know they are not differentially expressed
    n_genes_tested = nrow(dgelist_filtered)
    #print(n_genes_tested)
    
    pct_filtered = round((n_genes_tested / n_genes_initial) * 100, 2)
    
    #print(n_genes_initial)
    #print(pct_filtered)
    
  cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
              n_genes_tested, n_genes_initial, pct_filtered))
  
  # Step 3: Normalization #
  cat(sprintf(" Starting Step 3: Normalization: Normalizing with %s... \n", edgeR_norm))
  
  dgelist_norm = calcNormFactors(dgelist_filtered, method = edgeR_norm)
  #head(dgelist_norm)
  
  # Step 4: Design matrix and model fitting
  cat("Starting Step 4: Fitting glmQLFit model... \n")
  
  # Create valid R names for the design matrix groups
  group = factor(group)
  levels(group) = make.names(levels(group))
  #print(group)
  
  # Design matrix (no intercept was included since it is a simple two group comparison)
  design_mtx = model.matrix(~ 0 + group)
  colnames(design_mtx) = levels(group)
  #print(design_mtx)
  
  # Fit the quasi-likelihood model
  fit = glmQLFit(dgelist_norm, design = design_mtx, robust = TRUE)
  #dim(dgelist_norm) %>% print()
  
  # Step 5: Set the contrast matrix and DE analysis
  cat(" Starting Step 5: Testing contrasts... \n")
  
  # Create contrast: cancer vs healthy controls
  contrast_formula = paste(levels(group)[1], "-", levels(group)[2])
  #print(contrast_formula)
  contrasts = makeContrasts(contrasts = contrast_formula, levels = design_mtx)
  #print(contrasts)
  
  # Quasi-likelihood F-test
  qlf = glmQLFTest(fit, contrast = contrasts)
  #print(qlf)
  
  # Extract the full results
  results = as.data.frame(topTags(qlf, n = Inf))
  
  # Step 6: Calculate the metrics at each Pvalue threshold #
  cat(" Starting Step 6: Calculating metrics at each Pvalue threshold... \n")
  
  metrics = data.frame(rep_id = rep_id)
  metrics$n_genes_initial = n_genes_initial
  metrics$n_genes_tested = n_genes_tested
  metrics$pct_filtered = pct_filtered
  
  # For each FDR threshold
  for (pval in pvalue_thresholds) {
    n_sig = sum(results$PValue < pval, na.rm = TRUE)
    #n_sig = sum(results$PValue < pval, na.rm = TRUE)
    # false positive rate: number of genes that were found to be significant over the number of genes following the filtering step (genes tested for DE) 
    fpr = n_sig / n_genes_tested

    # Create column names based on the threshold
    pval_names = gsub("\\.", "", sprintf("%.3f", pval))
    print(pval_names)
    
    metrics[[paste0("n_sig_", pval_names)]] <- n_sig
    metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
  }
  
  # P-value distribution metrics (OUTSIDE the loop)
  metrics$mean_pval = round(mean(results$PValue, na.rm = TRUE), 4)
  metrics$median_pval = round(median(results$PValue, na.rm = TRUE), 4)
  metrics$min_pval = round(min(results$PValue, na.rm = TRUE), 6)
  
  # Runtime (OUTSIDE the loop)
  processing_time <- Sys.time()
  metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
  
  cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))
  
  return(metrics)  # OUTSIDE the loop
  }, error = function(e) {
    # If error occurs, return a row with NA values and error message
    cat(sprintf("  ERROR: %s\n", e$message))
    
    error_row <- data.frame(
      rep_id = rep_id,
      n_genes_initial = NA,
      n_genes_tested = NA,
      pct_filtered = NA
    )
    
    # Add NA columns for pval thresholds
    for (pval in pvalue_thresholds) {
      pval_label <- gsub("\\.", "", sprintf("%.3f", pval))
      error_row[[paste0("n_sig_", pval_label)]] <- NA
      error_row[[paste0("fpr_", pval_label)]] <- NA
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(error_row)
  })
}

# run function
test_result <- run_edgeR_DE0(
  bsj_counts_path  = "/Users/alexakouri/Downloads/Test_data4Alexa/BSJ_only/DE0_test-set/HCC-tissue_autofilter_counts_DE0_rep5_counts.csv",
  metadata_path = "/Users/alexakouri/Downloads/Test_data4Alexa/BSJ_only/DE0_test-set/HCC-tissue_autofilter_counts_DE0_rep5_design.csv",
  edgeR_norm  = "TMM"
)

print(test_result)

