#### limma-voom Pipeline for Null Simulations (DE=0) ####
#### Purpose: Flexible limma-voom analysis for multiple datasets with null hypothesis
#### To run the lmFit in the robust mode please uncomment line 106

library(dplyr)
library(edgeR)
library(limma)

#' Run limma-voom pipeline on a single simulation replicate
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
run_limma_default_DE0 <- function(bsj_counts_path,
                                    metadata_path,
                                    edgeR_norm = "TMM", # switch to TMMwsp and quantile
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
    
    # print the initial gene number (should be the same across all the simulated datasets)
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    bsj_fc_md$sample_id = colnames(bsj_fc)
    
    # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
    stopifnot("sample_id" %in% colnames(bsj_fc_md))
    
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), ]
    
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in the BSJ counts matrix and the metadata file do not match 1:1.")
    }
    
    # Step 2: Create DGEList and filter #
    cat(" Starting Step 2: Creating the DGEList and filtering... \n")
    
    # first define groups for the dgelist
    fc_sim1_md$Group <- make.names(as.character(fc_sim1_md$Group))
    group <- factor(fc_sim1_md$Group)
    group = as.factor(fc_sim1_md$Group)
    
    print(group)
    
    # Check that we have exactly 2 conditions
    if (length(levels(group)) != 2) {
      stop(paste("Expected 2 groups, found:", length(levels(group))))
    }
    
    # Create DGEList
    dgelist = DGEList(counts = bsj_fc, group = group)
    
    # Apply automated filtering
    keep = filterByExpr(dgelist, group = group)
    dgelist_filtered = dgelist[keep, , keep.lib.sizes = FALSE]
    
    # these are our truth set that are being tested in the DE for which we know they are not differentially expressed
    n_genes_tested = nrow(dgelist_filtered)
    
    pct_filtered = round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 3: Normalization #
    cat(sprintf(" Starting Step 3: Normalization: Normalizing with %s... \n", edgeR_norm))
    
    dgelist_norm = calcNormFactors(dgelist_filtered, method = edgeR_norm)
    
    # Step 4: Design matrix and model fitting
    cat("Starting Step 4: Fitting voomLmFit model... \n")
    
    # Design matrix (no intercept was included since it is a simple two group comparison)
    design_mtx = model.matrix(~ 0 + group)
    colnames(design_mtx) <- levels(group)
    
    # Create the design matrix
    colnames_cn <- colnames(design_mtx)
    contrasts <- makeContrasts(contrasts = paste0(colnames_cn[1], "-", colnames_cn[2]), levels = design_mtx)
    print(contrasts)
    
    #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
    # it replaces all the commands with one function
    dgelist_v <- voom(dgelist_norm, design_mtx, plot = F)
    dgelist_vfit <- lmFit(dgelist_v, design_mtx)
    # dgelist_vfit <- lmFit(dgelist_v, design_mtx, robust = TRUE)

    # Step 5: Set the contrast matrix and DE analysis
    cat(" Starting Step 5: Testing contrasts... \n")
    
    # Get the group names from design matrix columns
    cat(sprintf(" Comparing: %s vs %s\n", colnames_cn[1], colnames_cn[2]))
    
    dgelist_vfit <- contrasts.fit(dgelist_vfit, contrasts = contrasts)
    
    dgelist_efit <- eBayes(dgelist_vfit)
    
    # Extract the full results
    summary(decideTests(dgelist_efit))
    de_results <- topTable(dgelist_efit, number=Inf, adjust.method="none", sort.by="none")
    
    # Step 6: Calculate the metrics at each Pvalue threshold #
    cat(" Starting Step 6: Calculating metrics at each Pvalue threshold... \n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    
    # For each FDR threshold
    for (pval in pvalue_thresholds) {
      n_sig = sum(de_results$P.Value < pval, na.rm = TRUE)
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
    metrics$mean_pval = round(mean(de_results$P.Value, na.rm = TRUE), 4)
    metrics$median_pval = round(median(de_results$P.Value, na.rm = TRUE), 4)
    metrics$min_pval = round(min(de_results$P.Value, na.rm = TRUE), 6)
    
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
test_result <- run_limma_default_DE0(
  bsj_counts_path  = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/BC_autofilter_counts_DE0_rep1_counts.csv",
  metadata_path = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/BC_autofilter_counts_DE0_rep1_design.csv",
  edgeR_norm  = "TMM"
)

print(test_result_robust)

