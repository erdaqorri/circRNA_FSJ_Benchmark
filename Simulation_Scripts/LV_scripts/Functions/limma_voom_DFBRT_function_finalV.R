#### limma-voom Pipeline for Null Simulations (DE=0) ####
#### Purpose: Flexible limma-voom analysis for multiple datasets with null hypothesis

# Author: Erda Qorri
# Date: 01-01-2026

##############################################################################################
                                         # limma-voom default with ROBUST=TRUE
##############################################################################################

library(dplyr)
library(edgeR)
library(limma)

#' Run limma-voom pipeline on a single simulation replicate
#'
#' @param bsj_counts_path Path to counts CSV file
#' @param metadata_path Path to design/metadata CSV file
#' @param contrast_formula Specify a string "Group1-Group2" for the contrast matrix
#' @param edgeR_norm Normalization method (default: "TMM", options: "TMM", "TMMwsp)
#' @param pvalue_thresholds Vector of pvalue thresholds to evaluate, the default is set to 3 but you can change it to as many as you want to test (default: c(0.01, 0.05, 0.10))
#' @param rep_id Replicate identifier (for tracking the simulated datasets in the results table)
#' 
#' @return Single-row dataframe with summary metrics
#' 
#' 

###### Automated function ###### 
run_limma_DFRBT_DE0 <- function(bsj_counts_path,
                                    metadata_path,
                                    contrast_formula = NULL,
                                    edgeR_norm = "TMM", # switch to TMMwsp
                                    pvalue_thresholds = c(0.01, 0.05, 0.10),
                                    rep_id = NA) {
  # Start timer
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata #
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    
    # Not an mtx -> converts it to one by removing the X col which is the sample id col
    rownames(bsj_fc) = bsj_fc$X
    bsj_fc$X = NULL
    
    # print the initial gene number (should be the same across all the simulated datasets)
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    
    # Step 3: Add the sample id information to the bsj mtx for tracking purposes
    bsj_fc_md$sample_id = colnames(bsj_fc)
    
    # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
    
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
    }
    
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Step 2: Create DGEList and filter #
    cat(" Starting Step 2: Creating the DGEList and filtering... \n")
    
    # first define groups for the dgelist
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    group <- factor(bsj_fc_md$Group)
    
    print(group)
    
    # Check that we have exactly 2 conditions
    if (length(levels(group)) != 2) {
      stop(paste("Expected 2 groups, found:", length(levels(group))))
    }
    
    cat(sprintf("Hey there! I detected the groups: %s\n", paste(levels(group), collapse = " vs ")))
    cat(sprintf("And the number of samples per group is: %s\n", paste(table(group), collapse = " and ")))
    
    # Create DGEList
    dgelist = DGEList(counts = bsj_fc, group = group)
    
    # Apply automated filtering
    keep = filterByExpr(dgelist, group = group)
    dgelist_filtered = dgelist[keep, , keep.lib.sizes = FALSE]
    
    # these are our truth set that are being tested in the DE for which we know they are not differentially expressed
    n_genes_tested = nrow(dgelist_filtered)
    
    # pct of genes that passed the automatic filtering
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
    
    # Assigns the name of the levels to the design mtx colnames
    colnames(design_mtx) <- levels(group)
    
    # Sanity check: the user provides the contrast formula so there should not be any issues here anyways
    
    if (is.null(contrast_formula)) {
      stop("Provide the contrast formula as a string like 'Breast.Cancer.Tissue-Normal.Adjacent.Tissue' to control the direction of the comparison!")
    }
    
    # the provided colnames must reference the existing design columns
    cn <- colnames(design_mtx)
    
    parts <-  strsplit(gsub("\\s+", "", contrast_formula), "-", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("`contrast` must look like 'Group1-Group2' (exactly one '-').")
    
    if (!all(parts %in% cn)) {
      stop(sprintf(
        "Contrast groups not found in design. Provided: %s. Available: %s",
        paste(parts, collapse = ", "),
        paste(cn, collapse = ", ")
      ))
    }
    
    contrasts <- makeContrasts(contrasts = contrast_formula, levels = design_mtx)
    print(contrasts)
    
    #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
    # it replaces all the commands with one function
    
    # Fit the default limma pipeline but using quantile normalization
    dgelist_v <- voom(dgelist_norm, design_mtx)
    
    # Fit the original linear model with lmFit just set robust to TRUE
    dgelist_vfit <- lmFit(dgelist_v, design_mtx, robust = TRUE)
    
    # Step 5: Set the contrast matrix and DE analysis
    cat(" Starting Step 5: Testing contrasts... \n")
    
    # Get the group names from design matrix columns
    cat(sprintf(" Testing contrast: %s\n", contrast_formula))
    
    dgelist_vfit <- contrasts.fit(dgelist_vfit, contrasts = contrasts)
    
    dgelist_efit <- eBayes(dgelist_vfit)
    
    # Extract the full results
    summary(decideTests(dgelist_efit))
    de_results <- topTable(dgelist_efit, number = Inf, sort.by="none")
    
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
    
    return(list(metrics = metrics, de_results = de_results))
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
lima_DFRBT_results <- run_limma_DFRBT_DE0(
  bsj_counts_path  = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/BC_autofilter_counts_DE0_rep30_counts.csv",
  metadata_path = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/BC_autofilter_counts_DE0_rep30_design.csv",
  edgeR_norm  = "TMM",
  contrast_formula = "Breast.Cancer.Tissue-Normal.Adjacent.Tissue"
)

print(lima_DFRBT_results$metrics)

