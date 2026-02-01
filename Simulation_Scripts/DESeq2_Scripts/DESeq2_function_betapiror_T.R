#### DESeq2 Pipeline for Null Simulations (DE=0) ####
#### Purpose: Flexible DESeq2 analysis for multiple datasets with null hypothesis

# Author: Erda Qorri
# Date: 01-02-2026

library(dplyr)
library(DESeq2)

#' Run DESeq2 pipeline on a single simulation replicate
#'
#' @param bsj_counts_path Path to counts CSV file
#' @param metadata_path Path to design/metadata CSV file
#' @param pvalue_thresholds Vector of pvalue thresholds to evaluate (default: c(0.01, 0.05, 0.10))
#' @param rep_id Replicate identifier (for tracking the simulated datasets in the results table)
#' 
#' @return Single-row dataframe with summary metrics

run_DESeq2_DE0 <- function(bsj_counts_path,
                           metadata_path,
                           pvalue_thresholds = c(0.01, 0.05, 0.10),
                           rep_id = NA) {
  # Start timer
  start_time = Sys.time()
  
  tryCatch({
    
    # Step 1: Load the data and the metadata
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    
    if("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # Print the initial gene number
    n_genes_initial = nrow(bsj_fc)
    
    # Step 2: Load the metadata
    bsj_fc_md = read.csv(metadata_path)
    
    # Make DESeq2 compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # Step 3: Create the DESeq2 object
    cat("Starting Step 3: Creating DESeqDataSet...\n")
    dds <- DESeqDataSetFromMatrix(countData = bsj_fc,
                                  colData = bsj_fc_md,
                                  design = ~Group)
    
    # Step 4: Relevel to compare cancer vs healthy
   # dds$Group <- relevel(dds$Group, ref = "Breast.Cancer.Tissue")
    
    # Step 5: Run DESeq2
    cat("Starting Step 5: Running DESeq2...\n")
    dds <- DESeq(dds,
                 betaPrior = TRUE)
    
    
    # Extract number of genes after DESeq2 filtering**
    # DESeq2 automatically filters out genes with all zeros or very low counts
    # Genes that pass filtering are those with results
    res_cancer_healthy <- results(dds, 
                                  contrast = c("Group", 
                                               "Breast.Cancer.Tissue",
                                               "Normal.Adjacent.Tissue"))
    
    res_cancer_healthy.df <- as.data.frame(res_cancer_healthy)
    print(res_cancer_healthy.df)
    
    # Remove genes with NA pvalues (these were filtered by DESeq2)
    n_genes_tested <- sum(!is.na(res_cancer_healthy.df$pvalue))
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 6: Calculate the metrics at each Pvalue threshold
    cat("Starting Step 6: Calculating metrics at each Pvalue threshold...\n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    
    # For each p-value threshold
    for (pval in pvalue_thresholds) {
      
      # Calculate the sum of rows that fulfill the pval condition
      n_sig = sum(res_cancer_healthy.df$pvalue < pval, na.rm = TRUE)
      
      # False positive rate
      fpr = n_sig / n_genes_tested
      
      # Create column names based on the threshold
      pval_names = gsub("\\.", "", sprintf("%.3f", pval))
      
      metrics[[paste0("n_sig_", pval_names)]] <- n_sig
      metrics[[paste0("fpr_", pval_names)]] <- round(fpr, 6)
    }
    
    # P-value distribution metrics (OUTSIDE the loop)
    metrics$mean_pval = round(mean(res_cancer_healthy.df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval = round(median(res_cancer_healthy.df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval = round(min(res_cancer_healthy.df$pvalue, na.rm = TRUE), 6)
    
    # Runtime (OUTSIDE the loop)
    processing_time <- Sys.time()
    metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
    
    cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))
    
    return(metrics)
    
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

# Run function
DESeq2_result <- run_DESeq2_DE0(
  bsj_counts_path  = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/BC_autofilter_counts_DE0_rep11_counts.csv",
  metadata_path = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/BC_autofilter_counts_DE0_rep11_design.csv"
)

print(DESeq2_result)
