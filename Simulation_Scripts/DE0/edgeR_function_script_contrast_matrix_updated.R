#### edgeR Pipeline for Null Simulations (DE=0) ####
#### Purpose: Flexible edgeRanalysis for multiple datasets with null hypothesis

library(dplyr)
library(edgeR)
library(limma)

#' Run edgeRpipeline on a single simulation replicate
#'
#' @param bsj_counts_path Path to counts CSV file
#' @param metadata_path Path to design/metadata CSV file
#' @param tp_path Path to the simulated baseline (TP) CSV file with columns `X` and `DE.ind`
#' @param contrast_formula Specify a string "Group1-Group2" for the contrast matrix
#' @param edgeR_norm Normalization method (default: "TMM", options: "TMM", "TMMwsp)
#' @param pvalue_thresholds Vector of thresholds to evaluate (default: c(0.01, 0.05, 0.10))
#' @param threshold_column Which statistics column to threshold on (`PValue` or `FDR`, default: `PValue`)
#' @param rep_id Replicate identifier (for tracking the simulated datasets in the results table)
#' 
#' @return List with `metrics` (single-row dataframe) and `de_results` (full DE table)
#' 
#' 

###### Automated function ###### 
run_edgeR_DE0 <- function(bsj_counts_path,
                          metadata_path,
                          tp_path = NULL,
                          contrast_formula = NULL,
                          edgeR_norm = "TMM", # switch to TMMwsp
                          pvalue_thresholds = c(0.01, 0.05, 0.10),
                          threshold_column = c("PValue", "FDR"),
                          rep_id = NA) {
  # Start timer
  start_time = Sys.time()
  
  tryCatch({
    threshold_column <- match.arg(threshold_column)
    
    # Step 1: Load the data and the metadata #
    cat("Starting Step 1: Loading the BSJ matrix ...\n")
    
    bsj_fc = read.csv(bsj_counts_path)
    
    # Not an mtx -> converts it to one by removing the X col which is the sample id col
    rownames(bsj_fc) = bsj_fc$X
    bsj_fc$X = NULL
    
    # print the initial gene number (should be the same across all the simulated datasets)
    n_genes_initial_counts = nrow(bsj_fc)
    
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
    tested_ids <- rownames(dgelist_filtered$counts)
    n_genes_tested = length(tested_ids)
    
    # Optional: Load ground-truth baseline (TP) file
    if (!is.null(tp_path)) {
      simulated_baseline <- read.csv(tp_path)
      required_cols <- c("X", "DE.ind")
      missing_cols <- setdiff(required_cols, colnames(simulated_baseline))
      if (length(missing_cols) > 0) {
        stop(sprintf("TP file is missing required columns: %s", paste(missing_cols, collapse = ", ")))
      }
      
      if (anyDuplicated(simulated_baseline$X)) {
        stop("TP file has duplicated circRNA IDs in column `X`.")
      }
      
      # Coerce DE.ind to logical if needed
      if (!is.logical(simulated_baseline$DE.ind)) {
        if (is.numeric(simulated_baseline$DE.ind)) {
          simulated_baseline$DE.ind <- simulated_baseline$DE.ind != 0
        } else if (is.character(simulated_baseline$DE.ind)) {
          simulated_baseline$DE.ind <- tolower(simulated_baseline$DE.ind) %in% c("true", "t", "1", "yes", "y")
        } else {
          simulated_baseline$DE.ind <- as.logical(simulated_baseline$DE.ind)
        }
      }
      
      n_genes_initial_truth <- nrow(simulated_baseline)
      
      if (!all(tested_ids %in% simulated_baseline$X)) {
        missing_truth <- setdiff(tested_ids, simulated_baseline$X)
        stop(sprintf(
          "TP file is missing %d tested circRNAs (e.g. %s).",
          length(missing_truth),
          paste(head(missing_truth, 5), collapse = ", ")
        ))
      }
      
      ground_truth_set_tested <- simulated_baseline %>% dplyr::filter(X %in% tested_ids)
      DE_truth <- ground_truth_set_tested %>% dplyr::filter(DE.ind == TRUE)
      non_DE_truth <- ground_truth_set_tested %>% dplyr::filter(DE.ind == FALSE)
      
      n_genes_initial <- n_genes_initial_truth
    } else {
      n_genes_initial <- n_genes_initial_counts
      ground_truth_set_tested <- NULL
      DE_truth <- NULL
      non_DE_truth <- NULL
    }
    
    # pct of genes that passed the automatic filtering
    pct_filtered = round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf(" Retained %d/%d genes (%.1f%%)\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Step 3: Normalization #
    cat(sprintf(" Starting Step 3: Normalization: Normalizing with %s... \n", edgeR_norm))
    
    dgelist_norm = calcNormFactors(dgelist_filtered, method = edgeR_norm)
    
    # Step 4: Design matrix and model fitting
    cat("Starting Step 4: Fitting QGLM model... \n")
    
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
    
    #### Step 5: Fit glmQLFit with glmQLFTest
    
    # Fit the quasi-likelihood model
    fit = glmQLFit(dgelist_norm, design = design_mtx, robust = TRUE)
    
    # Quasi-likelihood F-test
    qlf = glmQLFTest(fit, contrast = contrasts)
    #print(qlf)
    
    # Extract the full results
    de_results = as.data.frame(topTags(qlf, n = Inf))
    de_results$X <- rownames(de_results)
    rownames(de_results) <- NULL
    
    # Step 6: Calculate the metrics at each Pvalue threshold #
    cat(" Starting Step 6: Calculating metrics at each threshold... \n")
    
    metrics = data.frame(rep_id = rep_id)
    metrics$n_genes_initial = n_genes_initial
    metrics$n_genes_tested = n_genes_tested
    metrics$pct_filtered = pct_filtered
    
    safe_divide <- function(num, den) ifelse(den == 0, NA_real_, num / den)
    
    # For each threshold
    for (thr in pvalue_thresholds) {
      stat_vals <- de_results[[threshold_column]]
      n_sig <- sum(!is.na(stat_vals) & stat_vals <= thr)
      
      # Create column names based on the threshold
      thr_name = gsub("\\.", "", sprintf("%.3f", thr))
      
      metrics[[paste0("n_sig_", thr_name)]] <- n_sig
      
      # If we have a TP file, compute confusion matrix + metrics
      if (!is.null(tp_path)) {
        pred_pos_ids <- de_results$X[!is.na(stat_vals) & stat_vals <= thr]
        pred_neg_ids <- setdiff(tested_ids, pred_pos_ids)
        
        TP <- nrow(dplyr::inner_join(data.frame(X = pred_pos_ids), DE_truth, by = "X"))
        FP <- nrow(dplyr::inner_join(data.frame(X = pred_pos_ids), non_DE_truth, by = "X"))
        FN <- nrow(dplyr::inner_join(data.frame(X = pred_neg_ids), DE_truth, by = "X"))
        TN <- nrow(dplyr::inner_join(data.frame(X = pred_neg_ids), non_DE_truth, by = "X"))
        
        sen <- safe_divide(TP, TP + FN)
        spe <- safe_divide(TN, TN + FP)
        precision <- safe_divide(TP, TP + FP)
        recall <- sen
        f1 <- safe_divide(2 * precision * recall, precision + recall)
        fdr <- safe_divide(FP, FP + TP)
        fpr <- safe_divide(FP, FP + TN)
        
        metrics[[paste0("TP_", thr_name)]] <- TP
        metrics[[paste0("FP_", thr_name)]] <- FP
        metrics[[paste0("FN_", thr_name)]] <- FN
        metrics[[paste0("TN_", thr_name)]] <- TN
        
        metrics[[paste0("sen_", thr_name)]] <- round(sen, 6)
        metrics[[paste0("spe_", thr_name)]] <- round(spe, 6)
        metrics[[paste0("precision_", thr_name)]] <- round(precision, 6)
        metrics[[paste0("recall_", thr_name)]] <- round(recall, 6)
        metrics[[paste0("f1_", thr_name)]] <- round(f1, 6)
        metrics[[paste0("fdr_", thr_name)]] <- round(fdr, 6)
        metrics[[paste0("fpr_", thr_name)]] <- round(fpr, 6)
      } else {
        # Null-only fallback: approximate FPR as significant/tested
        fpr <- safe_divide(n_sig, n_genes_tested)
        metrics[[paste0("fpr_", thr_name)]] <- round(fpr, 6)
      }
    }
    
    # P-value distribution metrics (OUTSIDE the loop)
    metrics$mean_pval = round(mean(de_results$PValue, na.rm = TRUE), 4)
    metrics$median_pval = round(median(de_results$PValue, na.rm = TRUE), 4)
    metrics$min_pval = round(min(de_results$PValue, na.rm = TRUE), 6)
    
    # Runtime (OUTSIDE the loop)
    processing_time <- Sys.time()
    metrics$runtime_sec = round(as.numeric(difftime(processing_time, start_time, units = "secs")), 2)
    
    cat(sprintf(" Completed in %.2f seconds\n", metrics$runtime_sec))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    
    return(list(metrics = metrics, de_results = de_results, ground_truth_set_tested = ground_truth_set_tested))
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
    for (thr in pvalue_thresholds) {
      thr_label <- gsub("\\.", "", sprintf("%.3f", thr))
      error_row[[paste0("n_sig_", thr_label)]] <- NA
      
      if (!is.null(tp_path)) {
        error_row[[paste0("TP_", thr_label)]] <- NA
        error_row[[paste0("FP_", thr_label)]] <- NA
        error_row[[paste0("FN_", thr_label)]] <- NA
        error_row[[paste0("TN_", thr_label)]] <- NA
        error_row[[paste0("sen_", thr_label)]] <- NA
        error_row[[paste0("spe_", thr_label)]] <- NA
        error_row[[paste0("precision_", thr_label)]] <- NA
        error_row[[paste0("recall_", thr_label)]] <- NA
        error_row[[paste0("f1_", thr_label)]] <- NA
        error_row[[paste0("fdr_", thr_label)]] <- NA
        error_row[[paste0("fpr_", thr_label)]] <- NA
      } else {
        error_row[[paste0("fpr_", thr_label)]] <- NA
      }
    }
    
    error_row$mean_pval <- NA
    error_row$median_pval <- NA
    error_row$min_pval <- NA
    error_row$runtime_sec <- NA
    error_row$error_message <- e$message
    
    return(list(metrics = error_row, de_results = NULL, ground_truth_set_tested = NULL))
  })
}

# run function
edgeR_default_TMM_DE0 <- run_edgeR_DE0(
  bsj_counts_path  = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/BC_autofilter_counts_DE0_rep5_counts.csv",
  metadata_path = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0/BC_autofilter_counts_DE0_rep5_design.csv",
  edgeR_norm  = "TMM",
  contrast_formula = "Breast.Cancer.Tissue-Normal.Adjacent.Tissue",
  threshold_column = "PValue"
)

print(edgeR_default_TMM_DE0$metrics)
