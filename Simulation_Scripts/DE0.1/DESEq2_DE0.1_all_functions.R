################### DESeq2 functions for the three evaluated pipelines ################### 
###### Author: Erda Alexandria Qorri ###### 
###### Date: 06-02-2026

set.seed(42)

##### 1. DESeq2 with Wald Test (default pipeline) #####
library(dplyr)
library(DESeq2)
library(tibble)

run_DESeq2_WaldTest_with_truth <- function(bsj_counts_path,
                                           metadata_path,
                                           truth_set_path,
                                           pvalue_thresholds = c(0.01, 0.05, 0.10),
                                           threshold_col = c("padj"),
                                           rep_id = NA,
                                           contrast = c("Group", "Breast.Cancer.Tissue", "Normal.Adjacent.Tissue"),
                                           quiet = TRUE) {
  
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  tryCatch({
    
    # Step 1: Load counts
    cat("Step 1: Loading counts...\n")
    bsj_fc <- read.csv(bsj_counts_path)
    
    if ("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # DESeq2 wants integer counts
    bsj_fc <- as.matrix(bsj_fc)
    storage.mode(bsj_fc) <- "integer"
    
    n_genes_initial <- nrow(bsj_fc)
    
    # Step 2: Load metadata
    cat("Step 2: Loading metadata...\n")
    bsj_fc_md <- read.csv(metadata_path)
    
    # Ensure sample_id exists (many of your design files may not include it)
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      bsj_fc_md$sample_id <- colnames(bsj_fc)
    }
    
    # Reorder metadata to match counts
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Make Group DESeq2-compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # Step 3: Load truth set
    cat("Step 3: Loading truth set...\n")
    truth_set <- read.csv(truth_set_path)
    if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
      stop("Truth file must contain columns: X, DE.ind")
    }
    
    truth_set$DE.ind <- as.logical(truth_set$DE.ind)
    
    # Important: removes duplicates
    truth_set <- truth_set %>% distinct(X, .keep_all = TRUE)
    
    n_truth_initial <- length(unique(truth_set$X))
    
    # Step 4: DESeq2 object
    cat("Step 4: Creating DESeqDataSet...\n")
    
    dds <- DESeqDataSetFromMatrix(
      countData = bsj_fc,
      colData = bsj_fc_md,
      design = ~ Group
    )
    
    
    # Step 5: Run DESeq2
    cat("Step 5: Running DESeq2...\n")
    
    dds <- DESeq(dds)
    
    # Step 6: Results
    cat("Step 6: Extracting results...\n")
    res <- results(dds,
                   contrast = c("Group",
                                "Breast.Cancer.Tissue",
                                "Normal.Adjacent.Tissue"))
    
    res_df <- as.data.frame(res)
    
    # Tested universe: non-NA in the chosen threshold column
    res_df <- as.data.frame(res) %>% rownames_to_column("X")
    tested_df <- res_df %>% filter(!is.na(.data[[threshold_col]]))
    
    # these are the tested circRNAs for DE
    tested_ids <- tested_df$X
    
    # Genes that were tested in the DE analysis by DESeq2
    n_genes_tested <- length(tested_ids)
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf("Retained %d/%d genes (%.1f%%) in tested set\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Restrict truth to tested universe
    truth_tested <- truth_set %>% filter(X %in% tested_ids)
    
    DE_truth     <- truth_tested %>% filter(DE.ind == TRUE)
    nonDE_truth  <- truth_tested %>% filter(DE.ind == FALSE)
    
    # Helpers
    safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
    
    # Store metrics per threshold (tidy format: one row per threshold)
    metrics <- lapply(pvalue_thresholds, function(th) {
      
      called <- tested_df %>% filter(.data[[threshold_col]] < th)
      not_called <- tested_df %>% filter(!(X %in% called$X))
      
      TP <- nrow(inner_join(called,     DE_truth,    by = "X"))
      FP <- nrow(inner_join(called,     nonDE_truth, by = "X"))
      TN <- nrow(inner_join(not_called, nonDE_truth, by = "X"))
      FN <- nrow(inner_join(not_called, DE_truth,    by = "X"))
      
      precision <- safe_div(TP, TP + FP)
      recall    <- safe_div(TP, TP + FN)
      specificity <- safe_div(TN, TN + FP)
      fpr <- safe_div(FP, FP + TN)
      fdr <- safe_div(FP, TP + FP)
      f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                   NA_real_,
                   2 * precision * recall / (precision + recall))
      
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        
        n_genes_initial = n_genes_initial,
        n_truth_initial = n_truth_initial,
        n_genes_tested = n_genes_tested,
        pct_filtered = pct_filtered,
        
        TP = TP, FP = FP, TN = TN, FN = FN,
        
        sensitivity = recall * 100,
        specificity = specificity * 100,
        precision = precision,
        recall = recall,
        F1 = f1,
        FDR = fdr * 100,
        FPR = fpr * 100,
        
        n_called = nrow(called)
      )
    }) %>% bind_rows()
    
    metrics$mean_pval <- round(mean(res_df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval <- round(median(res_df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval <- round(min(res_df$pvalue, na.rm = TRUE), 6)
    
    metrics$runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
    
    return(list(metrics = metrics, de_results = res_df))
    
  }, error = function(e) {
    
    cat(sprintf("ERROR: %s\n", e$message))
    
    # Return a tidy error table with one row per threshold (helps batch binding)
    error_rows <- lapply(pvalue_thresholds, function(th) {
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        n_genes_initial = NA,
        n_truth_initial = NA,
        n_genes_tested = NA,
        pct_filtered = NA,
        TP = NA, FP = NA, TN = NA, FN = NA,
        sensitivity = NA, specificity = NA,
        precision = NA, recall = NA, F1 = NA,
        FDR = NA, FPR = NA,
        n_called = NA,
        mean_pval = NA, median_pval = NA, min_pval = NA,
        runtime_sec = NA,
        error_message = e$message
      )
    }) %>% bind_rows()
    
    return(list(metrics = error_rows))
  })
}

set.seed(42)

res <- run_DESeq2_WaldTest_with_truth(
  bsj_counts_path  = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep2_counts.csv",
  metadata_path    = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep2_design.csv",
  truth_set_path   = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep2_TP.csv",
  pvalue_thresholds = c(0.01, 0.05, 0.10),
  threshold_col    = "padj",
  rep_id           = 11,
  contrast         = c("Group", "Breast.Cancer.Tissue", "Normal.Adjacent.Tissue")
)

print(res$metrics)


##### 2. DESeq2 with  LTR and reduce 1 #####
library(dplyr)
library(DESeq2)
library(tibble)

set.seed(42)

run_DESeq2_LTR_with_truth <- function(bsj_counts_path,
                                      metadata_path,
                                      truth_set_path,
                                      pvalue_thresholds = c(0.01, 0.05, 0.10),
                                      threshold_col = c("padj"),
                                      rep_id = NA,
                                      contrast = c("Group", "Breast.Cancer.Tissue", "Normal.Adjacent.Tissue"),
                                      quiet = TRUE) {
  
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  tryCatch({
    
    # Step 1: Load counts
    cat("Step 1: Loading counts...\n")
    bsj_fc <- read.csv(bsj_counts_path)
    
    if ("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # DESeq2 wants integer counts
    bsj_fc <- as.matrix(bsj_fc)
    storage.mode(bsj_fc) <- "integer"
    
    n_genes_initial <- nrow(bsj_fc)
    
    # Step 2: Load metadata
    cat("Step 2: Loading metadata...\n")
    bsj_fc_md <- read.csv(metadata_path)
    
    # Ensure sample_id exists (many of your design files may not include it)
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      bsj_fc_md$sample_id <- colnames(bsj_fc)
    }
    
    # Reorder metadata to match counts
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Make Group DESeq2-compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # Step 3: Load truth set
    cat("Step 3: Loading truth set...\n")
    truth_set <- read.csv(truth_set_path)
    if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
      stop("Truth file must contain columns: X, DE.ind")
    }
    
    truth_set$DE.ind <- as.logical(truth_set$DE.ind)
    
    # Important: removes duplicates
    truth_set <- truth_set %>% distinct(X, .keep_all = TRUE)
    
    n_truth_initial <- length(unique(truth_set$X))
    
    # Step 4: DESeq2 object
    cat("Step 4: Creating DESeqDataSet...\n")
    
    dds <- DESeqDataSetFromMatrix(
      countData = bsj_fc,
      colData = bsj_fc_md,
      design = ~ Group
    )
    
    
    # Step 5: Run DESeq2
    cat("Step 5: Running DESeq2...\n")
    
    dds <- DESeq(dds,
                 test = "LRT",
                 reduced = ~ 1)
    
    
    # Step 6: Results
    cat("Step 6: Extracting results...\n")
    res <- results(dds,
                   contrast = c("Group",
                                "Breast.Cancer.Tissue",
                                "Normal.Adjacent.Tissue"))
    
    res_df <- as.data.frame(res)
    
    # Tested universe: non-NA in the chosen threshold column
    res_df <- as.data.frame(res) %>% rownames_to_column("X")
    tested_df <- res_df %>% filter(!is.na(.data[[threshold_col]]))
    
    # these are the tested circRNAs for DE
    tested_ids <- tested_df$X
    
    # Genes that were tested in the DE analysis by DESeq2
    n_genes_tested <- length(tested_ids)
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf("Retained %d/%d genes (%.1f%%) in tested set\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Restrict truth to tested universe
    truth_tested <- truth_set %>% filter(X %in% tested_ids)
    
    DE_truth     <- truth_tested %>% filter(DE.ind == TRUE)
    nonDE_truth  <- truth_tested %>% filter(DE.ind == FALSE)
    
    # Helpers
    safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
    
    # Store metrics per threshold (tidy format: one row per threshold)
    metrics <- lapply(pvalue_thresholds, function(th) {
      
      called <- tested_df %>% filter(.data[[threshold_col]] < th)
      not_called <- tested_df %>% filter(!(X %in% called$X))
      
      TP <- nrow(inner_join(called,     DE_truth,    by = "X"))
      FP <- nrow(inner_join(called,     nonDE_truth, by = "X"))
      TN <- nrow(inner_join(not_called, nonDE_truth, by = "X"))
      FN <- nrow(inner_join(not_called, DE_truth,    by = "X"))
      
      precision <- safe_div(TP, TP + FP)
      recall    <- safe_div(TP, TP + FN)
      specificity <- safe_div(TN, TN + FP)
      fpr <- safe_div(FP, FP + TN)
      fdr <- safe_div(FP, TP + FP)
      f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                   NA_real_,
                   2 * precision * recall / (precision + recall))
      
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        
        n_genes_initial = n_genes_initial,
        n_truth_initial = n_truth_initial,
        n_genes_tested = n_genes_tested,
        pct_filtered = pct_filtered,
        
        TP = TP, FP = FP, TN = TN, FN = FN,
        
        sensitivity = recall * 100,
        specificity = specificity * 100,
        precision = precision,
        recall = recall,
        F1 = f1,
        FDR = fdr * 100,
        FPR = fpr * 100,
        
        n_called = nrow(called)
      )
    }) %>% bind_rows()
    
    metrics$mean_pval <- round(mean(res_df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval <- round(median(res_df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval <- round(min(res_df$pvalue, na.rm = TRUE), 6)
    
    metrics$runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
    
    return(list(metrics = metrics, de_results = res_df))
    
  }, error = function(e) {
    
    cat(sprintf("ERROR: %s\n", e$message))
    
    # Return a tidy error table with one row per threshold (helps batch binding)
    error_rows <- lapply(pvalue_thresholds, function(th) {
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        n_genes_initial = NA,
        n_truth_initial = NA,
        n_genes_tested = NA,
        pct_filtered = NA,
        TP = NA, FP = NA, TN = NA, FN = NA,
        sensitivity = NA, specificity = NA,
        precision = NA, recall = NA, F1 = NA,
        FDR = NA, FPR = NA,
        n_called = NA,
        mean_pval = NA, median_pval = NA, min_pval = NA,
        runtime_sec = NA,
        error_message = e$message
      )
    }) %>% bind_rows()
    
    return(list(metrics = error_rows))
  })
}

set.seed(42)

res <- run_DESeq2_LTR_with_truth(
  bsj_counts_path  = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep11_counts.csv",
  metadata_path    = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep11_design.csv",
  truth_set_path   = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep11_TP.csv",
  pvalue_thresholds = c(0.01, 0.05, 0.10),
  threshold_col    = "padj",
  rep_id           = 11,
  contrast         = c("Group", "Breast.Cancer.Tissue", "Normal.Adjacent.Tissue")
)

print(res$metrics)


##### 3. DESeq2 default pipeline with beta.prior = TRUE #####
library(dplyr)
library(DESeq2)
library(tibble)

run_DESeq2_BetaPriorT_with_truth <- function(bsj_counts_path,
                                             metadata_path,
                                             truth_set_path,
                                             pvalue_thresholds = c(0.01, 0.05, 0.10),
                                             threshold_col = c("padj"),
                                             rep_id = NA,
                                             contrast = c("Group", "Breast.Cancer.Tissue", "Normal.Adjacent.Tissue"),
                                             quiet = TRUE) {
  
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  tryCatch({
    
    # Step 1: Load counts
    cat("Step 1: Loading counts...\n")
    bsj_fc <- read.csv(bsj_counts_path)
    
    if ("X" %in% colnames(bsj_fc)) {
      rownames(bsj_fc) <- bsj_fc$X
      bsj_fc$X <- NULL
    }
    
    # DESeq2 wants integer counts
    bsj_fc <- as.matrix(bsj_fc)
    storage.mode(bsj_fc) <- "integer"
    
    n_genes_initial <- nrow(bsj_fc)
    
    # Step 2: Load metadata
    cat("Step 2: Loading metadata...\n")
    bsj_fc_md <- read.csv(metadata_path)
    
    # Ensure sample_id exists (many of your design files may not include it)
    if (!"sample_id" %in% colnames(bsj_fc_md)) {
      bsj_fc_md$sample_id <- colnames(bsj_fc)
    }
    
    # Reorder metadata to match counts
    bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
    if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
      stop("Sample IDs in counts and metadata do not match 1:1.")
    }
    
    # Make Group DESeq2-compatible
    bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
    bsj_fc_md$Group <- factor(bsj_fc_md$Group)
    
    # Step 3: Load truth set
    cat("Step 3: Loading truth set...\n")
    truth_set <- read.csv(truth_set_path)
    if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
      stop("Truth file must contain columns: X, DE.ind")
    }
    
    truth_set$DE.ind <- as.logical(truth_set$DE.ind)
    
    # Important: removes duplicates
    truth_set <- truth_set %>% distinct(X, .keep_all = TRUE)
    
    n_truth_initial <- length(unique(truth_set$X))
    
    # Step 4: DESeq2 object
    cat("Step 4: Creating DESeqDataSet...\n")
    
    dds <- DESeqDataSetFromMatrix(
      countData = bsj_fc,
      colData = bsj_fc_md,
      design = ~ Group
    )
    
    
    # Step 5: Run DESeq2
    cat("Step 5: Running DESeq2...\n")
    
    dds <- DESeq(dds,
                 betaPrior = TRUE)
    
    # Step 6: Results
    cat("Step 6: Extracting results...\n")
    res <- results(dds,
                   contrast = c("Group",
                                "Breast.Cancer.Tissue",
                                "Normal.Adjacent.Tissue"))
    
    res_df <- as.data.frame(res)
    
    # Tested universe: non-NA in the chosen threshold column
    res_df <- as.data.frame(res) %>% rownames_to_column("X")
    tested_df <- res_df %>% filter(!is.na(.data[[threshold_col]]))
    
    # these are the tested circRNAs for DE
    tested_ids <- tested_df$X
    
    # Genes that were tested in the DE analysis by DESeq2
    n_genes_tested <- length(tested_ids)
    pct_filtered <- round((n_genes_tested / n_genes_initial) * 100, 2)
    
    cat(sprintf("Retained %d/%d genes (%.1f%%) in tested set\n",
                n_genes_tested, n_genes_initial, pct_filtered))
    
    # Restrict truth to tested universe
    truth_tested <- truth_set %>% filter(X %in% tested_ids)
    
    DE_truth     <- truth_tested %>% filter(DE.ind == TRUE)
    nonDE_truth  <- truth_tested %>% filter(DE.ind == FALSE)
    
    # Helpers
    safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
    
    # Store metrics per threshold (tidy format: one row per threshold)
    metrics <- lapply(pvalue_thresholds, function(th) {
      
      called <- tested_df %>% filter(.data[[threshold_col]] < th)
      not_called <- tested_df %>% filter(!(X %in% called$X))
      
      TP <- nrow(inner_join(called,     DE_truth,    by = "X"))
      FP <- nrow(inner_join(called,     nonDE_truth, by = "X"))
      TN <- nrow(inner_join(not_called, nonDE_truth, by = "X"))
      FN <- nrow(inner_join(not_called, DE_truth,    by = "X"))
      
      precision <- safe_div(TP, TP + FP)
      recall    <- safe_div(TP, TP + FN)
      specificity <- safe_div(TN, TN + FP)
      fpr <- safe_div(FP, FP + TN)
      fdr <- safe_div(FP, TP + FP)
      f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                   NA_real_,
                   2 * precision * recall / (precision + recall))
      
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        
        n_genes_initial = n_genes_initial,
        n_truth_initial = n_truth_initial,
        n_genes_tested = n_genes_tested,
        pct_filtered = pct_filtered,
        
        TP = TP, FP = FP, TN = TN, FN = FN,
        
        sensitivity = recall * 100,
        specificity = specificity * 100,
        precision = precision,
        recall = recall,
        F1 = f1,
        FDR = fdr * 100,
        FPR = fpr * 100,
        
        n_called = nrow(called)
      )
    }) %>% bind_rows()
    
    metrics$mean_pval <- round(mean(res_df$pvalue, na.rm = TRUE), 4)
    metrics$median_pval <- round(median(res_df$pvalue, na.rm = TRUE), 4)
    metrics$min_pval <- round(min(res_df$pvalue, na.rm = TRUE), 6)
    
    metrics$runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
    
    return(list(metrics = metrics, de_results = res_df))
    
  }, error = function(e) {
    
    cat(sprintf("ERROR: %s\n", e$message))
    
    # Return a tidy error table with one row per threshold (helps batch binding)
    error_rows <- lapply(pvalue_thresholds, function(th) {
      data.frame(
        rep_id = rep_id,
        threshold_col = threshold_col,
        threshold = th,
        n_genes_initial = NA,
        n_truth_initial = NA,
        n_genes_tested = NA,
        pct_filtered = NA,
        TP = NA, FP = NA, TN = NA, FN = NA,
        sensitivity = NA, specificity = NA,
        precision = NA, recall = NA, F1 = NA,
        FDR = NA, FPR = NA,
        n_called = NA,
        mean_pval = NA, median_pval = NA, min_pval = NA,
        runtime_sec = NA,
        error_message = e$message
      )
    }) %>% bind_rows()
    
    return(list(metrics = error_rows))
  })
}

set.seed(42)

res <- run_DESeq2_BetaPriorT_with_truth(
  bsj_counts_path  = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep2_counts.csv",
  metadata_path    = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep2_design.csv",
  truth_set_path   = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep2_TP.csv",
  pvalue_thresholds = c(0.01, 0.05, 0.10),
  threshold_col    = "padj",
  rep_id           = 11,
  contrast         = c("Group", "Breast.Cancer.Tissue", "Normal.Adjacent.Tissue")
)

print(res$metrics)

