library(dplyr)
library(edgeR)
library(limma)
library(tibble)

run_limma_vDFBRT_with_truth <- function(bsj_counts_path,
                                 metadata_path,
                                 truth_set_path, # path to the TP file from SPSimSeq
                                 contrast_formula, # unchanged
                                 edgeR_norm = "TMM", # optional -> change to TMM
                                 thresholds = c(0.01, 0.05, 0.10),
                                 threshold_col = c("adj.P.Val"),
                                 rep_id = NA) {
  threshold_col <- match.arg(threshold_col)
  start_time <- Sys.time()
  
  # Step 1: Load the data and the metadata #
  cat("Starting Step 1: Loading the BSJ matrix ...\n")
  
  bsj_fc <- read.csv(bsj_counts_path)
  rownames(bsj_fc) <- bsj_fc$X
  bsj_fc$X <- NULL
  
  # circRNAs before filtering with filterByExpr
  n_genes_initial_counts <- nrow(bsj_fc)
  
  # Load the metadata 
  bsj_fc_md <- read.csv(metadata_path)
  bsj_fc_md$sample_id <- colnames(bsj_fc)
  
  # Sanity Check 1: Check if the sample order matches the metadata order to avoid issues with the dgelist
  
  if (!"sample_id" %in% colnames(bsj_fc_md)) {
    stop("metadata must contain a `sample_id` column matching colnames(bsj_counts).")
  }
  
  bsj_fc_md <- bsj_fc_md[match(colnames(bsj_fc), bsj_fc_md$sample_id), , drop = FALSE]
  
  if (any(is.na(bsj_fc_md$sample_id)) || !identical(colnames(bsj_fc), bsj_fc_md$sample_id)) {
    stop("Sample IDs in counts and metadata do not match 1:1.")
  }
  
  ##### Truth set loading #####
  cat("Starting Step 2: Loading the Truth Set ...\n")
  
  truth_set <- read.csv(truth_set_path)
  
  # Check that the colnames X and DE.ind are present in the truth set
  if (!all(c("X", "DE.ind") %in% colnames(truth_set))) {
    stop("Truth file must contain columns: X, DE.ind")
  }
  
  # to ensure this is logical (True and False values)
  truth_set$DE.ind <- as.logical(truth_set$DE.ind)
  
  # set the number of circRNAs that were in the truth set
  n_genes_initial_truth_set <- length(unique(truth_set$X))
  
  ##### DGEList and filtering #####
  cat("Starting Step 3: Creating the DGEList and starting the automated filtering ...\n")
  
  bsj_fc_md$Group <- make.names(as.character(bsj_fc_md$Group))
  group <- factor(bsj_fc_md$Group)
  
  # must always be only.2 groups for comparison in our case
  if (length(levels(group)) != 2) stop("Expected exactly 2 groups.")
  
  dgelist <- DGEList(counts = bsj_fc, group = group)
  keep <- filterByExpr(dgelist, group = group)
  dgelist_filtered <- dgelist[keep, , keep.lib.sizes = FALSE]
  
  # extract the universe size -> circRNA ids that will be tested by the DE methods
  tested_circRNAs <- rownames(dgelist_filtered)
  n_circRNAs_tested <- length(tested_circRNAs)
  
  # calculate the pct of retained circRNAs post filtering
  pct_retained <- round(n_circRNAs_tested / n_genes_initial_truth_set * 100, 2)
  
  # Restrict truth to TESTED universe
  truth_set_tested <- truth_set %>% filter(X %in% tested_circRNAs)
  
  # If truth has any duplicates, remove them to avoid join inflation
  truth_set_tested <- truth_set_tested %>% distinct(X, .keep_all = TRUE)
  
  ##### Normalizing and setting the contrasts #####
  cat("Starting Step 4: Normalizing the DGEList and setting the contrasts as defined by Alexa ...\n")  
  
  dgelist_filtered_norm <- calcNormFactors(dgelist_filtered, method = edgeR_norm)
  
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
  
  ##### Fit glmQLFit model and test with glmQLFTest the contrasts #####
  cat("Starting Step 5: Fitting voom and lmFit and testing the contrasts ...\n")  
  #### Step 4: Fit voomLmFit with sample.weights set to TRUE voom to the model ####
  # it replaces all the commands with one function
  dgelist_v <- voom(dgelist_filtered_norm, design_mtx)
  
  dgelist_vfit <- lmFit(dgelist_v, design_mtx, robust = TRUE)
  
  # Step 5: Set the contrast matrix and DE analysis
  cat(" Starting Step 5: Testing contrasts... \n")
  
  # Get the group names from design matrix columns
  cat(sprintf(" Testing contrast: %s\n", contrast_formula))
  
  dgelist_vfit <- contrasts.fit(dgelist_vfit, contrasts = contrasts)
  
  dgelist_efit <- eBayes(dgelist_vfit)
  
  # Extract the full results
  summary(decideTests(dgelist_efit))
  de_results <- topTable(dgelist_efit, number = Inf)
  
  ##### Create the results table #####
  cat("Starting Step 6: Calculating the metrics and generating the results table ...\n")  
  de_results <- as.data.frame(de_results)
  de_results$X <- rownames(de_results)
  rownames(de_results) <- NULL
  
  # Ensure we only evaluate on tested IDs that appear in results
  # (should match, but we keep this safe)
  tested_ids_try2 <- de_results$X
  truth_tested2 <- truth_set_tested %>% dplyr::filter(X %in% tested_ids_try2)
  truth_tested2 <- truth_tested2 %>% dplyr::distinct(X, .keep_all = TRUE)
  
  # Extract the true DE (TP) and non DE (TN)
  DE_truth <- truth_tested2 %>% dplyr::filter(DE.ind == TRUE)
  DE_notTruth <- truth_tested2 %>% dplyr::filter(DE.ind == FALSE)
  
  # Helper for safe division
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
  
  # Calculate metrics for every specified threshold
  out <- lapply(thresholds, function(th) {
    
    called <- de_results %>% filter(.data[[threshold_col]] <= th)
    not_called <- de_results %>% filter(!(X %in% called$X))
    
    TP <- nrow(inner_join(called,     DE_truth,  by = "X"))
    FP <- nrow(inner_join(called,     DE_notTruth, by = "X"))
    TN <- nrow(inner_join(not_called, DE_notTruth, by = "X"))
    FN <- nrow(inner_join(not_called, DE_truth,  by = "X"))
    
    precision <- safe_div(TP, TP + FP)
    recall    <- safe_div(TP, TP + FN)  # = sensitivity = TPR
    specificity <- safe_div(TN, TN + FP)
    fpr <- safe_div(FP, FP + TN)
    fdr <- safe_div(FP, TP + FP)
    f1  <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                  NA_real_,
                  2 * precision * recall / (precision + recall))
    
    data.frame(
      rep_id = rep_id,
      threshold_type = threshold_col,
      threshold = th,
      
      initial_genes_truth = n_genes_initial_truth_set,
      initial_genes_counts = n_genes_initial_counts,
      tested_genes = n_circRNAs_tested,
      pct_retained = pct_retained,
      
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
  
  runtime_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  out$runtime_sec <- runtime_sec
  
  list(metrics = out, de_results = de_results)
}

res <- run_limma_vDFBRT_with_truth(
  bsj_counts_path  = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep4_counts.csv",
  metadata_path    = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep4_design.csv",
  truth_set_path   = "/Users/alexakouri/Documents/DE_benchmarking_paper/simulations/BC/autofilter_counts/DE_0.1/BC_autofilter_counts_DE0.1_rep4_TP.csv",
  edgeR_norm       = "TMM",
  contrast_formula = "Breast.Cancer.Tissue-Normal.Adjacent.Tissue",
  thresholds       = c(0.01, 0.05, 0.10),
  threshold_col    = "adj.P.Val",   # recommend FDR for consistency
  rep_id           = 30
)

print(res$metrics)
