# ==============================================================================
# Beyond Varimax: Algorithmic Stability and Variance Redistribution in Orthogonal Component Rotation
# File: empirical_validation.R
# Purpose: Empirical application on contrasting dataset topologies (Holzinger-Swineford & IPIP Big-Five)
# ==============================================================================

# Required packages
required_packages <- c("psych", "lavaan", "GPArotation")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(psych)
library(lavaan)
library(GPArotation)

# Helper function to run 10 rotations and compute variance redistribution
evaluate_empirical_rotations <- function(data_matrix, n_factors = NULL, dataset_name = "") {
  # Handle complete cases
  clean_data <- na.omit(data_matrix)
  p <- ncol(clean_data)
  if (is.null(n_factors)) n_factors <- p # Full component model (k = p)
  
  # 1. Component Extraction (Unrotated PCA)
  pca_fit <- principal(clean_data, nfactors = n_factors, rotate = "none")
  L_unrotated <- as.matrix(unclass(pca_fit$loadings))
  
  # Unrotated variance distribution (%)
  unrot_vars <- colSums(L_unrotated^2) / p * 100
  
  # 2. Define 10 Orthogonal Rotations
  rotations <- list(
    "Varimax"    = function(L) Varimax(L),
    "Quartimax"  = function(L) Quartimax(L),
    "Equamax"    = function(L) Equamax(L),
    "CF-Varimax" = function(L) cfT(L, kappa = 1/p),
    "GeominT"    = function(L) geominT(L),
    "EntropyT"   = function(L) entropyT(L),
    "InfomaxT"   = function(L) infomaxT(L),
    "BentlerT"   = function(L) bentlerT(L),
    "BifactorT"  = function(L) bifactorT(L),
    "Varimin"    = function(L) Orthomax(L, gamma = 0)
  )
  
  results <- data.frame(
    Dataset = dataset_name,
    Method = character(),
    Converged = logical(),
    PC1_Var_Pct = numeric(),
    PC2_Var_Pct = numeric(),
    PC3_Var_Pct = numeric(),
    First3_Total_Pct = numeric(),
    Residual_Var_Pct = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Unrotated Baseline
  results <- rbind(results, data.frame(
    Dataset = dataset_name,
    Method = "Unrotated",
    Converged = TRUE,
    PC1_Var_Pct = round(unrot_vars[1], 2),
    PC2_Var_Pct = round(unrot_vars[2], 2),
    PC3_Var_Pct = round(unrot_vars[3], 2),
    First3_Total_Pct = round(sum(unrot_vars[1:min(3, n_factors)]), 2),
    Residual_Var_Pct = round(ifelse(n_factors > 3, sum(unrot_vars[4:n_factors]), 0), 2)
  ))
  
  # Apply Rotations
  for (name in names(rotations)) {
    rot_fit <- tryCatch({
      rotations[[name]](L_unrotated)
    }, error = function(e) NULL)
    
    if (!is.null(rot_fit)) {
      L_rot <- as.matrix(unclass(rot_fit$loadings))
      var_exp <- colSums(L_rot^2) / p * 100
      
      results <- rbind(results, data.frame(
        Dataset = dataset_name,
        Method = name,
        Converged = TRUE,
        PC1_Var_Pct = round(var_exp[1], 2),
        PC2_Var_Pct = round(var_exp[2], 2),
        PC3_Var_Pct = round(var_exp[3], 2),
        First3_Total_Pct = round(sum(var_exp[1:min(3, n_factors)]), 2),
        Residual_Var_Pct = round(ifelse(n_factors > 3, sum(var_exp[4:n_factors]), 0), 2)
      ))
    } else {
      results <- rbind(results, data.frame(
        Dataset = dataset_name, Method = name, Converged = FALSE,
        PC1_Var_Pct = NA, PC2_Var_Pct = NA, PC3_Var_Pct = NA,
        First3_Total_Pct = NA, Residual_Var_Pct = NA
      ))
    }
  }
  
  return(results)
}

# ==============================================================================
# 1. Dataset 1: Holzinger-Swineford 1939 (General Factor / Correlated Topology)
# ==============================================================================
cat("--- Evaluating Holzinger-Swineford 1939 (9 Cognitive Tests) ---\n")
hs_data <- HolzingerSwineford1939[, paste0("x", 1:9)]
hs_results <- evaluate_empirical_rotations(hs_data, dataset_name = "Holzinger-Swineford (k=9)")
print(hs_results)

# ==============================================================================
# 2. Dataset 2: IPIP Big-Five Factor Markers (Independent Silos Topology)
# ==============================================================================
cat("\n--- Evaluating IPIP Big-Five (25 Personality Items) ---\n")
data(bfi, package = "psych")
bfi_data <- bfi[, 1:25]
bfi_results <- evaluate_empirical_rotations(bfi_data, dataset_name = "IPIP Big-Five (k=25)")
print(bfi_results)
