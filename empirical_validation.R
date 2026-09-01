# ==============================================================================
# Beyond Varimax: Algorithmic Stability and Variance Redistribution in Orthogonal Component Rotation
# File: empirical_validation.R
# Purpose: Empirical application on contrasting dataset topologies (Matching Table 3 of the Manuscript)
# ==============================================================================

if(!require(psych)) install.packages("psych")
if(!require(GPArotation)) install.packages("GPArotation")
if(!require(lavaan)) install.packages("lavaan")

library(psych)
library(GPArotation)
library(lavaan)

# 1. Custom Rotation Wrapper (Bypasses psych's string limitations)
apply_rotation <- function(loadings, method, k_val) {
  res <- tryCatch({
    suppressWarnings({
      if (method == "VARIMAX")          GPArotation::Varimax(loadings)$loadings
      else if (method == "QUARTIMAX")   GPArotation::quartimax(loadings)$loadings
      else if (method == "EQUAMAX")     GPArotation::equamax(loadings)$loadings
      else if (method == "VARIMIN")     GPArotation::varimin(loadings)$loadings
      else if (method == "GEOMINT")     GPArotation::geominT(loadings)$loadings
      else if (method == "ENTROPY")     GPArotation::entropy(loadings)$loadings
      else if (method == "INFOMAXT")    GPArotation::infomaxT(loadings)$loadings
      else if (method == "BIFACTORT")   GPArotation::bifactorT(loadings)$loadings
      else if (method == "BENTLERT")    GPArotation::bentlerT(loadings)$loadings
      else if (method == "CF-VARIMAX")  GPArotation::cfT(loadings, kappa = 1/k_val)$loadings
      else if (method == "NONE")        loadings
    })
  }, error = function(e) NULL)
  return(res)
}

# 2. Performance Metrics Calculators (Top 3 VAF & Mean Loading Gap Diff12)
calc_vaf_3 <- function(rot_loadings, p_val) {
  sum(rot_loadings[, 1:3]^2) / p_val * 100
}

calc_diff12 <- function(rot_loadings) {
  mean(apply(abs(rot_loadings), 1, function(x) {
    sorted_x <- sort(x, decreasing = TRUE)
    sorted_x[5] - sorted_x[6]
  }))
}

# 3. List of All 10 Orthogonal Rotations + Unrotated
rotations <- c("NONE", "VARIMAX", "QUARTIMAX", "EQUAMAX", "CF-VARIMAX", 
               "VARIMIN", "GEOMINT", "ENTROPY", "INFOMAXT", "BIFACTORT", "BENTLERT")

# 4. Master Analysis Engine
analyze_topology <- function(df, dataset_name) {
  p_val <- ncol(df)
  cat("\n======================================================================\n")
  cat(sprintf("DATASET TOPOLOGY: %s (p = %d)", dataset_name, p_val), "\n")
  cat("======================================================================\n\n")
  
  # Unrotated PCA (k = p = 9)
  pca_unrot <- principal(df, nfactors = p_val, rotate = "none")
  unrot_loadings <- unclass(pca_unrot$loadings)
  
  results <- data.frame(
    Method = character(),
    Top3_VAF_Pct = numeric(),
    Mean_Loading_Gap_Diff12 = numeric(),
    Status = character(),
    stringsAsFactors = FALSE
  )
  
  for (rot in rotations) {
    rot_rot <- apply_rotation(unrot_loadings, rot, p_val)
    if (!is.null(rot_rot)) {
      vaf3 <- round(calc_vaf_3(rot_rot, p_val), 2)
      d12  <- round(calc_diff12(rot_rot), 4)
      results <- rbind(results, data.frame(
        Method = rot,
        Top3_VAF_Pct = vaf3,
        Mean_Loading_Gap_Diff12 = d12,
        Status = "Converged"
      ))
    } else {
      results <- rbind(results, data.frame(
        Method = rot,
        Top3_VAF_Pct = NA,
        Mean_Loading_Gap_Diff12 = NA,
        Status = "Convergence Failed"
      ))
    }
  }
  
  print(results)
  return(results)
}

# ==============================================================================
# Execution on Real Dataset Topologies (Matching Table 3 of the Manuscript)
# ==============================================================================

# Topology 1: Holzinger-Swineford 1939 (General Factor / Highly Correlated)
data(HolzingerSwineford1939, package = "lavaan")
hs_data <- HolzingerSwineford1939[, 7:15] # x1 to x9 (p = 9)
analyze_topology(hs_data, "HOLZINGER-SWINEFORD (Highly Correlated / General Factor)")

# Topology 2: Big Five Inventory (Orthogonal Silos / Independent Dimensions)
data(bfi, package = "psych")
bfi_data <- na.omit(bfi[, c("N1", "N2", "N3", "A1", "A2", "A3", "C1", "C2", "C3")]) # 9 items (p = 9)
analyze_topology(bfi_data, "BIG FIVE INVENTORY - BFI (Orthogonal Silos / No General Factor)")
