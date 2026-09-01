# ==============================================================================
# Beyond Varimax: Algorithmic Stability and Variance Redistribution in Orthogonal Component Rotation
# File: demo_local.R
# Purpose: Local test script reflecting the EXACT data generation mechanism of the dissertation
# ==============================================================================

# Required packages
required_packages <- c("psych", "GPArotation", "MASS", "corpcor")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(psych)
library(GPArotation)
library(MASS)
library(corpcor)

# 1. Simulation Setup (Matching Dissertation & TRUBA Pipeline)
set.seed(20250821)
p <- 10        # Number of variables (ds)
n <- 200       # Sample size (os = n/p * p)
k <- p         # Full component model (k = p)

# 2. Construct Population Correlation Matrix (R)
# Example Scenario: Dengesiz Tip-2 / Heterogeneous Pattern (0.2, 0.4, 0.6, 0.8)
R <- matrix(1, nrow = p, ncol = p)
pattern <- c(0.2, 0.4, 0.6, 0.8)
idx <- 1
for (i in 2:p) {
  for (j in 1:(i-1)) {
    val <- pattern[(idx - 1) %% length(pattern) + 1]
    R[i, j] <- val
    R[j, i] <- val
    idx <- idx + 1
  }
}

# Ensure positive definiteness (corpcor)
R <- make.positive.definite(R)

# 3. Generate Synthetic Data via Multivariate Normal Distribution
X <- mvrnorm(n = n, mu = rep(0, p), Sigma = R)

# 4. Unrotated Component Extraction (PCA, k = p)
pca_res <- principal(X, nfactors = k, rotate = "none")
L_unrotated <- as.matrix(unclass(pca_res$loadings))

# 5. Apply 10 Orthogonal Rotation Criteria
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
  Method = character(),
  Converged = logical(),
  PC1_Var_Pct = numeric(),
  PC2_Var_Pct = numeric(),
  PC3_Var_Pct = numeric(),
  Total_Var_Pct = numeric(),
  stringsAsFactors = FALSE
)

cat("--- Running 10 Orthogonal Rotations on Generated Correlation Structure ---\n\n")

for (name in names(rotations)) {
  rot_fit <- tryCatch({
    rotations[[name]](L_unrotated)
  }, error = function(e) NULL)
  
  if (!is.null(rot_fit)) {
    L_rot <- as.matrix(unclass(rot_fit$loadings))
    var_exp <- colSums(L_rot^2) / p * 100
    results <- rbind(results, data.frame(
      Method = name,
      Converged = TRUE,
      PC1_Var_Pct = round(var_exp[1], 2),
      PC2_Var_Pct = round(var_exp[2], 2),
      PC3_Var_Pct = round(var_exp[3], 2),
      Total_Var_Pct = round(sum(var_exp), 2)
    ))
  } else {
    results <- rbind(results, data.frame(
      Method = name,
      Converged = FALSE,
      PC1_Var_Pct = NA, PC2_Var_Pct = NA, PC3_Var_Pct = NA, Total_Var_Pct = NA
    ))
  }
}

print(results)
cat("\nDemo run matching dissertation data-generation architecture completed!\n")
