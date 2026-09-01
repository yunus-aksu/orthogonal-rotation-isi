# ==============================================================================
# Beyond Varimax: Algorithmic Stability and Variance Redistribution in Orthogonal Component Rotation
# File: isi_calculator.R
# Purpose: Calculates the Interpretability Stability Index (ISI) and performs sensitivity analysis
# ==============================================================================

library(dplyr)

#' Calculate Composite Interpretability Stability Index (ISI)
#' 
#' @param freq_data Dataframe containing scenario frequency wins (ds, np, kk, kk_tip, measure, method, freq)
#' @param w Numeric weight assigned to Global Dominance (default = 0.50 for balanced ISI)
#' @return A list containing GD, SD, and composite ISI scores
calculate_isi <- function(freq_data, w = 0.50) {
  
  # Total number of variable-level competitions across the dataset
  grand_total_races <- sum(freq_data$freq, na.rm = TRUE)
  
  # 1. Global Dominance (GD)
  method_proportions <- freq_data %>%
    group_by(method) %>%
    summarise(total_wins = sum(freq, na.rm = TRUE), .groups = "drop") %>%
    mutate(proportion = total_wins / grand_total_races)
  
  gd_val <- max(method_proportions$proportion, na.rm = TRUE)
  dominant_method <- method_proportions$method[which.max(method_proportions$proportion)]
  
  # 2. Selection Decisiveness (SD)
  scenario_decisiveness <- freq_data %>%
    group_by(ds, np, kk, kk_tip, measure) %>%
    summarise(
      max_freq = max(freq, na.rm = TRUE),
      current_ds = first(ds),
      .groups = "drop"
    ) %>%
    mutate(decisiveness_ratio = max_freq / (10000 * current_ds))
  
  sd_val <- mean(scenario_decisiveness$decisiveness_ratio, na.rm = TRUE)
  
  # 3. Composite ISI Score
  isi_score <- w * gd_val + (1 - w) * sd_val
  
  return(list(
    dominant_method = dominant_method,
    Global_Dominance = round(gd_val, 4),
    Selection_Decisiveness = round(sd_val, 4),
    ISI_Score = round(isi_score, 4)
  ))
}

# Example execution if MASTER_scenario_freq.rds is loaded
if (file.exists("MASTER_scenario_freq.rds")) {
  df_master <- readRDS("MASTER_scenario_freq.rds")
  results <- calculate_isi(df_master)
  
  cat("=== SYSTEM-LEVEL ISI METRICS ===\n")
  cat(sprintf("Dominant Method: %s\n", results$dominant_method))
  cat(sprintf("Global Dominance (GD): %.4f\n", results$Global_Dominance))
  cat(sprintf("Selection Decisiveness (SD): %.4f\n", results$Selection_Decisiveness))
  cat(sprintf("Composite ISI Score: %.4f\n", results$ISI_Score))
}
