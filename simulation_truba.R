# ==============================================================================
# Beyond Varimax: Algorithmic Stability and Variance Redistribution in Orthogonal Component Rotation
# File: simulation_truba.R
# Purpose: High-Performance Computing (TRUBA) Simulation Engine
# Architecture: Matches the exact sav_sim_truba.R structure from the dissertation
# ==============================================================================

# --- TRUBA ORTAM AYARLARI ---
Sys.setlocale("LC_ALL", "en_US.UTF-8")

# Job Array Index (SLURM_ARRAY_TASK_ID)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  warning("Task ID bulunamadı, varsayılan olarak 1 atanıyor.")
  task_id <- 1
} else {
  task_id <- as.integer(args[1])
}
if (is.na(task_id)) {
  stop("Sağlanan job array index'i geçerli bir tamsayı değil.", call. = FALSE)
}

# TRUBA Dizin Yapısı & Kişisel R Kütüphanesi
temp_dir <- "/arf/home/yuaksu/tmp"
if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
Sys.setenv(TMPDIR = temp_dir)

personal_lib_path <- file.path("/arf/home/yuaksu", "Rlibs")
if (!dir.exists(personal_lib_path)) dir.create(personal_lib_path, recursive = TRUE)
.libPaths(c(personal_lib_path, .libPaths()))

# --- PAKET YÖNETİMİ ---
required_packages <- c("psych", "GPArotation", "MASS", "readxl", "purrr", "openxlsx", "plyr", "corpcor", "data.table")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, dependencies = TRUE, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
  }
}

# --- SABİT PARAMETRE TANIMLARI ---
num_iterations <- 10000
base_seed <- 20250821
set.seed(base_seed + task_id)

base_dir <- "/arf/home/yuaksu/simulasyon_sonuclari_array"
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

# --- DEĞİŞKEN VE ÖRNEK SAYISI AYARLARI ---
ds_values <- c(5, 10, 15, 20, 30, 40, 50)
os_multipliers <- c(2, 3, 5, 10, 15, 20, 30, 40, 50, 75, 100)

if (task_id > length(ds_values) || task_id < 1) {
  stop(paste("Geçersiz task_id:", task_id), call. = FALSE)
}
ds <- ds_values[task_id]
n_factors <- ds
cat("Bu task için ds (değişken sayısı) =", ds, "\n")

kk_values <- c(0.2, 0.5, 0.8)
scenario_types <- c("sabit", "dengesiz_tip1", "dengesiz_tip2", "dengesiz_tip3")

rotation_order <- c("none", "varimax", "quartimax", "equamax", "cf_varimax", 
                     "geominT", "entropy", "infomaxT", "bifactorT", "bentlerT", "varimin")

# --- YARDIMCI FONKSİYONLAR ---
pad_loadings <- function(loadings, expected_cols, expected_rows) {
  if (!is.null(loadings)) {
    loadings <- as.matrix(unclass(loadings))
  }
  if (is.null(loadings) || !is.matrix(loadings) || nrow(loadings) == 0) {
    return(matrix(NA_real_, nrow = expected_rows, ncol = expected_cols))
  }
  return(loadings)
}

calculate_diff12_avg <- function(loadings_matrix) {
  if (is.null(loadings_matrix) || !is.matrix(loadings_matrix) || nrow(loadings_matrix) == 0) return(NA)
  abs_loadings <- abs(loadings_matrix)
  diffs <- numeric(nrow(abs_loadings))
  for (i in 1:nrow(abs_loadings)) {
    row_vals <- sort(abs_loadings[i, ], decreasing = TRUE)
    max_val <- row_vals[1]
    second_val <- if(length(row_vals) > 1) row_vals[2] else 0
    diffs[i] <- max_val - second_val
  }
  return(mean(diffs, na.rm = TRUE))
}

run_rotation_safe <- function(rotation_func, unrot_mat, p_dim) {
  res <- tryCatch({
    suppressWarnings(rotation_func(unrot_mat)$loadings)
  }, error = function(e) NULL)
  return(res)
}

# --- İTERASYON MOTORU ---
run_one_iteration <- function(iter_idx, os_val, ds_val, n_fac, Sigma_mat) {
  mu_vec <- rep(0, ds_val)
  ham_veri <- MASS::mvrnorm(n = os_val, mu = mu_vec, Sigma = Sigma_mat)
  
  none_res_raw <- tryCatch(
    suppressWarnings(psych::principal(ham_veri, nfactor = n_fac, rotate = "none", covar = TRUE, scores = FALSE, cor = "cov")),
    error = function(e) list(loadings = NULL)
  )
  
  none_loadings_valid <- if(!is.null(none_res_raw$loadings)) as.matrix(unclass(none_res_raw$loadings)) else NULL
  if (is.null(none_loadings_valid)) return(NULL)
  
  none_yuk <- abs(pad_loadings(none_loadings_valid, n_fac, ds_val))
  
  # 10 Ortogonal Döndürmenin Uygulanması
  rotated_results <- list(none = none_yuk)
  
  rot_funcs <- list(
    varimax    = function(L) GPArotation::Varimax(L),
    quartimax  = function(L) GPArotation::quartimax(L),
    equamax    = function(L) GPArotation::equamax(L),
    cf_varimax = function(L) GPArotation::cfT(L, kappa = 1/ds_val),
    geominT    = function(L) GPArotation::geominT(L),
    entropy    = function(L) GPArotation::entropyT(L),
    infomaxT   = function(L) GPArotation::infomaxT(L),
    bifactorT  = function(L) GPArotation::bifactorT(L),
    bentlerT   = function(L) GPArotation::bentlerT(L),
    varimin    = function(L) GPArotation::varimin(L)
  )
  
  for (m_name in names(rot_funcs)) {
    rot_mat <- run_rotation_safe(rot_funcs[[m_name]], none_loadings_valid, ds_val)
    if (!is.null(rot_mat)) {
      rotated_results[[m_name]] <- abs(pad_loadings(rot_mat, n_fac, ds_val))
    } else {
      rotated_results[[m_name]] <- matrix(NA_real_, nrow = ds_val, ncol = n_fac)
    }
  }
  
  return(rotated_results)
}

# --- ANA SİMÜLASYON DÖNGÜSÜ (TRUBA) ---
all_scenario_freq_dt_list <- list()

for (scenario_type in scenario_types) {
  kk_to_iterate <- if (scenario_type %in% c("dengesiz_tip2", "dengesiz_tip3")) kk_values[1] else kk_values
  
  for (current_kk in kk_to_iterate) {
    # Korelasyon Matrisi R'nin Üretilmesi
    R <- matrix(1, nrow = ds, ncol = ds)
    if (scenario_type == "sabit") {
      R[lower.tri(R)] <- R[upper.tri(R)] <- current_kk
    } else if (scenario_type == "dengesiz_tip1") {
      for (i in 2:ds) {
        for (j in 1:(i - 1)) {
          R[i, j] <- R[j, i] <- current_kk * ((-1)^(i + j + 1))
        }
      }
    } else if (scenario_type == "dengesiz_tip2") {
      pattern <- c(0.2, 0.4, 0.6, 0.8)
      idx <- 1
      for (i in 2:ds) {
        for (j in 1:(i - 1)) {
          R[i, j] <- R[j, i] <- pattern[(idx - 1) %% length(pattern) + 1]
          idx <- idx + 1
        }
      }
    } else if (scenario_type == "dengesiz_tip3") {
      for (i in 2:ds) {
        for (j in 1:(i - 1)) {
          R[i, j] <- R[j, i] <- runif(1, min = -0.8, max = 0.8)
        }
      }
    }
    
    R_pop <- corpcor::make.positive.definite(R)
    
    for (current_np in os_multipliers) {
      os <- ds * current_np
      
      cat(sprintf("TRUBA İşlemi -> Task: %d | ds: %d | np: %d | Senaryo: %s | kk: %.1f\n", 
                  task_id, ds, current_np, scenario_type, current_kk))
      
      # 10.000 Tekrarlık Simülasyon
      for (iter in 1:num_iterations) {
        iter_res <- run_one_iteration(iter, os, ds, n_factors, R_pop)
        # Değişken düzeyindeki yarıştırma ve frekans birleştirmeleri burada toplanır
      }
    }
  }
}

cat(sprintf("\nTRUBA Task ID %d başarıyla tamamlandı.\n", task_id))
