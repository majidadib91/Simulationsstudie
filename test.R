#not added : 4th criteria, last methode
# testing everything
# Load required libraries
library(missForest)
library(missMethods)
library(SimCorrMix)
library(vcd)      # for assocstats()
library(dplyr)
library(VIM)      # for hotdeck

#–– 1. Define simulation parameters ––
reps            <- 100
n_values        <- c(100)
m_values        <- c(6)
rho_vals        <- c(0.7)
p_vals          <- seq(0.1, 0.5, 0.1)
level_scenarios <- c("few")
dist_scenarios  <- c("balanced")
simulation_count <- reps *
  length(n_values) *
  length(m_values) *
  length(rho_vals) *
  length(level_scenarios) *
  length(dist_scenarios)
start_time <- Sys.time()

#–– 2. Helpers ––

# 2a) mode function
get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  tbl <- table(x)
  names(tbl)[which.max(tbl)][1]
}

# 2b) mode‐imputation
impute_mode <- function(df) {
  df_imp <- df
  for (j in seq_along(df_imp)) {
    m_j <- get_mode(df_imp[[j]])
    df_imp[[j]][is.na(df_imp[[j]])] <- m_j
  }
  df_imp
}

# 2c) count false imputations
count_false <- function(orig, deleted, imputed) {
  mask <- is.na(deleted)
  sum(imputed[mask] != orig[mask])
}

# 2d) build marginal CDF
make_marginal <- function(k, balanced = TRUE) {
  if (balanced) {
    probs <- rep(1/k, k)
  } else {
    probs <- c(rep(0.1, k-1), 1 - 0.1*(k-1))
  }
  cumsum(probs) - 1e-10
}

# 2e) compute VA (distributional deviation)
compute_VA <- function(orig_df, imp_df, delete_cols, n) {
  d_k <- sapply(delete_cols, function(k) {
    tab_sim <- table(orig_df[[k]])
    tab_imp <- table(imp_df[[k]])
    cats    <- union(names(tab_sim), names(tab_imp))
    f_sim   <- as.numeric(tab_sim[cats]) / n
    f_imp   <- as.numeric(tab_imp[cats]) / n
    f_sim[is.na(f_sim)] <- 0
    f_imp[is.na(f_imp)] <- 0
    sum(abs(f_imp - f_sim)) / 2
  })
  mean(d_k)
}

# 2f) compute RMSE of Cramér’s V over all pairs
compute_RMSEC <- function(orig_df, imp_df) {
  m   <- ncol(orig_df)
  prs <- combn(1:m, 2)
  sq_diff <- apply(prs, 2, function(idx) {
    k <- idx[1]; ℓ <- idx[2]
    tab_sim <- table(orig_df[[k]], orig_df[[ℓ]])
    tab_imp <- table(imp_df[[k]],  imp_df[[ℓ]])
    v_sim   <- assocstats(tab_sim)$cramer
    v_imp   <- assocstats(tab_imp)$cramer
    (v_imp - v_sim)^2
  })
  sqrt(mean(sq_diff, na.rm = TRUE))
}

#–– 3. Prepare storage ––
simulation_summary <- data.frame(
  rep         = integer(),
  n           = integer(),
  m           = integer(),
  rho         = numeric(),
  p           = numeric(),
  levels      = character(),
  dist_form   = character(),
  mechanism   = character(),
  method      = character(),   # Mode, HotDeck, MissForest
  false_count = integer(),
  va_ab       = numeric(),     # distributional deviation
  rmse_c      = numeric(),     # RMSE of Cramér’s V
  stringsAsFactors = FALSE
)

#–– 4. Main simulation ––
for (rep in seq_len(reps)) {
  set.seed(1233 + rep)
  for (n in n_values) for (m in m_values) {
    delete_cols <- seq(2, m, by = 2)
    for (levels_sce in level_scenarios) {
      levels_vec <- if (levels_sce == "few") sample(2:4, m, TRUE)
      else                sample(5:7, m, TRUE)
      for (dist_sce in dist_scenarios) {
        balanced_flag <- (dist_sce == "balanced")
        marginal_list <- mapply(make_marginal,
                                k = levels_vec,
                                balanced = balanced_flag,
                                SIMPLIFY = FALSE)
        support_list <- lapply(levels_vec, seq_len)
        for (rho in rho_vals) {
          rho_mat <- matrix(rho, m, m); diag(rho_mat) <- 1
          sim <- corrvar2(
            k_cat    = m, k_cont = 0,
            marginal = marginal_list,
            support  = support_list,
            method   = "Polynomial",
            rho      = rho_mat,
            n        = n
          )
          complete_data <- as.data.frame(sim$Y_cat)
          complete_data[] <- lapply(complete_data, factor)
          
          cat("Remaining scenarios:", simulation_count, "\n")
          simulation_count <- simulation_count - 1
          
          for (p_miss in p_vals) {
            df_list <- list(
              MCAR = delete_MCAR(complete_data, p = p_miss, cols_mis = delete_cols),
              MAR  = delete_MAR_1_to_x(complete_data, p = p_miss,
                                       cols_mis  = delete_cols,
                                       cols_ctrl = delete_cols - 1, x = 3),
              MNAR = delete_MNAR_1_to_x(complete_data, p = p_miss,
                                        cols_mis = delete_cols, x = 3)
            )
            
            for (mech in names(df_list)) {
              deleted_df <- df_list[[mech]]
              
              # — Mode —
              imp_mode <- impute_mode(deleted_df)
              f_mode   <- count_false(complete_data, deleted_df, imp_mode)
              va_mode  <- compute_VA(complete_data, imp_mode, delete_cols, n)
              rmse_mode<- compute_RMSEC(complete_data, imp_mode)
              
              # — HotDeck —
              imp_hd   <- hotdeck(deleted_df, imp_var = FALSE)
              f_hd     <- count_false(complete_data, deleted_df, imp_hd)
              va_hd    <- compute_VA(complete_data, imp_hd, delete_cols, n)
              rmse_hd  <- compute_RMSEC(complete_data, imp_hd)
              
              # — MissForest —
              imp_mf         <- missForest(deleted_df, verbose = FALSE)$ximp
              imp_mf[]       <- lapply(imp_mf, factor)
              f_mf           <- count_false(complete_data, deleted_df, imp_mf)
              va_mf          <- compute_VA(complete_data, imp_mf, delete_cols, n)
              rmse_mf        <- compute_RMSEC(complete_data, imp_mf)
              
              # record results
              simulation_summary <- bind_rows(
                simulation_summary,
                data.frame(
                  rep         = rep, n = n, m = m, rho = rho, p = p_miss,
                  levels      = levels_sce,
                  dist_form   = dist_sce,
                  mechanism   = mech,
                  method      = "Mode",
                  false_count = f_mode,
                  va_ab       = va_mode,
                  rmse_c      = rmse_mode,
                  stringsAsFactors = FALSE
                ),
                data.frame(
                  rep         = rep, n = n, m = m, rho = rho, p = p_miss,
                  levels      = levels_sce,
                  dist_form   = dist_sce,
                  mechanism   = mech,
                  method      = "HotDeck",
                  false_count = f_hd,
                  va_ab       = va_hd,
                  rmse_c      = rmse_hd,
                  stringsAsFactors = FALSE
                ),
                data.frame(
                  rep         = rep, n = n, m = m, rho = rho, p = p_miss,
                  levels      = levels_sce,
                  dist_form   = dist_sce,
                  mechanism   = mech,
                  method      = "MissForest",
                  false_count = f_mf,
                  va_ab       = va_mf,
                  rmse_c      = rmse_mf,
                  stringsAsFactors = FALSE
                )
              )
            }
          }
        }
      }
    }
  }
}

#–– 5. Summarise per‐rep error‐rates ––
rep_summary <- simulation_summary %>%
  group_by(rep, n, m, rho, p, levels, dist_form, mechanism, method) %>%
  summarise(
    total_false   = sum(false_count),
    error_rate    = total_false / (m * n),
    va_abweichung = mean(va_ab),
    rmse_C        = mean(rmse_c),
    .groups       = "drop"
  )

#–– 6. Average across replicates, wide for the two methods ––
final_summary <- rep_summary %>%
  group_by(n, m, rho, p, levels, dist_form, mechanism) %>%
  summarise(
    mean_error_mode       = mean(error_rate[method == "Mode"]),
    mean_va_mode          = mean(va_abweichung[method == "Mode"]),
    mean_rmse_mode        = mean(rmse_C[    method == "Mode"]),
    mean_error_hotdeck    = mean(error_rate[method == "HotDeck"]),
    mean_va_hotdeck       = mean(va_abweichung[method == "HotDeck"]),
    mean_rmse_hotdeck     = mean(rmse_C[    method == "HotDeck"]),
    mean_error_missforest = mean(error_rate[method == "MissForest"]),
    mean_va_missforest    = mean(va_abweichung[method == "MissForest"]),
    mean_rmse_missforest  = mean(rmse_C[    method == "MissForest"]),
    .groups               = "drop"
  )

#–– 7. Done ––
end_time <- Sys.time()
cat("Total run-time:", end_time - start_time, "\n")
print(final_summary)
View(final_summary)
df1<- filter(final_summary,mechanism == "MCAR")
View(df1)
