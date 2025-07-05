# chat gpt with all variables
#Mar not precise
#missForest and Hot deck to be added
# Load required libraries
library(missMethods)
library(SimCorrMix)
library(vcd)
library(dplyr)

#–– 1. Define simulation parameters ––
reps              <- 1
n_values          <- c(100,500)          # e.g. c(100, 500)
m_values          <- c(6,30)           # e.g. c(6, 30)
rho_vals          <- c(0.1,0.4,0.7)         # e.g. c(0.1, 0.4, 0.7)
p_vals            <- seq(0.1,0.5,0.1)         # e.g. seq(0.1, 0.5, 0.1)
level_scenarios   <- c("few", "many")
dist_scenarios    <- c("balanced", "unbalanced")
simulation_count  <-reps*length(n_values)*length(m_values)*length(rho_vals)*length(level_scenarios)*length(dist_scenarios)
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

# 2d) build marginal CDF for one variable
make_marginal <- function(k, balanced = TRUE) {
  if (balanced) {
    probs <- rep(1/k, k)
  } else {
    # first k−1 levels get 0.1 each, last level takes the rest
    probs <- c(rep(0.1, k-1), 1 - 0.1*(k-1))
  }
  cumsum(probs) - 1e-10
}

#–– 3. Prepare storage for raw false‐counts ––
simulation_summary <- data.frame(
  rep         = integer(),
  n           = integer(),
  m           = integer(),
  rho         = numeric(),
  p           = numeric(),
  levels      = character(),   # "few" or "many"
  dist_form   = character(),   # "balanced" or "unbalanced"
  mechanism   = character(),   # "MCAR", "MAR", "MNAR"
  false_count = integer(),
  stringsAsFactors = FALSE
)

#–– 4. Main simulation loops ––
for (rep in seq_len(reps)) {
  set.seed(1233 + rep)
  
  for (n in n_values) for (m in m_values) {
    
    # which columns to delete
    delete_cols <- seq(2, m, by = 2)
    
    for (levels_sce in level_scenarios) {
      # choose how many levels each variable has
      if (levels_sce == "few") {
        levels_vec <- sample(2:4, m, replace = TRUE)
      } else {
        levels_vec <- sample(5:7, m, replace = TRUE)
      }
      
      for (dist_sce in dist_scenarios) {
        # build marginal CDFs for each variable
        balanced_flag <- (dist_sce == "balanced")
        marginal_list <- mapply(
          make_marginal,
          k        = levels_vec,
          balanced = balanced_flag,
          SIMPLIFY = FALSE
        )
        support_list <- lapply(levels_vec, seq_len)
        
        for (rho in rho_vals) {
          # build correlation matrix
          rho_mat <- matrix(rho, nrow = m, ncol = m)
          diag(rho_mat) <- 1
          
          # generate complete categorical data
          sim <- corrvar2(
            k_cat    = m,
            k_cont   = 0,
            marginal = marginal_list,
            support  = support_list,
            method   = "Polynomial",
            rho      = rho_mat,
            n        = n
          )
          complete_data <- as.data.frame(sim$Y_cat)
          complete_data[] <- lapply(complete_data, factor)
          
          #print remaining datasets count
          cat("Remaining Datasets:")
          simulation_count <- simulation_count-1
          print(simulation_count)
          
          for (p_miss in p_vals) {
            # apply missingness
            df_mcar <- delete_MCAR  (complete_data, p = p_miss, cols_mis = delete_cols)
            df_mar  <- delete_MAR_1_to_x(complete_data, p = p_miss,
                                         cols_mis = delete_cols,
                                         cols_ctrl = delete_cols - 1, x = 3)
            df_mnar <- delete_MNAR_1_to_x(complete_data, p = p_miss,
                                          cols_mis = delete_cols, x = 3)
            
            # impute
            imp_mcar <- impute_mode(df_mcar)
            imp_mar  <- impute_mode(df_mar)
            imp_mnar <- impute_mode(df_mnar)
            
            # count errors
            f_mcar <- count_false(complete_data, df_mcar, imp_mcar)
            f_mar  <- count_false(complete_data, df_mar,  imp_mar)
            f_mnar <- count_false(complete_data, df_mnar, imp_mnar)
            
            # record
            simulation_summary <- bind_rows(
              simulation_summary,
              data.frame(
                rep         = rep, n = n, m = m, rho = rho, p = p_miss,
                levels      = levels_sce,
                dist_form   = dist_sce,
                mechanism   = "MCAR",
                false_count = f_mcar,
                stringsAsFactors = FALSE
              ),
              data.frame(
                rep         = rep, n = n, m = m, rho = rho, p = p_miss,
                levels      = levels_sce,
                dist_form   = dist_sce,
                mechanism   = "MAR",
                false_count = f_mar,
                stringsAsFactors = FALSE
              ),
              data.frame(
                rep         = rep, n = n, m = m, rho = rho, p = p_miss,
                levels      = levels_sce,
                dist_form   = dist_sce,
                mechanism   = "MNAR",
                false_count = f_mnar,
                stringsAsFactors = FALSE
              )
            )
          }
        }
      }
    }
  }
}

#–– 5. Compute per‐rep error rates ––
rep_summary <- simulation_summary %>%
  group_by(rep, n, m, rho, p, levels, dist_form, mechanism) %>%
  summarise(
    total_false = sum(false_count),
    error_rate  = total_false / (m * n),
    .groups = "drop"
  )

#–– 6. Average over replicates ––
final_summary <- rep_summary %>%
  group_by(n, m, rho, p, levels, dist_form, mechanism) %>%
  summarise(
    mean_error_rate = mean(error_rate),
    sd_error_rate   = sd(error_rate),
    .groups = "drop"
  )

#–– 7. Inspect results ––
end_time <- Sys.time()
run_time <- end_time-start_time
print(run_time)

print(final_summary)