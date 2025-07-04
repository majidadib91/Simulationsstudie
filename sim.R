# simple and fast
# Load libraries
library(missMethods)
library(SimCorrMix)
library(vcd)

# Define parameters

n <- (100)
m_values <- c(6)
reps <- 100
rho_1 <- c(0.7)
p <-c(0.5)
mechanisms <- c("mcar")

# Initialize result
simulation_summary <- data.frame()

# Mode function
get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  tbl <- table(x)
  names(tbl)[which.max(tbl)][1]
}

# Start simulation loop
for (m in m_values) {
  
  for (rho_val in rho_1) {
    
    rho_matrix <- matrix(rho_val, m, m)
    diag(rho_matrix) <- 1
    
    for (p_miss in p) {
      
      sim_data_list <- vector("list", reps)
      delete_cols <- seq(2, m, 2)
      
      # Step 1: Generate complete data
      for (i in 1:reps) {
        set.seed(i+100)
        levels <- sample(2:4, m, replace = TRUE)
        support <- lapply(levels, function(k) 1:k)
        marginal <- lapply(levels, function(k) {
          probs <- rep(1/k, k)
          cumsum(probs) - 1e-10
        })
        sim <- corrvar2(
          k_cat = m,
          k_cont = 0,
          marginal = marginal,
          support = support,
          method = "Polynomial",
          rho = rho_matrix,
          n = n
        )
        nominal_data <- as.data.frame(sim$Y_cat)
        nominal_data[] <- lapply(nominal_data, factor)
        sim_data_list[[i]] <- list(data = nominal_data)
      }
      
      # Step 2: Add 3 types of missingness
      for (i in 1:reps) {
        df0 <- sim_data_list[[i]]$data
        
        df_mcar <- delete_MCAR(df0, p = p_miss, cols_mis = delete_cols)
        df_mar  <- delete_MAR_1_to_x(df0, p = p_miss, cols_mis = delete_cols,
                                     cols_ctrl = delete_cols - 1, x = 3)
        df_mnar <- delete_MNAR_1_to_x(df0, p = p_miss, cols_mis = delete_cols, x = 3)
        
        sim_data_list[[i]]$mcar <- df_mcar
        sim_data_list[[i]]$mar  <- df_mar
        sim_data_list[[i]]$mnar <- df_mnar
      }
      
      # Step 3: Impute and evaluate each mechanism
      for (mech in mechanisms) {
        falsely_imputed_total <- 0
        rmse_values <- numeric(reps)
        
        for (i in 1:reps) {
          df0 <- sim_data_list[[i]]$data
          df_mis <- sim_data_list[[i]][[mech]]
          
          # Mode imputation
          imputed_data <- df_mis
          for (j in 1:m) {
            na_indices <- is.na(imputed_data[[j]])
            if (any(na_indices)) {
              mode_val <- get_mode(imputed_data[[j]])
              imputed_data[[j]][na_indices] <- mode_val
            }
          }
          imputed_data[] <- lapply(imputed_data, as.factor)
          sim_data_list[[i]][[paste0(mech, "_imp")]] <- imputed_data
          
          # PFC
          falsely_imputed <- sum(imputed_data != df0, na.rm = TRUE)
          falsely_imputed_total <- falsely_imputed_total + falsely_imputed
          
          # RMSE using Cramér’s V
          v_sim <- numeric()
          v_imp <- numeric()
          for (j in 1:(m - 1)) {
            for (k in (j + 1):m) {
              sim_table <- table(df0[[j]], df0[[k]])
              imp_table <- table(imputed_data[[j]], imputed_data[[k]])
              v_sim <- c(v_sim, assocstats(sim_table)$cramer)
              v_imp <- c(v_imp, assocstats(imp_table)$cramer)
            }
          }
          rmse <- sqrt(2 / (m * (m - 1)) * sum((v_imp - v_sim)^2, na.rm = TRUE))
          rmse_values[i] <- rmse
        }
        
        # Step 4: Store results
        pfc <- round(falsely_imputed_total / (n * m * reps), 4)
        mean_rmse <- round(mean(rmse_values, na.rm = TRUE), 4)
        
        new_row <- data.frame(
          mechanism = toupper(mech),
          wiederholungen = reps,
          n = n,
          m = m,
          rho_1 = rho_val,
          verteilung = "gleich",
          Ausprägungen = "wenig",
          p_miss = p_miss,
          pfc = pfc,
          mean_rmse = mean_rmse
        )
        simulation_summary <- rbind(simulation_summary, new_row)
        cat("Done: m =", m, "rho =", rho_val, "p =", p_miss, "mech =", toupper(mech), "\n")
      }
    }
  }
}

# Final result
print(simulation_summary)

#rm(list = setdiff(ls(), c("test_old","test_even_col","test_eve_col_clevels","test_even_cols_clevels_seed1234","simulation_summary_m_p_mech_rho")))

df3 
