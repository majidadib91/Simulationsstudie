#Mode Imputation MCAR
#Loading libraries
library(missMethods)
library(SimCorrMix)
library(vcd)


# Definig parameters
n      <- 100
m      <- 6
reps   <- 100
rho_1  <- c(0.1,0.4,0.7)
p <- seq(0.1,0.5,0.1)

sim_data_list <- vector("list", reps)
simulation_summary_mcar_mode <-vector()


# Defining rho matrix
for (rho_1 in rho_1){
  rho_matirx  = matrix(rho_1,m,m) ; diag(rho_matirx) <- 1
  
  for (p_miss in p){
    
    
    
    
    
    
    
    
    
    # Generating the datasets
    for (i in 1:reps){
      set.seed(i)
      levels <- sample(2:4, m, replace = TRUE)
      # Build support & marginal lists
      support  <- lapply(levels, function(k) 1:k)
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
        means = NULL,
        vars = NULL,
        rho = rho_matirx,
        n = n
      )
      nominal_data <- as.data.frame(sim$Y_cat)
      nominal_data[] <- lapply(nominal_data, factor)
      
      # Save  the data used
      sim_data_list[[i]] <- list(
        data = nominal_data)
    }
    
    
    
    # deleting data
    
    
    
    
    sim_mcar_list <- vector("list", reps)
    
    for (i in 1:reps) {
      df0 <- sim_data_list[[i]]$data
      df_mcar <- delete_MCAR(df0, p = p_miss)
      sim_mcar_list[[i]] <- list(
        data   = df_mcar
      )
    }
    
    # mode imputation
    
    
    get_mode <- function(x) {
      # Exclude NA values
      x <- x[!is.na(x)]
      if (length(x) == 0) return(NA)  # If all NA, return NA
      # Get the most frequent value
      tbl <- table(x)
      names(tbl)[which.max(tbl)][1]  # Pick first mode in case of ties
    }
    
    # Impute missing values using mode for each dataset
    for (i in 1:reps) {
      # Extract the data component
      data <- sim_mcar_list[[i]]$data
      # Create a copy for imputation
      imputed_data <- data
      # Apply mode imputation to each column
      for (j in 1:ncol(data)) {
        # Find indices of NA values
        na_indices <- is.na(data[[j]])
        if (any(na_indices)) {
          # Impute NA with mode of the column
          mode_val <- get_mode(data[[j]])
          imputed_data[[j]][na_indices] <- mode_val
        }
      }
      # Ensure imputed column remains a factor with original levels
      imputed_data[] <- lapply(imputed_data, as.factor)
      # Store imputed data in a new component
      sim_data_list[[i]]$imputed_data <- imputed_data
    }
    
    
    
    
    # Step 3a:  PFC
    total_falsely_imputed <- 0
    for (i in 1:reps){
      
      falsely_imputed <- sum(sim_data_list[[i]]$imputed_data != sim_data_list[[i]]$data)
      total_falsely_imputed <- total_falsely_imputed + falsely_imputed
    }
    
    pfc <- round(total_falsely_imputed/(n*m*reps),4)
    
    
    
    
    
    # Step 3b: Cramér's V and RMSE calculation
    rmse_values <- numeric()
    
    for (i in 1:reps) {
      sim_data <- sim_data_list[[i]]$data
      imp_data <- sim_data_list[[i]]$imputed_data
      
      # Initialize vectors to store Cramér's V for all column pairs
      v_sim <- numeric()
      v_imp <- numeric()
      
      # Calculate Cramér's V for each pair of columns
      for (j in 1:(m-1)) {
        for (k in (j+1):m) {
          sim_table <- table(sim_data[[j]], sim_data[[k]])
          imp_table <- table(imp_data[[j]], imp_data[[k]])
          v_sim <- c(v_sim, assocstats(sim_table)$cramer)
          v_imp <- c(v_imp, assocstats(imp_table)$cramer)
        }
      }
      
      # Calculate RMSE for this dataset
      diff_squared <- (v_imp - v_sim)^2
      rmse <- sqrt(2 / (m * (m-1)) * sum(diff_squared, na.rm = TRUE))
      rmse_values[i] <- rmse
      
      # Store average Cramér's V for this dataset (optional)
      sim_data_list[[i]]$mean_cramers_v_sim <- mean(v_sim, na.rm = TRUE)
      sim_data_list[[i]]$mean_cramers_v_imp <- mean(v_imp, na.rm = TRUE)
    }
    
    # Calculate average RMSE across all repetitions
    mean_rmse <- round(mean(rmse_values, na.rm = TRUE),4)
    
    
    
    
    
    
    
    
    
    
    new_row <- data.frame(
      wiederholungen = reps,
      n = n,
      m = m,
      rho_1 = rho_1,
      verteilung = "gleich",
      Ausprägungen = "wenig",
      p_miss = p_miss,
      pfc = pfc,
      mean_rmse = mean_rmse
    )
    simulation_summary_mcar_mode<- rbind(simulation_summary_mcar_mode,new_row)
    
    
    
  }
}

