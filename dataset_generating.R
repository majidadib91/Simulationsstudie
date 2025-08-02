library(SimCorrMix) # For generating correlated mixed data
library(dplyr)      # For data manipulation (though not strictly necessary for this part)

#–– 1. Define simulation parameters ––
reps            <- 10 # Now, this 'reps' refers to the number of *complete datasets* to generate for each scenario
#n_values        <- c(100, 500)
n_values        <- c(10)
#m_values        <- c(6, 12, 30)
m_values        <- c(6)
rho_vals        <- c(0.1, 0.4, 0.7)
#level_scenarios <- c("few", "many")
level_scenarios <- c("few")
#dist_scenarios  <- c("balanced", "unbalanced")
dist_scenarios  <- c("balanced")

# Calculate total number of unique scenarios (parameter combinations *within* one rep)
num_scenarios_per_rep <- length(n_values) *
  length(m_values) *
  length(rho_vals) *
  length(level_scenarios) *
  length(dist_scenarios)

# Total number of datasets to generate across all reps
total_datasets_to_generate <- num_scenarios_per_rep * reps

cat("Total repetitions (outermost loop):", reps, "\n")
cat("Total unique scenarios per repetition:", num_scenarios_per_rep, "\n")
cat("Total datasets to generate overall:", total_datasets_to_generate, "\n")

start_time <- Sys.time()

#–– 2. Helpers ––

# 2d) build marginal CDF
make_marginal <- function(k, balanced = TRUE) {
  if (balanced) {
    probs <- rep(1/k, k)
  } else {
    probs <- c(rep(0.1, k-1), 1 - sum(rep(0.1, k-1)))
  }
  cumsum(probs) - 1e-10
}

#–– 3. Prepare storage for generated complete datasets ––
# Using a list of lists:
# Outer list: indexed by rep_index (e.g., "rep_1", "rep_2")
# Next level list: indexed by scenario_id
# Innermost level: the actual dataframe
generated_complete_datasets <- list()
current_generated_count <- 0 # To track overall progress

cat("\n--- Starting complete dataset generation ---\n")

# -- Outermost loop: The 'reps' loop --
for (rep_idx in seq_len(reps)) {
  # Initialize a list for this specific repetition
  rep_label <- paste0("rep_", rep_idx)
  generated_complete_datasets[[rep_label]] <- list()
  
  # Set a seed for this entire repetition block to ensure all datasets
  # within this rep_idx are generated based on a specific, reproducible
  # random state that differs from other rep_idx's.
  set.seed(1233 + rep_idx) # Use your original seed logic for reps here
  
  cat("\n--- Generating datasets for Repetition", rep_idx, "/", reps, "---\n")
  
  for (n in n_values) {
    for (m in m_values) {
      for (levels_sce in level_scenarios) {
        levels_vec <- if (levels_sce == "few") sample(2:4, m, TRUE)
        else                       sample(5:7, m, TRUE)
        for (dist_sce in dist_scenarios) {
          balanced_flag <- (dist_sce == "balanced")
          marginal_list <- mapply(make_marginal,
                                  k = levels_vec,
                                  balanced = balanced_flag,
                                  SIMPLIFY = FALSE)
          support_list <- lapply(levels_vec, seq_len)
          for (rho in rho_vals) {
            
            # Define a unique identifier for this parameter combination scenario within this rep
            scenario_id <- paste0("n", n, "_m", m, "_levels", levels_sce, "_dist", dist_sce, "_rho", rho)
            
            rho_mat <- matrix(rho, m, m); diag(rho_mat) <- 1
            
            # Generate the dataset
            sim <- corrvar2(
              k_cat    = m, k_cont = 0,
              marginal = marginal_list,
              support  = support_list,
              method   = "Polynomial",
              rho      = rho_mat,
              n        = n
            )
            current_data <- as.data.frame(sim$Y_cat)
            current_data[] <- lapply(current_data, factor) # Ensure columns are factors
            
            # Store the generated dataset in the nested list: rep -> scenario -> dataframe
            generated_complete_datasets[[rep_label]][[scenario_id]] <- current_data
            
            current_generated_count <- current_generated_count + 1
            if (current_generated_count %% (total_datasets_to_generate / 10) == 0 || current_generated_count == total_datasets_to_generate) { # Update progress less frequently
              cat("Overall progress:", current_generated_count, "/", total_datasets_to_generate,
                  " (", round(current_generated_count / total_datasets_to_generate * 100, 2), "%)\n")
            }
          }
        }
      }
    }
  }
}

end_time <- Sys.time()
cat("\n--- All complete datasets generated and stored ---\n")
cat("Total data generation time:", format(end_time - start_time), "\n")

# --- OPTIONAL: Save all generated data to disk ---
# This is highly recommended given the computational cost
saveRDS(generated_complete_datasets, file = "all_simulated_complete_datasets_outer_reps.rds")
cat("\nAll generated datasets saved to 'all_simulated_complete_datasets_outer_reps.rds'\n")

# To load them back in a new R session:
# loaded_datasets <- readRDS("all_simulated_complete_datasets_outer_reps.rds")

# Example of how to access a dataset:
# Access the dataset for rep_1, scenario n=100, m=6, levels="few", dist="balanced", rho=0.1
# example_df <- loaded_datasets[["rep_1"]][["n100_m6_levelsfew_distbalanced_rho0.1"]]
# print(head(example_df))