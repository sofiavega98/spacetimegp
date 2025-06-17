# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)
#install.packages("robustbase")
library(robustbase)

#######################################################
## Step 1: Load Simulated Data and Predicted Results ##
#######################################################

# Define working directory 
wd <- "/n/holylfs05/LABS/nethery_lab/Users/svega/gp_nonsep/"

# Define kernel choices and seed values
kernel_choices <- c("ICM","sep","gneiting_v3")
#model_choices <- c("ICM","nonsep","sep") 
model_choices <- c("ICM","nonsep","sep") 
seed_values <- 1:100

# Initialize empty lists to store loaded data
simulated_data_list <- list()
stan_fit_list_sep <- list()
stan_fit_list_ICM <- list()
stan_fit_list_nonsep <- list()
#stan_fit_list_MC <- list()

# Loop over kernel choices and seed values to load all simulated data and results
for (kernel_choice in kernel_choices) {
  for (seed_value in seed_values) {
    
    # Construct file paths for simulated data and Stan model results
    sim_data_file <- paste0(wd, "code/simulations/data/simulated_normal_data_", kernel_choice, "_seed_", seed_value, ".rds")
    stan_fit_file_ICM <- paste0(wd, "code/simulations/results/normal/ICM/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds")
    stan_fit_file_sep <- paste0(wd, "code/simulations/results/normal/sep/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds")
    stan_fit_file_nonsep <- paste0(wd, "code/simulations/results/normal/nonsep/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds")
    #stan_fit_file_MC <- paste0(wd, "code/simulations/results/MC/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds")
    
    
    # Load simulated data and append to list
    if (file.exists(sim_data_file)) {
      sim_data <- readRDS(sim_data_file)
      simulated_data_list[[paste0(kernel_choice, "_", seed_value)]] <- sim_data
    } else {
      warning(paste("Simulated data file not found:", sim_data_file))
    }
    

    # Load Stan fit results for sep model and append to list
    if (file.exists(stan_fit_file_sep)) {
      stan_fit_sep <- readRDS(stan_fit_file_sep)
      stan_fit_list_sep[[paste0(kernel_choice, "_", seed_value)]] <- stan_fit_sep
    } else {
      warning(paste("Stan fit file not found:", stan_fit_file_sep))
      # REMOVE `break` TO PREVENT DROPPING ALL sep DGP ENTRIES
    }
    
    # Load Stan fit results for ICM model and append to list
    if (file.exists(stan_fit_file_ICM)) {
      stan_fit_ICM <- readRDS(stan_fit_file_ICM)
      stan_fit_list_ICM[[paste0(kernel_choice, "_", seed_value)]] <- stan_fit_ICM
    } else {
      warning(paste("Stan fit file not found:", stan_fit_file_ICM))
    }
    
    # Load Stan fit results for nonsep model and append to list
    if (file.exists(stan_fit_file_nonsep)) {
      stan_fit_nonsep <- readRDS(stan_fit_file_nonsep)
      stan_fit_list_nonsep[[paste0(kernel_choice, "_", seed_value)]] <- stan_fit_nonsep
    } else {
      warning(paste("Stan fit file not found:", stan_fit_file_nonsep))
    }
    
    # Load Stan fit results for nonsep model and append to list
    #if (file.exists(stan_fit_file_gneiting)) {
    #  stan_fit_gneiting <- readRDS(stan_fit_file_gneiting)
    #  stan_fit_list_gneiting[[paste0(kernel_choice, "_", seed_value)]] <- stan_fit_gneiting
    #} else {
   #   warning(paste("Stan fit file not found:", stan_fit_file_nonsep))
    #}
    
    
    # Load Stan fit results for nonsep model and append to list
    #if (file.exists(stan_fit_file_MC)) {
    #  stan_fit_MC <- readRDS(stan_fit_file_MC)
    #  stan_fit_list_MC[[paste0(kernel_choice, "_", seed_value)]] <- stan_fit_MC
   # } else {
    #  warning(paste("Stan fit file not found:", stan_fit_file_MC))
    #}
    
  }
}

##############################################################
## Step 2: Calculate Absolute Percent Bias for Treated Counties ##
##############################################################

# Define a function to calculate absolute percent bias
calc_abs_percent_bias <- function(true_values, predicted_values) {
  abs((true_values - predicted_values) / true_values) * 100
}

medianAbsPercBiasY0_year <- function(Y_pred_all, Y0_matrix) {
  # Y_pred_all: 500 by 40 matrix
  # Y0_matrix: 5 by 8 matrix
  
  # Extract true counterfactual values for treated counties
  true_Y0_mean <- colMeans(Y0_matrix, na.rm = TRUE) # Take mean across years for each county
  
  # Extract predictions: Median across 500 MCMC samples per time point
  Y_pred_treated <- colMedians(Y_pred_all, na.rm = TRUE)
  
  # Compute county-level averages across time
  county_means <- colMeans(matrix(Y_pred_treated, ncol = 8, byrow = TRUE), na.rm = TRUE)
  
  # Count the number of counties where the average is exactly 0
  zero_count <- sum(true_Y0_mean == 0, na.rm = TRUE)
  
  # Compute Absolute Percent Bias for each time point, ensuring no division by zero
  absPercBias <- ifelse(true_Y0_mean == 0, NA, abs(((true_Y0_mean - county_means) / true_Y0_mean) * 100))
  
  # Store the median bias of the time points
  medianAbsPercBias <- median(absPercBias, na.rm = TRUE)
  
  # Compute summary statistics
  LB <- quantile(absPercBias, probs = 0.25, na.rm = TRUE)
  UB <- quantile(absPercBias, probs = 0.75, na.rm = TRUE)
  IQR <- IQR(absPercBias, na.rm = TRUE)
  
  return(list(
    median = medianAbsPercBias, 
    LB = LB, 
    UB = UB, 
    IQR = IQR, 
    zero_count = zero_count  # Count of counties with zero averages
  ))
}

calc_MSE <- function(true_values, predicted_values) {
  mean((true_values - predicted_values) ^ 2, na.rm = TRUE)
}

calc_abs_difference <- function(true_values, predicted_values) {
  abs(true_values - predicted_values)
}

calc_coverage <- function(true_values, predicted_intervals) {
  mean(true_values >= predicted_intervals[, 1] & true_values <= predicted_intervals[, 2], na.rm = TRUE)
}



# Initialize an empty list to store bias results
bias_results <- list()
true_Y0 <- list()



# Loop through each simulated dataset and calculate bias for treated counties only
# Loop through each model, kernel, and seed value to calculate bias
time_ind <- c()
for(i in 0:9){
  time_ind <- append(time_ind,seq(8,15) + (i*15))
}

for (model_choice in model_choices) {
  for (kernel_choice in kernel_choices) {
    for (seed_value in seed_values) {
      
      sim_key <- paste0(kernel_choice, "_", seed_value)
      
      if (is.null(simulated_data_list[[sim_key]])) next  
      
      sim_data <- simulated_data_list[[sim_key]] %>% arrange(county_id)
      
      set.seed(501)
      trt_counties <- sample(1:49, 4, replace = FALSE)
      
      treated_data <- sim_data %>%
        filter(county_id %in% trt_counties & time > 7)
      
      true_counterfactuals <- treated_data$mean
      
      if (model_choice == "sep") {
        if (is.null(stan_fit_list_sep[[sim_key]])) next
        post_samp <- rstan::extract(stan_fit_list_sep[[sim_key]])
        Y_pred_all <- post_samp$f_samples
        
        county_index <- unique(treated_data$county_id)
        time_indices <- 8:15
        
        
        reordered <- Y_pred_all[, time_indices, county_index]
        flattened_matrix <- matrix(reordered, nrow = 500, byrow = F)
        
        if (!exists("flattened_matrix") || is.null(flattened_matrix)) {
          warning(paste("Skipping", model_choice, "for", kernel_choice, "seed", seed_value, "- No valid predictions"))
          next  # Skip this model-seed pair but keep others
        }
        
      } else if (model_choice == "ICM") {
        if (is.null(stan_fit_list_ICM[[sim_key]])) next
        post_samp <- rstan::extract(stan_fit_list_ICM[[sim_key]])
        Y_pred_all <- post_samp$f_samples
        
        county_index <- unique(treated_data$county_id)
        time_indices <- 8:15
        
        
        reordered <- Y_pred_all[, time_indices, county_index]
        flattened_matrix <- matrix(reordered, nrow = 500, byrow = F)
        
      } else if (model_choice == "nonsep") {
        if (is.null(stan_fit_list_nonsep[[sim_key]])) next
        post_samp <- rstan::extract(stan_fit_list_nonsep[[sim_key]])
        Y_pred_all <- post_samp$f_samples
        
        county_index <- unique(treated_data$county_id)
        time_indices <- 8:15
        
        
        reordered <- Y_pred_all[, time_indices, county_index]
        flattened_matrix <- matrix(reordered, nrow = 500, byrow = F)
      } else if (model_choice == "gneiting") {
        if (is.null(stan_fit_list_gneiting[[sim_key]])) next
        post_samp <- rstan::extract(stan_fit_list_gneiting[[sim_key]])
        Y_pred_all <- post_samp$f_samples
        
        county_index <- unique(treated_data$county_id)
        time_indices <- 8:15
        
        
        reordered <- Y_pred_all[, time_indices, county_index]
        flattened_matrix <- matrix(reordered, nrow = 500, byrow = F)
        
        }else if (model_choice == "MC") {
        if (is.null(stan_fit_list_MC[[sim_key]])) next
        post_samp <- rstan::extract(stan_fit_list_MC[[sim_key]])
        Y_pred_all <- post_samp$Y_pred
        
        reordered <- Y_pred_all[,time_ind]
      }

      Y_pred_treated <- colMedians(flattened_matrix)
      
      bias <- median(calc_abs_percent_bias(true_counterfactuals, Y_pred_treated))
      
      Y_pred_treated_mat <- matrix(Y_pred_treated, nrow = length(trt_counties), byrow = T)
      true_Y0_mat <- matrix(true_counterfactuals, nrow = length(trt_counties), byrow = T)
      
      att_results <- medianAbsPercBiasY0_year(Y_pred_treated_mat, true_Y0_mat)
      
      att_median <- att_results$median
      att_IQR <- att_results$IQR
      zero_count <- att_results$zero_count  
      
      MSE <- calc_MSE(true_counterfactuals, Y_pred_treated)
      abs_difference <- calc_abs_difference(true_counterfactuals, Y_pred_treated)
      
      lower_bound <- apply(flattened_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
      upper_bound <- apply(flattened_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
      
      coverage <- calc_coverage(true_counterfactuals, cbind(lower_bound, upper_bound))
      
      bias_results[[paste0(model_choice, "_", kernel_choice, "_", seed_value)]] <- data.frame(
        Bias = bias,
        Bias_year = att_median,
        Model = model_choice,
        Kernel = kernel_choice,
        Zero_Count = zero_count,
        rhat = NA,
        MSE = MSE,
        Mean_Abs_Difference = mean(abs_difference, na.rm = TRUE),
        Coverage = coverage
      )
      
      true_Y0[[paste0(model_choice, "_", kernel_choice, "_", seed_value)]] <- data.frame(
        Kernel = kernel_choice,
        Seed = seed_value,
        true_counterfactuals
      )
    }
  }
}


# Convert bias_results from list to a single data frame, ensuring missing models don't remove entire DGPs
bias_results_df <- bind_rows(bias_results, .id = "id") %>%
  complete(Kernel, Model, fill = list(Mean_Bias = NA, Mean_MSE = NA, Mean_Abs_Difference = NA, Mean_Coverage = NA))


# Combine all the individual result data frames into one large data frame
bias_df <- do.call(rbind, bias_results)
true_Y0_df <- do.call(rbind, true_Y0)



##########################################
## Step 4: Summarize Results in a table ##
##########################################

# Load necessary library
library(gt)
library(dplyr)

# Convert bias_results from list to a single data frame
bias_results_df <- bind_rows(bias_results, .id = "id")

# Summarize to calculate the mean of each metric (excluding Bias_year)
summary_table_pois <- bias_results_df %>%
  group_by(Kernel, Model) %>%
  summarise(
    Mean_Bias = mean(Bias, na.rm = TRUE),
    Mean_MSE = mean(MSE, na.rm = TRUE),
    Mean_Abs_Difference = mean(Mean_Abs_Difference, na.rm = TRUE),
    Mean_Coverage = mean(Coverage, na.rm = TRUE),
    .groups = "drop"  # This prevents grouping issues
  ) %>%
  # Rename Kernel values using case_when()
  mutate(
    Kernel = case_when(
      Kernel == "ICM" ~ "Separable ICM",
      Kernel == "sep" ~ "Separable Spatio-Temporal",
      Kernel == "nonsep" ~ "Nonseparable Spatio-Temporal",
      Kernel == "gneiting" ~ "Nonseparable Gneiting",
      #Kernel == "MC" ~ "Spatio-Temporal MC",
      TRUE ~ Kernel  # Keep original name if not matched
    ))%>%
  # Rename Model values using case_when()
  mutate(
    Model = case_when(
      Model == "ICM" ~ "Separable ICM",
      Model == "sep" ~ "Separable Spatio-Temporal",
      Model == "nonsep" ~ "Nonseparable Spatio-Temporal",
      #Model == "MC" ~ "Spatio-Temporal MC",
      TRUE ~ Model  # Keep original name if not matched
    )
  )

# Load necessary libraries
library(knitr)
library(kableExtra)

# Save the LaTeX table to a .tex file
# Prepare the table: drop Mean_Abs_Difference and round values
table_to_save <- summary_table_pois %>%
  select(-Mean_Abs_Difference) %>%       # remove this column
  mutate(across(where(is.numeric), ~ round(.x, 2)))  # round numeric columns

# Save as LaTeX
table_to_save %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    caption = "Summary statistics for posterior bias, mean squared error (MSE), and coverage across different kernel and model combinations. Metrics are averaged over simulated datasets."
  ) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  cat(file = paste0(wd, "/code/simulations/results/normal/summary_table_mean_0529.tex"))


# Save results
saveRDS(summary_table_pois, file=paste0(wd, "/code/simulations/results/normal/summary_table_mean.rds"))

# Create the table using gt
gt_table <- summary_table_pois %>%
  gt() %>%
  tab_header(
    title = "Summary of Model Performance by DGP Kernel"
  ) %>%
  fmt_number(
    columns = c(Mean_Bias, Mean_MSE, Mean_Abs_Difference, Mean_Coverage),
    decimals = 3
  ) %>%
  cols_label(
    Kernel = "DGP Kernel",
    Model = "Model",
    Mean_Bias = "Mean Bias",
    Mean_MSE = "Mean MSE",
    Mean_Abs_Difference = "Mean Abs Difference",
    Mean_Coverage = "Mean Coverage"
  ) %>%
  tab_options(
    table.font.size = px(14),
    heading.title.font.size = px(18),
    heading.title.font.weight = "bold"
  )

# Print the table
gt_table

