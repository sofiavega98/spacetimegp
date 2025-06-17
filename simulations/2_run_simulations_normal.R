# Script to run simulations using continuous data

#############################################################################################################
## 0. Read command line arguments and load packages and functions                                          ##
##    Arguments should include:                                                                            ##
##    - kernel_dgp: kernel used for the DGP                                                                ##
##    - model: model used to estimate Y0                                                                   ##
##    - simnum:
#############################################################################################################

# For testing
#kernel = "sep"
#model = "sep"
#simnum = 1

# Read command line arguments 
args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }

start_time <- Sys.time()

# Load packages
library(tidyverse)
library(netdiffuseR)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)
library(mvtnorm)

wd = "/n/holylfs05/LABS/nethery_lab/Users/svega/gp_nonsep/"

# Source functions for simulations
source(paste0(wd,"code/simulations/code/0_functions.R"))

simulated_data <- readRDS(paste0("/n/holylfs05/LABS/nethery_lab/Users/svega/gp_nonsep/code/simulations/data/simulated_normal_data_",kernel,"_seed_",simnum,".rds"))

# Randomly sample 10 counties to be considered treated
set.seed(501)
trt_counties <- sample(1:49, 10, replace = FALSE)  # Without replacement (unique numbers)

load("/n/holylfs05/LABS/nethery_lab/Users/svega/gp_nonsep/code/simulations/data/adj_matrix.RData")


if(model == "ICM"){
  
  run_simulation_and_stan_model_ICM_normal(seed_value = simnum, simulated_data = simulated_data, kernel_choice = kernel, trt_counties)
  
} else if(model == "sep"){

  run_simulation_and_stan_model_sep_normal(seed_value = simnum, simulated_data = simulated_data, kernel_choice = kernel, trt_counties)

} else if(model == "nonsep"){
  
  run_simulation_and_stan_model_nonsep_normal(seed_value = simnum, simulated_data = simulated_data, kernel_choice = kernel, trt_counties)
  
}else if(model == "nonsep_v2"){
  
  run_simulation_and_stan_model_nonsep_v2_normal(seed_value = simnum, simulated_data = simulated_data, kernel_choice = kernel, trt_counties)
  
} else if(model == "MC"){
  
  run_simulation_and_stan_model_MC_normal(seed_value = simnum, simulated_data = simulated_data, kernel_choice = kernel, trt_counties, adj_matrix)
  
}else if(model == "sep_NNGP"){
  
  run_simulation_and_stan_model_sep_nngp_normal(seed_value = simnum, simulated_data = simulated_data, kernel_choice = kernel, trt_counties, M = 10)
  
}else if(model == "ICM_NNGP"){
  
  run_simulation_and_stan_model_ICM_nngp_normal(seed_value = simnum, simulated_data = simulated_data, kernel_choice = kernel, trt_counties, M = 10)
  
}else if(model == "nonsep_NNGP"){
  
  run_simulation_and_stan_model_nonsep_nngp_normal(seed_value = simnum, simulated_data = simulated_data, kernel_choice = kernel, trt_counties, M = 10)
  
}

end_time <- Sys.time()
run_time <- end_time - start_time

# Calculate run time in minutes
run_time_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))

# Create a data frame
run_time_df <- tibble(
  kernel = kernel,
  model = model,
  simnum = simnum,
  run_time_minutes = run_time_minutes
)

# Define path
log_path <- paste0(wd, "code/simulations/runtime_logs/run_times_normal.csv")

# Append to CSV (create if it doesn’t exist)
if (!file.exists(log_path)) {
  write_csv(run_time_df, log_path)
} else {
  write_csv(run_time_df, log_path, append = TRUE)
}

