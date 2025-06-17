# Script to get simulated data


# Load packages
library(tidyverse)
library(netdiffuseR)
library(mvtnorm)

wd = "/n/holylfs05/LABS/nethery_lab/Users/svega/gp_nonsep/"

# Source functions for simulations
source(paste0(wd,"code/simulations/code/0_functions.R"))


for(i in 1:100){
  
  kernel_choice <- "ICM"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                  n_timepoints = 15, distribution = "poisson")
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_pois_data_", kernel_choice, "_seed_", i, ".rds"))
  
  kernel_choice <- "sep"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, distribution = "poisson")
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_pois_data_", kernel_choice, "_seed_", i, ".rds"))
  
  kernel_choice <- "nonsep"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, distribution = "poisson", temporal_lengthscale = .95)
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_pois_data_", kernel_choice, "_seed_", i, ".rds"))
  
  kernel_choice <- "gneiting"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, distribution = "poisson")
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_pois_data_", kernel_choice, "_seed_", i, "_v2.rds"))
  
}

for(i in 1:100){
  
  kernel_choice <- "ICM"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                  n_timepoints = 15, distribution = "normal")
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_normal_data_", kernel_choice, "_seed_", i, ".rds"))
  
  kernel_choice <- "sep"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, distribution = "normal")
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_normal_data_", kernel_choice, "_seed_", i, ".rds"))
  
  kernel_choice <- "nonsep"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, distribution = "normal", temporal_lengthscale = .95)
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_normal_data_", kernel_choice, "_seed_", i, ".rds"))
  
  kernel_choice <- "gneiting"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, distribution = "normal")
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_normal_data_", kernel_choice, "_seed_", i, "_v2.rds"))
  
}

###### temp

for (i in 84:100) {
  
  kernel_choice <- "gneiting_v2"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, distribution = "poisson")
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_pois_data_", kernel_choice, "_seed_", i, ".rds"))
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, distribution = "normal")
  
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_normal_data_", kernel_choice, "_seed_", i, ".rds"))
  
}

for (i in 1:100) {
  
  kernel_choice <- "gneiting_v3"
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, spatial_lengthscale = 8, temporal_lengthscale = 1.75, distribution = "poisson")
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_pois_data_", kernel_choice, "_seed_", i, ".rds"))
  
  # Simulate data using chosen kernel and seed
  simulated_data <- simulate_grid_data(kernel = kernel_choice, seed = i, n_counties = 49, 
                                       n_timepoints = 15, spatial_lengthscale = 8, temporal_lengthscale = 1.75, distribution = "normal")
  
  
  # Save results to files
  saveRDS(simulated_data, file=paste0(wd, "code/simulations/data/simulated_normal_data_", kernel_choice, "_seed_", i, ".rds"))
  
}

