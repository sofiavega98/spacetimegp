# Script: 1_create_sim_data.R
# Purpose: Generate simulated spatiotemporal data using different kernel types and parameters
# This script creates the initial dataset for the kernel demonstration analysis

# Set working directory to the kernel demo folder
wd <- "~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/GPs and CI/code/kernel_demo/"

# Source the functions file containing simulation and analysis functions
source(paste0(wd,"/code/0_functions.R"))

# Part 1: Generate ICM (Intrinsic Coregionalization Model) Data
# This creates a baseline dataset using the ICM kernel
kernel_choice <- "ICM"

# Simulate data using ICM kernel with specified parameters:
# - 49 counties in a 7x7 grid
# - 15 time points
# - Temporal lengthscale of 0.9
# - Normal distribution for observations
simulated_data <- simulate_grid_data(
  kernel = kernel_choice, 
  seed = 1, 
  n_counties = 49, 
  n_timepoints = 15, 
  temporal_lengthscale = .9, 
  distribution = "normal"
)

# Save the ICM simulated data to an RDS file
saveRDS(simulated_data, file=paste0(wd, "/data/simulated_data_", kernel_choice, ".rds"))

# Part 2: Generate Data with Different Kernel Types and Lengthscales
# Define sets of spatial lengthscales to test
# These values control the spatial correlation range in the simulated data
sep_scales    <- c(0.3, 0.7, 0.9)  # Lengthscales for separable kernel
nonsep_scales <- c(0.3, 0.7, 0.9)  # Lengthscales for nonseparable kernels

# Generate data using separable kernel with different spatial lengthscales
# The separable kernel assumes independence between space and time
for (l in sep_scales) {
  cat("Separable: spatial_lengthscale =", l, "\n")
  sim <- simulate_grid_data(
    kernel               = "sep",
    seed                 = 1,
    n_counties           = 49,
    n_timepoints         = 15,
    spatial_lengthscale  = l,
    temporal_lengthscale = 0.9,
    distribution         = "normal"
  )
  # Save each simulation with a unique filename indicating the lengthscale
  saveRDS(
    sim,
    file = file.path(
      wd,
      "data",
      sprintf("simulated_data_sep_sls_%0.1f.rds", l)
    )
  )
}

# Generate data using nonseparable kernels (nonsep and Gneiting) with different lengthscales
# These kernels allow for space-time interaction in the correlation structure
for (l in nonsep_scales) {
  for (kernel_choice in c("nonsep", "gneiting")) {
    cat(kernel_choice, ": spatial_lengthscale =", l, "\n")
    sim <- simulate_grid_data(
      kernel               = kernel_choice,
      seed                 = 1,
      n_counties           = 49,
      n_timepoints         = 15,
      spatial_lengthscale  = l,
      temporal_lengthscale = 0.9,
      distribution         = "normal"
    )
    # Save each simulation with a unique filename indicating the kernel type and lengthscale
    saveRDS(
      sim,
      file = file.path(
        wd,
        "data",
        sprintf("simulated_data_%s_sls_%0.1f.rds", kernel_choice, l)
      )
    )
  }
}

