# Script: 3_analyze_unit6_weights.R
# Purpose: Analyze and visualize the weights for treated unit 6 across different kernel types
# This script focuses on understanding how the weights for unit 6 change over time
# and compares different kernel specifications
# Date created: 5/15/2025
# Author: Sofia Vega

# Load required packages for data manipulation and visualization
library(tidyverse)
library(ggplot2)
library(patchwork)  # For combining multiple plots

# Set working directory and load functions
wd <- "/Users/sofiavega/Library/Mobile Documents/com~apple~CloudDocs/Harvard/GPs and CI/code/kernel_demo/"
source(paste0(wd,"/code/0_functions.R"))

# Set plot directory for saving output figures
plot_dir <- paste0(wd,"figures/")

# Load simulated data for different kernel types
# Using lengthscale of 0.9 for all kernels
simulated_data_ICM <- readRDS(paste0(wd,"/data/simulated_data_ICM.rds"))
simulated_data_sep <- readRDS(paste0(wd,"/data/simulated_data_sep_sls_0.9.rds"))
simulated_data_nonsep <- readRDS(paste0(wd,"/data/simulated_data_gneiting_sls_0.9.rds"))

# Set analysis parameters
trt_counties <- c(6, 24, 38, 49)  # All treated counties
target_unit <- 6  # The unit we want to analyze in detail
n_counties <- 49  # Total number of counties in the grid
n_timepoints <- 15  # Total number of time points
trt_year <- 8  # Treatment start time

# Compute true weights for each kernel type
# ICM (Intrinsic Coregionalization Model)
true_w_ICM <- compute_true_weights(
  simulated_data = simulated_data_ICM, 
  kernel = "ICM", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .9, 
  t_lengthscale = .9
)

# Separable kernel
true_w_sep <- compute_true_weights(
  simulated_data = simulated_data_sep, 
  kernel = "sep", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .9, 
  t_lengthscale = .9
)

# Gneiting (nonseparable) kernel
true_w_nonsep <- compute_true_weights(
  simulated_data = simulated_data_nonsep, 
  kernel = "gneiting", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .9, 
  t_lengthscale = .9
)

#' Create County Grid Coordinates
#' 
#' @param n_counties Number of counties in the grid
#' @return Data frame with county coordinates and IDs
create_county_grid <- function(n_counties) {
  # Create a regular grid of coordinates
  county_coords <- expand.grid(
    lon = seq(0, 1, length.out = sqrt(n_counties)),
    lat = seq(0, 1, length.out = sqrt(n_counties))
  )
  county_coords <- county_coords[1:n_counties, ]
  county_coords$control_county <- 1:n_counties
  return(county_coords)
}

#' Create Pre/Post Treatment Plot
#' 
#' @param weight_df Data frame containing weight information
#' @param kernel_name Name of the kernel type for plot title
#' @return ggplot object showing weights before and after treatment
create_pre_post_plot <- function(weight_df, kernel_name) {
  # Create county grid for spatial coordinates
  county_coords <- create_county_grid(n_counties)
  
  # Filter for target unit and aggregate weights by period
  aggregated_data <- weight_df %>%
    filter(treated_county == target_unit) %>%
    mutate(period = ifelse(control_time < trt_year, "Pre-treatment", "Post-treatment")) %>%
    group_by(control_county, period) %>%
    summarize(avg_weight = mean(weights, na.rm = TRUE), .groups = "drop") %>%
    mutate(period = factor(period, levels = c("Pre-treatment", "Post-treatment")))
  
  # Merge weights with spatial coordinates
  merged_weights <- merge(
    aggregated_data, 
    county_coords,
    by.x = "control_county",
    by.y = "control_county",
    all.x = TRUE
  )
  
  # Get treated unit coordinates for both periods
  treated_coords <- county_coords %>%
    filter(control_county == target_unit) %>%
    tidyr::expand_grid(
      period = factor(
        c("Pre-treatment", "Post-treatment"),
        levels = c("Pre-treatment", "Post-treatment")
      )
    )
  
  # Create faceted plot showing pre and post treatment weights
  p <- ggplot() +
    geom_tile(
      data = merged_weights,
      aes(x = lon, y = lat, fill = avg_weight)
    ) +
    # Add black border for treated unit
    geom_tile(
      data = treated_coords,
      aes(x = lon, y = lat),
      fill = NA, 
      color = "black", 
      size = 0.5, 
      width = 0.5, 
      height = 0.5
    ) +
    facet_wrap(~ period) +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 0, 
      name = "Avg. Weight"
    ) +
    labs(
      title = paste("Spatial Weights for Unit", target_unit, "-", kernel_name),
      x = "Longitude", 
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(size = 12))
  
  return(p)
}

#' Create Time-Specific Plot (Version 1)
#' 
#' @param weight_df Data frame containing weight information
#' @param kernel_name Name of the kernel type for plot title
#' @return ggplot object showing weights at specific time points
create_time_specific_plot_v1 <- function(weight_df, kernel_name) {
  # Create county grid
  county_coords <- create_county_grid(n_counties)
  
  # Filter for target unit at treated time 11
  filtered_data <- weight_df %>%
    filter(treated_county == target_unit, treated_time == 11)
  
  # Merge with coordinates
  merged_weights <- merge(
    filtered_data, 
    county_coords,
    by.x = "control_county",
    by.y = "control_county",
    all.x = TRUE
  )
  
  # Get treated unit coordinates
  treated_coords <- county_coords %>%
    filter(control_county == target_unit)
  
  # Create faceted plot showing weights at each time point
  p <- ggplot() +
    geom_tile(
      data = merged_weights,
      aes(x = lon, y = lat, fill = weights)
    ) +
    geom_tile(
      data = treated_coords,
      aes(x = lon, y = lat),
      fill = NA, 
      color = "black", 
      size = 0.5, 
      width = 0.165, 
      height = 0.165
    ) +
    facet_wrap(~ control_time) +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 0, 
      name = "Weight"
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12))
  
  return(p)
}

#' Create Time-Specific Plot (Current Version)
#' 
#' @param weight_df Data frame containing weight information
#' @param kernel_name Name of the kernel type for plot title
#' @return ggplot object showing weights at specific time points with treated units grayed out
create_time_specific_plot <- function(weight_df, kernel_name) {
  # Create county grid
  county_coords <- create_county_grid(n_counties)
  
  # Filter for target unit at treated time 11
  filtered_data <- weight_df %>%
    filter(treated_county == target_unit, treated_time == 11)
  
  # Merge with coordinates
  merged_weights <- merge(
    filtered_data, 
    county_coords,
    by.x = "control_county",
    by.y = "control_county",
    all.x = TRUE
  )
  
  # Get treated unit coordinates
  treated_coords <- county_coords %>%
    filter(control_county == target_unit)
  
  # Get unique time points
  time_vals <- sort(unique(merged_weights$control_time))
  
  # Create data frame for grayed out treated units
  gray_coords <- county_coords %>%
    filter(control_county %in% trt_counties & control_county != target_unit)
  
  # Expand gray coordinates for post-treatment time points
  gray_df <- expand.grid(
    control_time = time_vals, 
    control_county = gray_coords$control_county
  ) %>%
    filter(control_time >= trt_year) %>%
    left_join(county_coords, by = "control_county")
  
  # Create faceted plot with grayed out treated units
  p <- ggplot() +
    # Main weight tiles
    geom_tile(
      data = merged_weights,
      aes(x = lon, y = lat, fill = weights)
    ) +
    # Gray out treated units (except target) after treatment
    geom_tile(
      data = gray_df,
      aes(x = lon, y = lat),
      fill = "grey80", 
      color = NA, 
      width = 0.165, 
      height = 0.165
    ) +
    # Black border for treated unit
    geom_tile(
      data = treated_coords,
      aes(x = lon, y = lat),
      fill = NA, 
      color = "black", 
      size = 0.5, 
      width = 0.165, 
      height = 0.165
    ) +
    facet_wrap(~ control_time) +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red",
      midpoint = 0, 
      name = "Weight"
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12))
  
  return(p)
}

# Generate and save plots for each kernel type
# Pre/post treatment plots
pre_post_ICM <- create_pre_post_plot(true_w_ICM$weight_df, "ICM")
pre_post_sep <- create_pre_post_plot(true_w_sep$weight_df, "Separable")
pre_post_nonsep <- create_pre_post_plot(true_w_nonsep$weight_df, "Gneiting")

# Time-specific plots
time_specific_ICM <- create_time_specific_plot(true_w_ICM$weight_df, "ICM")
time_specific_sep <- create_time_specific_plot(true_w_sep$weight_df, "Separable")
time_specific_nonsep <- create_time_specific_plot(true_w_nonsep$weight_df, "Gneiting")

# Save pre/post treatment plots
ggsave(
  paste0(plot_dir, "pre_post_ICM.pdf"), 
  plot = pre_post_ICM,
  width = 10, 
  height = 6
)
ggsave(
  paste0(plot_dir, "pre_post_sep.pdf"), 
  plot = pre_post_sep,
  width = 10, 
  height = 6
)
ggsave(
  paste0(plot_dir, "pre_post_nonsep.pdf"), 
  plot = pre_post_nonsep,
  width = 10, 
  height = 6
)

# Save time-specific plots
ggsave(
  paste0(plot_dir, "time_specific_ICM.pdf"), 
  plot = time_specific_ICM,
  width = 12, 
  height = 8
)
ggsave(
  paste0(plot_dir, "time_specific_sep.pdf"), 
  plot = time_specific_sep,
  width = 12, 
  height = 8
)
ggsave(
  paste0(plot_dir, "time_specific_nonsep.pdf"), 
  plot = time_specific_nonsep,
  width = 12, 
  height = 8
) 
