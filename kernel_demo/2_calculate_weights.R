# Script: 2_calculate_weights.R
# Purpose: Calculate and visualize treated unit weights for different kernel types and parameters
# This script analyzes the spatial distribution of weights across different kernel specifications
# Date created: 5/15/25
# Author: Sofia Vega

##########################
## Load and Set Up Data ##
##########################

# Load required libraries for data manipulation and visualization
library(tidyverse)
library(ggplot2)

# Set working directory and source functions
wd <- "~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/GPs and CI/code/kernel_demo/"
source(paste0(wd,"/code/0_functions.R"))

# Set path for saving output plots
plot_dir <- paste0(wd,"figures/")

# Load simulated data for different kernel types and lengthscales
# ICM (Intrinsic Coregionalization Model) data
simulated_data_ICM <- readRDS(paste0(wd,"data/simulated_data_ICM.rds"))

# Separable kernel data with different spatial lengthscales
simulated_data_sep_.3 <- readRDS(paste0(wd,"data/simulated_data_sep_sls_0.3.rds"))
simulated_data_sep_.7 <- readRDS(paste0(wd,"data/simulated_data_sep_sls_0.7.rds"))
simulated_data_sep_.9 <- readRDS(paste0(wd,"data/simulated_data_sep_sls_0.9.rds"))

# Gneiting (nonseparable) kernel data with different spatial lengthscales
simulated_data_gneiting_.3 <- readRDS(paste0(wd,"data/simulated_data_gneiting_sls_0.3.rds"))
simulated_data_gneiting_.7 <- readRDS(paste0(wd,"data/simulated_data_gneiting_sls_0.7.rds"))
simulated_data_gneiting_.9 <- readRDS(paste0(wd,"data/simulated_data_gneiting_sls_0.9.rds"))

# Set up treatment assignment
# Randomly select 4 counties for treatment, starting at time point 8
set.seed(501)
n_counties <- 49
n_timepoints <- 15
trt_year <- 8
trt_counties <- sample(1:n_counties, 4, replace = FALSE)  # Without replacement (unique numbers)

#################################################
## Obtain true weights from the simulated data ##
#################################################

# Calculate true weights for ICM kernel
true_w_ICM <- compute_true_weights(
  simulated_data = simulated_data_ICM, 
  kernel = "ICM", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .3, 
  t_lengthscale = .9
)

# Calculate true weights for separable kernel with different lengthscales
true_w_sep_.3 <- compute_true_weights(
  simulated_data = simulated_data_sep_.3, 
  kernel = "sep", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .3, 
  t_lengthscale = .9
)

true_w_sep_.7 <- compute_true_weights(
  simulated_data = simulated_data_sep_.7, 
  kernel = "sep", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .7, 
  t_lengthscale = .9
)

true_w_sep_.9 <- compute_true_weights(
  simulated_data = simulated_data_sep_.9, 
  kernel = "sep", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .9, 
  t_lengthscale = .9
)

# Calculate true weights for Gneiting kernel with different lengthscales
true_w_gneiting_.3 <- compute_true_weights(
  simulated_data = simulated_data_gneiting_.3, 
  kernel = "gneiting", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .3, 
  t_lengthscale = .9
)

true_w_gneiting_.7 <- compute_true_weights(
  simulated_data = simulated_data_gneiting_.7, 
  kernel = "gneiting", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .7, 
  t_lengthscale = .9
)

true_w_gneiting_.9 <- compute_true_weights(
  simulated_data = simulated_data_gneiting_.9, 
  kernel = "gneiting", 
  trt_counties = trt_counties, 
  seed_for_icm = 1, 
  nugget = 1e-6,
  s_lengthscale = .9, 
  t_lengthscale = .9
)

##############################################################
## Plot the true spatial weights averaged over time periods ##
##############################################################

# Aggregate weights for separable kernel across different lengthscales
# This computes the mean weight for each treated-control county pair
aggregated_data_.3 <- true_w_sep_.3$weight_df %>%
  group_by(treated_county, control_county) %>%
  summarize(avg_weight = mean(weights, na.rm = TRUE), .groups = 'drop')
aggregated_data_.7 <- true_w_sep_.7$weight_df %>%
  group_by(treated_county, control_county) %>%
  summarize(avg_weight = mean(weights, na.rm = TRUE), .groups = 'drop')
aggregated_data_.9 <- true_w_sep_.9$weight_df %>%
  group_by(treated_county, control_county) %>%
  summarize(avg_weight = mean(weights, na.rm = TRUE), .groups = 'drop')

# Combine all separable kernel weights for consistent color scaling
sep_weights <- bind_rows(
  aggregated_data_.3,
  aggregated_data_.7,
  aggregated_data_.9
)
sep_global_min <- min(sep_weights$avg_weight, na.rm = TRUE)
sep_global_max <- max(sep_weights$avg_weight, na.rm = TRUE)

# Aggregate weights for Gneiting kernel across different lengthscales
aggregated_data_.3 <- true_w_gneiting_.3$weight_df %>%
  group_by(treated_county, control_county) %>%
  summarize(avg_weight = mean(weights, na.rm = TRUE), .groups = 'drop')
aggregated_data_.7 <- true_w_gneiting_.7$weight_df %>%
  group_by(treated_county, control_county) %>%
  summarize(avg_weight = mean(weights, na.rm = TRUE), .groups = 'drop')
aggregated_data_.9 <- true_w_gneiting_.9$weight_df %>%
  group_by(treated_county, control_county) %>%
  summarize(avg_weight = mean(weights, na.rm = TRUE), .groups = 'drop')

# Combine all Gneiting kernel weights for consistent color scaling
nonsep_weights <- bind_rows(
  aggregated_data_.3,
  aggregated_data_.7,
  aggregated_data_.9
)
nonsep_global_min <- min(nonsep_weights$avg_weight, na.rm = TRUE)
nonsep_global_max <- max(nonsep_weights$avg_weight, na.rm = TRUE)

# Generate heatmap plots for each kernel type and lengthscale
# Using consistent color scales within each kernel type
spatial_w_heatmap_ICM <- spatial_w_heatmap(true_w_ICM$weight_df)
spatial_w_heatmap_sep_.3 <- spatial_w_heatmap(true_w_sep_.3$weight_df, fill_limits = c(sep_global_min, sep_global_max))
spatial_w_heatmap_gneiting_.3 <- spatial_w_heatmap(true_w_gneiting_.3$weight_df, fill_limits = c(nonsep_global_min, nonsep_global_max))

spatial_w_heatmap_sep_.7 <- spatial_w_heatmap(true_w_sep_.7$weight_df, fill_limits = c(sep_global_min, sep_global_max))
spatial_w_heatmap_gneiting_.7 <- spatial_w_heatmap(true_w_gneiting_.7$weight_df, fill_limits = c(nonsep_global_min, nonsep_global_max))

spatial_w_heatmap_sep_.9 <- spatial_w_heatmap(true_w_sep_.9$weight_df, fill_limits = c(sep_global_min, sep_global_max))
spatial_w_heatmap_gneiting_.9 <- spatial_w_heatmap(true_w_gneiting_.9$weight_df, fill_limits = c(nonsep_global_min, nonsep_global_max))

# Create a combined plot layout using patchwork
library(patchwork)
spatial_kernels_plot <- (
  (spatial_w_heatmap_sep_.3 | spatial_w_heatmap_gneiting_.3) /
    (spatial_w_heatmap_sep_.7 | spatial_w_heatmap_gneiting_.7) /
    (spatial_w_heatmap_sep_.9 | spatial_w_heatmap_gneiting_.9)
)

# Create individual rows for each lengthscale
row_03 <- spatial_w_heatmap_sep_.3 + spatial_w_heatmap_gneiting_.3
row_07 <- spatial_w_heatmap_sep_.7 + spatial_w_heatmap_gneiting_.7
row_09 <- spatial_w_heatmap_sep_.9 + spatial_w_heatmap_gneiting_.9

# Create title plots for the columns
separable_title_plot <- ggplot() + 
  ggtitle("Separable RBF-RBF DGP") + 
  theme_void()

nonseparable_title_plot <- ggplot() + 
  ggtitle("Nonseparable Gneiting DGP") + 
  theme_void()

# Create top row with column labels
top_labels <- plot_spacer() + separable_title_plot + nonseparable_title_plot
top_labels <- top_labels + plot_layout(ncol = 3, widths = c(1.2, 1.7, 1))

# Create row labels for each lengthscale
row_label_03 <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "Spatial lengthscale = 0.3", size = 5, angle = 90) +
  theme_void()

row_label_07 <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "Spatial lengthscale = 0.7", size = 5, angle = 90) +
  theme_void()

row_label_09 <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "Spatial lengthscale = 0.9", size = 5, angle = 90) +
  theme_void()

# Assemble final plot with all components
final_plot <- (
  top_labels /
    (row_label_03 + spatial_w_heatmap_sep_.3 + spatial_w_heatmap_gneiting_.3) /
    (row_label_07 + spatial_w_heatmap_sep_.7 + spatial_w_heatmap_gneiting_.7) /
    (row_label_09 + spatial_w_heatmap_sep_.9 + spatial_w_heatmap_gneiting_.9)
) +
  plot_layout(
    widths = c(-1, 10, 10),  # Keep the left column narrower for labels
    heights = c(0.1, 1, 1, 1),
    guides = "keep"
  )

# Display the final plot
print(final_plot)

# Save individual plots
ggsave(paste0(plot_dir, "spatial_kernels_plot_0616.pdf"), plot = final_plot,
       width = 12, height = 12)

# Save individual heatmap plots
ggsave(paste0(plot_dir, "spatial_w_heatmap_ICM.pdf"), plot = spatial_w_heatmap_ICM,
       width = 10, height = 6)
ggsave(paste0(plot_dir, "spatial_w_heatmap_sep_.3.pdf"), plot = spatial_w_heatmap_sep_.3,
       width = 10, height = 6)
ggsave(paste0(plot_dir, "spatial_w_heatmap_nonsep_.3.pdf"), plot = spatial_w_heatmap_gneiting_.3,
       width = 10, height = 6)

ggsave(paste0(plot_dir, "spatial_w_heatmap_sep_.7.pdf"), plot = spatial_w_heatmap_sep_.7,
       width = 10, height = 6)
ggsave(paste0(plot_dir, "spatial_w_heatmap_nonsep_.7.pdf"), plot = spatial_w_heatmap_gneiting_.7,
       width = 10, height = 6)

ggsave(paste0(plot_dir, "spatial_w_heatmap_sep_.9.pdf"), plot = spatial_w_heatmap_sep_.9,
       width = 10, height = 6)
ggsave(paste0(plot_dir, "spatial_w_heatmap_nonsep_.9.pdf"), plot = spatial_w_heatmap_gneiting_.9,
       width = 10, height = 6)


