# Script to plot means of the simulated dataset for paper

# Load required libraries
library(MASS)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)


wd = "/n/holylfs05/LABS/nethery_lab/Users/svega/gp_nonsep/"

# Source functions for simulations
source(paste0(wd,"code/simulations/code/0_functions.R"))

# Simulate data for each kernel
data_nonsep <- simulate_grid_data("gneiting_v2", seed = 1, n_counties = 49, n_timepoints = 15, distribution = "poisson",
                                  spatial_lengthscale = 8, temporal_lengthscale = 1.75)
data_sep <- simulate_grid_data("sep", seed = 1, n_counties = 49, n_timepoints = 15, distribution = "poisson")
data_icm <- simulate_grid_data("ICM", seed = 1, n_counties = 49, n_timepoints = 15, distribution = "poisson")

# Filter for three consecutive time points
selected_timepoints <- c(1,5,10) # Adjust these indices if needed
data_nonsep_filtered <- data_nonsep %>% filter(time %in% selected_timepoints)
data_sep_filtered <- data_sep %>% filter(time %in% selected_timepoints)
data_icm_filtered <- data_icm %>% filter(time %in% selected_timepoints)

# Function to create a count matrix for a single time point
create_count_matrix <- function(data, timepoint) {
  subset_data <- data %>% filter(time == timepoint)
  matrix(subset_data$gp_sample, nrow = sqrt(nrow(subset_data)), byrow = FALSE) # Ensure correct dimensions
}

# Generate count matrices for each kernel and time point
matrices_nonsep <- lapply(selected_timepoints, function(t) create_count_matrix(data_nonsep_filtered, t))
matrices_sep <- lapply(selected_timepoints, function(t) create_count_matrix(data_sep_filtered, t))
matrices_icm <- lapply(selected_timepoints, function(t) create_count_matrix(data_icm_filtered, t))

# Function to plot a single heatmap with consistent color scale and axis ticks
plot_heatmap_from_matrix <- function(matrix_data, fill_limits, grid_coords) {
  melted_data <- melt(matrix_data)
  
  # Update melted data with actual grid coordinates
  melted_data$Var1 <- rep(grid_coords$x, each = length(grid_coords$y))
  melted_data$Var2 <- rep(grid_coords$y, times = length(grid_coords$x))
  
  # Set tick marks every 0.2
  ticks <- seq(0, 1, by = 0.2)
  
  ggplot(melted_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 5)),
      axis.title.y = element_text(size = 10, margin = margin(r = 5)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line()
    ) +
    scale_x_continuous(
      name = "Longitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    scale_y_continuous(
      name = "Latitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    labs(fill = "GP Sample")
}


# Determine consistent color scales for each kernel
limits_nonsep <- range(unlist(matrices_nonsep))
limits_sep <- range(unlist(matrices_sep))
limits_icm <- range(unlist(matrices_icm))

# Generate spatial grid coordinates (assuming n_counties is a perfect square)
n_grid <- sqrt(49) # Adjust based on n_counties
grid_coords <- list(
  x = seq(0, 1, length.out = n_grid),
  y = seq(0, 1, length.out = n_grid)
)

# Create heatmaps for each kernel and time point
plots_nonsep <- lapply(matrices_nonsep, function(mat) plot_heatmap_from_matrix(mat, limits_nonsep, grid_coords))
plots_sep <- lapply(matrices_sep, function(mat) plot_heatmap_from_matrix(mat, limits_sep, grid_coords))
plots_icm <- lapply(matrices_icm, function(mat) plot_heatmap_from_matrix(mat, limits_icm, grid_coords))

# Combine heatmaps into rows for each kernel
row_nonsep <- plot_grid(plotlist = plots_nonsep, ncol = length(selected_timepoints), align = "h")
row_sep <- plot_grid(plotlist = plots_sep, ncol = length(selected_timepoints), align = "h")
row_icm <- plot_grid(plotlist = plots_icm, ncol = length(selected_timepoints), align = "h")

# Add row labels for kernel types
final_plot <- plot_grid(
  ggdraw() + draw_label("Separable \n ICM-RBF DGP", angle = 90),
  row_icm,
  ggdraw() + draw_label("Separable \n RBF-RBF DGP", angle = 90),
  row_sep,
  ggdraw() + draw_label("Nonseparable \n Gneiting DGP", angle = 90),
  row_nonsep,
  
  
  ncol = 2,
  rel_widths = c(0.1, 1)
)


# Add column labels for time points
time_labels <- plot_grid(
  ggdraw() + draw_label("Time 1"),
  ggdraw() + draw_label("Time 5"),
  ggdraw() + draw_label("Time 10"),
  
  ncol = length(selected_timepoints),
  rel_widths = c(1 / length(selected_timepoints), 
                 rep(1 / length(selected_timepoints), length(selected_timepoints) - 1))
)


final_plot_with_labels <- plot_grid(
  plot_grid(NULL, time_labels, NULL, ncol = 3, rel_widths = c(.06, .8, 0.03)),  # Time labels only over heatmaps
  final_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)


# Save plot
pdf(paste0(wd,"code/simulations/figures/GP_presentation_sim_gp_new_nonsep_v3_0617.pdf"),
    paper = 'special', width = 10, height = 7)


# Display the final plot
print(final_plot_with_labels)

dev.off()

# Function to create a count matrix for a single time point
create_count_matrix <- function(data, timepoint) {
  subset_data <- data %>% filter(time == timepoint)
  matrix(subset_data$mean, nrow = sqrt(nrow(subset_data)), byrow = FALSE) # Ensure correct dimensions
}

# Generate count matrices for each kernel and time point
matrices_nonsep <- lapply(selected_timepoints, function(t) create_count_matrix(data_nonsep_filtered, t))
matrices_sep <- lapply(selected_timepoints, function(t) create_count_matrix(data_sep_filtered, t))
matrices_icm <- lapply(selected_timepoints, function(t) create_count_matrix(data_icm_filtered, t))

# Function to plot a single heatmap with consistent color scale and axis ticks
plot_heatmap_from_matrix <- function(matrix_data, fill_limits, grid_coords) {
  melted_data <- melt(matrix_data)
  
  # Update melted data with actual grid coordinates
  melted_data$Var1 <- rep(grid_coords$x, each = length(grid_coords$y))
  melted_data$Var2 <- rep(grid_coords$y, times = length(grid_coords$x))
  
  # Set tick marks every 0.2
  ticks <- seq(0, 1, by = 0.2)
  
  ggplot(melted_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 5)),
      axis.title.y = element_text(size = 10, margin = margin(r = 5)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line()
    ) +
    scale_x_continuous(
      name = "Longitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    scale_y_continuous(
      name = "Latitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    labs(fill = "Poisson \nMean")
}


# Determine consistent color scales for each kernel
limits_nonsep <- range(unlist(matrices_nonsep))
limits_sep <- range(unlist(matrices_sep))
limits_icm <- range(unlist(matrices_icm))

# Generate spatial grid coordinates (assuming n_counties is a perfect square)
n_grid <- sqrt(49) # Adjust based on n_counties
grid_coords <- list(
  x = seq(0, 1, length.out = n_grid),
  y = seq(0, 1, length.out = n_grid)
)

# Create heatmaps for each kernel and time point
plots_nonsep <- lapply(matrices_nonsep, function(mat) plot_heatmap_from_matrix(mat, limits_nonsep, grid_coords))
plots_sep <- lapply(matrices_sep, function(mat) plot_heatmap_from_matrix(mat, limits_sep, grid_coords))
plots_icm <- lapply(matrices_icm, function(mat) plot_heatmap_from_matrix(mat, limits_icm, grid_coords))

# Combine heatmaps into rows for each kernel
row_nonsep <- plot_grid(plotlist = plots_nonsep, ncol = length(selected_timepoints), align = "h")
row_sep <- plot_grid(plotlist = plots_sep, ncol = length(selected_timepoints), align = "h")
row_icm <- plot_grid(plotlist = plots_icm, ncol = length(selected_timepoints), align = "h")

# Add row labels for kernel types
final_plot <- plot_grid(
  ggdraw() + draw_label("Non-Separable \n Gneiting DGP", angle = 90),
  row_nonsep,
  ggdraw() + draw_label("Separable \n RBF-RBF DGP", angle = 90),
  row_sep,
  ggdraw() + draw_label("Separable \n ICM-RBF DGP", angle = 90),
  row_icm,
  ncol = 2,
  rel_widths = c(0.1, 1)
)


# Add column labels for time points
time_labels <- plot_grid(
  ggdraw() + draw_label("Time 1"),
  ggdraw() + draw_label("Time 2"),
  ggdraw() + draw_label("Time 3"),
  
  ncol = length(selected_timepoints),
  rel_widths = c(1 / length(selected_timepoints), 
                 rep(1 / length(selected_timepoints), length(selected_timepoints) - 1))
)


final_plot_with_labels <- plot_grid(
  plot_grid(NULL, time_labels, NULL, ncol = 3, rel_widths = c(.06, .8, 0.03)),  # Time labels only over heatmaps
  final_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)


# Save plot
pdf(paste0(wd,"code/simulations/figures/GP_presentation_sim_mean_poisson_v2.pdf"),
    paper = 'special', width = 10, height = 7)


# Display the final plot
print(final_plot_with_labels)

dev.off()

# Function to create a count matrix for a single time point
create_count_matrix <- function(data, timepoint) {
  subset_data <- data %>% filter(time == timepoint)
  matrix(subset_data$count, nrow = sqrt(nrow(subset_data)), byrow = FALSE) # Ensure correct dimensions
}

# Generate count matrices for each kernel and time point
matrices_nonsep <- lapply(selected_timepoints, function(t) create_count_matrix(data_nonsep_filtered, t))
matrices_sep <- lapply(selected_timepoints, function(t) create_count_matrix(data_sep_filtered, t))
matrices_icm <- lapply(selected_timepoints, function(t) create_count_matrix(data_icm_filtered, t))

# Function to plot a single heatmap with consistent color scale and axis ticks
plot_heatmap_from_matrix <- function(matrix_data, fill_limits, grid_coords) {
  melted_data <- melt(matrix_data)
  
  # Update melted data with actual grid coordinates
  melted_data$Var1 <- rep(grid_coords$x, each = length(grid_coords$y))
  melted_data$Var2 <- rep(grid_coords$y, times = length(grid_coords$x))
  
  # Set tick marks every 0.2
  ticks <- seq(0, 1, by = 0.2)
  
  ggplot(melted_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 5)),
      axis.title.y = element_text(size = 10, margin = margin(r = 5)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line()
    ) +
    scale_x_continuous(
      name = "Longitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    scale_y_continuous(
      name = "Latitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    labs(fill = "Count")
}


# Determine consistent color scales for each kernel
limits_nonsep <- range(unlist(matrices_nonsep))
limits_sep <- range(unlist(matrices_sep))
limits_icm <- range(unlist(matrices_icm))

# Generate spatial grid coordinates (assuming n_counties is a perfect square)
n_grid <- sqrt(49) # Adjust based on n_counties
grid_coords <- list(
  x = seq(0, 1, length.out = n_grid),
  y = seq(0, 1, length.out = n_grid)
)

# Create heatmaps for each kernel and time point
plots_nonsep <- lapply(matrices_nonsep, function(mat) plot_heatmap_from_matrix(mat, limits_nonsep, grid_coords))
plots_sep <- lapply(matrices_sep, function(mat) plot_heatmap_from_matrix(mat, limits_sep, grid_coords))
plots_icm <- lapply(matrices_icm, function(mat) plot_heatmap_from_matrix(mat, limits_icm, grid_coords))

# Combine heatmaps into rows for each kernel
row_nonsep <- plot_grid(plotlist = plots_nonsep, ncol = length(selected_timepoints), align = "h")
row_sep <- plot_grid(plotlist = plots_sep, ncol = length(selected_timepoints), align = "h")
row_icm <- plot_grid(plotlist = plots_icm, ncol = length(selected_timepoints), align = "h")

# Add row labels for kernel types
final_plot <- plot_grid(
  ggdraw() + draw_label("Non-Separable \n Gneiting DGP", angle = 90),
  row_nonsep,
  ggdraw() + draw_label("Separable \n RBF-RBF DGP", angle = 90),
  row_sep,
  ggdraw() + draw_label("Separable \n ICM-RBF DGP", angle = 90),
  row_icm,
  ncol = 2,
  rel_widths = c(0.1, 1)
)


# Add column labels for time points
time_labels <- plot_grid(
  ggdraw() + draw_label("Time 1"),
  ggdraw() + draw_label("Time 2"),
  ggdraw() + draw_label("Time 3"),
  
  ncol = length(selected_timepoints),
  rel_widths = c(1 / length(selected_timepoints), 
                 rep(1 / length(selected_timepoints), length(selected_timepoints) - 1))
)


final_plot_with_labels <- plot_grid(
  plot_grid(NULL, time_labels, NULL, ncol = 3, rel_widths = c(.06, .8, 0.03)),  # Time labels only over heatmaps
  final_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)


# Save plot
pdf(paste0(wd,"code/simulations/figures/GP_presentation_sim_count_poisson.pdf"),
    paper = 'special', width = 10, height = 7)


# Display the final plot
print(final_plot_with_labels)

dev.off()

############
## Normal ##
############

# Simulate data for each kernel
data_nonsep <- simulate_grid_data("nonsep", seed = 1, n_counties = 49, n_timepoints = 15, distribution = "normal")
data_sep <- simulate_grid_data("sep", seed = 1, n_counties = 49, n_timepoints = 15, distribution = "normal")
data_icm <- simulate_grid_data("ICM", seed = 1, n_counties = 49, n_timepoints = 15, distribution = "normal")

# Filter for three consecutive time points
selected_timepoints <- c(1:3) # Adjust these indices if needed
data_nonsep_filtered <- data_nonsep %>% filter(time %in% selected_timepoints)
data_sep_filtered <- data_sep %>% filter(time %in% selected_timepoints)
data_icm_filtered <- data_icm %>% filter(time %in% selected_timepoints)

# Function to create a count matrix for a single time point
create_count_matrix <- function(data, timepoint) {
  subset_data <- data %>% filter(time == timepoint)
  matrix(subset_data$mean, nrow = sqrt(nrow(subset_data)), byrow = FALSE) # Ensure correct dimensions
}

# Generate count matrices for each kernel and time point
matrices_nonsep <- lapply(selected_timepoints, function(t) create_count_matrix(data_nonsep_filtered, t))
matrices_sep <- lapply(selected_timepoints, function(t) create_count_matrix(data_sep_filtered, t))
matrices_icm <- lapply(selected_timepoints, function(t) create_count_matrix(data_icm_filtered, t))

# Function to plot a single heatmap with consistent color scale and axis ticks
plot_heatmap_from_matrix <- function(matrix_data, fill_limits, grid_coords) {
  melted_data <- melt(matrix_data)
  
  # Update melted data with actual grid coordinates
  melted_data$Var1 <- rep(grid_coords$x, each = length(grid_coords$y))
  melted_data$Var2 <- rep(grid_coords$y, times = length(grid_coords$x))
  
  # Set tick marks every 0.2
  ticks <- seq(0, 1, by = 0.2)
  
  ggplot(melted_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 5)),
      axis.title.y = element_text(size = 10, margin = margin(r = 5)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line()
    ) +
    scale_x_continuous(
      name = "Longitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    scale_y_continuous(
      name = "Latitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    labs(fill = "Normal \nMean")
}


# Determine consistent color scales for each kernel
limits_nonsep <- range(unlist(matrices_nonsep))
limits_sep <- range(unlist(matrices_sep))
limits_icm <- range(unlist(matrices_icm))

# Generate spatial grid coordinates (assuming n_counties is a perfect square)
n_grid <- sqrt(49) # Adjust based on n_counties
grid_coords <- list(
  x = seq(0, 1, length.out = n_grid),
  y = seq(0, 1, length.out = n_grid)
)

# Create heatmaps for each kernel and time point
plots_nonsep <- lapply(matrices_nonsep, function(mat) plot_heatmap_from_matrix(mat, limits_nonsep, grid_coords))
plots_sep <- lapply(matrices_sep, function(mat) plot_heatmap_from_matrix(mat, limits_sep, grid_coords))
plots_icm <- lapply(matrices_icm, function(mat) plot_heatmap_from_matrix(mat, limits_icm, grid_coords))

# Combine heatmaps into rows for each kernel
row_nonsep <- plot_grid(plotlist = plots_nonsep, ncol = length(selected_timepoints), align = "h")
row_sep <- plot_grid(plotlist = plots_sep, ncol = length(selected_timepoints), align = "h")
row_icm <- plot_grid(plotlist = plots_icm, ncol = length(selected_timepoints), align = "h")

# Add row labels for kernel types
final_plot <- plot_grid(
  ggdraw() + draw_label("Non-Separable \n Gneiting DGP", angle = 90),
  row_nonsep,
  ggdraw() + draw_label("Separable \n RBF-RBF DGP", angle = 90),
  row_sep,
  ggdraw() + draw_label("Separable \n ICM-RBF DGP", angle = 90),
  row_icm,
  ncol = 2,
  rel_widths = c(0.1, 1)
)


# Add column labels for time points
time_labels <- plot_grid(
  ggdraw() + draw_label("Time 1"),
  ggdraw() + draw_label("Time 2"),
  ggdraw() + draw_label("Time 3"),
  
  ncol = length(selected_timepoints),
  rel_widths = c(1 / length(selected_timepoints), 
                 rep(1 / length(selected_timepoints), length(selected_timepoints) - 1))
)


final_plot_with_labels <- plot_grid(
  plot_grid(NULL, time_labels, NULL, ncol = 3, rel_widths = c(.06, .8, 0.03)),  # Time labels only over heatmaps
  final_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)


# Save plot
pdf(paste0(wd,"code/simulations/figures/GP_presentation_sim_mean_normal.pdf"),
    paper = 'special', width = 10, height = 7)


# Display the final plot
print(final_plot_with_labels)

dev.off()

# Function to create a count matrix for a single time point
create_count_matrix <- function(data, timepoint) {
  subset_data <- data %>% filter(time == timepoint)
  matrix(subset_data$count, nrow = sqrt(nrow(subset_data)), byrow = FALSE) # Ensure correct dimensions
}

# Generate count matrices for each kernel and time point
matrices_nonsep <- lapply(selected_timepoints, function(t) create_count_matrix(data_nonsep_filtered, t))
matrices_sep <- lapply(selected_timepoints, function(t) create_count_matrix(data_sep_filtered, t))
matrices_icm <- lapply(selected_timepoints, function(t) create_count_matrix(data_icm_filtered, t))

# Function to plot a single heatmap with consistent color scale and axis ticks
plot_heatmap_from_matrix <- function(matrix_data, fill_limits, grid_coords) {
  melted_data <- melt(matrix_data)
  
  # Update melted data with actual grid coordinates
  melted_data$Var1 <- rep(grid_coords$x, each = length(grid_coords$y))
  melted_data$Var2 <- rep(grid_coords$y, times = length(grid_coords$x))
  
  # Set tick marks every 0.2
  ticks <- seq(0, 1, by = 0.2)
  
  ggplot(melted_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 5)),
      axis.title.y = element_text(size = 10, margin = margin(r = 5)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line()
    ) +
    scale_x_continuous(
      name = "Longitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    scale_y_continuous(
      name = "Latitude",
      expand = c(0, 0),
      breaks = ticks,
      labels = ticks
    ) +
    labs(fill = "Count")
}


# Determine consistent color scales for each kernel
limits_nonsep <- range(unlist(matrices_nonsep))
limits_sep <- range(unlist(matrices_sep))
limits_icm <- range(unlist(matrices_icm))

# Generate spatial grid coordinates (assuming n_counties is a perfect square)
n_grid <- sqrt(49) # Adjust based on n_counties
grid_coords <- list(
  x = seq(0, 1, length.out = n_grid),
  y = seq(0, 1, length.out = n_grid)
)

# Create heatmaps for each kernel and time point
plots_nonsep <- lapply(matrices_nonsep, function(mat) plot_heatmap_from_matrix(mat, limits_nonsep, grid_coords))
plots_sep <- lapply(matrices_sep, function(mat) plot_heatmap_from_matrix(mat, limits_sep, grid_coords))
plots_icm <- lapply(matrices_icm, function(mat) plot_heatmap_from_matrix(mat, limits_icm, grid_coords))

# Combine heatmaps into rows for each kernel
row_nonsep <- plot_grid(plotlist = plots_nonsep, ncol = length(selected_timepoints), align = "h")
row_sep <- plot_grid(plotlist = plots_sep, ncol = length(selected_timepoints), align = "h")
row_icm <- plot_grid(plotlist = plots_icm, ncol = length(selected_timepoints), align = "h")

# Add row labels for kernel types
final_plot <- plot_grid(
  ggdraw() + draw_label("Non-Separable \n Gneiting DGP", angle = 90),
  row_nonsep,
  ggdraw() + draw_label("Separable \n RBF-RBF DGP", angle = 90),
  row_sep,
  ggdraw() + draw_label("Separable \n ICM-RBF DGP", angle = 90),
  row_icm,
  ncol = 2,
  rel_widths = c(0.1, 1)
)


# Add column labels for time points
time_labels <- plot_grid(
  ggdraw() + draw_label("Time 1"),
  ggdraw() + draw_label("Time 2"),
  ggdraw() + draw_label("Time 3"),
  
  ncol = length(selected_timepoints),
  rel_widths = c(1 / length(selected_timepoints), 
                 rep(1 / length(selected_timepoints), length(selected_timepoints) - 1))
)


final_plot_with_labels <- plot_grid(
  plot_grid(NULL, time_labels, NULL, ncol = 3, rel_widths = c(.06, .8, 0.03)),  # Time labels only over heatmaps
  final_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)


# Save plot
pdf(paste0(wd,"code/simulations/figures/GP_presentation_sim_count_normal.pdf"),
    paper = 'special', width = 10, height = 7)


# Display the final plot
print(final_plot_with_labels)

dev.off()

