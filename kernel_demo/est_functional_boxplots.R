# ---------------------------------------------------------------------
# Script: est_functional_boxplots.R
# Purpose: Estimate and visualize functional boxplots of f̂ₕ(u) for different kernels
#          using Stan posterior parameters and simulated data
#          Also summarizes separability and f̂ₕ(u) statistics
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Load Required Libraries
# ---------------------------------------------------------------------
library(dplyr)      # Data manipulation
library(tidyr)      # Data tidying
library(fda)        # For fbplot (functional boxplot)
library(geosphere)  # For distance calculations

# ---------------------------------------------------------------------
# Define Gneiting Covariance Function
# ---------------------------------------------------------------------
#' Gneiting Covariance Function
#' 
#' @param h Spatial distance
#' @param u Temporal lag
#' @param sigma2 Variance parameter
#' @param beta, alpha, gamma, a, c Kernel parameters
#' @return Covariance value
#' @examples
#' gneiting_cov(0.1, 2, 1, 1, 1, 1, 0.9, 3)
gneiting_cov <- function(h, u, sigma2, beta, alpha, gamma, a, c) {
  denom <- 1/a * abs(u)^(2 * alpha) + 1
  cov <- sigma2 / denom^(beta / 2) * exp(-1/c * abs(h)^(2 * gamma) / denom^(beta * gamma))
  return(cov)
}

# ---------------------------------------------------------------------
# Parametric f̂ₕ(u) Estimator Using Stan Parameters
# ---------------------------------------------------------------------
#' Estimate f̂ₕ(u) from Stan posterior parameters
#' 
#' @param data Data frame with county_id
#' @param n_timepoints Number of time points
#' @param H Matrix of spatial distances
#' @param sigma2, beta, alpha, gamma, a, c Kernel parameters
#' @param max_u Maximum temporal lag to compute
#' @return Matrix of f̂ₕ(u) values (columns: county pairs, rows: lags)
#' @details Computes f̂ₕ(u) = C(h,u)/C(h,0) - C(0,u)/C(0,0) for each county pair
estimate_fhat_from_data_parametric <- function(data, n_timepoints, H,
                                               sigma2, beta, alpha, gamma, a, c,
                                               max_u = NULL) {
  if (is.null(max_u)) max_u <- n_timepoints - 1
  counties <- sort(unique(data$county_id))
  n_counties <- length(counties)
  fhat_list <- list()
  
  for (i in 1:(n_counties - 1)) {
    for (j in (i + 1):n_counties) {
      f_hat <- numeric(max_u + 1)
      h <- H[i, j]
      C_h0 <- gneiting_cov(h, 0, sigma2, beta, alpha, gamma, a, c)
      C_00 <- gneiting_cov(0, 0, sigma2, beta, alpha, gamma, a, c)
      
      if (abs(C_h0) < 1e-9 || abs(C_00) < 1e-9) {
        f_hat[] <- NA
      } else {
        for (u in 0:max_u) {
          C_hu <- gneiting_cov(h, u, sigma2, beta, alpha, gamma, a, c)
          C_0u <- gneiting_cov(0, u, sigma2, beta, alpha, gamma, a, c)
          f_hat[u + 1] <- (C_hu / C_h0) - (C_0u / C_00)
        }
      }
      
      fhat_list[[paste0("(", counties[i], ",", counties[j], ")")]] <- f_hat
    }
  }
  
  fhat_matrix <- do.call(cbind, fhat_list)
  fhat_matrix <- fhat_matrix[, colSums(is.na(fhat_matrix)) < nrow(fhat_matrix), drop = FALSE]
  if (ncol(fhat_matrix) == 0) return(NULL)
  return(fhat_matrix)
}

# ---------------------------------------------------------------------
# Side-by-Side f̂ₕ(u) Plots for Different Kernels
# ---------------------------------------------------------------------
# Settings for analysis
kernel_choices <- c("ICM", "sep", "gneiting_v2")  # Kernel types to compare
titles <- c("Separable ICM-RBF", "Separable RBF-RBF", "Nonseparable Gneiting")
wd <- "/n/holylfs05/LABS/nethery_lab/Users/svega/gp_nonsep/"  # Working directory for data/results
seed_value <- 1
n_counties <- 49
n_timepoints <- 15
max_u <- 14

# Create spatial grid and compute distance matrix
county_grid <- expand.grid(lon = seq(0, 1, length.out = 7),
                           lat = seq(0, 1, length.out = 7))
county_grid <- county_grid[1:n_counties, ]
H <- as.matrix(dist(county_grid))

# Open PNG device for output
png(paste0(wd,"code/simulations/figures/est_functional_boxplots_0617.png"), width = 1200, height = 600, res = 200)

# Set up plot layout for three panels
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
y_range <- NULL
fhat_all <- list()

# Loop through kernels and estimate f̂ₕ(u) for each
for (k in seq_along(kernel_choices)) {
  kernel <- kernel_choices[k]
  
  # Load Stan fit and simulated data
  stan_fit_file <- paste0(wd, "code/simulations/results/normal/nonsep_oldprior/stan_fit_", kernel, "_seed_", seed_value, ".rds")
  sim_data_file <- paste0(wd, "code/simulations/data/simulated_normal_data_", kernel, "_seed_", seed_value, ".rds")
  stan_fit <- readRDS(stan_fit_file)
  posterior <- rstan::extract(stan_fit)
  data <- readRDS(sim_data_file)
  data <- data %>% mutate(county_id = as.integer(factor(county_id)))
  
  # Extract posterior medians for kernel parameters
  sigma2 <- median(posterior$export_sigma_global^2)
  beta   <- median(posterior$export_beta)
  alpha  <- median(posterior$export_alpha)
  gamma  <- median(posterior$export_gamma)
  a      <- median(posterior$export_lengthscale_time)
  c      <- median(posterior$export_lengthscale_space)
  
  # Estimate f̂ₕ(u) matrix for this kernel
  fhat_matrix <- estimate_fhat_from_data_parametric(data,
                                                    n_timepoints,
                                                    H,
                                                    sigma2, beta, alpha, gamma, a, c,
                                                    max_u = max_u)
  fhat_all[[k]] <- fhat_matrix
  if (!is.null(fhat_matrix)) {
    y_range <- range(c(y_range, range(fhat_matrix, na.rm = TRUE)))
  }
}

# Final plotting loop: plot functional boxplots for each kernel
for (k in seq_along(kernel_choices)) {
  fhat_matrix <- fhat_all[[k]]
  if (!is.null(fhat_matrix)) {
    fbplot(fit = fhat_matrix,
           x = 0:max_u,
           xlab = "Temporal Lag (u)",
           ylab = expression(hat(f)[h](u)),
           main = titles[k],
           ylim = y_range)
    abline(h = 0, col = "red", lty = 2)
  } else {
    plot.new()
    title(main = paste(titles[k], "\n(no curves)"))
  }
}

dev.off()  # Close the PNG device

# ----------------------------------------------------
# Summarize separability for each kernel using posterior draws
# ----------------------------------------------------

# Pre-allocate storage for summary statistics
prop_exact_zero <- numeric(length(kernel_choices))
prop_near_zero  <- numeric(length(kernel_choices))

# For each kernel, compute proportion of posterior samples with beta=0 or near zero
for (k in seq_along(kernel_choices)) {
  # Reload posterior draws
  kernel <- kernel_choices[k]
  stan_fit_file <- paste0(wd, "code/simulations/results/normal/nonsep_oldprior/stan_fit_", kernel, "_seed_", seed_value, ".rds")
  stan_fit <- readRDS(stan_fit_file)
  posterior <- rstan::extract(stan_fit)
  beta_samples <- posterior$export_beta
  
  # Proportion exactly zero and near zero
  prop_exact_zero[k] <- mean(beta_samples == 0)
  prop_near_zero[k]  <- mean(beta_samples < 1e-1)
}

# Create and print summary table
separability_tbl <- data.frame(
  Kernel              = titles,
  Prop_Exact_Zero     = round(prop_exact_zero, 3),
  Prop_Effectively_0  = round(prop_near_zero,  3),
  stringsAsFactors    = FALSE
)
print(separability_tbl)

# Optionally save as CSV
write.csv(separability_tbl,
          file = paste0(wd, "code/simulations/figures/separability_summary.csv"),
          row.names = FALSE)

# ----------------------------------------------------
# Summarize max and mean |fhat| by kernel
# ----------------------------------------------------

# Compute max absolute fhat for each kernel
max_abs_fhat <- sapply(fhat_all, function(mat) {
  if (is.null(mat)) return(NA_real_)
  max(abs(mat), na.rm = TRUE)
})

# Assemble into a data.frame and print
summary_tbl <- data.frame(
  Kernel             = titles,
  Max_Abs_fhat       = round(max_abs_fhat, 3),
  stringsAsFactors   = FALSE
)
print(summary_tbl)

# Optionally save to CSV
write.csv(summary_tbl,
          file = file.path(wd, "code/simulations/figures/max_abs_fhat_summary.csv"),
          row.names = FALSE)

# Compute mean absolute fhat for each kernel
mean_abs_fhat <- sapply(fhat_all, function(mat) {
  if (is.null(mat)) return(NA_real_)
  mean(abs(mat), na.rm = TRUE)
})

# Assemble into a data.frame and print
summary_tbl <- data.frame(
  Kernel             = titles,
  Mean_Abs_fhat       = round(mean_abs_fhat, 3),
  stringsAsFactors   = FALSE
)
print(summary_tbl)

# Optionally save to CSV
write.csv(summary_tbl,
          file = file.path(wd, "code/simulations/figures/mean_abs_fhat_summary.csv"),
          row.names = FALSE)
