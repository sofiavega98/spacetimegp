# Functions for the kernel demo in paper
# This file contains core functions for simulating and analyzing spatiotemporal data
# with different kernel types and distributions.

#' Simulate Spatiotemporal Data
#'
#' @param kernel Type of kernel to use ("nonsep", "sep", "ICM", or "gneiting")
#' @param seed Random seed for reproducibility
#' @param n_counties Number of spatial locations (default: 49)
#' @param n_timepoints Number of time points (default: 15)
#' @param distribution Type of distribution ("poisson" or "normal")
#' @param spatial_lengthscale Lengthscale determining the spread of association in the spatial dimension
#' @param temporal_lengthscale Lengthscale determining the spread of association in the temporal dimension
#' @param J Number of latent factors in the ICM kernel (default: 5)
#' 
#' @return A data frame containing the simulated data with columns:
#'         - kernel: Type of kernel used
#'         - county_id: Spatial location identifier
#'         - time: Time point
#'         - lon: Longitude coordinate
#'         - lat: Latitude coordinate
#'         - gp_sample: Gaussian process sample
#'         - mean: Mean of the distribution
#'         - count: Simulated count/continuous value
#' @export
#'
#' @examples
#' # Simulate Poisson data with separable kernel
#' data <- simulate_grid_data(kernel = "sep", seed = 1, distribution = "poisson")
#' 
#' # Simulate Normal data with nonseparable kernel
#' data <- simulate_grid_data(kernel = "nonsep", seed = 1, distribution = "normal")
simulate_grid_data <- function(kernel, seed, n_counties = 49, n_timepoints = 15, distribution = "poisson",
                               spatial_lengthscale = .3, temporal_lengthscale = .9, J = 5) {
  # Load required libraries
  library(MASS)
  
  # Create a grid of counties in a 7x7 grid (49 locations)
  # This creates a regular spatial grid for simulation
  county_grid <- expand.grid(lon = seq(0, 1, length.out = 7), lat = seq(0, 1, length.out = 7))
  county_grid <- county_grid[1:n_counties, ]
  
  # Calculate spatial distance matrix between all pairs of locations
  # This will be used to compute spatial correlations
  dist_matrix <- as.matrix(dist(county_grid))
  
  # Create normalized timepoints (0 to 1)
  # This ensures consistent temporal scaling across simulations
  timepoints <- seq(0, 1, length.out = n_timepoints)
  
  # Calculate temporal distance matrix
  # Used to compute temporal correlations
  temporal_dist_matrix <- as.matrix(dist(timepoints))
  
  # Initialize covariance matrix based on kernel type
  # Each kernel type has different correlation structures
  
  if (kernel == "nonsep") {
    # Non-separable kernel implementation
    # This kernel allows for space-time interaction
    beta <- 1
    sigma <- 1
    alpha <- 1
    a <- spatial_lengthscale
    b <- temporal_lengthscale
    K <- matrix(0, nrow = n_counties * n_timepoints, ncol = n_counties * n_timepoints)
    
    # Compute covariance matrix for non-separable kernel
    # The kernel function is: K(s,t) = σ² * exp(-u/b) * exp(-exp(-u/b) * (h/a)²)
    for (i in seq_len(nrow(K))) {
      for (j in seq_len(ncol(K))) {
        county_i <- (i - 1) %/% n_timepoints + 1
        time_i <- (i - 1) %% n_timepoints + 1
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j <- (j - 1) %% n_timepoints + 1
        
        h <- dist_matrix[county_i, county_j]
        u <- temporal_dist_matrix[time_i, time_j]
        
        K[i, j] <- sigma^2 * exp(-u / b) * exp(-exp(-u / b) * (h / a)^2)
      }
    }
    
    # Generate Gaussian process sample
    set.seed(1)
    gp_sample <- mvrnorm(1, mu = rep(0, nrow(K)), Sigma = K)
    
    # Generate observations based on distribution type
    if (distribution == "poisson") {
      # For Poisson, use exponential link function
      mean <- exp(4 + gp_sample)
      set.seed(seed)
      count <- rpois(length(mean), lambda = mean)
    } else if (distribution == "normal") {
      # For Normal, use identity link function
      mean <- 4 + gp_sample
      set.seed(seed)
      count <- rnorm(length(gp_sample), mean = mean, sd = 1)
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
    
    # Return data frame with all relevant information
    return(data.frame(
      kernel = kernel,
      county_id = rep(seq_len(n_counties), each = n_timepoints),
      time = rep(seq_len(n_timepoints), times = n_counties),
      lon = rep(county_grid$lon, each = n_timepoints),
      lat = rep(county_grid$lat, each = n_timepoints),
      gp_sample = gp_sample,
      mean = mean,
      count = count
    ))
    
  } else if (kernel == "sep") {
    # Separable kernel implementation
    # This kernel assumes independence between space and time
    spatial_lengthscale <- spatial_lengthscale
    temporal_lengthscale <- temporal_lengthscale
    sigma <- 1
    
    K <- matrix(0, nrow = n_counties * n_timepoints, ncol = n_counties * n_timepoints)
    
    # Compute covariance matrix for separable kernel
    # The kernel function is: K(s,t) = σ² * exp(-0.5*(h/a)²) * exp(-0.5*(u/b)²)
    for (i in seq_len(nrow(K))) {
      for (j in seq_len(ncol(K))) {
        county_i <- (i - 1) %/% n_timepoints + 1
        time_i <- (i - 1) %% n_timepoints + 1
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j <- (j - 1) %% n_timepoints + 1
        
        K[i, j] <- sigma^2 * exp(-0.5 * (dist_matrix[county_i, county_j] / spatial_lengthscale)^2) * 
          exp(-0.5 * (temporal_dist_matrix[time_i, time_j] / temporal_lengthscale)^2)
      }
    }
    
    # Generate Gaussian process sample
    set.seed(1)
    gp_sample <- mvrnorm(1, mu = rep(0, nrow(K)), Sigma = K)
    
    # Generate observations based on distribution type
    if (distribution == "poisson") {
      mean <- exp(4 + gp_sample)
      set.seed(seed)
      count <- rpois(length(mean), lambda = mean)
    } else if (distribution == "normal") {
      mean <- 4 + gp_sample
      set.seed(seed)
      count <- rnorm(length(gp_sample), mean = mean, sd = 1)
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
    
    return(data.frame(
      kernel = kernel,
      county_id = rep(seq_len(n_counties), each = n_timepoints),
      time = rep(seq_len(n_timepoints), times = n_counties),
      lon = rep(county_grid$lon, each = n_timepoints),
      lat = rep(county_grid$lat, each = n_timepoints),
      gp_sample = gp_sample,
      mean = mean,
      count = count
    ))
    
  } else if (kernel == "ICM") {
    # Intrinsic Coregionalization Model (ICM) implementation
    # This model uses latent factors to capture spatial dependence
    J <- J
    sigma <- .4
    temporal_lengthscale <- temporal_lengthscale
    
    # Generate random coefficients for latent processes
    set.seed(1)
    beta_matrix <- matrix(rnorm(n_counties * J), nrow = n_counties, ncol = J)
    
    # Compute spatial covariance matrix using latent factors
    K_unit <- beta_matrix %*% t(beta_matrix)
    
    K <- matrix(0, nrow = n_counties * n_timepoints, ncol = n_counties * n_timepoints)
    
    # Compute covariance matrix for ICM kernel
    # The kernel function combines spatial latent factors with temporal correlation
    for (i in seq_len(nrow(K))) {
      for (j in seq_len(ncol(K))) {
        county_i <- (i - 1) %/% n_timepoints + 1
        time_i <- (i - 1) %% n_timepoints + 1
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j <- (j - 1) %% n_timepoints + 1
        
        K[i, j] <- sigma^2 * K_unit[county_i, county_j] * 
          exp(-0.5 * (temporal_dist_matrix[time_i, time_j] / temporal_lengthscale)^2)
      }
    }
    
    # Generate Gaussian process sample
    set.seed(1)
    gp_sample <- mvrnorm(1, mu = rep(0, nrow(K)), Sigma = K)
    
    # Generate observations based on distribution type
    if (distribution == "poisson") {
      mean <- exp(4 + gp_sample)
      set.seed(seed)
      count <- rpois(length(mean), lambda = mean)
    } else if (distribution == "normal") {
      mean <- 4 + gp_sample
      set.seed(seed)
      count <- rnorm(length(gp_sample), mean = mean, sd = 1)
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
    
    return(data.frame(
      kernel = kernel,
      county_id = rep(seq_len(n_counties), each = n_timepoints),
      time = rep(seq_len(n_timepoints), times = n_counties),
      lon = rep(county_grid$lon, each = n_timepoints),
      lat = rep(county_grid$lat, each = n_timepoints),
      gp_sample = gp_sample,
      mean = mean,
      count = count
    ))
    
  } else if (kernel == "gneiting") {
    # Gneiting kernel implementation
    # This is a non-separable kernel with specific space-time interaction
    sigma <- 1
    a <- temporal_lengthscale
    c <- spatial_lengthscale  # space scaling parameter
    alpha <- 1
    beta <- 1
    gamma <- 1
    d <- 2
    
    K <- matrix(0, nrow = n_counties * n_timepoints, ncol = n_counties * n_timepoints)
    
    # Compute covariance matrix for Gneiting kernel
    # The kernel function is: K(s,t) = σ² * (1/(ψ(u)^(βd/2))) * exp(-h^(2γ)/(c*ψ(u)^(βγ)))
    # where ψ(u) = (1/a * u^(2α) + 1)
    for (i in seq_len(nrow(K))) {
      for (j in seq_len(ncol(K))) {
        county_i <- (i - 1) %/% n_timepoints + 1
        time_i <- (i - 1) %% n_timepoints + 1
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j <- (j - 1) %% n_timepoints + 1
        
        h <- dist_matrix[county_i, county_j]
        u <- abs(timepoints[time_i] - timepoints[time_j])
        
        psi_u <- (1/a * u^(2 * alpha) + 1)
        space_term <- exp(-1/c * h^(2 * gamma) / psi_u^(beta * gamma))
        time_term <- 1 / psi_u^(beta * d / 2)
        
        K[i, j] <- sigma^2 * time_term * space_term
      }
    }
    
    # Generate Gaussian process sample
    set.seed(1)
    gp_sample <- mvrnorm(1, mu = rep(0, nrow(K)), Sigma = K)
    
    # Generate observations based on distribution type
    if (distribution == "poisson") {
      mean <- exp(4 + gp_sample)
      set.seed(seed)
      count <- rpois(length(mean), lambda = mean)
    } else if (distribution == "normal") {
      set.seed(seed)
      mean <- 4 + gp_sample
      count <- rnorm(length(gp_sample), mean = mean, sd = 1)
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
    
    return(data.frame(
      kernel = kernel,
      county_id = rep(seq_len(n_counties), each = n_timepoints),
      time = rep(seq_len(n_timepoints), times = n_counties),
      lon = rep(county_grid$lon, each = n_timepoints),
      lat = rep(county_grid$lat, each = n_timepoints),
      gp_sample = gp_sample,
      mean = mean,
      count = count
    ))
  } else {
    stop("Invalid kernel type. Choose from 'nonsep', 'sep', 'ICM', or 'gneiting'.")
  }
}

#' Compute True Weights from Simulated Data
#'
#' @param simulated_data Data frame containing simulated spatiotemporal data
#' @param kernel Type of kernel used ("gneiting", "sep", or "ICM")
#' @param trt_counties Vector of treated county IDs
#' @param seed_for_icm Random seed for ICM kernel (optional)
#' @param nugget Small value for numerical stability (default: 1e-6)
#' @param s_lengthscale Spatial lengthscale parameter
#' @param t_lengthscale Temporal lengthscale parameter
#'
#' @return List containing:
#'         - weights: Matrix of computed weights
#'         - treated_idx: Indices of treated observations
#'         - control_idx: Indices of control observations
#'         - K: Full covariance matrix
#'         - weight_df: Data frame with weights and county/time information
#' @export
compute_true_weights <- function(simulated_data, kernel, trt_counties, seed_for_icm = NULL, nugget = 1e-6, s_lengthscale, t_lengthscale) {
  # Determine dimensions from the simulated data
  n_timepoints <- length(unique(simulated_data$time))
  counties <- sort(unique(simulated_data$county_id))
  n_counties <- length(counties)
  
  # Reconstruct county grid following the simulation convention:
  # A 7x7 grid with only the first 49 counties used
  county_grid <- expand.grid(lon = seq(0, 1, length.out = 7),
                             lat = seq(0, 1, length.out = 7))
  county_grid <- county_grid[1:n_counties, ]
  
  # Compute spatial and temporal distance matrices
  # These matrices are used to calculate correlations between locations and time points
  dist_matrix <- as.matrix(dist(county_grid))
  timepoints_norm <- seq(0, 1, length.out = n_timepoints)
  temporal_dist_matrix <- as.matrix(dist(timepoints_norm))
  
  # Initialize full covariance matrix K (size: n_counties * n_timepoints)
  N_total <- n_counties * n_timepoints
  K <- matrix(0, nrow = N_total, ncol = N_total)
  
  # Construct the covariance matrix based on the specified kernel
  if (kernel == "gneiting") {
    # Gneiting kernel parameters
    sigma <- 1
    c <- s_lengthscale
    a <- t_lengthscale
    alpha <- 1
    beta <- 1
    gamma <- 1
    d <- 2
    
    # Compute covariance matrix using Gneiting kernel function
    # K(s,t) = σ² * (1/(ψ(u)^(βd/2))) * exp(-h^(2γ)/(c*ψ(u)^(βγ)))
    # where ψ(u) = (1/a * u^(2α) + 1)
    for (i in 1:N_total) {
      county_i <- (i - 1) %/% n_timepoints + 1
      time_i   <- (i - 1) %% n_timepoints + 1
      for (j in 1:N_total) {
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j   <- (j - 1) %% n_timepoints + 1
        h <- dist_matrix[county_i, county_j]
        u <- temporal_dist_matrix[time_i, time_j]
        psi_u <- (1/a * u^(2 * alpha) + 1)
        space_term <- exp(-1/c * h^(2 * gamma) / psi_u^(beta * gamma))
        time_term <- 1 / psi_u^(beta * d / 2)
        
        K[i, j] <- sigma^2 * time_term * space_term
      }
    }
    
  } else if (kernel == "sep") {
    # Separable kernel parameters
    sigma <- 1.5
    spatial_lengthscale <- s_lengthscale
    temporal_lengthscale <- t_lengthscale
    
    # Compute covariance matrix using separable kernel function
    # K(s,t) = σ² * exp(-0.5*(h/a)²) * exp(-0.5*(u/b)²)
    for (i in 1:N_total) {
      county_i <- (i - 1) %/% n_timepoints + 1
      time_i   <- (i - 1) %% n_timepoints + 1
      for (j in 1:N_total) {
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j   <- (j - 1) %% n_timepoints + 1
        h <- dist_matrix[county_i, county_j]
        u <- temporal_dist_matrix[time_i, time_j]
        K[i, j] <- sigma^2 * exp(-0.5 * (h / spatial_lengthscale)^2) *
          exp(-0.5 * (u / temporal_lengthscale)^2)
      }
    }
    
  } else if (kernel == "ICM") {
    # ICM kernel parameters
    sigma <- 0.4
    temporal_lengthscale <- t_lengthscale
    J <- 5  # number of latent processes
    
    # Generate random coefficients for latent processes if seed provided
    if (!is.null(seed_for_icm)) {
      set.seed(seed_for_icm)
    }
    beta_matrix <- matrix(rnorm(n_counties * J), nrow = n_counties, ncol = J)
    K_unit <- beta_matrix %*% t(beta_matrix)
    
    # Compute covariance matrix using ICM kernel function
    # Combines spatial latent factors with temporal correlation
    for (i in 1:N_total) {
      county_i <- (i - 1) %/% n_timepoints + 1
      time_i   <- (i - 1) %% n_timepoints + 1
      for (j in 1:N_total) {
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j   <- (j - 1) %% n_timepoints + 1
        u <- temporal_dist_matrix[time_i, time_j]
        K[i, j] <- sigma^2 * K_unit[county_i, county_j] *
          exp(-0.5 * (u / temporal_lengthscale)^2)
      }
    }
    
  } else {
    stop("Invalid kernel type. Please choose 'gneiting', 'sep', or 'ICM'.")
  }
  
  # Identify treated and control indices
  # Treated observations are those with counties in trt_counties and time > 7
  treated_idx <- which(simulated_data$county_id %in% trt_counties & simulated_data$time > 7)
  control_idx <- setdiff(1:N_total, treated_idx)
  
  # Extract features for treated and control observations
  treated_counties <- simulated_data$county_id[treated_idx]
  treated_time <- simulated_data$time[treated_idx]
  control_counties <- simulated_data$county_id[control_idx]
  control_time <- simulated_data$time[control_idx]
  
  # Partition the covariance matrix into treated-control and control-control blocks
  K_tc <- K[treated_idx, control_idx, drop = FALSE]  # treated-control covariance
  K_cc <- K[control_idx, control_idx, drop = FALSE]  # control-control covariance
  
  # Compute weights using the formula: w = K_tc * (K_cc + σ²I)^(-1)
  # Add nugget to diagonal for numerical stability
  sigma <- 1
  weights <- K_tc %*% solve(K_cc + diag(sigma^2, nrow(K_cc)))
  
  # Create a tidy data frame from the weight matrix
  weight_df <- data.frame(
    weights = as.vector(t(weights)),
    treated_county = rep(treated_counties, each = ncol(weights)),
    treated_time = rep(treated_time, each = ncol(weights)),
    control_county = rep(control_counties, times = nrow(weights)),
    control_time = rep(control_time, times = nrow(weights))
  )
  
  # Return list with computed weights and related information
  list(
    weights = weights,
    treated_idx = treated_idx,
    control_idx = control_idx,
    K = K,
    weight_df = weight_df
  )
}

#' Create Spatial Weight Heatmap (Old Version)
#'
#' @param weight_df Data frame containing weight information
#'
#' @return ggplot object showing spatial distribution of weights
#' @export
spatial_w_heatmap_old <- function(weight_df) {
  # Create county grid coordinates
  county_coords <- expand.grid(
    lon = seq(0, 1, length.out = sqrt(n_counties)),
    lat = seq(0, 1, length.out = sqrt(n_counties))
  )
  county_coords <- county_coords[1:n_counties, ]
  county_coords$control_county <- 1:n_counties
  
  # Aggregate weights by treated and control county
  aggregated_data <- weight_df %>%
    group_by(treated_county, control_county) %>%
    summarize(avg_weight = mean(weights, na.rm = TRUE)) %>%
    ungroup()
  
  # Merge weights with county coordinates
  merged_true_weights <- merge(
    aggregated_data, county_coords,
    by.x = "control_county",
    by.y = "control_county",
    all.x = TRUE
  )
  
  # Define treated units for border overlay
  treated_units <- c(24, 38, 6, 49)
  treated_coords <- county_coords %>%
    filter(control_county %in% treated_units) %>%
    mutate(treated_county = control_county)
  
  # Create faceted heatmap plot
  ggplot() +
    geom_tile(
      data = merged_true_weights,
      aes(x = lon, y = lat, fill = avg_weight)
    ) +
    geom_tile(
      data = treated_coords,
      aes(x = lon, y = lat),
      fill = NA, color = "black", size = 0.5
    ) +
    facet_wrap(~ treated_county) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = "Avg. True Weight"
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12))
}

#' Create Spatial Weight Heatmap (New Version)
#'
#' @param weight_df Data frame containing weight information
#' @param fill_limits Optional vector of length 2 specifying the limits for the fill scale
#'
#' @return ggplot object showing spatial distribution of weights
#' @export
spatial_w_heatmap <- function(weight_df, fill_limits = NULL) {
  # Create county grid coordinates
  county_coords <- expand.grid(
    lon = seq(0, 1, length.out = sqrt(n_counties)),
    lat = seq(0, 1, length.out = sqrt(n_counties))
  )
  county_coords <- county_coords[1:n_counties, ]
  county_coords$control_county <- 1:n_counties
  
  # Aggregate weights by treated and control county
  aggregated_data <- weight_df %>%
    group_by(treated_county, control_county) %>%
    summarize(avg_weight = mean(weights, na.rm = TRUE), .groups = 'drop')
  
  # Merge weights with county coordinates
  merged_true_weights <- merge(
    aggregated_data, county_coords,
    by.x = "control_county", by.y = "control_county",
    all.x = TRUE
  )
  
  # Define treated units for border overlay
  treated_units <- c(24, 38, 6, 49)
  treated_coords <- county_coords %>%
    filter(control_county %in% treated_units) %>%
    mutate(treated_county = control_county)
  
  # Create faceted heatmap plot with optional fill limits
  ggplot() +
    geom_tile(
      data = merged_true_weights,
      aes(x = lon, y = lat, fill = avg_weight)
    ) +
    geom_tile(
      data = treated_coords,
      aes(x = lon, y = lat),
      fill = NA, color = "black", size = 0.5
    ) +
    facet_wrap(~ treated_county) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = "Avg. True Weight",
      limits = fill_limits
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12))
}
