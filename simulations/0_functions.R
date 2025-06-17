## Script to create data generating processes for simulations_grid_small_RCN_newDGP_upint_v3
## Author: Sofia Vega
## Date Created: May 11, 2021

#' Simulate Spatiotemporal Data with Various Kernel Structures
#'
#' This function generates simulated spatiotemporal data using different kernel structures
#' (nonseparable, separable, ICM, or Gneiting) and distributions (Poisson or Normal).
#' The data is generated on a grid of spatial locations over multiple time points.
#'
#' @param kernel Character. Type of kernel to use:
#'   - "nonsep": Non-separable spatiotemporal kernel
#'   - "sep": Separable spatiotemporal kernel
#'   - "ICM": Intrinsic Coregionalization Model
#'   - "gneiting": Gneiting's non-separable covariance function
#' @param seed Integer. Random seed for reproducibility
#' @param n_counties Integer. Number of spatial locations (default: 49)
#' @param n_timepoints Integer. Number of time points (default: 15)
#' @param distribution Character. Type of distribution:
#'   - "poisson": Poisson distribution with log link
#'   - "normal": Normal distribution with identity link
#' @param spatial_lengthscale Numeric. Lengthscale parameter for spatial correlation (default: 0.3)
#' @param temporal_lengthscale Numeric. Lengthscale parameter for temporal correlation (default: 0.9)
#' @param J Integer. Number of latent factors for ICM kernel (default: 5)
#' 
#' @return A data frame containing:
#'   - kernel: Type of kernel used
#'   - county_id: Spatial location identifier
#'   - time: Time point identifier
#'   - lon: Longitude coordinate
#'   - lat: Latitude coordinate
#'   - gp_sample: Generated Gaussian process sample
#'   - mean: Mean of the distribution
#'   - count: Simulated response variable
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
  
  # Step 1: Create spatial grid
  # Generate a 7x7 grid of locations and subset to desired number of counties
  county_grid <- expand.grid(lon = seq(0, 1, length.out = 7), lat = seq(0, 1, length.out = 7))
  county_grid <- county_grid[1:n_counties, ]
  
  # Step 2: Calculate distance matrices
  # Spatial distance matrix between all pairs of locations
  dist_matrix <- as.matrix(dist(county_grid))
  
  # Create normalized timepoints and their distance matrix
  timepoints <- seq(0, 1, length.out = n_timepoints)
  temporal_dist_matrix <- as.matrix(dist(timepoints))
  
  # Step 3: Initialize and compute covariance matrix based on kernel type
  if (kernel == "nonsep") {
    # Non-separable kernel parameters
    beta <- 1      # Interaction parameter
    sigma <- 1     # Overall variance
    alpha <- 1     # Smoothness parameter
    a <- spatial_lengthscale    # Spatial lengthscale
    b <- temporal_lengthscale   # Temporal lengthscale
    
    # Initialize covariance matrix
    K <- matrix(0, nrow = n_counties * n_timepoints, ncol = n_counties * n_timepoints)
    
    # Compute covariance for each pair of space-time points
    for (i in seq_len(nrow(K))) {
      for (j in seq_len(ncol(K))) {
        # Extract spatial and temporal indices
        county_i <- (i - 1) %/% n_timepoints + 1
        time_i <- (i - 1) %% n_timepoints + 1
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j <- (j - 1) %% n_timepoints + 1
        
        # Calculate spatial and temporal distances
        h <- dist_matrix[county_i, county_j]
        u <- temporal_dist_matrix[time_i, time_j]
        
        # Non-separable kernel formula:
        # K(s,t) = σ² * exp(-u/b) * exp(-exp(-u/b) * (h/a)²)
        K[i, j] <- sigma^2 * exp(-u / b) * exp(-exp(-u / b) * (h / a)^2)
      }
    }
    
    # Step 4: Generate GP sample
    set.seed(1)  # Fixed seed for consistent mean across simulations
    gp_sample <- mvrnorm(1, mu = rep(0, nrow(K)), Sigma = K)
    
    # Step 5: Generate response variable based on distribution
    if (distribution == "poisson") {
      # Poisson: mean = exp(intercept + GP)
      mean <- exp(4 + gp_sample)
      set.seed(seed)  # User-specified seed for response generation
      count <- rpois(length(mean), lambda = mean)
    } else if (distribution == "normal") {
      # Normal: mean = intercept + GP, with unit variance
      mean <- 4 + gp_sample
      set.seed(seed)
      count <- rnorm(length(gp_sample), mean = mean, sd = 1)
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
    
    # Step 6: Return formatted data frame
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
    # Separable kernel parameters
    spatial_lengthscale <- spatial_lengthscale
    temporal_lengthscale <- temporal_lengthscale
    sigma <- 1
    
    # Initialize covariance matrix
    K <- matrix(0, nrow = n_counties * n_timepoints, ncol = n_counties * n_timepoints)
    
    # Compute covariance for each pair of space-time points
    for (i in seq_len(nrow(K))) {
      for (j in seq_len(ncol(K))) {
        # Extract spatial and temporal indices
        county_i <- (i - 1) %/% n_timepoints + 1
        time_i <- (i - 1) %% n_timepoints + 1
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j <- (j - 1) %% n_timepoints + 1
        
        # Separable kernel formula:
        # K(s,t) = σ² * exp(-0.5(h/ℓₛ)²) * exp(-0.5(u/ℓₜ)²)
        K[i, j] <- sigma^2 * exp(-0.5 * (dist_matrix[county_i, county_j] / spatial_lengthscale)^2) * 
          exp(-0.5 * (temporal_dist_matrix[time_i, time_j] / temporal_lengthscale)^2)
      }
    }
    
    # Generate GP sample and response variable (same as nonsep)
    set.seed(1)
    gp_sample <- mvrnorm(1, mu = rep(0, nrow(K)), Sigma = K)
    
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
    # ICM kernel parameters
    J <- J                    # Number of latent factors
    sigma <- .4              # Overall variance
    temporal_lengthscale <- temporal_lengthscale
    
    # Generate random coefficients for latent processes
    set.seed(1)
    beta_matrix <- matrix(rnorm(n_counties * J), nrow = n_counties, ncol = J)
    
    # Compute spatial covariance matrix using latent factors
    K_unit <- beta_matrix %*% t(beta_matrix)
    
    # Initialize full covariance matrix
    K <- matrix(0, nrow = n_counties * n_timepoints, ncol = n_counties * n_timepoints)
    
    # Compute covariance for each pair of space-time points
    for (i in seq_len(nrow(K))) {
      for (j in seq_len(ncol(K))) {
        county_i <- (i - 1) %/% n_timepoints + 1
        time_i <- (i - 1) %% n_timepoints + 1
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j <- (j - 1) %% n_timepoints + 1
        
        # ICM kernel formula:
        # K(s,t) = σ² * K_unit(s,s') * exp(-0.5(u/ℓₜ)²)
        K[i, j] <- sigma^2 * K_unit[county_i, county_j] * 
          exp(-0.5 * (temporal_dist_matrix[time_i, time_j] / temporal_lengthscale)^2)
      }
    }
    
    # Generate GP sample and response variable (same as other kernels)
    set.seed(1)
    gp_sample <- mvrnorm(1, mu = rep(0, nrow(K)), Sigma = K)
    
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
    # Gneiting kernel parameters
    sigma <- 1
    a <- temporal_lengthscale
    c <- spatial_lengthscale
    alpha <- 1
    beta <- 1
    gamma <- 1
    d <- 2
    
    # Initialize covariance matrix
    K <- matrix(0, nrow = n_counties * n_timepoints, ncol = n_counties * n_timepoints)
    
    # Compute covariance for each pair of space-time points
    for (i in seq_len(nrow(K))) {
      for (j in seq_len(ncol(K))) {
        county_i <- (i - 1) %/% n_timepoints + 1
        time_i <- (i - 1) %% n_timepoints + 1
        county_j <- (j - 1) %/% n_timepoints + 1
        time_j <- (j - 1) %% n_timepoints + 1
        
        # Calculate distances
        h <- dist_matrix[county_i, county_j]
        u <- abs(timepoints[time_i] - timepoints[time_j])
        
        # Gneiting kernel formula:
        # K(s,t) = σ² * (1/ψ(u)^(βd/2)) * exp(-c*h^(2γ)/ψ(u)^(βγ))
        # where ψ(u) = (a|u|^(2α) + 1)
        psi_u <- (a * u^(2 * alpha) + 1)
        space_term <- exp(-c * h^(2 * gamma) / psi_u^(beta * gamma))
        time_term <- 1 / psi_u^(beta * d / 2)
        
        K[i, j] <- sigma^2 * time_term * space_term
      }
    }
    
    # Generate GP sample and response variable (same as other kernels)
    set.seed(1)
    gp_sample <- mvrnorm(1, mu = rep(0, nrow(K)), Sigma = K)
    
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
    stop("Invalid kernel type. Choose from 'nonsep', 'sep', or 'ICM'.")
  }
}

#' Format Data for Non-separable Gaussian Process Stan Model
#'
#' This function prepares data for fitting a non-separable Gaussian Process model using Stan.
#' It handles the creation of distance matrices and indices for spatial and temporal components.
#'
#' @param simulated_data Data frame containing the simulated data with columns:
#'   - county_id: Spatial location identifier
#'   - time: Time point identifier
#'   - lon: Longitude coordinate
#'   - lat: Latitude coordinate
#'   - count: Response variable
#' @param trt_counties Vector of county IDs that received treatment
#'
#' @return A list containing:
#'   - N: Number of time points
#'   - D: Number of spatial locations
#'   - n_k_f: Number of knots for basis functions
#'   - time_points: Vector of time points
#'   - y: Response variable vector
#'   - num_treated: Number of treated observations
#'   - control_idx: Indices of control observations
#'   - H: Vector of spatial distances
#'   - U: Vector of temporal distances
#'   - H_ind: Matrix of spatial distance indices
#'   - U_ind: Matrix of temporal distance indices
#'
#' @examples
#' # Format data for non-separable GP model
#' stan_data <- format_for_stan_nonsep(sim_data, trt_counties = c(1, 2, 3))
format_for_stan_nonsep <- function(simulated_data, trt_counties) {
  # Define treated counties
  trt_counties <- trt_counties
  
  # Order simulated_data by county_id
  simulated_data <- simulated_data[order(simulated_data$county_id), ]
  
  # Identify indices for Connecticut counties (treated) and only include times after time = 7
  treated_indices <- which(simulated_data$county_id %in% trt_counties & simulated_data$time > 7)
  
  # Total number of observations
  N <- length(unique(simulated_data$time))
  D <- length(unique(simulated_data$county_id))
  
  # Covariate and target variable
  x_covariate <- unique(simulated_data$time)
  y_target_variable <- simulated_data$count
  
  # Determine control indices (non-Connecticut counties)
  # Identifies the control group indices by excluding the indices of the treated group
  control_indices <- setdiff(seq_along(y_target_variable), treated_indices)
  
  # Compute spatial distance matrix H and its index matrix H_ind
  coords <- unique(simulated_data[, c("county_id", "lon", "lat")])
  H_matrix <- as.matrix(dist(coords[, c("lon", "lat")]))
  H <- as.vector(H_matrix)
  
  # Convert H and U into long vectors
  H_vec = as.vector(H)
  DN <- D * N
  H_ind <- matrix(0, nrow = D * N, ncol = D * N)
  
  for (i in 1:(D*N)) {
    for (j in 1:(D*N)) {
      unit_i <- ((i - 1) %/% N) + 1  # outer index is COUNTY
      unit_j <- ((j - 1) %/% N) + 1
      H_ind[i, j] <- (unit_i - 1) * D + unit_j
    }
  }
  
  
  # Make matrix of the H indices 
  #H_seq <- seq_along(H_vec)
  #H_temp <- c()
  
  #for(i in seq_len(length(H_seq) / D)) {
  # H_temp <- rbind(H_temp, matrix(H_seq[i], N, N))
  #}
  
  # H_ind <- H_temp
  # for(i in seq_len(D - 1)) {
  #  H_ind <- cbind(H_ind, H_temp + i * D)
  #}
  
  # Compute temporal distance matrix U and its index matrix U_ind
  # Temporal coordinates (time points)
  temporal_coords <- unique(simulated_data$time)
  U_matrix <- as.matrix(dist(temporal_coords))
  
  U <- as.vector(U_matrix)
  
  U_ind <- matrix(0, nrow = D * N, ncol = D * N)
  
  for (i in 1:(D*N)) {
    for (j in 1:(D*N)) {
      time_i <- ((i - 1) %% N) + 1
      time_j <- ((j - 1) %% N) + 1
      U_ind[i, j] <- (time_i - 1) * N + time_j
    }
  }
  
  # make matrix of the U indices
  # this should be an m*nxm*n matrix
  #U_seq <- 1:length(U)
  #U_seq_mat <- matrix(U_seq,N,N)
  #U_temp <- U_seq_mat
  #for(i in 1:(D-1)){
  #  U_temp <- rbind(U_temp,U_seq_mat)
  #}
  
  #U_ind <- c()
  # for(i in 1:(D)){ 
  #  U_ind <- cbind(U_ind,U_temp)
  #}
  
  #print(summary(H_ind))  # Should be in 1:(D*D)
  #print(summary(U_ind))  # Should be in 1:(N*N)
  
  list(
    N = N,
    D = D,
    n_k_f = 15,
    time_points = x_covariate,
    y = y_target_variable,
    num_treated = length(treated_indices),
    control_idx = control_indices,
    H = H,
    U = U,
    H_ind = H_ind,
    U_ind = U_ind
  )
}

#' Format Data for Separable Gaussian Process Stan Model
#'
#' This function prepares data for fitting a separable Gaussian Process model using Stan.
#' It handles the creation of spatial distance matrices and indices for the separable kernel structure.
#'
#' @param simulated_data Data frame containing the simulated data with columns:
#'   - county_id: Spatial location identifier
#'   - time: Time point identifier
#'   - lon: Longitude coordinate
#'   - lat: Latitude coordinate
#'   - count: Response variable
#' @param trt_counties Vector of county IDs that received treatment
#'
#' @return A list containing:
#'   - N: Number of time points
#'   - D: Number of spatial locations
#'   - n_k_f: Number of knots for basis functions
#'   - time_points: Vector of time points
#'   - y: Response variable vector
#'   - num_treated: Number of treated observations
#'   - control_idx: Indices of control observations
#'   - H: Vector of spatial distances
#'   - H_ind: Matrix of spatial distance indices
#'   - distance_matrix: Full spatial distance matrix
#'
#' @examples
#' # Format data for separable GP model
#' stan_data <- format_for_stan_sep(sim_data, trt_counties = c(1, 2, 3))
format_for_stan_sep <- function(simulated_data, trt_counties) {
  # Define treated counties
  trt_counties <- trt_counties
  
  # Order simulated_data by county_id
  simulated_data <- simulated_data[order(simulated_data$county_id), ]
  
  # Identify indices for Connecticut counties (treated) and only include times after time = 7
  treated_indices <- which(simulated_data$county_id %in% trt_counties & simulated_data$time > 7)
  
  # Total number of observations
  N <- length(unique(simulated_data$time))
  D <- length(unique(simulated_data$county_id))
  
  # Covariate and target variable
  x_covariate <- unique(simulated_data$time)
  y_target_variable <- simulated_data$count
  
  # Determine control indices (non-Connecticut counties)
  control_indices <- setdiff(seq_along(y_target_variable), treated_indices)
  
  # Compute spatial distance matrix H and its index matrix H_ind
  coords <- unique(simulated_data[, c("county_id", "lon", "lat")])
  H_matrix <- as.matrix(dist(coords[, c("lon", "lat")]))
  H <- as.vector(H_matrix)
  
  # Convert H and U into long vectors
  H_vec = as.vector(H)
  DN <- D * N
  H_ind <- matrix(0, nrow = D * N, ncol = D * N)
  
  for (i in 1:(D*N)) {
    for (j in 1:(D*N)) {
      unit_i <- ((i - 1) %/% N) + 1  # outer index is COUNTY
      unit_j <- ((j - 1) %/% N) + 1
      H_ind[i, j] <- (unit_i - 1) * D + unit_j
    }
  }
  
  
  # Make matrix of the H indices 
  #H_seq <- seq_along(H_vec)
  #H_temp <- c()
  
  #for(i in seq_len(length(H_seq) / D)) {
  # H_temp <- rbind(H_temp, matrix(H_seq[i], N, N))
  #}
  
  # H_ind <- H_temp
  # for(i in seq_len(D - 1)) {
  #  H_ind <- cbind(H_ind, H_temp + i * D)
  #}
  
  # Compute temporal distance matrix U and its index matrix U_ind
  # Temporal coordinates (time points)
  
  
  # make matrix of the U indices
  # this should be an m*nxm*n matrix
  #U_seq <- 1:length(U)
  #U_seq_mat <- matrix(U_seq,N,N)
  #U_temp <- U_seq_mat
  #for(i in 1:(D-1)){
  #  U_temp <- rbind(U_temp,U_seq_mat)
  #}
  
  #U_ind <- c()
  # for(i in 1:(D)){ 
  #  U_ind <- cbind(U_ind,U_temp)
  #}
  
  #print(summary(H_ind))  # Should be in 1:(D*D)
  #print(summary(U_ind))  # Should be in 1:(N*N)
  
  list(
    N = N,
    D = D,
    n_k_f = 15,
    time_points = x_covariate,
    y = y_target_variable,
    num_treated = length(treated_indices),
    control_idx = control_indices,
    H = H,
    H_ind = H_ind,
    distance_matrix = H_matrix 
  )
}

# Function to format data for poisson_RBF.stan
format_for_stan_sep <- function(simulated_data, trt_counties) {
  # Define treated counties
  trt_counties <- trt_counties
  
  # Order simulated_data by county_id
  simulated_data <- simulated_data[order(simulated_data$county_id), ]
  
  # Identify indices for Connecticut counties (treated) and only include times after time = 7
  treated_indices <- which(simulated_data$county_id %in% trt_counties & simulated_data$time > 7)
  
  # Total number of observations
  N <- length(unique(simulated_data$time))
  D <- length(unique(simulated_data$county_id))
  
  # Covariate and target variable
  x_covariate <- unique(simulated_data$time)
  y_target_variable <- simulated_data$count
  
  # Determine control indices (non-Connecticut counties)
  control_indices <- setdiff(seq_along(y_target_variable), treated_indices)
  
  # Compute spatial distance matrix H and its index matrix H_ind
  coords <- unique(simulated_data[, c("county_id", "lon", "lat")])
  H_matrix <- as.matrix(dist(coords[, c("lon", "lat")]))
  H <- as.vector(H_matrix)
  
  # Convert H and U into long vectors
  H_vec = as.vector(H)
  DN <- D * N
  H_ind <- matrix(0, nrow = D * N, ncol = D * N)
  
  for (i in 1:(D*N)) {
    for (j in 1:(D*N)) {
      unit_i <- ((i - 1) %/% N) + 1  # outer index is COUNTY
      unit_j <- ((j - 1) %/% N) + 1
      H_ind[i, j] <- (unit_i - 1) * D + unit_j
    }
  }
  
  
  # Make matrix of the H indices 
  #H_seq <- seq_along(H_vec)
  #H_temp <- c()
  
  #for(i in seq_len(length(H_seq) / D)) {
  # H_temp <- rbind(H_temp, matrix(H_seq[i], N, N))
  #}
  
  # H_ind <- H_temp
  # for(i in seq_len(D - 1)) {
  #  H_ind <- cbind(H_ind, H_temp + i * D)
  #}
  
  # Compute temporal distance matrix U and its index matrix U_ind
  # Temporal coordinates (time points)
  
  
  # make matrix of the U indices
  # this should be an m*nxm*n matrix
  #U_seq <- 1:length(U)
  #U_seq_mat <- matrix(U_seq,N,N)
  #U_temp <- U_seq_mat
  #for(i in 1:(D-1)){
  #  U_temp <- rbind(U_temp,U_seq_mat)
  #}
  
  #U_ind <- c()
  # for(i in 1:(D)){ 
  #  U_ind <- cbind(U_ind,U_temp)
  #}
  
  #print(summary(H_ind))  # Should be in 1:(D*D)
  #print(summary(U_ind))  # Should be in 1:(N*N)
  
  list(
    N = N,
    D = D,
    n_k_f = 15,
    time_points = x_covariate,
    y = y_target_variable,
    num_treated = length(treated_indices),
    control_idx = control_indices,
    H = H,
    H_ind = H_ind,
    distance_matrix = H_matrix 
  )
}

# Function to format for poisson_RBF_NNGP.stan
format_for_stan_sep_nngp <- function(simulated_data, trt_counties, M = 10) {
  library(FNN) # For fast nearest neighbors
  
  # Order simulated_data by county_id and time
  simulated_data <- simulated_data[order(simulated_data$county_id, simulated_data$time), ]
  
  # Define treated counties
  trt_counties <- trt_counties
  
  # Identify indices for treated counties at times > 7
  treated_indices <- which(simulated_data$county_id %in% trt_counties & simulated_data$time > 7)
  
  # Determine dimensions
  N <- length(unique(simulated_data$time))
  D <- length(unique(simulated_data$county_id))
  
  # Covariates and target variable
  time_points <- sort(unique(simulated_data$time))
  y <- simulated_data$count
  
  # Control indices
  control_indices <- setdiff(seq_along(y), treated_indices)
  
  # Compute spatial distance matrix
  coords <- unique(simulated_data[, c("county_id", "lon", "lat")])
  distance_matrix <- as.matrix(dist(coords[, c("lon", "lat")]))
  
  # Compute nearest neighbors based on spatio-temporal proximity
  # Create expanded coordinates (space-time)
  spatial_coords <- coords[, c("lon", "lat")]
  st_coords <- expand.grid(time = time_points, county = 1:D)
  st_coords$lon <- spatial_coords[st_coords$county, 1]
  st_coords$lat <- spatial_coords[st_coords$county, 2]
  st_coords_scaled <- scale(st_coords[, c("time", "lon", "lat")])
  
  # Find nearest neighbors
  nn_result <- get.knn(st_coords_scaled, k = M)
  
  NN_ind <- nn_result$nn.index
  NN_num <- rep(M, nrow(NN_ind))
  
  # Adjust NN_num for any cases with fewer neighbors (optional, depending on your scenario)
  
  list(
    N = N,
    D = D,
    time_points = time_points,
    y = y,
    num_treated = length(treated_indices),
    control_idx = control_indices,
    distance_matrix = distance_matrix,
    M = M,
    NN_ind = NN_ind,
    NN_num = NN_num
  )
}

# Function to format data for poisson.stan
#' Format Data for Intrinsic Coregionalization Model (ICM) Stan Model
#'
#' This function prepares data for fitting an Intrinsic Coregionalization Model using Stan.
#' The ICM model uses latent factors to capture spatial dependencies.
#'
#' @param simulated_data Data frame containing the simulated data with columns:
#'   - county_id: Spatial location identifier
#'   - time: Time point identifier
#'   - count: Response variable
#' @param trt_counties Vector of county IDs that received treatment
#'
#' @return A list containing:
#'   - N: Number of time points
#'   - D: Number of spatial locations
#'   - n_k_f: Number of knots for basis functions
#'   - x: Vector of time points
#'   - y: Response variable vector
#'   - num_treated: Number of treated observations
#'   - control_idx: Indices of control observations
#'
#' @examples
#' # Format data for ICM model
#' stan_data <- format_for_stan_ICM(sim_data, trt_counties = c(1, 2, 3))
format_for_stan_ICM <- function(simulated_data, trt_counties) {
  # Define treated counties
  trt_counties <- trt_counties
  
  # Order simulated_data by county_id
  simulated_data <- simulated_data[order(simulated_data$county_id), ]
  
  # Identify indices for Connecticut counties (treated) and only include times after time = 7
  treated_indices <- which(simulated_data$county_id %in% trt_counties & simulated_data$time > 7)
  
  # Total number of observations
  N <- length(unique(simulated_data$time))
  D <- length(unique(simulated_data$county_id))
  
  # Covariate and target variable
  x_covariate <- unique(simulated_data$time)
  y_target_variable <- simulated_data$count
  
  # Determine control indices (non-Connecticut counties)
  control_indices <- setdiff(seq_along(y_target_variable), treated_indices)
  
  list(
    N = N,
    D = D,
    n_k_f = 15,
    x = x_covariate,
    y = y_target_variable,
    num_treated = length(treated_indices),
    control_idx = control_indices
  )
}

#' Format Data for Normal Intrinsic Coregionalization Model (ICM) Stan Model
#'
#' This function prepares data for fitting a Normal Intrinsic Coregionalization Model using Stan.
#' The ICM model uses latent factors to capture spatial dependencies, with a Normal likelihood.
#'
#' @param simulated_data Data frame containing the simulated data with columns:
#'   - county_id: Spatial location identifier
#'   - time: Time point identifier
#'   - count: Response variable
#' @param trt_counties Vector of county IDs that received treatment
#'
#' @return A list containing:
#'   - N: Number of time points
#'   - D: Number of spatial locations
#'   - n_k_f: Number of knots for basis functions
#'   - x: Vector of time points
#'   - y: Matrix of response variables (time points × locations)
#'   - num_treated: Number of treated observations
#'   - control_idx: Indices of control observations
#'
#' @examples
#' # Format data for Normal ICM model
#' stan_data <- format_for_stan_ICM_normal(sim_data, trt_counties = c(1, 2, 3))
format_for_stan_ICM_normal <- function(simulated_data, trt_counties) {
  # Define treated counties
  trt_counties <- trt_counties
  
  # Order simulated_data by county_id
  simulated_data <- simulated_data[order(simulated_data$county_id), ]
  
  # Identify indices for Connecticut counties (treated) and only include times after time = 7
  treated_indices <- which(simulated_data$county_id %in% trt_counties & simulated_data$time > 7)
  
  # Total number of observations
  N <- length(unique(simulated_data$time))
  D <- length(unique(simulated_data$county_id))
  
  # Covariate and target variable
  x_covariate <- unique(simulated_data$time)
  y_target_variable <- simulated_data$count
  
  # Determine control indices (non-Connecticut counties)
  control_indices <- setdiff(seq_along(y_target_variable), treated_indices)
  
  list(
    N = N,
    D = D,
    n_k_f = 15,
    x = x_covariate,
    y = matrix(y_target_variable, nrow=N, ncol = D, byrow = T),
    num_treated = length(treated_indices),
    control_idx = control_indices
  )
}

#' Run Simulation and Fit ICM Poisson Model
#'
#' This function runs a simulation and fits an Intrinsic Coregionalization Model (ICM)
#' with Poisson likelihood using Stan. It handles data formatting, model compilation,
#' and saving of results.
#'
#' @param seed_value Integer. Random seed for reproducibility
#' @param simulated_data Data frame containing the simulated data
#' @param kernel_choice Character. Type of kernel used for simulation
#' @param trt_counties Vector of county IDs that received treatment
#'
#' @return The fitted Stan model object
#'
#' @examples
#' # Run simulation and fit ICM Poisson model
#' fit <- run_simulation_and_stan_model_ICM_pois(
#'   seed_value = 1,
#'   simulated_data = sim_data,
#'   kernel_choice = "ICM",
#'   trt_counties = c(1, 2, 3)
#' )
run_simulation_and_stan_model_ICM_pois <- function(seed_value, simulated_data, kernel_choice, trt_counties) {
  
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_ICM(simulated_data, trt_counties)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/poisson_ICM.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f','f_samples'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0("simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "/code/simulations/results/poisson/ICM/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}


#' Run Simulation and Fit Separable Poisson Model
#'
#' This function runs a simulation and fits a separable Gaussian Process model
#' with Poisson likelihood using Stan. It handles data formatting, model compilation,
#' and saving of results.
#'
#' @param seed_value Integer. Random seed for reproducibility
#' @param simulated_data Data frame containing the simulated data
#' @param kernel_choice Character. Type of kernel used for simulation
#' @param trt_counties Vector of county IDs that received treatment
#'
#' @return The fitted Stan model object
#'
#' @examples
#' # Run simulation and fit separable Poisson model
#' fit <- run_simulation_and_stan_model_sep_pois(
#'   seed_value = 1,
#'   simulated_data = sim_data,
#'   kernel_choice = "sep",
#'   trt_counties = c(1, 2, 3)
#' )
run_simulation_and_stan_model_sep_pois <- function(seed_value, simulated_data, kernel_choice, trt_counties) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_sep(simulated_data, trt_counties)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/poisson_RBF.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f_samples'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0(wd, "simulations_smaller/data/"simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "/code/simulations/results/poisson/sep/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}

# Main function to run simulation and Stan separable model with NNGP
run_simulation_and_stan_model_sep_nngp_pois <- function(seed_value, simulated_data, kernel_choice, trt_counties, M=10) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_sep_nngp(simulated_data, trt_counties, M)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/poisson_RBF_NNGP.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f_samples', 'f_draw'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0(wd, "simulations_smaller/data/"simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "code/simulations/results/poisson/sep_NNGP/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}



# Main function to run simulation and Stan nonseparable model
#' Run Simulation and Fit Non-separable Poisson Model
#'
#' This function runs a simulation and fits a non-separable Gaussian Process model
#' with Poisson likelihood using Stan. It handles data formatting, model compilation,
#' and saving of results.
#'
#' @param seed_value Integer. Random seed for reproducibility
#' @param simulated_data Data frame containing the simulated data
#' @param kernel_choice Character. Type of kernel used for simulation
#' @param trt_counties Vector of county IDs that received treatment
#'
#' @return The fitted Stan model object
#'
#' @examples
#' # Run simulation and fit non-separable Poisson model
#' fit <- run_simulation_and_stan_model_nonsep_pois(
#'   seed_value = 1,
#'   simulated_data = sim_data,
#'   kernel_choice = "nonsep",
#'   trt_counties = c(1, 2, 3)
#' )
run_simulation_and_stan_model_nonsep_pois <- function(seed_value, simulated_data, kernel_choice, trt_counties) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_nonsep(simulated_data, trt_counties)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/poisson_nonsep.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f_samples'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0(wd, "simulations_small/data/"simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "code/simulations/results/poisson/nonsep/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}

run_simulation_and_stan_model_nonsep_v2_pois <- function(seed_value, simulated_data, kernel_choice, trt_counties) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_nonsep(simulated_data, trt_counties)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/poisson_nonsep_v2.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f_samples'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0(wd, "simulations_small/data/"simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "code/simulations/results/poisson/nonsep_v2/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}

# Main function to run simulation and Stan ICM model
run_simulation_and_stan_model_ICM_normal <- function(seed_value, simulated_data, kernel_choice, trt_counties) {
  
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_ICM_normal(simulated_data, trt_counties)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/normal_ICM.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f','f_samples'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0("simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "/code/simulations/results/normal/ICM/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}


# Main function to run simulation and Stan normal separable model
run_simulation_and_stan_model_sep_normal <- function(seed_value, simulated_data, kernel_choice, trt_counties) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_sep(simulated_data, trt_counties)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/normal_RBF.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f_samples'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0(wd, "simulations_smaller/data/"simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "/code/simulations/results/normal/sep/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}

# Main function to run simulation and Stan normal separable model with NNGP
run_simulation_and_stan_model_sep_nngp_normal <- function(seed_value, simulated_data, kernel_choice, trt_counties, M=10) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_sep_nngp(simulated_data, trt_counties, M)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/normal_RBF_NNGP.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f_samples', 'f_draw'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0(wd, "simulations_smaller/data/"simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "code/simulations/results/normal/sep_NNGP/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}

# Main function to run simulation and Stan normal nonseparable model
run_simulation_and_stan_model_nonsep_normal <- function(seed_value, simulated_data, kernel_choice, trt_counties) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_nonsep(simulated_data, trt_counties)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/normal_nonsep.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f_samples', 'f_draw', 'export_sigma_global',
                                            'export_beta', 'export_alpha', 'export_gamma',
                                            'export_lengthscale_time', 'export_lengthscale_space'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0(wd, "simulations_small/data/"simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "code/simulations/results/normal/nonsep/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}

run_simulation_and_stan_model_nonsep_v2_normal <- function(seed_value, simulated_data, kernel_choice, trt_counties) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_nonsep(simulated_data, trt_counties)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/normal_nonsep_v2.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('f_samples', 'f_draw', 'export_sigma_global',
                                            'export_beta', 'export_alpha', 'export_gamma',
                                            'export_lengthscale_time', 'export_lengthscale_space'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0(wd, "simulations_small/data/"simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "code/simulations/results/normal/nonsep_v2/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}

## Abs Percent Bias Y0 Function GSC
#' Calculate Median Absolute Percent Bias for GSC Model
#'
#' This function calculates the median absolute percent bias of the counterfactual predictions
#' from a Gaussian Synthetic Control (GSC) model compared to the true values.
#'
#' @param fit_gsc List of fitted GSC model objects
#' @param Y0_matrix List of true counterfactual matrices
#' @param ind Vector of indices to evaluate
#' @param time_ind Vector of time indices to evaluate
#'
#' @return A list containing:
#'   - median: Median absolute percent bias
#'   - LB: Lower bound of the interquartile range
#'   - UB: Upper bound of the interquartile range
#'   - IQR: Interquartile range
#'
#' @examples
#' # Calculate bias for GSC model
#' bias <- medianAbsPercBiasY0_year_gsc(
#'   fit_gsc = gsc_fits,
#'   Y0_matrix = true_values,
#'   ind = 1:10,
#'   time_ind = 8:15
#' )
medianAbsPercBiasY0_year_gsc <- function(fit_gsc, Y0_matrix, ind, time_ind) {
  medianAbsPercBias <- c()
  #for(i in c(seq(1:79),seq(81,84),seq(86,100))){
  #print(i)
  for(i in 1:length(fit_gsc)){
    # Find true treatment effect rate
    true_Y0_rate <- (Y0_matrix[[i]][ind])# Find Rate
    
    # Take average of counties for each time point
    true_Y0_mean_rate <- colMeans(matrix(true_Y0_rate , ncol = 8, byrow = T)) # take median across counties
    
    # extract the treated indices
    # Convert to rate
    mean_vec_rate <- as.vector(fit_gsc[[i]]$Y.ct)[time_ind]
    
    # For each simulation, take average counterfactual for each year over all counties 
    # should be a vector of length 8
    est_Y0 <- colMeans(matrix(mean_vec_rate, ncol = 8, byrow = T))
    
    # find absolute percent bias Y0 
    absPercBias <- abs(((est_Y0  - true_Y0_mean_rate)/true_Y0_mean_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  print(summary(medianAbsPercBias))
  # median of the simulations
  median <- median(medianAbsPercBias)
  LB <- quantile(medianAbsPercBias, probs = 0.25, na.rm = T)
  UB <- quantile(medianAbsPercBias, probs = 0.75, na.rm = T)
  IQR <- IQR(medianAbsPercBias, na.rm = T)
  
  return(list(median = median, LB = LB, UB = UB, IQR = IQR))
}

## Abs Percent Bias Y0 ASC
#' Calculate Median Absolute Percent Bias for ASC Model
#'
#' This function calculates the median absolute percent bias of the counterfactual predictions
#' from an Augmented Synthetic Control (ASC) model compared to the true values.
#'
#' @param fit_asc List of fitted ASC model objects
#' @param Y0_matrix List of true counterfactual matrices
#' @param ind Vector of indices to evaluate
#'
#' @return A list containing:
#'   - median: Median absolute percent bias
#'   - LB: Lower bound of the interquartile range
#'   - UB: Upper bound of the interquartile range
#'   - IQR: Interquartile range
#'
#' @examples
#' # Calculate bias for ASC model
#' bias <- medianAbsPercBiasY0_year_asc(
#'   fit_asc = asc_fits,
#'   Y0_matrix = true_values,
#'   ind = 1:10
#' )
medianAbsPercBiasY0_year_asc <- function(fit_asc, Y0_matrix, ind) {
  medianAbsPercBias <- c()
  #for(i in c(seq(1:79),seq(81,84),seq(86,100))){
  #print(i)
  for(i in 1:length(fit_asc)){
    # Find true treatment effect rate
    true_Y0_rate <- (Y0_matrix[[i]][ind])# Find Rate
    
    # Take average of counties for each time point
    true_Y0_mean_rate <- colMeans(matrix(true_Y0_rate , ncol = 8, byrow = T)) # take median across counties
    
    # extract the treated indices
    # Convert to rate
    mean_vec_rate <- predict(fit_asc[[i]], att = F)[8:15]
    
    # For each simulation, take average counterfactual for each year over all counties 
    # should be a vector of length 8
    est_Y0 <- colMeans(matrix(mean_vec_rate, ncol = 8, byrow = T))
    
    # find absolute percent bias Y0 
    absPercBias <- abs(((est_Y0  - true_Y0_mean_rate)/true_Y0_mean_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  
  # median of the simulations
  median <- median(medianAbsPercBias)
  LB <- quantile(medianAbsPercBias, probs = 0.25, na.rm = T)
  UB <- quantile(medianAbsPercBias, probs = 0.75, na.rm = T)
  IQR <- IQR(medianAbsPercBias)
  
  return(list(median = median, LB = LB, UB = UB, IQR = IQR))
}

# Function to format data for SpaceTimeAR.stan
#' Format Data for Markov Chain Stan Model
#'
#' This function prepares data for fitting a Markov Chain model using Stan.
#' It handles the creation of spatial and temporal adjacency matrices and
#' missing value indicators.
#'
#' @param sim_data Data frame containing the simulated data with columns:
#'   - county_id: Spatial location identifier
#'   - time: Time point identifier
#'   - count: Response variable
#' @param trt_counties Vector of county IDs that received treatment
#' @param adj_matrix Matrix indicating spatial adjacency between locations
#'
#' @return A list containing:
#'   - m: Number of time points
#'   - k: Number of latent factors
#'   - n: Number of spatial locations
#'   - nmiss: Number of missing values
#'   - n_exp: Number of exposed (treated) locations
#'   - y_miss: Matrix of missing value indicators
#'   - y: Matrix of response variables
#'   - N: Number of spatial locations
#'   - N_edges: Number of spatial edges
#'   - s_node1: First node of each spatial edge
#'   - s_node2: Second node of each spatial edge
#'   - M: Number of temporal edges
#'   - M_edges: Number of temporal edges
#'   - t_node1: First node of each temporal edge
#'   - t_node2: Second node of each temporal edge
#'
#' @examples
#' # Format data for Markov Chain model
#' stan_data <- format_for_stan_MC(
#'   sim_data = sim_data,
#'   trt_counties = c(1, 2, 3),
#'   adj_matrix = adjacency_matrix
#' )
format_for_stan_MC <- function(sim_data, trt_counties, adj_matrix) {
  library(dplyr)
  library(tidyr)
  
  # Extract dimensions
  n <- length(unique(sim_data$county_id))  # Number of counties
  m <- length(unique(sim_data$time))  # Number of time points
  
  # Reshape data into a matrix (counties in rows, time points in columns)
  y_matrix <- sim_data %>%
    select(county_id, time, count) %>%
    pivot_wider(names_from = time, values_from = count) %>%
    select(-county_id) %>%
    as.matrix()
  
  # Order simulated_data by county_id
  sim_data <- sim_data[order(sim_data$county_id), ]
  
  # Identify indices for Connecticut counties (treated) and only include times after time = 7
  treated_indices <- which(sim_data$county_id %in% trt_counties & sim_data$time > 7)
  
  # Determine control indices 
  control_indices <- setdiff(seq_along(y_matrix), treated_indices)
  
  # Create missingness indicator matrix
  y_miss <- sim_data$count
  y_miss[treated_indices] <- 1
  y_miss[control_indices] <- 0
  y_miss <- matrix(y_miss, nrow = 49, byrow = T)
  
  # Count total number of missing values
  nmiss <- sum(y_miss)
  
  # Identify exposed (treated) counties
  n_exp <- length(trt_counties)
  
  # Extract adjacency list (edges) from input adjacency matrix
  edges <- which(adj_matrix == 1, arr.ind = TRUE)
  
  # Convert to s_node format (ensure unique pairs)
  edges <- edges[edges[,1] < edges[,2], , drop = FALSE]  # Remove duplicates
  N_edges <- nrow(edges)
  s_node1 <- edges[, 1]
  s_node2 <- edges[, 2]
  
  # Temporal adjacency (each time point is adjacent to the previous one)
  M_edges <- m - 1
  t_node1 <- 1:M_edges
  t_node2 <- 2:(M_edges + 1)
  
  # Compile data list for Stan
  stan_data <- list(
    m = m,
    k = 3,  # Assuming 5 latent factors (can be adjusted)
    n = n,
    nmiss = nmiss,
    n_exp = n_exp,
    y_miss = y_miss,
    y = y_matrix,
    N = n,
    N_edges = N_edges,
    s_node1 = s_node1,
    s_node2 = s_node2,
    M = M_edges,
    M_edges = M_edges,
    t_node1 = t_node1,
    t_node2 = t_node2
  )
  
  print(stan_data)
  
  return(stan_data)
}



run_simulation_and_stan_model_MC <- function(seed_value, simulated_data, kernel_choice, trt_counties, adj_matrix) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  # Simulate data using chosen kernel and seed
  #simulated_data <- simulate_dataset(kernel_choice, seed_value)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_MC(simulated_data, trt_counties, adj_matrix)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/SpaceTimeAR.stan"),auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, data=stan_data_list, chains = 1, iter = 1000, 
                       warmup=500, pars = c('Y_pred','Mu_trt'), init_r = .1, seed = 300, save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  #saveRDS(simulated_data, file=paste0(wd, "simulations_smaller/data/"simulated_data_", kernel_choice, "_seed_", seed_value, ".rds"))
  saveRDS(stan_fit, file=paste0(wd, "simulations_grid_small_RCN_newDGP_upint_v3/results/MC/stan_fit_", kernel_choice, "_seed_", seed_value, ".rds"))
  
}

# Function to format data for nonseparable NNGP Stan model: poisson_nonsep_NNGP.stan and normal_nonsep_NNGP.stan
format_for_stan_nonsep_nngp <- function(simulated_data, trt_counties, M = 10) {
  library(FNN) # For fast nearest neighbors
  
  # Order simulated_data by county_id and time
  simulated_data <- simulated_data[order(simulated_data$county_id, simulated_data$time), ]
  
  # Define treated counties
  trt_counties <- trt_counties
  
  # Identify indices for treated counties at times > 7
  treated_indices <- which(simulated_data$county_id %in% trt_counties & simulated_data$time > 7)
  
  # Determine dimensions
  N <- length(unique(simulated_data$time))
  D <- length(unique(simulated_data$county_id))
  
  # Covariates and target variable
  time_points <- sort(unique(simulated_data$time))
  y <- simulated_data$count
  
  # Control indices
  control_indices <- setdiff(seq_along(y), treated_indices)
  
  # Compute spatial distance matrix
  coords <- unique(simulated_data[, c("county_id", "lon", "lat")])
  distance_matrix <- as.matrix(dist(coords[, c("lon", "lat")]))
  
  # Create expanded coordinates (space-time)
  spatial_coords <- coords[, c("lon", "lat")]
  st_coords <- expand.grid(time = time_points, county = 1:D)
  st_coords$lon <- spatial_coords[st_coords$county, 1]
  st_coords$lat <- spatial_coords[st_coords$county, 2]
  
  # Scale coordinates for better neighbor search
  st_coords_scaled <- scale(st_coords[, c("time", "lon", "lat")])
  
  # Find nearest neighbors using scaled coordinates
  nn_result <- get.knn(st_coords_scaled, k = M)
  
  NN_ind <- nn_result$nn.index
  NN_num <- rep(M, nrow(NN_ind))
  
  list(
    N = N,
    D = D,
    time_points = time_points,
    y = y,
    num_treated = length(treated_indices),
    control_idx = control_indices,
    distance_matrix = distance_matrix,
    M = M,
    NN_ind = NN_ind,
    NN_num = NN_num
  )
}

# Main function to run simulation with poisson nonseparable NNGP model
run_simulation_and_stan_model_nonsep_nngp_pois <- function(seed_value, simulated_data, kernel_choice, trt_counties, M=10) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_nonsep_nngp(simulated_data, trt_counties, M)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/poisson_nonsep_NNGP.stan"), auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, 
                       data = stan_data_list, 
                       chains = 1, 
                       iter = 1000, 
                       warmup = 500, 
                       pars = c('f_samples'), 
                       init_r = .1, 
                       seed = 300, 
                       save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  saveRDS(stan_fit, 
          file = paste0(wd, "code/simulations/results/poisson/nonsep_NNGP/stan_fit_", 
                        kernel_choice, "_seed_", seed_value, ".rds"))
  
  return(stan_fit)
}

# Main function to run simulation with poisson normal nonseparable NNGP model
run_simulation_and_stan_model_nonsep_nngp_normal <- function(seed_value, simulated_data, kernel_choice, trt_counties, M=10) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_nonsep_nngp(simulated_data, trt_counties, M)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/normal_nonsep_NNGP.stan"), auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, 
                       data = stan_data_list, 
                       chains = 1, 
                       iter = 1000, 
                       warmup = 500, 
                       pars = c('f_samples'), 
                       init_r = .1, 
                       seed = 300, 
                       save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  saveRDS(stan_fit, 
          file = paste0(wd, "code/simulations/results/normal/nonsep_NNGP/stan_fit_", 
                        kernel_choice, "_seed_", seed_value, ".rds"))
  
  return(stan_fit)
}



# Function to format data for ICM NNGP model
format_for_stan_ICM_nngp <- function(simulated_data, trt_counties, M = 10) {
  library(FNN) # For fast nearest neighbors
  
  # Order simulated_data by county_id and time
  simulated_data <- simulated_data[order(simulated_data$county_id, simulated_data$time), ]
  
  # Define treated counties
  trt_counties <- trt_counties
  
  # Identify indices for treated counties at times > 7
  treated_indices <- which(simulated_data$county_id %in% trt_counties & simulated_data$time > 7)
  
  # Determine dimensions
  N <- length(unique(simulated_data$time))
  D <- length(unique(simulated_data$county_id))
  J <- 5  # Number of latent processes, same as in ICM
  
  # Covariates and target variable
  time_points <- sort(unique(simulated_data$time))
  y <- simulated_data$count
  
  # Control indices
  control_indices <- setdiff(seq_along(y), treated_indices)
  
  # Compute spatial coordinates
  coords <- unique(simulated_data[, c("county_id", "lon", "lat")])
  distance_matrix <- as.matrix(dist(coords[, c("lon", "lat")]))
  
  # Create expanded coordinates (space-time)
  spatial_coords <- coords[, c("lon", "lat")]
  st_coords <- expand.grid(time = time_points, county = 1:D)
  st_coords$lon <- spatial_coords[st_coords$county, 1]
  st_coords$lat <- spatial_coords[st_coords$county, 2]
  
  # Scale coordinates for better neighbor search
  st_coords_scaled <- scale(st_coords[, c("time", "lon", "lat")])
  
  # Find nearest neighbors using scaled coordinates
  nn_result <- get.knn(st_coords_scaled, k = M)
  
  NN_ind <- nn_result$nn.index
  NN_num <- rep(M, nrow(NN_ind))
  
  list(
    N = N,
    D = D,
    J = J,
    time_points = time_points,
    y = y,
    num_treated = length(treated_indices),
    control_idx = control_indices,
    distance_matrix = distance_matrix,
    M = M,
    NN_ind = NN_ind,
    NN_num = NN_num
  )
}

# Main function to run simulation and Stan ICM NNGP model
run_simulation_and_stan_model_ICM_nngp_pois <- function(seed_value, simulated_data, kernel_choice, trt_counties, M=10) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_ICM_nngp(simulated_data, trt_counties, M)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/poisson_ICM_NNGP.stan"), auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, 
                       data = stan_data_list, 
                       chains = 1, 
                       iter = 1000, 
                       warmup = 500, 
                       pars = c('f_samples'), 
                       init_r = .1, 
                       seed = 300, 
                       save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  saveRDS(stan_fit, 
          file = paste0(wd, "/code/simulations/results/poisson/ICM_NNGP/stan_fit_", 
                        kernel_choice, "_seed_", seed_value, ".rds"))
  
  return(stan_fit)
}


run_simulation_and_stan_model_ICM_nngp_normal <- function(seed_value, simulated_data, kernel_choice, trt_counties, M=10) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = FALSE)
  
  # Format data for Stan model
  stan_data_list <- format_for_stan_ICM_nngp(simulated_data, trt_counties, M)
  
  # Load Stan model
  stan_model <- stan_model(paste0(wd,"/code/STAN/normal_ICM_NNGP.stan"), auto_write=F)
  
  # Run Stan model 
  print("Begin Sampling")
  stan_fit <- sampling(stan_model, 
                       data = stan_data_list, 
                       chains = 1, 
                       iter = 1000, 
                       warmup = 500, 
                       pars = c('f_samples'), 
                       init_r = .1, 
                       seed = 300, 
                       save_warmup = FALSE)
  print("Finished Sampling")
  
  # Save results to files
  saveRDS(stan_fit, 
          file = paste0(wd, "/code/simulations/results/normal/ICM_NNGP/stan_fit_", 
                        kernel_choice, "_seed_", seed_value, ".rds"))
  
  return(stan_fit)
}

#' Fit Gaussian Process Model with Stan
#'
#' This function fits a Gaussian Process model using Stan, supporting both Normal and Poisson likelihoods.
#' The model can use different kernel structures (nonseparable, separable, or ICM) and includes
#' spatial and temporal lengthscale parameters.
#'
#' @param data Data frame containing the response variable and covariates
#' @param kernel Character. Type of kernel to use:
#'   - "nonsep": Non-separable spatiotemporal kernel
#'   - "sep": Separable spatiotemporal kernel
#'   - "ICM": Intrinsic Coregionalization Model
#' @param distribution Character. Type of distribution:
#'   - "poisson": Poisson distribution with log link
#'   - "normal": Normal distribution with identity link
#' @param spatial_lengthscale Numeric. Initial value for spatial lengthscale parameter
#' @param temporal_lengthscale Numeric. Initial value for temporal lengthscale parameter
#' @param J Integer. Number of latent factors for ICM kernel (default: 5)
#' @param chains Integer. Number of MCMC chains (default: 4)
#' @param iter Integer. Number of iterations per chain (default: 2000)
#' @param warmup Integer. Number of warmup iterations per chain (default: 1000)
#' @param cores Integer. Number of CPU cores to use (default: 4)
#' @param adapt_delta Numeric. Target acceptance rate for HMC (default: 0.8)
#' @param max_treedepth Integer. Maximum tree depth for NUTS (default: 10)
#' 
#' @return A list containing:
#'   - model: The fitted Stan model
#'   - data: The data used for fitting
#'   - kernel: The kernel type used
#'   - distribution: The distribution type used
#'   - parameters: A list of model parameters
#'
#' @examples
#' # Fit GP model with nonseparable kernel and Poisson likelihood
#' fit <- fit_gp(data = sim_data, kernel = "nonsep", distribution = "poisson")
#' 
#' # Fit GP model with separable kernel and Normal likelihood
#' fit <- fit_gp(data = sim_data, kernel = "sep", distribution = "normal")
fit_gp <- function(data, kernel, distribution, spatial_lengthscale = .3, temporal_lengthscale = .9, J = 5,
                   chains = 4, iter = 2000, warmup = 1000, cores = 4, adapt_delta = 0.8, max_treedepth = 10) {
  # Load required libraries
  library(rstan)
  library(MASS)
  
  # Step 1: Data preparation
  # Extract unique spatial locations and create distance matrix
  unique_locations <- unique(data[, c("lon", "lat")])
  dist_matrix <- as.matrix(dist(unique_locations))
  
  # Create temporal distance matrix
  timepoints <- unique(data$time)
  temporal_dist_matrix <- as.matrix(dist(timepoints))
  
  # Step 2: Prepare Stan data
  stan_data <- list(
    N = nrow(data),
    N_unit = nrow(unique_locations),
    N_time = length(timepoints),
    y = data$count,
    unit_id = data$county_id,
    time_id = data$time,
    dist_matrix = dist_matrix,
    temporal_dist_matrix = temporal_dist_matrix,
    spatial_lengthscale = spatial_lengthscale,
    temporal_lengthscale = temporal_lengthscale,
    J = J
  )
  
  # Step 3: Select appropriate Stan model based on kernel and distribution
  if (kernel == "nonsep") {
    if (distribution == "poisson") {
      model_file <- "models/gp_nonsep_poisson.stan"
    } else if (distribution == "normal") {
      model_file <- "models/gp_nonsep_normal.stan"
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
  } else if (kernel == "sep") {
    if (distribution == "poisson") {
      model_file <- "models/gp_sep_poisson.stan"
    } else if (distribution == "normal") {
      model_file <- "models/gp_sep_normal.stan"
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
  } else if (kernel == "ICM") {
    if (distribution == "poisson") {
      model_file <- "models/gp_icm_poisson.stan"
    } else if (distribution == "normal") {
      model_file <- "models/gp_icm_normal.stan"
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
  } else {
    stop("Invalid kernel type. Choose from 'nonsep', 'sep', or 'ICM'.")
  }
  
  # Step 4: Compile and fit Stan model
  # Set Stan options for better sampling
  options(mc.cores = cores)
  rstan_options(auto_write = TRUE)
  
  # Compile model
  model <- stan_model(model_file)
  
  # Fit model with specified parameters
  fit <- sampling(
    model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  )
  
  # Step 5: Return results
  return(list(
    model = fit,
    data = data,
    kernel = kernel,
    distribution = distribution,
    parameters = list(
      spatial_lengthscale = spatial_lengthscale,
      temporal_lengthscale = temporal_lengthscale,
      J = J
    )
  ))
}

#' Fit Nearest Neighbor Gaussian Process Model with Stan
#'
#' This function fits a Nearest Neighbor Gaussian Process (NNGP) model using Stan.
#' NNGP is a computationally efficient approximation of full GP that uses a sparse
#' precision matrix based on nearest neighbors. The model supports both Normal and
#' Poisson likelihoods and can use different kernel structures.
#'
#' @param data Data frame containing the response variable and covariates
#' @param kernel Character. Type of kernel to use:
#'   - "nonsep": Non-separable spatiotemporal kernel
#'   - "sep": Separable spatiotemporal kernel
#'   - "ICM": Intrinsic Coregionalization Model
#' @param distribution Character. Type of distribution:
#'   - "poisson": Poisson distribution with log link
#'   - "normal": Normal distribution with identity link
#' @param spatial_lengthscale Numeric. Initial value for spatial lengthscale parameter
#' @param temporal_lengthscale Numeric. Initial value for temporal lengthscale parameter
#' @param J Integer. Number of latent factors for ICM kernel (default: 5)
#' @param m Integer. Number of nearest neighbors to use (default: 10)
#' @param chains Integer. Number of MCMC chains (default: 4)
#' @param iter Integer. Number of iterations per chain (default: 2000)
#' @param warmup Integer. Number of warmup iterations per chain (default: 1000)
#' @param cores Integer. Number of CPU cores to use (default: 4)
#' @param adapt_delta Numeric. Target acceptance rate for HMC (default: 0.8)
#' @param max_treedepth Integer. Maximum tree depth for NUTS (default: 10)
#' 
#' @return A list containing:
#'   - model: The fitted Stan model
#'   - data: The data used for fitting
#'   - kernel: The kernel type used
#'   - distribution: The distribution type used
#'   - parameters: A list of model parameters
#'   - neighbors: The nearest neighbor structure used
#'
#' @examples
#' # Fit NNGP model with nonseparable kernel and Poisson likelihood
#' fit <- fit_nngp(data = sim_data, kernel = "nonsep", distribution = "poisson")
#' 
#' # Fit NNGP model with separable kernel and Normal likelihood
#' fit <- fit_nngp(data = sim_data, kernel = "sep", distribution = "normal")
fit_nngp <- function(data, kernel, distribution, spatial_lengthscale = .3, temporal_lengthscale = .9, J = 5, m = 10,
                     chains = 4, iter = 2000, warmup = 1000, cores = 4, adapt_delta = 0.8, max_treedepth = 10) {
  # Load required libraries
  library(rstan)
  library(MASS)
  library(FNN)
  
  # Step 1: Data preparation
  # Extract unique spatial locations and create distance matrix
  unique_locations <- unique(data[, c("lon", "lat")])
  dist_matrix <- as.matrix(dist(unique_locations))
  
  # Create temporal distance matrix
  timepoints <- unique(data$time)
  temporal_dist_matrix <- as.matrix(dist(timepoints))
  
  # Step 2: Find nearest neighbors for each location
  # For each location, find m nearest neighbors based on spatial distance
  nn_indices <- matrix(0, nrow = nrow(unique_locations), ncol = m)
  nn_distances <- matrix(0, nrow = nrow(unique_locations), ncol = m)
  
  for (i in 1:nrow(unique_locations)) {
    # Get distances to all other locations
    distances <- dist_matrix[i, ]
    # Find m nearest neighbors (excluding self)
    nn <- order(distances)[2:(m + 1)]
    nn_indices[i, ] <- nn
    nn_distances[i, ] <- distances[nn]
  }
  
  # Step 3: Prepare Stan data
  stan_data <- list(
    N = nrow(data),
    N_unit = nrow(unique_locations),
    N_time = length(timepoints),
    y = data$count,
    unit_id = data$county_id,
    time_id = data$time,
    dist_matrix = dist_matrix,
    temporal_dist_matrix = temporal_dist_matrix,
    spatial_lengthscale = spatial_lengthscale,
    temporal_lengthscale = temporal_lengthscale,
    J = J,
    m = m,
    nn_indices = nn_indices,
    nn_distances = nn_distances
  )
  
  # Step 4: Select appropriate Stan model based on kernel and distribution
  if (kernel == "nonsep") {
    if (distribution == "poisson") {
      model_file <- "models/nngp_nonsep_poisson.stan"
    } else if (distribution == "normal") {
      model_file <- "models/nngp_nonsep_normal.stan"
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
  } else if (kernel == "sep") {
    if (distribution == "poisson") {
      model_file <- "models/nngp_sep_poisson.stan"
    } else if (distribution == "normal") {
      model_file <- "models/nngp_sep_normal.stan"
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
  } else if (kernel == "ICM") {
    if (distribution == "poisson") {
      model_file <- "models/nngp_icm_poisson.stan"
    } else if (distribution == "normal") {
      model_file <- "models/nngp_icm_normal.stan"
    } else {
      stop("Invalid distribution. Choose from 'poisson' or 'normal'.")
    }
  } else {
    stop("Invalid kernel type. Choose from 'nonsep', 'sep', or 'ICM'.")
  }
  
  # Step 5: Compile and fit Stan model
  # Set Stan options for better sampling
  options(mc.cores = cores)
  rstan_options(auto_write = TRUE)
  
  # Compile model
  model <- stan_model(model_file)
  
  # Fit model with specified parameters
  fit <- sampling(
    model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  )
  
  # Step 6: Return results
  return(list(
    model = fit,
    data = data,
    kernel = kernel,
    distribution = distribution,
    parameters = list(
      spatial_lengthscale = spatial_lengthscale,
      temporal_lengthscale = temporal_lengthscale,
      J = J,
      m = m
    ),
    neighbors = list(
      indices = nn_indices,
      distances = nn_distances
    )
  ))
}






