/**
 * Poisson RBF (Radial Basis Function) Spatiotemporal Model
 * 
 * This Stan model implements a spatiotemporal Gaussian Process model using separable
 * RBF (squared exponential) kernels for analyzing count data with treatment effects.
 * The model uses Kronecker product structure for computational efficiency in handling
 * spatiotemporal covariance.
 * 
 * Model Structure:
 * - Separable RBF kernels for spatial and temporal components
 * - Kronecker product structure for efficient computation
 * - Poisson likelihood for count data
 * - Treatment effect analysis capability
 * - State-specific offsets
 * 
 * Parameters:
 * - lengthscale_time: Temporal correlation length
 * - lengthscale_space: Spatial correlation length
 * - sigma_global: Overall variance parameter
 * - state_offset: Unit-specific offsets
 * - sigma_f: Scale parameter for latent functions
 * - intercept: Global intercept
 * - f_vec: Latent GP function values
 * 
 * Data Requirements:
 * - N: Number of time points
 * - D: Number of spatial units
 * - n_k_f: Number of latent functions
 * - time_points: Vector of time points
 * - y: Count observations
 * - num_treated: Number of treated units
 * - control_idx: Indices of control observations
 * - H: Precomputed distance matrix (vectorized)
 * - H_ind: Matrix of indices for H
 * - distance_matrix: Spatial distances
 */

functions {
  /**
   * Compute Kronecker product of two matrices
   * 
   * This function implements the Kronecker product operation,
   * which is essential for efficient spatiotemporal covariance
   * computation.
   * 
   * @param A First matrix
   * @param B Second matrix
   * @return Kronecker product matrix
   */
  matrix kronecker_product(matrix A, matrix B) {
    matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
    int m = rows(A);
    int n = cols(A);
    int p = rows(B);
    int q = cols(B);
    
    for (i in 1:m) {
      for (j in 1:n) {
        // Define row and column ranges
        int row_start = (i - 1) * p + 1;
        int row_end = (i - 1) * p + p;
        int col_start = (j - 1) * q + 1;
        int col_end = (j - 1) * q + q;
        
        // Assign submatrix
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
    return C;
  }
}

data {
  int<lower=1> N;                     // Number of time points
  int<lower=1> D;                     // Number of spatial units
  int<lower=1> n_k_f;                 // Number of latent functions
  vector[N] time_points;              // Time points
  int<lower=0> y[N * D];              // Count observations
  int num_treated;                    // Number of treated units
  int control_idx[N * D - num_treated]; // Control indices
  
  // Spatial components
  real H[D * D];                      // Vectorized distance matrix
  int H_ind[D*N, D*N];                // Distance matrix indices
  matrix[D,D] distance_matrix;         // Spatial distances
}

transformed data {
  real sigma_intercept = 0.1;         // Prior scale for intercept
}

parameters {
  real<lower=0> lengthscale_time;     // Temporal correlation
  real<lower=0> lengthscale_space;    // Spatial correlation
  real<lower=0> sigma_global;         // Overall variance
  row_vector[D] state_offset;         // Unit-specific offsets
  real<lower=0> sigma_f;              // Scale for latent functions
  real intercept;                     // Global intercept
  vector[N * D] f_vec;                // GP function values
}

transformed parameters {
  matrix[N * D, N * D] K_global;      // Global covariance
  matrix[N * D, N * D] L_global;      // Cholesky decomposition
  
  // Compute temporal covariance
  real time_points_vec[N] = to_array_1d(time_points);
  matrix[N, N] K_time = gp_exp_quad_cov(time_points_vec, 
                                       sigma_global, 
                                       lengthscale_time);
  
  // Compute spatial covariance
  matrix[D, D] K_space = sigma_global^2 * 
                        exp((-.5 * distance_matrix) / 
                            square(lengthscale_space));
}

model {
  // Prior distributions
  lengthscale_time ~ inv_gamma(5,5);
  lengthscale_space ~ inv_gamma(5,5);
  sigma_f ~ std_normal();
  sigma_global ~ std_normal();
  state_offset ~ std_normal();
  
  // GP prior using Kronecker precision
  f_vec ~ multi_normal_prec(rep_vector(0, N * D),
                           kronecker_product(inverse(K_space),
                                          inverse(K_time)));
  
  // Likelihood for control observations
  y[control_idx] ~ poisson_log(intercept + f_vec[control_idx]);
}

generated quantities {
  matrix[N, D] f_samples;             // Posterior predictive samples
  vector[N * D] log_lik;              // Log-likelihood values
  matrix[N, D] f_draw;                // Reshaped GP values
  
  // Export parameters for post-processing
  real export_sigma_global = sigma_global;
  real export_lengthscale_time = lengthscale_time;
  real export_lengthscale_space = lengthscale_space;
  
  // Generate posterior samples and reshape GP values
  for (d in 1:D) {
    for (n in 1:N) {
      int idx = (d - 1) * N + n;
      f_samples[n, d] = poisson_log_rng(intercept + f_vec[idx]);
    }
  }
  
  // Reshape f_vec to [N x D] matrix
  f_draw = to_matrix(f_vec, D, N)';
}

