/**
 * Normal RBF (Radial Basis Function) Spatiotemporal Gaussian Process Model
 * 
 * This Stan model implements a spatiotemporal Gaussian Process model using separable
 * RBF (squared exponential) kernels for spatial and temporal components. The model 
 * is designed for analyzing continuous spatiotemporal data with treatment effects.
 * 
 * Model Structure:
 * - Separable spatiotemporal covariance using RBF kernels
 * - Normal likelihood for continuous observations
 * - Unit-specific offsets
 * - Control-only fitting for causal inference
 * - Kronecker product structure for computational efficiency
 * 
 * Parameters:
 * - lengthscale_time: Temporal correlation length
 * - lengthscale_space: Spatial correlation length
 * - sigma_global: Overall variance parameter
 * - sigma_f: Scale parameter for latent functions
 * - state_offset: Unit-specific offsets
 * - sigma_obs: Observation noise standard deviation
 * - intercept: Global intercept
 * - f_vec: GP function values
 * 
 * Data Requirements:
 * - N: Number of time points
 * - D: Number of spatial units
 * - n_k_f: Number of latent functions
 * - time_points: Vector of time points
 * - y: Response vector (continuous)
 * - num_treated: Number of treated units
 * - control_idx: Indices of control units
 * - H: Spatial distance matrix (vectorized)
 * - distance_matrix: Precomputed spatial distances
 */

functions {
  /**
   * Compute Kronecker product of two matrices
   * 
   * This function computes the Kronecker product of matrices A and B,
   * which is essential for constructing separable spatiotemporal covariance.
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
        int row_start = (i - 1) * p + 1;
        int row_end = (i - 1) * p + p;
        int col_start = (j - 1) * q + 1;
        int col_end = (j - 1) * q + q;
        
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
    
    return C;
  }
}

data {
  int<lower=1> N;        // Number of time points
  int<lower=1> D;        // Number of spatial units
  int<lower=1> n_k_f;    // Number of latent functions
  vector[N] time_points; // Time points
  
  real y[N * D];        // Response vector (continuous)
  int num_treated;      // Number of treated units
  int control_idx[N * D - num_treated];  // Control unit indices
  
  // Spatial distances
  real H[D * D];        // Spatial distance matrix (vectorized)
  int H_ind[D*N, D*N];  // Spatial distance indices
  matrix[D,D] distance_matrix;  // Precomputed spatial distances
}

transformed data {
  real sigma_intercept = 0.1;  // Prior scale for intercept
}

parameters {
  real<lower=0> lengthscale_time;    // Temporal correlation length
  real<lower=0> lengthscale_space;   // Spatial correlation length
  real<lower=0> sigma_global;        // Overall variance
  row_vector[D] state_offset;        // Unit-specific offsets
  real<lower=0> sigma_f;             // Scale of latent functions
  real intercept;                    // Global intercept
  vector[N * D] f_vec;               // GP function values
  real<lower=0> sigma_obs;           // Observation noise std
}

transformed parameters {
  matrix[N * D, N * D] K_global;  // Global covariance matrix
  matrix[N * D, N * D] L_global;  // Cholesky factor
  
  // Compute temporal covariance (RBF kernel)
  real time_points_vec[N] = to_array_1d(time_points);
  matrix[N, N] K_time = gp_exp_quad_cov(
    time_points_vec, 
    sigma_global, 
    lengthscale_time
  );
  
  // Compute spatial covariance (RBF kernel)
  matrix[D, D] K_space = sigma_global^2 * 
    exp((-.5 * distance_matrix) / square(lengthscale_space));
}

model {
  // Prior distributions
  lengthscale_time ~ inv_gamma(5,5);
  lengthscale_space ~ inv_gamma(5,5);
  sigma_f ~ std_normal();
  sigma_global ~ std_normal();
  state_offset ~ std_normal();
  sigma_obs ~ normal(0, 1);
  
  // GP prior using Kronecker structure
  f_vec ~ multi_normal_prec(
    rep_vector(0, N * D),
    kronecker_product(inverse(K_space), inverse(K_time))
  );
  
  // Likelihood (control units only)
  y[control_idx] ~ normal(intercept + f_vec[control_idx], sigma_obs);
}

generated quantities {
  matrix[N, D] f_samples;   // Posterior predictive samples
  vector[N * D] log_lik;    // Log-likelihood values
  matrix[N, D] f_draw;      // Reshaped GP values
  
  // Export parameter values for post-processing
  real export_sigma_global = sigma_global;
  real export_lengthscale_time = lengthscale_time;
  real export_lengthscale_space = lengthscale_space;
  
  {
    // Generate posterior predictive samples
    for (d in 1:D) {
      for (n in 1:N) {
        int idx = (d - 1) * N + n;
        f_samples[n, d] = normal_rng(intercept + f_vec[idx], sigma_obs);
      }
    }
    
    // Reshape GP values for export
    f_draw = to_matrix(f_vec, D, N)';
  }
}

