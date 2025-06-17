/**
 * Normal Non-Separable Spatiotemporal Gaussian Process Model
 * 
 * This Stan model implements a spatiotemporal Gaussian Process model using Gneiting's 
 * non-separable kernel function for continuous data with a Normal likelihood. The model 
 * is designed for analyzing spatiotemporal data with treatment effects.
 * 
 * Model Structure:
 * - Non-separable spatiotemporal covariance using Gneiting's kernel
 * - Normal likelihood for continuous observations
 * - Unit-specific offsets
 * - Control-only fitting for causal inference
 * 
 * Parameters:
 * - lengthscale_time: Temporal correlation length
 * - lengthscale_space: Spatial correlation length
 * - sigma_global: Overall variance parameter
 * - beta: Interaction parameter in Gneiting's kernel
 * - alpha: Temporal smoothness parameter (0,1)
 * - gamma: Spatial smoothness parameter (0,1)
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
 * - U: Temporal distance matrix (vectorized)
 * - H_ind, U_ind: Index matrices for H and U
 */

functions {
  /**
   * Compute Gneiting's non-separable spatiotemporal kernel
   * 
   * This function implements the covariance function from:
   * Gneiting (2002) "Nonseparable, Stationary Covariance Functions for Space-Time Data"
   * 
   * @param H Vector of spatial distances
   * @param U Vector of temporal distances
   * @param H_ind Matrix of spatial distance indices
   * @param U_ind Matrix of temporal distance indices
   * @param sigma2 Overall variance
   * @param beta Space-time interaction
   * @param alpha Temporal smoothness
   * @param gamma Spatial smoothness
   * @param a Temporal scale
   * @param c Spatial scale
   * @param delta Numerical stability constant
   * @param n Number of spatial locations
   * @param m Number of time points
   * @return Covariance matrix
   */
  matrix spacetime_nonsep_kernel(real[] H, real[] U, int[,] H_ind, int[,] U_ind, 
                               real sigma2, real beta, real alpha, 
                               real gamma, real a, real c, real delta, int n, int m) {
    real denom;
    matrix[n*m,n*m] K;
    
    // Build covariance matrix
    for(i in 1:((n*m)-1)) {
      K[i,i] = sigma2 + delta;  // Add jitter to diagonal
      for(j in (i+1):(n*m)) {
        real h;  // Spatial distance
        real u;  // Temporal distance
        int h_ind;
        int u_ind;
        
        // Get distances from index matrices
        h_ind = H_ind[i,j];
        u_ind = U_ind[i,j];
        h = H[h_ind];
        u = U[u_ind];
        
        // Compute Gneiting kernel
        denom = a * pow(abs(u),(2*alpha))+1;
        K[i, j] = sigma2/pow(denom,(beta/2)) * 
                  exp(- c * pow(abs(h),2*gamma)/pow(denom,(beta*gamma)));
        K[j, i] = K[i, j];  // Symmetric
      }
    }
    K[n*m, n*m] = sigma2 + delta;
    
    return K;
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
  
  // Temporal distances
  real U[N*N];          // Temporal distance matrix (vectorized)
  int U_ind[D*N, D*N];  // Temporal distance indices
}

transformed data {
  real sigma_intercept = 0.1;  // Prior scale for intercept
  real delta = 1e-9;          // Numerical stability constant
}

parameters {
  real<lower=0> lengthscale_time;    // Temporal correlation length
  real<lower=0> lengthscale_space;   // Spatial correlation length
  real<lower=0> sigma_global;        // Overall variance
  real<lower=0> beta;                // Space-time interaction
  real<lower=0,upper=1> alpha;       // Temporal smoothness
  real<lower=0,upper=1> gamma;       // Spatial smoothness
  row_vector[D] state_offset;        // Unit-specific offsets
  real<lower=0> sigma_obs;           // Observation noise std
  real intercept;                    // Global intercept
  vector[N * D] f_vec;               // GP function values
}

transformed parameters {
  matrix[N * D, N * D] K_global;  // Global covariance matrix
  matrix[N * D, N * D] L_global;  // Cholesky factor
  matrix[D,N] f_mat;              // Reshaped GP values
  
  // Compute covariance using Gneiting kernel
  K_global = spacetime_nonsep_kernel(
    to_array_1d(H), U, H_ind, U_ind,
    sigma_global, beta, alpha, gamma,
    lengthscale_time, lengthscale_space,
    delta, D, N
  );
  
  // Reshape GP values to matrix form
  f_mat = to_matrix(f_vec, D,N);
}

model {
  // Prior distributions
  lengthscale_time ~ inv_gamma(5,5);
  lengthscale_space ~ inv_gamma(5,5);
  sigma_global ~ std_normal();
  sigma_obs ~ normal(0, 1);
  beta ~ beta(1, 1);
  state_offset ~ std_normal();
  
  // GP prior
  f_vec ~ multi_normal(rep_vector(0, N * D), K_global);
  
  // Likelihood (control units only)
  y[control_idx] ~ normal(intercept + f_vec[control_idx], sigma_obs);
}

generated quantities {
  matrix[N, D] f_samples;   // Posterior predictive samples
  vector[N * D] log_lik;    // Log-likelihood values
  matrix[N, D] f_draw;      // Reshaped GP values
  
  // Export parameter values for post-processing
  real<lower=0> export_sigma_global = sigma_global;
  real<lower=0> export_beta = beta;
  real<lower=0, upper=1> export_alpha = alpha;
  real<lower=0, upper=1> export_gamma = gamma;
  real<lower=0> export_lengthscale_time = lengthscale_time;
  real<lower=0> export_lengthscale_space = lengthscale_space;
  
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


