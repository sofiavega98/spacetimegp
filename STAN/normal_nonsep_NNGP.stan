/**
 * Normal Non-Separable NNGP (Nearest Neighbor Gaussian Process) Model
 * 
 * This Stan model implements a spatiotemporal Gaussian Process model using Gneiting's
 * non-separable covariance function combined with a Nearest Neighbor Gaussian Process
 * (NNGP) approximation for scalability. The model is designed for analyzing continuous
 * spatiotemporal data with treatment effects.
 * 
 * Model Structure:
 * - Non-separable spatiotemporal covariance using Gneiting's kernel
 * - NNGP approximation for computational efficiency
 * - Normal likelihood for continuous observations
 * - Treatment effect analysis capability
 * - Sparse neighbor-based structure
 * 
 * Parameters:
 * - sigma_global: Overall standard deviation
 * - beta: Space-time interaction parameter
 * - alpha: Temporal smoothness parameter
 * - gamma: Spatial smoothness parameter
 * - a: Temporal scale parameter
 * - c: Spatial scale parameter
 * - intercept: Global intercept
 * - f_vec: Latent GP function values
 * - sigma_obs: Observation noise standard deviation
 * 
 * Data Requirements:
 * - N: Number of time points
 * - D: Number of spatial units
 * - time_points: Vector of time points
 * - y: Response vector (continuous)
 * - num_treated: Number of treated units
 * - control_idx: Indices of control observations
 * - distance_matrix: Spatial distances
 * - M: Number of nearest neighbors
 * - NN_ind: Nearest neighbor indices
 * - NN_num: Number of neighbors per point
 */

functions {
  /**
   * Compute Gneiting's non-separable spatiotemporal covariance
   * 
   * This function implements the covariance function from:
   * Gneiting (2002) "Nonseparable, Stationary Covariance Functions for Space-Time Data"
   * 
   * @param h Spatial distance
   * @param u Temporal distance
   * @param sigma2 Overall variance
   * @param beta Space-time interaction
   * @param alpha Temporal smoothness
   * @param gamma Spatial smoothness
   * @param a Temporal scale
   * @param c Spatial scale
   * @return Covariance value
   */
  real st_cov_nonsep(real h, real u, real sigma2, real beta, real alpha,
                     real gamma, real a, real c) {
    real denom = a * pow(fabs(u), 2*alpha) + 1;
    real cov = sigma2 / pow(denom, beta/2) * 
               exp(-c * pow(fabs(h), 2*gamma) / pow(denom, beta*gamma));
    return cov;
  }
}

data {
  int<lower=1> N;                     // Number of time points
  int<lower=1> D;                     // Number of spatial units
  vector[N] time_points;              // Time points
  real y[N*D];                        // Observations (continuous)
  int<lower=1> num_treated;           // Number of treated units
  int<lower=1> control_idx[N*D-num_treated]; // Control indices
  matrix[D,D] distance_matrix;        // Spatial distances
  int<lower=1> M;                     // Number of neighbors
  int<lower=1> NN_ind[N*D, M];       // Neighbor indices
  int<lower=1> NN_num[N*D];          // Neighbors per point
}

parameters {
  real<lower=0.1,upper=10> sigma_global; // Overall variance
  real<lower=0.1,upper=5> beta;          // Space-time interaction
  real<lower=0.1,upper=5> alpha;         // Temporal smoothness
  real<lower=0.1,upper=5> gamma;         // Spatial smoothness
  real<lower=0.1,upper=10> a;            // Temporal scale
  real<lower=0.1,upper=10> c;            // Spatial scale
  real intercept;                        // Global intercept
  vector[N*D] f_vec;                     // GP function values
  real<lower=0> sigma_obs;               // Observation noise std
}

model {
  // Prior distributions
  sigma_global ~ inv_gamma(5, 5);
  beta ~ inv_gamma(5, 5);
  alpha ~ inv_gamma(5, 5);
  gamma ~ inv_gamma(5, 5);
  a ~ inv_gamma(5, 5);
  c ~ inv_gamma(5, 5);
  intercept ~ normal(0, 1);
  sigma_obs ~ normal(0, 1);
  
  // NNGP construction
  for (i in 1:(N*D)) {
    if (NN_num[i] == 0) {
      // No neighbors case: Independent GP
      f_vec[i] ~ normal(0, sigma_global);
    } else {
      // Initialize neighbor computations
      vector[NN_num[i]] cov_neighbors;
      matrix[NN_num[i], NN_num[i]] C_NN;
      vector[NN_num[i]] f_neighbors;
      row_vector[NN_num[i]] weights;
      real mean_cond;
      real var_cond;
      
      // Extract temporal and spatial indices
      int ti = 1 + ((i - 1) % N);
      int di = 1 + ((i - 1) / N);
      
      // Compute covariances with neighbors
      for (j in 1:NN_num[i]) {
        int neighbor_idx = NN_ind[i, j];
        int tj = 1 + ((neighbor_idx - 1) % N);
        int dj = 1 + ((neighbor_idx - 1) / N);
        
        // Compute distances
        real h = distance_matrix[di, dj];
        real u = fabs(time_points[ti] - time_points[tj]);
        
        // Compute covariance with neighbor j
        cov_neighbors[j] = st_cov_nonsep(h, u, square(sigma_global), 
                                       beta, alpha, gamma, a, c);
        
        // Compute neighbor-to-neighbor covariances
        for (k in 1:NN_num[i]) {
          int neighbor_idx_k = NN_ind[i, k];
          int tk = 1 + ((neighbor_idx_k - 1) % N);
          int dk = 1 + ((neighbor_idx_k - 1) / N);
          
          real h_kk = distance_matrix[dj, dk];
          real u_kk = fabs(time_points[tj] - time_points[tk]);
          
          C_NN[j, k] = st_cov_nonsep(h_kk, u_kk, square(sigma_global), 
                                    beta, alpha, gamma, a, c);
          if (j == k)
            C_NN[j, k] += 1e-8;  // Numerical stability
        }
        f_neighbors[j] = f_vec[neighbor_idx];
      }
      
      // Compute conditional distribution
      weights = cov_neighbors' / C_NN;
      mean_cond = dot_product(weights, f_neighbors);
      var_cond = square(sigma_global) - dot_product(weights, cov_neighbors);
      var_cond = var_cond <= 0 ? 1e-8 : var_cond;  // Ensure positive variance
      
      // Sample from conditional
      f_vec[i] ~ normal(mean_cond, sqrt(var_cond));
    }
  }
  
  // Likelihood for control observations
  y[control_idx] ~ normal(intercept + f_vec[control_idx], sigma_obs);
}

generated quantities {
  matrix[N, D] f_samples;  // Posterior predictive samples
  matrix[N, D] f_draw;     // Reshaped GP values
  
  // Generate posterior samples and reshape GP values
  for (d in 1:D) {
    for (n in 1:N) {
      int idx = (d - 1) * N + n;
      f_samples[n, d] = normal_rng(intercept + f_vec[idx], sigma_obs);
    }
  }
}
