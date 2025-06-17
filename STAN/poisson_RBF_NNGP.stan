/**
 * Poisson RBF-NNGP (Radial Basis Function with Nearest Neighbor Gaussian Process)
 * 
 * This Stan model implements a spatiotemporal Gaussian Process model using separable
 * RBF (squared exponential) kernels combined with a Nearest Neighbor Gaussian Process
 * (NNGP) approximation for scalability. The model is designed for analyzing count data
 * with spatiotemporal dependence and treatment effects.
 * 
 * Model Structure:
 * - Separable RBF kernels for spatial and temporal components
 * - NNGP approximation for computational efficiency
 * - Poisson likelihood for count data
 * - Treatment effect analysis capability
 * - Sparse neighbor-based structure
 * 
 * Parameters:
 * - lengthscale_time: Temporal correlation length
 * - lengthscale_space: Spatial correlation length
 * - sigma_global: Overall variance parameter
 * - intercept: Global intercept
 * - f_vec: Latent GP function values
 * 
 * Data Requirements:
 * - N: Number of time points
 * - D: Number of spatial units
 * - time_points: Vector of time points
 * - distance_matrix: Spatial distances
 * - y: Count observations
 * - num_treated: Number of treated units
 * - control_idx: Indices of control observations
 * - M: Number of nearest neighbors
 * - NN_ind: Nearest neighbor indices
 * - NN_num: Number of neighbors per point
 */

functions {
  /**
   * Compute separable spatiotemporal covariance using RBF kernels
   * 
   * This function implements a separable covariance function that
   * combines RBF kernels for both spatial and temporal components.
   * 
   * @param dist_space Spatial distance
   * @param dist_time Temporal distance
   * @param sigma_sq Overall variance
   * @param rho_space Spatial lengthscale
   * @param rho_time Temporal lengthscale
   * @return Covariance value
   */
  real st_cov(real dist_space, real dist_time, real sigma_sq, 
              real rho_space, real rho_time) {
    return sigma_sq * 
           exp(-0.5 * square(dist_space / rho_space)) * 
           exp(-0.5 * square(dist_time / rho_time));
  }
}

data {
  int<lower=1> N;            // Number of time points
  int<lower=1> D;            // Number of spatial units
  vector[N] time_points;     // Time points
  matrix[D, D] distance_matrix; // Spatial distances
  int<lower=0> y[N * D];     // Count observations
  int num_treated;           // Number of treated units
  int control_idx[N * D - num_treated]; // Control indices
  
  int<lower=1> M;           // Maximum neighbors
  int<lower=0> NN_ind[N*D, M]; // Neighbor indices
  int<lower=0> NN_num[N*D];    // Neighbors per point
}

parameters {
  real<lower=0> lengthscale_time;   // Temporal correlation
  real<lower=0> lengthscale_space;  // Spatial correlation
  real<lower=0> sigma_global;       // Overall variance
  real intercept;                   // Global intercept
  vector[N*D] f_vec;               // GP function values
}

model {
  // Prior distributions
  lengthscale_time ~ inv_gamma(5, 5);
  lengthscale_space ~ inv_gamma(5, 5);
  sigma_global ~ std_normal();
  
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
        real dist_space = distance_matrix[di, dj];
        real dist_time = fabs(time_points[ti] - time_points[tj]);
        
        // Compute covariance with neighbor j
        cov_neighbors[j] = st_cov(dist_space, dist_time, 
                                square(sigma_global), 
                                lengthscale_space, 
                                lengthscale_time);
        
        // Compute neighbor-to-neighbor covariances
        for (k in 1:NN_num[i]) {
          int neighbor_idx_k = NN_ind[i, k];
          int tk = 1 + ((neighbor_idx_k - 1) % N);
          int dk = 1 + ((neighbor_idx_k - 1) / N);
          
          real dist_space_kk = distance_matrix[dj, dk];
          real dist_time_kk = fabs(time_points[tj] - time_points[tk]);
          
          C_NN[j, k] = st_cov(dist_space_kk, dist_time_kk, 
                             square(sigma_global), 
                             lengthscale_space, 
                             lengthscale_time);
          if (j == k)
            C_NN[j, k] += 1e-8;  // Numerical stability
        }
        f_neighbors[j] = f_vec[neighbor_idx];
      }
      
      // Compute conditional distribution
      weights = cov_neighbors' / C_NN;
      mean_cond = dot_product(weights, f_neighbors);
      var_cond = square(sigma_global) - dot_product(weights, cov_neighbors);
      
      // Sample from conditional
      f_vec[i] ~ normal(mean_cond, sqrt(var_cond));
    }
  }
  
  // Likelihood for control observations
  y[control_idx] ~ poisson_log(intercept + f_vec[control_idx]);
}

generated quantities {
  matrix[N, D] f_samples;  // Posterior predictive samples
  matrix[N, D] f_draw;     // Reshaped GP values
  
  // Generate posterior samples and reshape GP values
  for (d in 1:D) {
    for (n in 1:N) {
      int idx = (d - 1) * N + n;
      f_samples[n, d] = poisson_log_rng(intercept + f_vec[idx]);
      f_draw[n, d] = f_vec[idx];
    }
  }
}
