/**
 * Poisson ICM-NNGP (Intrinsic Coregionalization Model with Nearest Neighbor Gaussian Process)
 * 
 * This Stan model implements a spatiotemporal Gaussian Process model using an ICM kernel
 * combined with a Nearest Neighbor Gaussian Process (NNGP) approximation for scalability.
 * The model is designed for analyzing count data with spatiotemporal dependence and
 * treatment effects.
 * 
 * Model Structure:
 * - ICM kernel for multivariate spatiotemporal modeling
 * - NNGP approximation for computational efficiency
 * - Poisson likelihood for count data
 * - Treatment effect analysis capability
 * - Sparse neighbor-based structure
 * 
 * Parameters:
 * - beta: Spatial random effects matrix (D x J)
 * - sigma: Overall variance parameter
 * - temporal_lengthscale: Temporal correlation length
 * - intercept: Global intercept
 * - f: Latent GP function values
 * 
 * Data Requirements:
 * - N: Number of time points
 * - D: Number of spatial units
 * - J: Number of latent factors
 * - time_points: Vector of time points
 * - y: Count observations
 * - num_treated: Number of treated units
 * - control_idx: Indices of control observations
 * - distance_matrix: Spatial distances
 * - M: Number of nearest neighbors
 * - NN_ind: Nearest neighbor indices
 * - NN_num: Number of neighbors per point
 */

functions {
  /**
   * Compute ICM (Intrinsic Coregionalization Model) kernel
   * 
   * This function implements a separable kernel that combines spatial
   * correlation through a linear model of coregionalization with
   * temporal correlation through an RBF kernel.
   * 
   * @param beta Matrix of spatial random effects
   * @param sigma Overall variance parameter
   * @param temporal_lengthscale Temporal correlation length
   * @param di First spatial index
   * @param dj Second spatial index
   * @param time_diff Time difference between points
   * @return Kernel value
   */
  real icm_kernel(matrix beta, real sigma, real temporal_lengthscale,
                 int di, int dj, real time_diff) {
    return sigma * dot_product(beta[di,], beta[dj,]) * 
           exp(-0.5 * square(time_diff / temporal_lengthscale));
  }
}

data {
  int<lower=1> N;                     // Number of time points
  int<lower=1> D;                     // Number of spatial units
  int<lower=1> J;                     // Number of latent factors
  vector[N] time_points;              // Time points
  int<lower=0> y[N*D];               // Count observations
  int<lower=1> num_treated;           // Number of treated units
  int<lower=1> control_idx[N*D-num_treated]; // Control indices
  matrix[D,D] distance_matrix;        // Spatial distances
  int<lower=1> M;                     // Number of neighbors
  int<lower=1> NN_ind[N*D, M];       // Neighbor indices
  int<lower=1> NN_num[N*D];          // Neighbors per point
}

parameters {
  matrix[D,J] beta;                   // Spatial random effects
  real<lower=0> sigma;                // Overall variance
  real<lower=0> temporal_lengthscale; // Temporal correlation
  real intercept;                     // Global intercept
  vector[N*D] f;                      // GP function values
}

model {
  // Prior distributions
  for (i in 1:D) {
    beta[i] ~ normal(0, 1);
  }
  sigma ~ inv_gamma(5, 5);
  temporal_lengthscale ~ inv_gamma(5, 5);
  intercept ~ normal(0, 1);
  
  // NNGP construction
  for (i in 1:(N*D)) {
    if (NN_num[i] == 0) {
      // No neighbors case: Independent GP
      f[i] ~ normal(0, sigma);
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
        
        // Compute time difference
        real time_diff = time_points[ti] - time_points[tj];
        
        // Compute covariance with neighbor j
        cov_neighbors[j] = icm_kernel(beta, sigma, temporal_lengthscale, 
                                    di, dj, time_diff);
        
        // Compute neighbor-to-neighbor covariances
        for (k in 1:NN_num[i]) {
          int neighbor_idx_k = NN_ind[i, k];
          int tk = 1 + ((neighbor_idx_k - 1) % N);
          int dk = 1 + ((neighbor_idx_k - 1) / N);
          
          real time_diff_kk = time_points[tj] - time_points[tk];
          
          C_NN[j, k] = icm_kernel(beta, sigma, temporal_lengthscale, 
                                 dj, dk, time_diff_kk);
          if (j == k)
            C_NN[j, k] += 1e-8;  // Numerical stability
        }
        f_neighbors[j] = f[neighbor_idx];
      }
      
      // Compute conditional distribution
      weights = cov_neighbors' / C_NN;
      mean_cond = dot_product(weights, f_neighbors);
      var_cond = sigma - dot_product(weights, cov_neighbors);
      
      // Sample from conditional
      f[i] ~ normal(mean_cond, sqrt(var_cond));
    }
  }
  
  // Likelihood for control observations
  y[control_idx] ~ poisson_log(intercept + f[control_idx]);
}

generated quantities {
  vector[num_treated] f_samples;    // Posterior samples for treated
  matrix[N,D] f_draw;              // Reshaped GP values
  
  // Export parameters for post-processing
  real export_sigma = sigma;
  real export_temporal_lengthscale = temporal_lengthscale;
  
  // Generate posterior samples for treated units
  for (i in 1:num_treated) {
    f_samples[i] = poisson_log_rng(intercept + f[i]);
  }
  
  // Reshape GP values to [N x D] matrix
  for (i in 1:N) {
    for (j in 1:D) {
      f_draw[i,j] = f[(j-1)*N + i];
    }
  }
} 

