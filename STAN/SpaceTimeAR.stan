/**
 * Spatiotemporal AR(1) Factor Model
 * 
 * This Stan model implements a spatiotemporal factor model with an AR(1) prior
 * on factor loadings and an ICAR (Intrinsic Conditional Autoregressive) prior
 * on factor scores. The model is designed for analyzing count data with missing
 * values and treatment effects.
 * 
 * Model Structure:
 * - Low-rank factor decomposition for spatiotemporal structure
 * - AR(1) prior on factor loadings for temporal smoothness
 * - ICAR prior on factor scores for spatial smoothness
 * - Negative binomial likelihood for overdispersed counts
 * - Treatment effect analysis capability
 * - Handles missing data
 * 
 * Parameters:
 * - alpha: Global intercept
 * - d0: County-specific deviations
 * - c0: Time-specific deviations
 * - FS: Factor scores matrix (spatial components)
 * - L: Factor loadings matrix (temporal components)
 * - log_phi: Negative binomial dispersion
 * - alpha_AR: AR(1) intercept
 * - beta_AR: AR(1) coefficient
 * - sigma_AR: AR(1) innovation variance
 * - tau_ICAR: ICAR precision parameter
 * 
 * Data Requirements:
 * - m: Number of time points
 * - k: Number of latent factors
 * - n: Number of spatial units
 * - nmiss: Total number of missing values
 * - n_exp: Number of exposed units
 * - y_miss: Missing data indicators
 * - y: Count observations
 * - s_node1, s_node2: Spatial adjacency edges
 * - t_node1, t_node2: Temporal adjacency edges
 */

data {
  int<lower = 0> m;                   // Number of time points
  int<lower = 0> k;                   // Number of latent factors
  int<lower = 0> n;                   // Number of spatial units
  int<lower = 0> nmiss;               // Number of missing values
  int<lower = 0> n_exp;               // Number of exposed units
  int<lower=0, upper=1> y_miss[n,m];  // Missing data indicators
  int<lower=0> y[n,m];                // Count observations
  
  // Spatial adjacency structure
  int<lower=0> N;                     // Number of spatial nodes
  int<lower=0> N_edges;               // Number of spatial edges
  int<lower=1,upper=N> s_node1[N_edges];  // Spatial edge pairs
  int<lower=1, upper=N> s_node2[N_edges];
  
  // Temporal adjacency structure
  int<lower=0> M;                     // Number of temporal nodes
  int<lower=0> M_edges;               // Number of temporal edges
  int<lower=1,upper=N> t_node1[M_edges];  // Temporal edge pairs
  int<lower=1, upper=N> t_node2[M_edges];
}

parameters {
  real alpha;                         // Global intercept
  vector[n] d0;                       // Spatial random effects
  vector[m] c0;                       // Temporal random effects
  matrix[n,k] FS;                     // Factor scores (spatial)
  matrix[m, k] L;                     // Factor loadings (temporal)
  real<lower=0> log_phi;              // NB dispersion
  real alpha_AR;                      // AR(1) intercept
  real beta_AR;                       // AR(1) coefficient
  real<lower=0> sigma_AR;             // AR(1) innovation std
  real<lower=0> tau_ICAR;             // ICAR precision
}

transformed parameters {
  matrix[n,m] Ups;                    // Linear predictor
  matrix<lower=0>[n,m] Mu;            // Expected counts
  
  // Compute linear predictor
  Ups = FS * L';
  
  // Transform to response scale
  for(j in 1:m) {
    for (i in 1:n) {
      Mu[i,j] = exp(alpha + d0[i] + c0[j] + Ups[i,j]);
    }
  }
}

model {
  // Prior for ICAR precision
  tau_ICAR ~ gamma(1, .01);
  
  // Prior for NB dispersion
  log_phi ~ normal(0, 1);
  
  // ICAR prior on factor scores
  for (i in 1:k) {
    target += (n * 0.5) * log(tau_ICAR) - 
              (0.5 * tau_ICAR) * 
              dot_self((FS[s_node1,i])' - (FS[s_node2,i])');
    sum(FS[,i]) ~ normal(0, 0.01 * N);  // Soft sum-to-zero constraint
  }
  
  // AR(1) prior on factor loadings
  for (i in 1:k) {
    L[1,i] ~ normal(0, 1);  // Initial state
    L[2:m,i] ~ normal(alpha_AR + beta_AR * L[1:(m-1),i], sigma_AR);
  }
  
  // Likelihood with missing data handling
  for(i in 1:n) {
    for (j in 1:m) {
      if (1-y_miss[i,j])  // Only include non-missing observations
        y[i,j] ~ neg_binomial_2(Mu[i,j], exp(log_phi));
    }
  }
}

generated quantities {
  int Y_pred[n_exp*m];              // Predictions for treated units
  real<lower=0> Mu_trt[n_exp*m];    // Expected values for treated
  
  {
    int idy = 0;
    int idz = 0;
    
    // Generate predictions for treated units
    for (i in 1:n) {
      if (sum(y_miss[i,]) > 0) {  // Check if unit is treated
        // Store expected values
        Mu_trt[(idz*m+1):(idz*m+m)] = to_array_1d(Mu[i,]);
        idz = idz + 1;
        
        // Generate predictions
        for (j in 1:m) {
          if (log(Mu[i,j]) > 20) {  // Numerical stability check
            Y_pred[idy+1] = -1;      // Flag extreme values
          } else {
            Y_pred[idy+1] = neg_binomial_2_rng(Mu[i,j], exp(log_phi));
          }
          idy = idy + 1;
        }
      }
    }
  }
}




