/**
 * Normal ICM (Intrinsic Coregionalization Model)
 * 
 * This Stan model implements a multivariate Gaussian Process model using an ICM
 * (Intrinsic Coregionalization Model) for analyzing continuous spatiotemporal data
 * with treatment effects. The model uses a combination of global and unit-specific
 * components with RBF kernels.
 * 
 * Model Structure:
 * - ICM kernel for multivariate spatiotemporal modeling
 * - Combination of global and unit-specific GPs
 * - Normal likelihood for continuous observations
 * - Treatment effect analysis capability
 * - State-specific offsets
 * 
 * Parameters:
 * - lengthscale_global: Global GP lengthscale
 * - sigma_global: Global GP variance
 * - lengthscale_f: Unit-specific GP lengthscale
 * - sigma_f: Unit-specific GP variance
 * - sigman: Observation noise standard deviation
 * - state_offset: Unit-specific offsets
 * - z_global, z_f: GP random effects
 * - k_f: ICM mixing coefficients
 * - global_offset: Global intercept
 * 
 * Data Requirements:
 * - N: Number of observations (time points)
 * - D: Number of units (spatial locations)
 * - n_k_f: Number of latent functions
 * - x: Univariate covariate (e.g., time)
 * - y: Response matrix [N x D]
 * - num_treated: Number of treated units
 * - control_idx: Indices of control observations
 */

data {
  int<lower=1> N;          // Number of observations
  int<lower=1> D;          // Number of units
  int<lower=1> n_k_f;      // Number of latent functions
  vector[N] x;             // Univariate covariate
  matrix[N, D] y;          // Response matrix
  int num_treated;         // Number of treated units
  int control_idx[N * D - num_treated]; // Control indices
}

transformed data {
  // Normalize covariate for numerical stability
  real xmean = mean(x);
  real xsd = sd(x);
  real xn[N] = to_array_1d((x - xmean)/xsd);
  real sigma_intercept = 0.1;
  vector[N] jitter = rep_vector(1e-9, N);  // Numerical stability term
}

parameters {
  real<lower=0> lengthscale_global;  // Global GP lengthscale
  real<lower=0> sigma_global;        // Global GP variance
  real<lower=0> lengthscale_f;       // Unit-specific GP lengthscale
  real<lower=0> sigma_f;             // Unit-specific GP variance
  real<lower=0> sigman;              // Observation noise std
  vector[D] state_offset;            // Unit-specific offsets
  vector[N] z_global;                // Global GP random effects
  matrix[N, n_k_f] z_f;             // Unit-specific GP random effects
  matrix[n_k_f, D] k_f;             // ICM mixing coefficients
  real global_offset;                // Global intercept
}

model {
  // Compute covariance matrices and their Cholesky decompositions
  matrix[N, N] K_f = gp_exp_quad_cov(xn, sigma_f, lengthscale_f);
  matrix[N, N] L_f = cholesky_decompose(add_diag(K_f, jitter));
  
  matrix[N, N] K_global = gp_exp_quad_cov(xn, sigma_global, lengthscale_global);
  matrix[N, N] L_global = cholesky_decompose(add_diag(K_global, jitter));
  
  // Prior distributions
  to_vector(z_f) ~ std_normal();
  to_vector(k_f) ~ std_normal();
  z_global ~ std_normal();
  lengthscale_f ~ inv_gamma(5,5);
  lengthscale_global ~ inv_gamma(5,5);
  sigma_f ~ normal(0, 1);
  sigma_global ~ normal(0, 1);
  sigman ~ normal(0, 1);
  state_offset ~ normal(0, 1);
  
  // Likelihood for control observations
  to_vector(y)[control_idx] ~ normal(
      to_vector(
        global_offset +                    // Global intercept
        rep_matrix(state_offset, N)' +     // Unit-specific offsets
        rep_matrix(L_global * z_global, D) + // Global GP component
        L_f * z_f * k_f                    // Unit-specific GP component
      )[control_idx],
      sigman
  );
}

generated quantities {
  matrix[N, D] f;              // Latent function values
  matrix[N, D] f_samples;      // Posterior predictive samples
  vector[N * D - num_treated] log_lik;  // Log-likelihood for WAIC/LOO
  
  // Recompute covariance matrices for predictions
  matrix[N, N] K_f = gp_exp_quad_cov(xn, sigma_f, lengthscale_f);
  matrix[N, N] L_f = cholesky_decompose(add_diag(K_f, jitter));
  
  matrix[N, N] K_global = gp_exp_quad_cov(xn, sigma_global, lengthscale_global);
  matrix[N, N] L_global = cholesky_decompose(add_diag(K_global, jitter));
  
  // Compute latent function values
  f = global_offset + 
      rep_matrix(state_offset, N)' + 
      rep_matrix(L_global * z_global, D) + 
      L_f * z_f * k_f;
  
  // Compute log-likelihood for model comparison
  for(n in 1:(N * D - num_treated)) {
    log_lik[n] = normal_lpdf(
      to_vector(y)[control_idx[n]] |
      to_vector(f)[control_idx[n]],
      sigman
    );
  }
  
  // Generate posterior predictive samples
  for(i in 1:N) {
    f_samples[i] = to_row_vector(normal_rng(f[i], sigman));
  }
}
