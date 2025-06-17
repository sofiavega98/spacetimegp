# Practical considerations for Gaussian Process modeling for causal inference quasi-experimental studies with panel data

This repository contains R code and analysis for research work on Gaussian Process modeling for causal inference quasi-experimental studies with panel data. The project includes various components for simulation, case studies, and applications using R and Stan.

## Project Structure

### `kernel_demo/`: Kernel Function Demonstrations
This directory contains code for demonstrating and analyzing kernel functions and their properties:
- `0_functions.R`: Core functions for kernel operations and calculations
- `1_create_sim_data.R`: Data generation for kernel demonstrations
- `2_calculate_weights.R`: Weight calculations for kernel functions
- `3_analyze_unit6_weights.R`: Analysis of unit weights
- `view_kernel_3d.R`: 3D visualization of kernel functions
- `est_functional_boxplots.R`: Functional boxplot estimation

### `simulations/`: Simulation Code
This directory contains code for running various simulations:
- `0_functions.R`: Core simulation functions
- `1_create_sim_data.R`: Data generation for simulations
- `2_run_simulations_normal.R`: Normal distribution simulations
- `2_run_simulations_poisson.R`: Poisson distribution simulations
- `3_export_results_normal.R`: Export normal simulation results
- `3_export_results_poisson.R`: Export Poisson simulation results
- `plot_counts_paper.R`: Visualization of simulation results
- `plot_DGP.R`: Data generating process visualization

### `STAN/`: Stan Model Implementations
This directory contains various Stan model implementations:
- Normal models:
  - `normal_nonsep.stan`: Non-separable normal model
  - `normal_RBF.stan`: Radial Basis Function normal model
  - `normal_ICM.stan`: Intrinsic Coregionalization Model
  - `normal_nonsep_NNGP.stan`: Non-separable NNGP normal model
  - `normal_RBF_NNGP.stan`: RBF NNGP normal model
  - `normal_ICM_NNGP.stan`: ICM NNGP normal model
- Poisson models:
  - `poisson_nonsep.stan`: Non-separable Poisson model
  - `poisson_RBF.stan`: RBF Poisson model
  - `poisson_ICM.stan`: ICM Poisson model
  - `poisson_nonsep_NNGP.stan`: Non-separable NNGP Poisson model
  - `poisson_RBF_NNGP.stan`: RBF NNGP Poisson model
  - `poisson_ICM_NNGP.stan`: ICM NNGP Poisson model
- Other:
  - `SpaceTimeAR.stan`: Space-time autoregressive model

## Setup Instructions

1. Clone this repository:
```bash
git clone [repository-url]
cd [repository-name]
```

2. Install R (version 4.0.0 or later) from [CRAN](https://cran.r-project.org/)

3. Install required R packages:
```R
install.packages(c("tidyverse", "rstan", "ggplot2", "dplyr", "tidyr", "purrr", "magrittr", "knitr", "rmarkdown"))
```

4. Install Stan:
```R
install.packages("rstan")
```

## Usage Instructions

### Kernel Demonstrations
1. Start with `0_functions.R` to understand the core kernel functions
2. Run `1_create_sim_data.R` to generate demonstration data
3. Use `2_calculate_weights.R` to compute kernel weights
4. Analyze results with `3_analyze_unit6_weights.R`
5. Visualize kernels in 3D using `view_kernel_3d.R`

### Running Simulations
1. Begin with `0_functions.R` to set up simulation functions
2. Generate simulation data using `1_create_sim_data.R`
3. Run simulations:
   - For normal distributions: `2_run_simulations_normal.R`
   - For Poisson distributions: `2_run_simulations_poisson.R`
4. Export results:
   - Normal results: `3_export_results_normal.R`
   - Poisson results: `3_export_results_poisson.R`
5. Visualize results using `plot_counts_paper.R`

### Using Stan Models
1. Choose the appropriate model from the STAN directory based on your needs:
   - Normal models for continuous data
   - Poisson models for count data
   - NNGP variants for large datasets
2. Compile the model in R:
```R
library(rstan)
model <- stan_model("STAN/your_model.stan")
```
3. Prepare your data according to the model's requirements
4. Run the model:
```R
fit <- sampling(model, data = your_data, chains = 4, iter = 2000)
```

## Dependencies

The project requires the following main dependencies:
- R 4.0.0+
- rstan
- tidyverse
- ggplot2
- dplyr
- tidyr
- purrr
- magrittr
- knitr
- rmarkdown

## Contributing

[Guidelines for contributing to the project]

## License

[License information]

## Contact

[Contact information] 