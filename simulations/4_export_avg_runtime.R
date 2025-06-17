library(dplyr)
library(readr)
library(tidyr)
library(knitr)
library(kableExtra)

# Read data
run_times_normal <- read_csv("code/simulations/runtime_logs/run_times_normal.csv", show_col_types = FALSE) %>%
  mutate(sim_type = "Normal")

run_times_poisson <- read_csv("code/simulations/runtime_logs/run_times_poisson.csv", show_col_types = FALSE) %>%
  mutate(sim_type = "Poisson")

# Combine and clean names
all_run_times <- bind_rows(run_times_normal, run_times_poisson) %>%
  mutate(
    kernel = case_when(
      kernel == "ICM" ~ "Separable ICM-RBF",
      kernel == "sep" ~ "Separable RBF-RBF",
      kernel == "nonsep" ~ "Nonseparable Gneiting",
      TRUE ~ kernel
    ),
    model = case_when(
      model == "ICM" ~ "Separable ICM-RBF",
      model == "sep" ~ "Separable RBF-RBF",
      model == "sep_NNGP" ~ "Separable RBF-RBF NNGP",
      model == "nonsep" ~ "Nonseparable Gneiting",
      model == "nonsep_NNGP" ~ "Nonseparable Gneiting NNGP",
      TRUE ~ model
    )
  )

# Summarize and pivot wider
runtime_table <- all_run_times %>%
  group_by(kernel, model, sim_type) %>%
  summarise(Average_Runtime = mean(run_time_minutes, na.rm = TRUE), .groups = "drop") %>%
  mutate(Average_Runtime = round(Average_Runtime, 2)) %>%
  pivot_wider(
    names_from = sim_type,
    values_from = Average_Runtime
  ) %>%
  arrange(kernel, model) %>%
  rename(
    `DGP Kernel` = kernel,
    `Model` = model
  )

# Save as LaTeX table
runtime_table %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    align = "l l r r",
    caption = "Average runtime (in minutes) for each kernel-model combination, for Normal and Poisson simulations."
  ) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  cat(file = "code/simulations/runtime_logs/average_runtime_table.tex")
