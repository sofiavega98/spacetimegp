# ---------------------------------------------------------------------
# Script: view_kernel_3d.R
# Purpose: Visualize 3D surfaces of different spatiotemporal kernels (Separable, Nonseparable, ICM)
#          for comparison and illustration in the Gaussian Process project
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Set Kernel Parameters
# ---------------------------------------------------------------------
s_lengthscale <- 0.7      # Spatial lengthscale for separable/nonseparable kernels
t_lengthscale <- 0.9      # Temporal lengthscale for all kernels
a <- .9                   # Gneiting kernel parameter
d <- 2                    # Spatial dimension
c <- 3                    # Gneiting kernel parameter
alpha <- 1                # Gneiting kernel parameter
beta <- 1                 # Gneiting kernel parameter
gamma <- 1                # Gneiting kernel parameter
sigma_sep <- 1            # Scale parameter for separable kernel
sigma_nonsep <- 1.0       # Scale parameter for nonseparable kernel
sigma_icm <- 1            # Scale parameter for ICM kernel
J <- 5                    # Number of latent processes for ICM

# Set seed for reproducibility
set.seed(1)

# ---------------------------------------------------------------------
# Generate ICM Kernel Matrix
# ---------------------------------------------------------------------
n_points <- 50  # Number of points in each dimension for the grid
beta_matrix <- matrix(rnorm(n_points * J), nrow = n_points, ncol = J)
K_unit <- beta_matrix %*% t(beta_matrix)

# ---------------------------------------------------------------------
# Create Distance Grids
# ---------------------------------------------------------------------
spatial_dist <- seq(0, 1, length.out = n_points)    # Normalized spatial distances
temporal_dist <- seq(0, 1, length.out = n_points)   # Normalized temporal distances

# Initialize kernel value matrices
K_sep <- K_nonsep <- K_icm <- matrix(0, nrow = length(spatial_dist), ncol = length(temporal_dist))

# ---------------------------------------------------------------------
# Compute Kernel Values for Each Grid Point
# ---------------------------------------------------------------------
for (i in seq_along(spatial_dist)) {
  for (j in seq_along(temporal_dist)) {
    h <- spatial_dist[i]
    u <- temporal_dist[j]
    
    # Separable kernel: product of spatial and temporal RBFs
    K_sep[i,j] <- sigma_sep^2 * exp(-0.5 * (h / s_lengthscale)^2) * 
      exp(-0.5 * (u / t_lengthscale)^2)
    
    # Nonseparable (Gneiting) kernel
    psi_u <- (a * u^(2 * alpha) + 1)
    space_term <- exp(-c * h^(2 * gamma) / psi_u^(beta * gamma))
    time_term <- 1 / psi_u^(beta * d / 2)
    K_nonsep[i, j] <- sigma_nonsep^2 * time_term * space_term
    
    # ICM kernel: latent process structure with temporal RBF
    K_icm[i,j] <- sigma_icm^2 * K_unit[i,j] * exp(-0.5 * (u / t_lengthscale)^2)
  }
}

# Normalize kernel values for visualization
K_sep <- K_sep / max(K_sep)
K_nonsep <- K_nonsep / max(K_nonsep)
K_icm <- K_icm / max(K_icm)

# ---------------------------------------------------------------------
# Plot 3D Kernel Surfaces and Save to PDF and PNG
# ---------------------------------------------------------------------
# Create figures directory if it doesn't exist
dir.create("figures", showWarnings = FALSE)

# Set up the plotting device for PDF
pdf("figures/kernel_comparison_3d_static.pdf", width = 15, height = 5)
par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))  # 1 row, 3 columns

# Common perspective parameters for all plots
theta <- 30  # Rotation angle around z
phi <- 30    # Viewing angle
expand <- 0.5 # Expansion factor for z-axis
ticktype <- "detailed"
nticks <- 5

# Plot ICM Kernel
persp(spatial_dist, temporal_dist, K_icm,
      main = "ICM Kernel",
      xlab = "Spatial Distance",
      ylab = "Temporal Distance",
      zlab = "Kernel Value",
      theta = theta, phi = phi,
      expand = expand,
      ticktype = ticktype,
      nticks = nticks,
      col = "lightpink",
      border = NA,
      shade = 0.5)

# Plot Separable Kernel
persp(spatial_dist, temporal_dist, K_sep,
      main = "Separable Kernel",
      xlab = "Spatial Distance",
      ylab = "Temporal Distance",
      zlab = "Kernel Value",
      theta = theta, phi = phi,
      expand = expand,
      ticktype = ticktype,
      nticks = nticks,
      col = "lightblue",
      border = NA,
      shade = 0.5)

# Plot Nonseparable Kernel
persp(spatial_dist, temporal_dist, K_nonsep,
      main = "Nonseparable Kernel",
      xlab = "Spatial Distance",
      ylab = "Temporal Distance",
      zlab = "Kernel Value",
      theta = theta, phi = phi,
      expand = expand,
      ticktype = ticktype,
      nticks = nticks,
      col = "lightgreen",
      border = NA,
      shade = 0.5)

dev.off()  # Close the PDF device

# Also save as PNG for easy viewing
png("figures/kernel_comparison_3d_static.png", width = 1500, height = 500, res = 100)
par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))

# Plot Separable Kernel
persp(spatial_dist, temporal_dist, K_sep,
      main = "Separable Kernel",
      xlab = "Spatial Distance",
      ylab = "Temporal Distance",
      zlab = "Kernel Value",
      theta = theta, phi = phi,
      expand = expand,
      ticktype = ticktype,
      nticks = nticks,
      col = "lightblue",
      border = NA,
      shade = 0.5)

# Plot Nonseparable Kernel
persp(spatial_dist, temporal_dist, K_nonsep,
      main = "Nonseparable Kernel",
      xlab = "Spatial Distance",
      ylab = "Temporal Distance",
      zlab = "Kernel Value",
      theta = theta, phi = phi,
      expand = expand,
      ticktype = ticktype,
      nticks = nticks,
      col = "lightgreen",
      border = NA,
      shade = 0.5)

# Plot ICM Kernel
persp(spatial_dist, temporal_dist, K_icm,
      main = "ICM Kernel",
      xlab = "Spatial Distance",
      ylab = "Temporal Distance",
      zlab = "Kernel Value",
      theta = theta, phi = phi,
      expand = expand,
      ticktype = ticktype,
      nticks = nticks,
      col = "lightpink",
      border = NA,
      shade = 0.5)

dev.off()  # Close the PNG device

# ---------------------------------------------------------------------
# High-Resolution PNG Export for Publication
# ---------------------------------------------------------------------
library(viridis)  # For color palettes
library(ragg)     # For high-res PNG output

#' Helper function to create a color matrix for surface plots
#' @param z Matrix of values
#' @param palette Color palette function
#' @return Matrix of colors
get_col_matrix <- function(z, palette = colorRampPalette(c("blue", "red"))) {
  z_scaled <- (z - min(z)) / (max(z) - min(z))
  color_index <- round(z_scaled * 99) + 1
  matrix(palette(100)[color_index], nrow = nrow(z), ncol = ncol(z))
}

# Common perspective parameters
theta <- 30
phi <- 30
expand <- 0.5
ticktype <- "detailed"
nticks <- 5

# Export to high-resolution PNG for publication
ragg::agg_png("~/Library/Mobile Documents/com~apple~CloudDocs/Harvard/GPs and CI/code/kernel_demo/figures/fig1_kernel_surfaces.png", width = 4100, height = 1200, res = 300)
par(mfrow = c(1, 3), bg = "white", mar = c(3, 3, 2,2))

# Plot ICM Kernel
persp(spatial_dist, temporal_dist, K_icm,
      col = "darkblue",
      border = "gray80",
      lwd = 0.5,
      shade = .0001,
      main = "Separable ICM-RBF Kernel",
      xlab = "Spatial Distance",
      ylab = "Temporal Distance",
      zlab = "Kernel Value",
      theta = theta, phi = phi,
      expand = expand,
      ticktype = ticktype,
      nticks = nticks,
      cex.lab = 1,
      cex.axis = .8,
      cex.main = 1.4)

# Plot Separable Kernel
persp(spatial_dist, temporal_dist, K_sep,
      col = "darkblue",
      border = "gray80",
      lwd = 0.5,   
      shade = .0001,
      main = "Separable RBF-RBF Kernel",
      xlab = "Spatial Distance",
      ylab = "Temporal Distance",
      zlab = "Kernel Value",
      theta = theta, phi = phi,
      expand = expand,
      ticktype = ticktype,
      nticks = nticks,
      cex.lab = 1,
      cex.axis = .8,
      cex.main = 1.4)

# Plot Nonseparable Kernel
persp(spatial_dist, temporal_dist, K_nonsep,
      col = "darkblue",
      border = "gray80",
      lwd = 0.5,
      shade = .0001,
      main = "Nonseparable Gneiting Kernel",
      xlab = "Spatial Distance",
      ylab = "Temporal Distance",
      zlab = "Kernel Value",
      theta = theta, phi = phi,
      expand = expand,
      ticktype = ticktype,
      nticks = nticks,
      cex.lab = 1,
      cex.axis = .8,
      cex.main = 1.4)

dev.off()
