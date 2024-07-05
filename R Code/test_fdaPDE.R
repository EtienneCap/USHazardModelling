# Install and load necessary packages
library(fdaPDE)
library(tidyverse)


# Generate synthetic spatiotemporal data
set.seed(123)
n <- 100
x <- runif(n, 0, 10)
y <- runif(n, 0, 10)
t <- runif(n, 0, 10)
data <- data.frame(x, y, t)

boundary_x = c(seq(0,10, length.out = 10), rep(10,10), seq(10, 0, length.out = 10), rep(0, 10))
boundary_y = c(rep(0, 10), seq(0,10, length.out = 10), rep(10,10), seq(10, 0, length.out = 10))

boundaries = data.frame(x = boundary_x, y = boundary_y) %>% distinct()
boundary_segments = matrix(c(seq(1, nrow(boundaries)), c(seq(2, nrow(boundaries)), 1)), ncol = 2)

# Create a mesh for the spatial domain
mesh <- create.mesh.2D(nodes = rbind(boundaries, data[c("x", "y")]), segments = boundary_segments)
plot(mesh)

# Create a finite element basis using the mesh
basisobj <- create.FEM.basis(mesh)

# Define time range and discretization
time_range <- range(data$t)
time_steps <- seq(time_range[1], time_range[2], length.out = 20)

# Define smoothing parameters
lambda_s <- 1e-4  # Spatial regularization parameter
lambda_t <- 1e-4  # Temporal regularization parameter

# Prepare spatiotemporal data
spatiotemporal_data <- data.frame(x = data$x, y = data$y, z = rep(1, n), t = data$t)  # z can be a constant if you don't have a third spatial dimension

observations = matrix(1, nrow = nrow(spatiotemporal_data), ncol = length(time_steps))

# Perform the spatiotemporal smoothing
result <- smooth.FEM.time(observations = observations, FEMbasis = basisobj, locations = as.matrix(spatiotemporal_data[,1:3]), time_mesh = time_steps, lambdaS = lambda_s, lambdaT = lambda_t)

# Print the result
print(result)


