# Install and load necessary packages
library(fdaPDE)
library(tidyverse)


domain_area <- 4*pi

# Spatial 2D Mesh over the Unitary Sphere for Estimation
vertices <- read.table("mesh/sphere.606.vertices.txt", quote = "\"", comment.char = "")
triangles <- read.table("mesh/sphere.606.triangles.txt", quote = "\"", comment.char = "")
mesh <- create.mesh.2.5D(nodes = vertices[,1:3], triangles = triangles[,1:3])
FEMbasis <- create.FEM.basis(mesh)

plot.mesh(mesh)
rgl.close()

# Fine Spatial Mesh over the Unitary Sphere for Evaluation
vertices.eval <- read.table("mesh/sphere.5016.vertices.txt", row.names = 1)
vertices.eval.proj <- projection.points.2.5D(mesh, vertices.eval)
triangles.eval <- read.table("mesh/sphere.5016.triangles.txt", row.names = 1)
mesh.eval <- create.mesh.2.5D(nodes = vertices.eval, triangles = triangles.eval)
FEMbasis.eval <- create.FEM.basis(mesh.eval)

plot.mesh(mesh.eval)
rgl.close()

# Fine Grid Size for Evaluation
n <- 32

# Temporal 1D Mesh Over [0,1]
mesh_time <- seq(from = 0, to = 1, by = 0.1)





# Time for One Process
t_proc <- proc.time()

### DATA -------------------------------------------------------------------
# Generate the Data

N = 10000

generate.data(N, proc) # run only once to generate the data

# Read the Data
data <- read.table(paste0("data/",N,"data_",proc,".txt"))

# Locations
locations <- data[,1:3]

# Times
times <- data[,4]

# Discrete Times
discrete_times <- (times-0.05) %/% 0.1
discrete_times[discrete_times == -1] <- 0
discrete_times[discrete_times == 9] <- 8

## STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. --------
t0 <- proc.time()

# Smoothing Parameters
lambda <- 0.01
# alternative: lambda <- 10^seq(from = -3, to = -1, by = 1)

lambda_time <- 0.001
# alternative: lambda_time <- 10^seq(from = -4, to = -2, by = 1)

# Solution
# [If lambda and/or lambda_time are vectors, to select the best proposals by
# 10-folds CV, please specify: preprocess_method = "RightCV"]
solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                mesh_time = mesh_time, lambda = lambda, tol1 = 1e-7,
                                lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                heatIter = 10, print = T, nfolds = 10, nsimulations = 10000,
                                step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                preprocess_method = "NoCrossValidation")

FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC=F)

# CPU Time
CPUtime <- proc.time() - t0
CPUtimes_STDEPDE[proc] <- CPUtime[3]




for(time_index in 1:length(t)) {
  # STDE-PDE
  evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = vertices.eval.proj, time.instants = t[time_index])
  evaluation_STDEPDE <- exp(evaluation_STDEPDE)
  evaluation_STDEPDE <- evaluation_STDEPDE / sum(evaluation_STDEPDE, na.rm = TRUE)
  
  mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
   
}

time_index = 5

plot.FEM(FEM(mean_sol_STDEPDE[,time_index], FEMbasis.eval), m = 0, M = M)