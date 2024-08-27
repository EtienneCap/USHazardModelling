# install.packages("fdaPDE")


### Simulation of a Space-Time non-homogeneous point process
# Radial-basis intenisty
intensity_pp = function(x, y, t, A = 5, T_tot = 10, d0 = 0.2){
  # Intensity in the [0,1] x [0,1] domain
  x1 = 0.8*cos(2*pi*t/T_tot)
  y1 = 0.8*sin(2*pi*t/T_tot)
  
  x2 = -0.8*cos(2*pi*t/T_tot)
  y2 = -0.8*sin(2*pi*t/T_tot)
  
  d1 = (x - x1)^2 + (y - y1)^2
  d2 = (x - x2)^2 + (y - y2)^2
  value = A * exp(-d1/d0^2) + A * exp(-d2/d0^2)
}

# Plot of the intensity function (contour plot)

x <- seq(-1, 1, length.out = 100)
y <- seq(-1, 1, length.out = 100)
# Create a matrix of z values
z <- outer(x, y, function(x, y) { intensity_pp(x, y, t = 7)})

contour(x, y, z, main = "Intensity plot", xlab = "X", ylab = "Y")

# Homogeneous poisson process

lambda = 75

t_mesh = seq(0, 10, length.out = 21)

points = rep(list(NULL), length(t_mesh))

for (j in seq_along(t_mesh)){

  N = rpois(1, lambda) # Random number of points
  
  point_pattern_x = runif(N, min = -1, max = 1)
  point_pattern_y = runif(N, min = -1, max = 1)
  
  #Thinning
  final_point_pattern_x = c()
  final_point_pattern_y = c()
  
  for (i in 1:N){
    if (runif(1) < intensity_pp(point_pattern_x[i], point_pattern_y[i], t = t_mesh[j], A = lambda)/lambda){
      final_point_pattern_x = c(final_point_pattern_x, point_pattern_x[i])
      final_point_pattern_y = c(final_point_pattern_y, point_pattern_y[i])
    }
  }
  
  points[[j]] = matrix(data = c(final_point_pattern_x, final_point_pattern_y), ncol = 2)
}

# Visualisation
J = 8
z <- outer(x, y, function(x, y) { intensity_pp(x, y, t = t_mesh[[J]], A = lambda)})
contour(x, y, z, main = "Intensity plot", xlab = "X", ylab = "Y")
points(points[[J]][,1], points[[J]][,2], col = 'red', pch = 16)

data = data.frame(x = points[[1]][,1], y = points[[1]][,2], t = t_mesh[1])
for (j in seq_along(t_mesh[-1])){
  add_data = data.frame(x = points[[j+1]][,1], y = points[[j+1]][,2], t = t_mesh[j+1])
  data = rbind(data, add_data)
}
dim(data)

########################################

### Intensity estimation ###

# STKDE

library(spatstat.geom)

window = spatstat.geom::owin(c(-1, 1), c(-1, 1))

ppp_data = spatstat.geom::ppp(data$x, data$y, window = window)

data_ppx = spatstat.geom::ppx(data, domain = list(spatial = window, temporal = c(0, 10)))

st_intensity = sparr::spattemp.density(ppp_data, h = 0.25, tt = data$t)

plot(st_intensity)


# STDE-PDE

library(fdaPDE)
library(ggplot2)

## Estimation mesh

x = c(seq(-1.3, 1.3, length.out = 20), rep(1.3, 20), seq(1.3, -1.3, length.out = 20), rep(-1.3, 20))
y = c(rep(-1.3, 20), seq(-1.3, 1.3, length.out = 20), rep(1.3, 20), seq(1.3, -1.3, length.out = 20))

boundary_nodes = matrix(data = c(x, y), ncol = 2)[-c(21, 41, 61, 80),]
boundary_segments = matrix(data = c(seq(1, dim(boundary_nodes)[1], 1), c(seq(2, dim(boundary_nodes)[1], 1)), 1), ncol = 2)
locations = matrix(data = c(data$x, data$y), ncol = 2)
time_locations = data$t


mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.15, minimum_angle = 30)


FEMbasis <- create.FEM.basis(mesh)

x11()
plot(mesh)
dev.off()


# Evaluation mesh
n = 60
X <- seq(-1.3, 1.3, length.out = n)
Y <- seq(-1.3, 1.3, length.out = n)
grid <- expand.grid(X, Y)
mesh.eval <- create.mesh.2D(grid)
FEMbasis.eval <- create.FEM.basis(mesh.eval)

x11()
plot(mesh.eval)
dev.off()

mesh_time = seq(range(time_locations)[1], range(time_locations)[2], length.out = 21)

lambda_space = 1e-3
lambda_time = 1e-2

solution_STDEPDE <- DE.FEM.time(data = locations, data_time = time_locations, FEMbasis = FEMbasis,
                                mesh_time = mesh_time, lambda = lambda_space, tol1 = 1e-5,
                                lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                heatIter = 10, print = T, nfolds = NULL, nsimulations = 1500,
                                step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                preprocess_method = "NoCrossValidation") 

FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC = F)

t <- mesh_time

mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))

for (time_index in 1:length(t)) {

  evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = mesh.eval$nodes, time.instants = t[time_index])
  evaluation_STDEPDE <- exp(evaluation_STDEPDE) * nrow(data)
  # evaluation_STDEPDE <- evaluation_STDEPDE/sum(evaluation_STDEPDE, na.rm = TRUE)
  
  mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
}

# M <- max(max(mean_sol_STDEPDE, na.rm = TRUE), na.rm = TRUE)

time_index = 20

z_contour <- outer(X,Y, function(x, y) { intensity_pp(x, y, t = mesh_time[time_index])})


ggplot(data.frame(x = mesh.eval$nodes[,1], y = mesh.eval$nodes[,2], z = mean_sol_STDEPDE[,time_index], z_contour = as.vector(z_contour)))+
  geom_tile(mapping = aes(x = x, y=y, fill = z))+
  geom_contour(mapping = aes(x = x, y = y, z = z_contour), col = "red") +
  theme_minimal()






