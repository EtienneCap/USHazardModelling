library(tidyverse)
library(fdaPDE)


### Handling the storm data
storm_data = read.csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

peril = "Tropical Storm"

storm_analysis = storm_data %>%
  filter(EVENT_CAT == peril) %>%
  filter(STATE != "HAWAII") %>%
  filter(STATE != "ALASKA") %>%
  filter(STATE != "PUERTO RICO") %>%
  filter(DAMAGE_PROPERTY > 0) %>%
  # filter(YEAR > 2004) %>%
  # filter(YEAR < 2013) %>%
  group_by(EPISODE_ID) %>%
  summarise(YEAR = first(YEAR),
            DAMAGE_PROPERTY = sum(DAMAGE_PROPERTY, na.rm = T),
            X_CART = mean(X_CART),
            Y_CART = mean(Y_CART),
            EPSIODE_NARRATIVE = first(EPISODE_NARRATIVE))

storm_analysis = storm_analysis[complete.cases(storm_analysis %>% select(c("X_CART", "Y_CART"))),]


# locations_mesh = storm_analysis %>% select(c('X_CART', "Y_CART")) %>% distinct()
# 
# tryCatch({
#   locations_mesh = locations_mesh[sample(1:nrow(locations_mesh), size = 400, replace = FALSE),]
# },
# error = function(e){
#   print(e)
#   print('Proceeding with the execution.')
# })
# 
# 
# locations_mesh = matrix(data = c(locations_mesh$X_CART, locations_mesh$Y_CART), ncol = 2)
locations = storm_analysis %>% select(c("X_CART", "Y_CART", "YEAR")) %>% 
  na.omit() %>%
  mutate(X_CART = X_CART + rnorm(nrow(storm_analysis), sd = 0.0003),
         Y_CART = Y_CART + rnorm(nrow(storm_analysis), sd = 0.0003))
 

begin_year = min(locations$YEAR)

time_locations = locations$YEAR-begin_year
locations = matrix(data = c(locations$X_CART, locations$Y_CART), ncol = 2)


### Creating the mesh

# Estimation mesh
boundaries_US = read.csv("./Data/US_map/boundaries_US_mesh.csv")
boundaries_US = boundaries_US[-nrow(boundaries_US),]

boundary_nodes = matrix(data = c(boundaries_US$X_CART, boundaries_US$Y_CART), ncol = 2)
boundary_segments = matrix(data = c(seq(1, dim(boundary_nodes)[1], 1), c(seq(2, dim(boundary_nodes)[1], 1)), 1), ncol = 2)

mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
mesh = refine.mesh.2D(mesh, maximum_area = 0.001, minimum_angle = 25)
FEMbasis <- create.FEM.basis(mesh)
plot(mesh)
saveRDS(mesh, "./R Code/R Data/mesh_estimation.rds")

# Evaluation mesh

data.mesh.eval = read.csv("./Data/US_map/evaluation_US_mesh.csv")
eval.nodes = matrix(data = c(data.mesh.eval$X_CART, data.mesh.eval$Y_CART), ncol = 2)
mesh.eval <- create.mesh.2D(nodes = rbind(boundary_nodes, eval.nodes), segments = boundary_segments)
# mesh.eval = refine.mesh.2D(mesh, maximum_area = 0.00005)
FEMbasis.eval <- create.FEM.basis(mesh.eval)
plot(mesh.eval)


mesh_time = seq(range(time_locations)[1], range(time_locations)[2], by = 1)

# lambda_space = 10^seq(from = -4, to = -2, by = 1)
# lambda_time = 10^seq(from = -4, to = -2, by = 1)

lambda_space = 1e-2
lambda_time = 1e-2

tryCatch({rm(storm_data)}, error = function(e){print(e)})

solution_STDEPDE <- DE.FEM.time(data = locations, data_time = time_locations, FEMbasis = FEMbasis,
                                mesh_time = mesh_time, lambda = lambda_space, tol1 = 1e-5, tol2 = 0,
                                lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                heatIter = 10, print = T, nfolds = NULL, nsimulations = 750,
                                step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                preprocess_method = "NoCrossValidation") 

saveRDS(solution_STDEPDE, paste0("./R Code/R Data/sol_", tolower(peril),".rds"))

FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC = F)


t <- mesh_time


mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
lim_max = -Inf
for (time_index in 1:length(t)) {
  
  evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = mesh.eval$nodes, time.instants = t[time_index])
  evaluation_STDEPDE <- exp(evaluation_STDEPDE)
  # evaluation_STDEPDE <- evaluation_STDEPDE/sum(evaluation_STDEPDE, na.rm = TRUE)
  
  mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE) * nrow(locations)/
  maxi = max(log10(mean_sol_STDEPDE[,time_index]), na.rm = T)
  if (maxi > lim_max){
    lim_max = maxi
  }
}
print(lim_max)

saveRDS(mean_sol_STDEPDE, paste0("./R Code/R Data/first_result_", tolower(peril), ".rds"))
custom_colors <- c("blue", "green", "yellow", "red")


t = 1:ncol(mean_sol_STDEPDE)

intensity_data = data.frame(x = mesh.eval$nodes[,1], y = mesh.eval$nodes[,2], z = mean_sol_STDEPDE[,3])


for (time_index in 1:length(t)){
  p = ggplot() +
    geom_point(data = intensity_data, mapping = aes(x = x, y = y, col = log10(z)), size = 1) +
    geom_path(data = boundaries_US, mapping = aes(x = X_CART, y = Y_CART)) +
    scale_color_gradientn(colors = custom_colors, limits = c(-5, lim_max))+
    labs(title = paste0("STKDE for ",  peril, " at year ", year))+
    theme_void()
  
  ggsave(paste0("./R Code/R Data/STDE-PDE/Plots_", peril,"/plot_", tolower(peril), "_", t[time_index]+begin_year, ".png"), plot = p, width = 10, height = 7, dpi = 300)
  
}


