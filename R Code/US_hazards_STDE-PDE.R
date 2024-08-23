rm(list = ls())
library(tidyverse)
library(fdaPDE)


### Handling the storm data
storm_data = read.csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

peril = "Hurricane"

storm_analysis = storm_data %>%
  filter(EVENT_CAT == peril) %>%
  filter(STATE != "HAWAII") %>%
  filter(STATE != "ALASKA") %>%
  filter(STATE != "PUERTO RICO") %>%
  filter(YEAR > 1995) %>%
  filter(YEAR < 2024) %>%
  group_by(EPISODE_ID) %>%
  summarise(YEAR = first(YEAR),
            EVENT_CAT = first(EVENT_CAT),
            STATE = first(STATE),
            STATE_FIPS = first(STATE_FIPS),
            DAMAGE = sum(DISCOUNT_DAMAGE_PROPERTY, na.rm = T),
            X_CART = mean(X_CART, na.rm = T),
            Y_CART = mean(Y_CART, na.rm = T),
            BEGIN_DATE_TIME = first(BEGIN_DATE_TIME),
            EPSIODE_NARRATIVE = first(EPISODE_NARRATIVE)) %>%
  filter(DAMAGE > 0)


storm_analysis = storm_analysis[complete.cases(storm_analysis %>% dplyr::select(c("X_CART", "Y_CART"))),]

locations = storm_analysis %>% dplyr::select(c("X_CART", "Y_CART", "YEAR")) %>% 
  na.omit()

# Jittering the duplicate coordinates

locations = locations %>% 
  mutate(dup = duplicated(locations %>% dplyr::select(c("X_CART", "Y_CART")))) %>%
  mutate(noise =  rnorm(n(), sd = 1e-4)) %>%
  mutate(X_CART = if_else(dup, X_CART + noise, X_CART)) %>%
  mutate(Y_CART = if_else(dup, Y_CART + noise, Y_CART))
 
locations = locations %>% dplyr::select(-c("noise", 'dup'))

begin_year = min(locations$YEAR)

time_locations = locations$YEAR #-begin_year
locations = matrix(data = c(locations$X_CART, locations$Y_CART), ncol = 2)

### Creating the mesh

# Estimation mesh
boundaries_US = read.csv("./Data/US_map/boundaries_US_mesh.csv")
boundaries_US = boundaries_US[-nrow(boundaries_US),]

# Filtering out the points that are outside the boundaries of the study region
bound = boundaries_US[nrow(boundaries_US):1,]
poly = list(x = bound$X_CART, y = bound$Y_CART)
window = spatstat.geom::owin(poly = poly)
inside_bool = spatstat.geom::inside.owin(x = locations[,1], y = locations[,2], w = window)


boundary_nodes = matrix(data = c(boundaries_US$X_CART, boundaries_US$Y_CART), ncol = 2)
boundary_segments = matrix(data = c(seq(1, dim(boundary_nodes)[1], 1), c(seq(2, dim(boundary_nodes)[1], 1)), 1), ncol = 2)

mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations[inside_bool,]), segments = boundary_segments)
mesh = refine.mesh.2D(mesh, maximum_area = 0.005)#, minimum_angle = 15)
FEMbasis <- create.FEM.basis(mesh)
plot(mesh)
saveRDS(mesh, "./R Code/R Data/mesh_estimation.rds")

# Evaluation mesh

data.mesh.eval = read.csv("./Data/US_map/evaluation_US_mesh.csv")
eval.nodes = matrix(data = c(data.mesh.eval$X_CART, data.mesh.eval$Y_CART), ncol = 2)
mesh.eval <- create.mesh.2D(nodes = rbind(boundary_nodes, eval.nodes), segments = boundary_segments)
FEMbasis.eval <- create.FEM.basis(mesh.eval)
plot(mesh.eval)
saveRDS(mesh.eval, "./R Code/R Data/mesh_evaluation.rds")

mesh_time = seq(range(time_locations)[1], range(time_locations)[2], by = 1)


lambda_space = 10^seq(from = -4, to = -2, by = 0.5)
lambda_time = 10^seq(from = -2, to = 0, by = 0.5)


tryCatch({rm(storm_data)}, error = function(e){print(e)}) # saving some memory space to prevent the R session to crash

solution_STDEPDE <- DE.FEM.time(data = locations, data_time = time_locations, FEMbasis = FEMbasis,
                                mesh_time = mesh_time, lambda = lambda_space, tol1 = 1e-6, tol2 = 0,
                                lambda_time = lambda_time, print = T, nfolds = 4, nsimulations = 500,
                                step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                preprocess_method = "RightCV", inference = T,
                                isTimeDiscrete = TRUE) 

saveRDS(solution_STDEPDE, paste0("./R Code/R Data/sol_", tolower(peril),".rds"))
solution_STDEPDE = readRDS(paste0("./R Code/R Data/sol_", tolower(peril),".rds"))
FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, solution_STDEPDE$FEMbasis, FLAG_PARABOLIC = F)


t <- mesh_time


mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
lim_max = -Inf
for (time_index in 1:length(t)) {
  
  evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = mesh.eval$nodes, time.instants = t[time_index])
  evaluation_STDEPDE <- exp(evaluation_STDEPDE)
  evaluation_STDEPDE <- evaluation_STDEPDE/sum(evaluation_STDEPDE, na.rm = TRUE)
  
  mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE) * nrow(locations)
  maxi = max(log10(mean_sol_STDEPDE[,time_index]), na.rm = T)
  if (maxi > lim_max){
    lim_max = maxi
  }
}
print(lim_max)

saveRDS(mean_sol_STDEPDE, paste0("./R Code/R Data/first_result_", tolower(peril), ".rds"))
custom_colors <- c("blue", "green", "yellow", "red", "purple")


t = 1:ncol(mean_sol_STDEPDE)




for (time_index in 1:length(t)){
  
  intensity_data = data.frame(x = mesh.eval$nodes[,1], y = mesh.eval$nodes[,2], z = mean_sol_STDEPDE[,time_index])
  p = ggplot(data = intensity_data, aes(x = x, y = y)) +
    geom_point(mapping = aes(col = log10(z)))+#, size = 1) +
    geom_path(data = boundaries_US, mapping = aes(x = X_CART, y = Y_CART)) +
    scale_color_gradientn(colors = custom_colors, na.value = "white", limits = c(-3,lim_max), name = "Intensity (log-scale)")+
    labs(title = paste0(t[time_index]+begin_year-1))+
    guides(color = guide_colorbar(title.position = "right", 
                                 title.theme = element_text(angle = -90),
                                 title.hjust = 0.5,
                                 barheight = 9))+
    theme_void()+
    theme(plot.title = element_text(size = 23))
  
  # ggsave(paste0("./R Code/R Data/STDE-PDE/Plots_", peril,"/plot_bis_", tolower(peril), "_", t[time_index]+begin_year-1, ".png"), plot = p, width = 10, height = 7, dpi = 600)
}


