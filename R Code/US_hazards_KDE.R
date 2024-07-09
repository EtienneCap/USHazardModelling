library(tidyverse)
library(fdaPDE)


### Handling the storm data
storm_data = read.csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

storm_analysis = storm_data %>%
  filter(EVENT_CAT == "Blizzard") %>%
  filter(STATE != "HAWAII") %>%
  filter(STATE != "ALASKA") %>%
  filter(STATE != "PUERTO RICO") %>%
  filter(DAMAGE_PROPERTY > 0) %>%
  group_by(EPISODE_ID) %>%
  summarise(YEAR = first(YEAR),
            DAMAGE_PROPERTY = sum(DAMAGE_PROPERTY),
            X_CART = mean(X_CART),
            Y_CART = mean(Y_CART),
            EPSIODE_NARRATIVE = first(EPISODE_NARRATIVE))

test = storm_analysis %>% select(c('X_CART', "Y_CART")) %>% distinct()

test = test[sample(1:nrow(test), size = 400, replace = FALSE),]

locations = matrix(data = c(test$X_CART, test$Y_CART), ncol = 2)
time_locations = storm_analysis$YEAR

### Creating the mesh

# Estimation mesh
boundaries_US = read.csv("./Data/US_map/boundaries_US_mesh.csv")
boundaries_US = boundaries_US[-nrow(boundaries_US),]

boundary_nodes = matrix(data = c(boundaries_US$X_CART, boundaries_US$Y_CART), ncol = 2)
boundary_segments = matrix(data = c(seq(1, dim(boundary_nodes)[1], 1), c(seq(2, dim(boundary_nodes)[1], 1)), 1), ncol = 2)


mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)

plot(mesh)

saveRDS(mesh, "US_mesh_estimation.rds")

# Evaluation mesh





