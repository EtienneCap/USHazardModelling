library(tidyverse)
library(spatstat.geom)


### Handling the storm data
storm_data = read.csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

peril = "Hurricane"

storm_analysis = storm_data %>%
  filter(EVENT_CAT == peril) %>%
  # filter(STATE != "HAWAII") %>%
  # filter(STATE != "ALASKA") %>%
  # filter(STATE != "PUERTO RICO") %>%
  filter(DAMAGE_PROPERTY > 0) %>%
  # filter(YEAR > 2004) %>%
  # filter(YEAR < 2013) %>%
  group_by(EPISODE_ID) %>%
  summarise(YEAR = first(YEAR),
            DAMAGE_PROPERTY = sum(DAMAGE_PROPERTY),
            X_CART = mean(X_CART),
            Y_CART = mean(Y_CART),
            EPSIODE_NARRATIVE = first(EPISODE_NARRATIVE))

storm_analysis = storm_analysis[complete.cases(storm_analysis %>% dplyr::select(c("X_CART", "Y_CART"))),] 
# Jittering
storm_analysis = storm_analysis %>%
  mutate(X_CART = X_CART + rnorm(nrow(storm_analysis), sd = 0.0003),
         Y_CART = Y_CART + rnorm(nrow(storm_analysis), sd = 0.0003))

years = unique(storm_analysis$YEAR) %>% sort()

# locations = storm_analysis %>% dplyr::select(c("X_CART", "Y_CART", "YEAR")) %>% na.omit()
# 
# begin_year = min(locations$YEAR)
# 
# time_locations = locations$YEAR-begin_year
# locations = matrix(data = c(locations$X_CART, locations$Y_CART), ncol = 2)


### Creating the mesh

boundaries_US = read.csv("./Data/US_map/boundaries_US_mesh.csv")
boundaries_US = boundaries_US[-nrow(boundaries_US),]
boundaries_US = boundaries_US[nrow(boundaries_US):1,]

poly = list(x = boundaries_US$X_CART, y = boundaries_US$Y_CART)

window = spatstat.geom::owin(poly = poly)

inside_bool = spatstat.geom::inside.owin(x = storm_analysis$X_CART, y = storm_analysis$Y_CART, w = window)
data = (storm_analysis %>% dplyr::select(c("X_CART", "Y_CART", "YEAR")))[inside_bool, ]

ppp_data = spatstat.geom::ppp(data$X_CART, data$Y_CART, window = window)

# data_ppx = spatstat.geom::ppx(data, domain = list(spatial = window, temporal = range(time_locations)))

st_intensity = sparr::spattemp.density(ppp_data, tt = data$YEAR, sres = 300)

# Transforming into an intensity function
lim_max = -Inf
for (year in unique(data$YEAR)){
  st_intensity$z[[as.character(year)]]$v = st_intensity$z[[as.character(year)]]$v * nrow(data)/ (max(years) - min(years))
  maxi = max(log10(st_intensity$z[[as.character(year)]]$v), na.rm = T)
  if (maxi > lim_max){
    lim_max = maxi
  }
}
print(lim_max)

x_coords = st_intensity$z[[as.character(unique(data$YEAR)[1])]]$xcol
y_coords = st_intensity$z[[as.character(unique(data$YEAR)[1])]]$yrow

intensity_data = expand.grid(x = x_coords, y = y_coords)
custom_colors <- c("blue", "green", "yellow", "red")

years = st_intensity$tgrid
# years = c(2015)

for (year in years){
  intensity_data$z = as.vector(st_intensity$z[[as.character(year)]]$v %>% t())
  
  p = ggplot() +
    geom_tile(data = intensity_data, mapping = aes(x = x, y = y, fill = log10(z)), linewidth = 0) +
    geom_path(data = boundaries_US, mapping = aes(x = X_CART, y = Y_CART)) +
    # scale_fill_viridis_c(option = "plasma", na.value = "white", limits = c(-5, 7)) +
    scale_fill_gradientn(colors = custom_colors, limits = c(-5, lim_max), na.value = "white")+
    # scale_fill_gradient(palette = 'plasma', na.value = "white")+
    labs(title = paste0("STKDE fo r",  peril, " at year ", year))+
    theme_void()
    
  # plot(p)
  ggsave(paste0("./R Code/R Data/STKDE/", peril,"/Plot_", tolower(peril),"_", year, ".png"), p, bg = "white", dpi = 600, width = 7, height = 5, unit = "in")
}
