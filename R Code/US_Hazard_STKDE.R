library(tidyverse)
library(cowplot)
library(spatstat.geom)

### Handling the storm data
storm_data = read.csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

peril = "Hail"

storm_analysis = storm_data %>%
  # filter(EVENT_CAT == "Tropical Storm" | EVENT_CAT == "Hurricane" | EVENT_CAT == "Tropical Depression") %>%
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

# Jittering the duplicated locations
storm_analysis = storm_analysis %>% 
  mutate(dup = duplicated(storm_analysis %>% dplyr::select(c("X_CART", "Y_CART")))) %>%
  mutate(noise =  rnorm(n(), sd = 1e-4)) %>%
  mutate(X_CART = if_else(dup, X_CART + noise, X_CART)) %>%
  mutate(Y_CART = if_else(dup, Y_CART + noise, Y_CART))

storm_analysis = storm_analysis %>% dplyr::select(-c("noise", 'dup'))
years = unique(storm_analysis$YEAR) %>% sort()

### Creating the study window

boundaries_US = read.csv("./Data/US_map/boundaries_US_mesh.csv")
boundaries_US = boundaries_US[-nrow(boundaries_US),]
boundaries_US = boundaries_US[nrow(boundaries_US):1,]

poly = list(x = boundaries_US$X_CART, y = boundaries_US$Y_CART)

window = spatstat.geom::owin(poly = poly)

inside_bool = spatstat.geom::inside.owin(x = storm_analysis$X_CART, y = storm_analysis$Y_CART, w = window)
data = (storm_analysis %>% dplyr::select(c("X_CART", "Y_CART", "YEAR")))[inside_bool, ]

ppp_data = spatstat.geom::ppp(data$X_CART, data$Y_CART, window = window)

st_intensity = sparr::spattemp.density(ppp_data, tt = data$YEAR, sres = 300)

# Transforming into an intensity function
lim_max = -Inf
for (year in unique(data$YEAR)){
  n = data %>% filter(YEAR == year) %>% nrow()
  
  st_intensity$z[[as.character(year)]]$v = st_intensity$z[[as.character(year)]]$v * n
  maxi = max(log10(st_intensity$z[[as.character(year)]]$v), na.rm = T)
  if (maxi > lim_max){
    lim_max = maxi
  }
}

x_coords = st_intensity$z[[as.character(unique(data$YEAR)[1])]]$xcol
y_coords = st_intensity$z[[as.character(unique(data$YEAR)[1])]]$yrow

intensity_data = expand.grid(x = x_coords, y = y_coords)
custom_colors <- c("blue", "green", "yellow", "red", "purple")

years = st_intensity$tgrid

bound = read.csv("./Data/US_map/boundaries_US_mesh.csv")

for (year in years){
  intensity_data$z = as.vector(st_intensity$z[[as.character(year)]]$v %>% t())
  
  p = ggplot() +
    geom_tile(data = intensity_data, mapping = aes(x = x, y = y, fill = log10(z)), linewidth = 0) +
    geom_path(data = bound, mapping = aes(x = X_CART, y = Y_CART)) +
    scale_fill_gradientn(colors = custom_colors, na.value = "white", limits = c(0,lim_max), name = "Intensity (log-scale)")+ # , limits = c(-5, lim_max)
    labs(title = paste0("STKDE for ",  peril, " at year ", year))+
    guides(fill = guide_colorbar(title.position = "right",
                                 title.theme = element_text(angle = -90),
                                 title.hjust = 0.5,
                                 barheight = 9))+
    theme_void() + 
    theme(plot.title = element_text(size = 15))
    
  # ggsave(paste0("./R Code/R Data/STKDE/", peril,"/Plot_", tolower(peril),"_", year, ".png"), p, bg = "white", dpi = 600, width = 7, height = 5, unit = "in")
  # ggsave(paste0("./R Code/R Data/STKDE/", peril,"/grouped_plots/Plot_", tolower(peril),"_", year, ".png"), p, bg = "white", dpi = 600, width = 7, height = 5, unit = "in")
}

# Grouped plot

plot_list = vector("list", length = length(years))
i=1
for (year in years){
  intensity_data$z = as.vector(st_intensity$z[[as.character(year)]]$v %>% t())
  
  p = ggplot() +
    geom_tile(data = intensity_data, mapping = aes(x = x, y = y, fill = log10(z)), linewidth = 0) +
    geom_path(data = bound, mapping = aes(x = X_CART, y = Y_CART)) +
    scale_fill_gradientn(colors = custom_colors, na.value = "white", limits = c(0,lim_max), name = "Intensity (log-scale)")+
    labs(title = paste0(peril, " (", year,")"))+
    guides(fill = guide_colorbar(title.position = "right", 
                                 title.theme = element_text(angle = -90),
                                 title.hjust = 0.5,
                                 barheight = 9))+
    theme_void()
  
  plot_list[[i]] = p
  i = i+1
  # ggsave(paste0("./R Code/R Data/STKDE/", peril,"/Plot_", tolower(peril),"_", year, ".png"), p, bg = "white", dpi = 600, width = 7, height = 5, unit = "in")
}

p = plot_grid(plotlist = plot_list, nrow = 6, ncol = 5, scale = 0.75) 
# ggsave(paste0("./R Code/R Data/STKDE/", peril,"/Plot_all_", tolower(peril), ".png"), p, bg = "white", dpi = 800, width = 21, height = 29.7, unit = "cm")

# Colorbar
dummy_data <- data.frame(x = 1, y = 1, fill = 1)

p <- ggplot(dummy_data, aes(x = x, y = y, fill = fill)) +
  geom_tile() +  
  scale_fill_gradientn(colors = custom_colors, na.value = "white", limits = c(-3,lim_max), name = "Intensity (log-scale)")+ # , limits = c(-5, lim_max) +  # Customize color bar
  guides(fill = guide_colorbar(direction = "horizontal",
                               title.position = "top", 
                               title.hjust = 0.5,
                               barheight = 2.5,
                               barwidth = 30))+ 
  theme_void() + 
  theme(legend.position = "top")
