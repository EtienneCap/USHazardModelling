library(tidyverse)
library(extRemes)

storm_data = read.csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

peril = "Blizzard"

storm_analysis_raw = storm_data %>%
  # filter(EVENT_CAT == "Tropical Storm" | EVENT_CAT == "Hurricane" | EVENT_CAT == "Tropical Depression") %>%
  filter(EVENT_CAT == peril) %>%
  filter(STATE != "HAWAII") %>%
  filter(STATE != "ALASKA") %>%
  filter(STATE != "PUERTO RICO") %>%
  filter(DAMAGE_PROPERTY > 0) %>%
  filter(YEAR > 1996) %>%
  filter(YEAR < 2024) %>%
  group_by(EPISODE_ID) %>%
  summarise(YEAR = first(YEAR),
            EVENT_CAT = first(EVENT_CAT),
            STATE = first(STATE),
            STATE_FIPS = first(STATE_FIPS),
            DAMAGE_PROPERTY = sum(DISCOUNT_DAMAGE_PROPERTY, na.rm = T),
            X_CART = mean(X_CART, na.rm = T),
            Y_CART = mean(Y_CART, na.rm = T),
            BEGIN_DATE_TIME = first(BEGIN_DATE_TIME),
            EPSIODE_NARRATIVE = first(EPISODE_NARRATIVE))

storm_analysis_raw = storm_analysis_raw[complete.cases(storm_analysis_raw %>% select(c("X_CART", "Y_CART"))),]


## Adding covariates

GDP_state = read.csv("./Data/mapping_tables/US_GDP_by_state.csv")
GDP_state = GDP_state %>% pivot_longer(cols = starts_with("X"), names_to = "YEAR", values_to = "GDP")
GDP_state = GDP_state %>%
  mutate(YEAR = as.numeric(substr(YEAR, 2,5))) %>%
  mutate(GeoFips = as.integer(GeoFips/1000))

HPI_state = read.csv("./Data/mapping_tables/US_HPI_by_state.csv")
HPI_state = HPI_state %>%
  filter(level == "State") %>%
  filter(hpi_type == "traditional") %>%
  filter(hpi_flavor == "all-transactions") %>%
  filter(frequency == "quarterly") %>%
  mutate(date = yr + 0.25*period - 0.01) %>%
  mutate(place_name = toupper(place_name)) %>%
  mutate(HPI = index_nsa) %>%
  select(-index_nsa)

mapping_month = data.frame(Month = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"), MONTH = seq(1, 12))

SOI_data = read.csv("./Data/mapping_tables/SOI_data.csv") %>%
  pivot_longer(cols = 2:13, names_to = "Month", values_to = "SOI") %>%
  left_join(mapping_month, by = "Month") %>%
  select(c("YEAR", "MONTH", "SOI"))

storm_analysis = storm_analysis_raw %>% left_join(GDP_state %>% select(c("GeoFips", "GDP", "YEAR")), by = c("YEAR" = "YEAR", "STATE_FIPS" = "GeoFips")) %>%
  mutate(quart = quarter(BEGIN_DATE_TIME)) %>%
  mutate(MONTH = month(BEGIN_DATE_TIME)) %>%
  left_join(SOI_data, by = c("YEAR", "MONTH")) %>%
  left_join(HPI_state %>% select(c("yr", "period", "place_name", "HPI")), by = c("YEAR" = "yr", "quart" = "period", "STATE" = "place_name")) %>%
  mutate(LOG_DAMAGE_PROPERTY = log10(DAMAGE_PROPERTY)) %>%
  mutate(LN_DAMAGE_PROPERTY = log(DAMAGE_PROPERTY)) %>%
  mutate(STATE_FIPS = factor(STATE_FIPS)) %>%
  group_by(STATE_FIPS, EVENT_CAT) %>%
  mutate(scaled_damages = scale(DAMAGE_PROPERTY, center = FALSE)[,1]) %>%
  mutate(thresh_log = quantile(LOG_DAMAGE_PROPERTY, 0.9),
         thresh = quantile(DAMAGE_PROPERTY, 0.9),
         thresh_ln = quantile(LN_DAMAGE_PROPERTY, 0.9),
         n = n()) %>%
  mutate(thresh_scaled = 5* sd(scaled_damages)) %>%
  filter(!is.nan(scaled_damages))

## Linear regression
simple_OLS = lm(LN_DAMAGE_PROPERTY ~ log(GDP), data = storm_analysis)
summary(simple_OLS)

## GLM
glm_fit = glm(LOG_DAMAGE_PROPERTY ~  STATE_FIPS, data = storm_analysis, family = inverse.gaussian())
summary(glm_fit)

## GLM with EVD

x = storm_analysis %>%
  group_by(YEAR, STATE_FIPS) %>%
  mutate(BOOL = DAMAGE_PROPERTY == max(DAMAGE_PROPERTY)) %>%
  filter(BOOL)
  # summarise(DAMAGE_PROPERTY = max(DAMAGE_PROPERTY), LOG_DAMAGE_PROPERTY = max(log(DAMAGE_PROPERTY)), GDP = first(GDP), HPI = mean(HPI)) #%>%
  # mutate(GDP = scale(GDP)) %>%
  # mutate(DAMAGE_PROPERTY = scale(DAMAGE_PROPERTY)) %>%
  # mutate(LOG_DAMAGE_PROPERTY= scale(LOG_DAMAGE_PROPERTY))

res_GEV_cov = fevd(x = LOG_DAMAGE_PROPERTY, data.frame(x), location.fun = ~HPI + log(GDP) + SOI , type = 'GEV', method = "MLE")
res_GEV_cov
ci.fevd(res_GEV_cov, type = "parameter")
plot(res_GEV_cov)

x$test = (x$LOG_DAMAGE_PROPERTY - res_GEV_cov$results$par[1] - res_GEV_cov$results$par[2]*log10(x$GDP)) / res_GEV_cov$results$par[3]

res_GEV_nocov = fevd(x = LOG_DAMAGE_PROPERTY, data.frame(x %>% group_by(YEAR) %>% summarise(LOG_DAMAGE_PROPERTY = max(LOG_DAMAGE_PROPERTY))), type = "GEV")
res_GEV_nocov
ci.fevd(res_GEV_nocov, type = "parameter")
plot(res_GEV_nocov)

## GLM with POT

storm_analysis %>% pull(scaled_damages) %>% plot()
storm_analysis %>% pull(thresh_scaled) %>% lines(type = 'l', col = 'red')

VGAM::meplot(storm_analysis %>% pull(scaled_damages))

# res_GPD_cov = fevd(x = LN_DAMAGE_PROPERTY, data = data.frame(storm_analysis), threshold = storm_analysis$thresh_ln, type = "GP", scale.fun = ~ EVENT_CAT, method = "GMLE")
res_GPD_cov = fevd(x = scaled_damages, data = data.frame(storm_analysis), threshold = 1, type = "GP", scale.fun = ~1, method = "MLE")
res_GPD_cov
ci.fevd(res_GPD_cov, type = "parameter")
plot(res_GPD_cov, type = 'probprob')
plot(res_GPD_cov, type = 'qq')
plot(res_GPD_cov, type = "qq2")
plot(res_GPD_cov, type = 'rl', rperiods = c(10))
plot(res_GPD_cov, type = 'hist', hist.args = c("breaks" = 15, "freq" = F), main = "", xlab= "x")
plot(res_GPD_cov, rperiods = c(10))
par(mfrow = c(1,1))

## Geoplotting the residuals
boundaries_US = read.csv("./Data/US_map/boundaries_US_mesh.csv")
# boundaries_US = boundaries_US[-nrow(boundaries_US),]
boundaries_US = boundaries_US[nrow(boundaries_US):1,]

model = simple_OLS

storm_data_regression = storm_analysis
storm_data_regression["residuals"] = model$residuals

custom_colors <- c("blue", "green", "red")

ggplot() +
  geom_path(data = boundaries_US, mapping = aes(x = X_CART, y = Y_CART))+
  geom_point(data = storm_analysis, mapping = aes(x = X_CART, y = Y_CART, col = log10(DAMAGE_PROPERTY)))+
  scale_color_gradientn(colors = custom_colors)

## Quantile regression

library(quantreg)

quant_reg = rq(log(DAMAGE_PROPERTY) ~ STATE_FIPS, data = storm_analysis, tau = 0.8)
summary(quant_reg)  

## SpatialExtremes

library(SpatialExtremes)
# 
# locations = x %>% ungroup() %>% select(c("X_CART", "Y_CART")) %>% as.matrix()
# 
# fitspatgev(x$DAMAGE_PROPERTY, locations, loc.form = loc ~ X_CART + Y_CART, scale.form = scale ~ X_CART + Y_CART, shape.form = shape ~1)

