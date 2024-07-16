library(tidyverse)

storm_data = read.csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

peril = "Blizzard"

storm_analysis = storm_data %>%
  filter(EVENT_CAT == peril) %>%
  filter(STATE != "HAWAII") %>%
  filter(STATE != "ALASKA") %>%
  filter(STATE != "PUERTO RICO") %>%
  filter(DAMAGE_PROPERTY > 0) %>%
  filter(YEAR > 1996) %>%
  filter(YEAR < 2024) %>%
  group_by(EPISODE_ID) %>%
  summarise(YEAR = first(YEAR),
            STATE_FIPS = first(STATE_FIPS),
            DAMAGE_PROPERTY = sum(DAMAGE_PROPERTY),
            X_CART = mean(X_CART),
            Y_CART = mean(Y_CART),
            BEGIN_DATE_TIME = first(BEGIN_DATE_TIME),
            EPSIODE_NARRATIVE = first(EPISODE_NARRATIVE))

storm_analysis = storm_analysis[complete.cases(storm_analysis %>% select(c("X_CART", "Y_CART"))),]


## Simple Generalised pareto distribution fit
data = storm_analysis %>% select(c("DAMAGE_PROPERTY", "YEAR"))
VGAM::meplot(data %>% pull(DAMAGE_PROPERTY) %>% log(), type ='p')
res_GPD = extRemes::fevd(x = log(data$DAMAGE_PROPERTY), type ='GP', threshold = 15, method = "Bayesian")
res_GPD
plot(res_GPD)

## Simple Generalised Extreme value distribution fit

x = data %>% group_by(YEAR) %>% summarise(max = max(DAMAGE_PROPERTY)) %>% pull(max)
res_EVD = extRemes::fevd(x = log(x), type = 'GEV', method = 'MLE')
res_EVD
plot(res_EVD)

## Adding covariates

GDP_state = read.csv("./Data/mapping_tables/US_GDP_by_state.csv")
GDP_state = GDP_state %>% pivot_longer(cols = starts_with("X"), names_to = "YEAR", values_to = "GDP")
GDP_state = GDP_state %>%
  mutate(YEAR = as.numeric(substr(YEAR, 2,5))) %>%
  mutate(GeoFips = as.integer(GeoFips/1000))


storm_analysis = storm_analysis %>% left_join(GDP_state %>% select(c("GeoFips", "GDP", "YEAR")), by = c("YEAR" = "YEAR", "STATE_FIPS" = "GeoFips")) %>%
  mutate(LOG_DAMAGE_PROPERTY = log(DAMAGE_PROPERTY)) 


simple_OLS = lm(LOG_DAMAGE_PROPERTY ~ YEAR, data = storm_analysis)
summary(simple_OLS)
