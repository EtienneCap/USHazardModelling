library(tidyverse)
library(extRemes)
library(quantreg)
library(evmix)

storm_data = read.csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

peril = "Hail"

storm_analysis_raw = storm_data %>%
  # filter(EVENT_CAT == "Tropical Storm" | EVENT_CAT == "Hurricane" | EVENT_CAT == "Tropical Depression") %>%
  filter(EVENT_CAT == peril) %>%
  filter(STATE != "HAWAII") %>%
  filter(STATE != "ALASKA") %>%
  filter(STATE != "PUERTO RICO") %>%
  filter(DAMAGE_PROPERTY > 0) %>%
  # filter(YEAR > 1996) %>%
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

storm_analysis_raw = storm_analysis_raw[complete.cases(storm_analysis_raw %>% dplyr::select(c("EPISODE_ID", "X_CART", "Y_CART"))),]


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
  dplyr::select(-index_nsa)

mapping_month = data.frame(Month = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"), MONTH = seq(1, 12))

SOI_data = read.csv("./Data/mapping_tables/SOI_data.csv") %>%
  pivot_longer(cols = 2:13, names_to = "Month", values_to = "SOI") %>%
  left_join(mapping_month, by = "Month") %>%
  dplyr::select(c("YEAR", "MONTH", "SOI"))

US_regions = read.csv("./Data/mapping_tables/US_STATE.csv") %>%
  mutate(State = toupper(State)) %>%
  rename(c("REGION" = "Region", "STATE" = "State", "CLIM_REGION" = "Clim_Region", "CUSTOM_REGION" = "New_reg"))

tau = 0.9
lambda_bxcx = 0

storm_analysis = storm_analysis_raw %>% 
  mutate(quart = quarter(BEGIN_DATE_TIME)) %>%
  mutate(MONTH = month(BEGIN_DATE_TIME)) %>%
  # left_join(GDP_state %>% dplyr::select(c("GeoFips", "GDP", "YEAR")), by = c("YEAR" = "YEAR", "STATE_FIPS" = "GeoFips")) %>%
  left_join(SOI_data, by = c("YEAR", "MONTH")) %>%
  left_join(HPI_state %>% dplyr::select(c("yr", "period", "place_name", "HPI")), by = c("YEAR" = "yr", "quart" = "period", "STATE" = "place_name")) %>%
  left_join(US_regions %>% dplyr::select(c("STATE", "REGION", "CLIM_REGION", "CUSTOM_REGION")), by = "STATE") %>%
  mutate(REGION = factor(REGION)) %>%
  mutate(DAMAGE_LOG = log10(DAMAGE_PROPERTY)) %>%
  mutate(DAMAGE_LN = log(DAMAGE_PROPERTY)) %>%
  mutate(STATE_FIPS = factor(STATE_FIPS)) %>%
  # group_by(REGION, EVENT_CAT) %>%
  mutate(DAMAGE_SCALED = scale(DAMAGE_PROPERTY, center = F)[,1]) %>%
  mutate(thresh_log = quantile(DAMAGE_LOG, tau),
         thresh = quantile(DAMAGE_PROPERTY, tau),
         thresh_ln = quantile(DAMAGE_LN, tau),
         n = n()) %>%
  filter(!is.nan(DAMAGE_SCALED)) %>%
  # ungroup() %>%
  mutate(DAMAGE_BOX_COX = sae::bxcx(DAMAGE_PROPERTY/sd(DAMAGE_PROPERTY), lambda_bxcx)) %>%
  mutate(min_BXCX = min(DAMAGE_BOX_COX)) %>%
  mutate(DAMAGE_BOX_COX = DAMAGE_BOX_COX - min(DAMAGE_BOX_COX)) %>%
  mutate(MEAN = mean(DAMAGE_PROPERTY)) %>%
  mutate(MED = median(DAMAGE_PROPERTY)) %>%
  mutate(BEGIN_DATE_TIME = as.Date(BEGIN_DATE_TIME)) %>%
  mutate(SCALER = sqrt(mean(DAMAGE_PROPERTY^2))) %>%
  filter(!is.na(REGION)) %>%
  filter(DAMAGE_PROPERTY > 500) %>%
  mutate(YEAR = YEAR - min(YEAR))

### ASSESSING TIME IMPACT ###

data_quantile = storm_analysis_raw %>%
  group_by(YEAR) %>%
  summarise(q = quantile(DAMAGE_PROPERTY, 0.9))

quant_reg = rq(DAMAGE_BOX_COX ~ YEAR + REGION , data = storm_analysis)
summary(quant_reg)


### DAMAGE MODELLING ###


par(mfrow = c(1,1))
storm_analysis %>% pull(DAMAGE_BOX_COX) %>% hist(breaks = 30)
storm_analysis %>% pull(DAMAGE_BOX_COX) %>% VGAM::meplot()

storm_analysis %>% group_by(CLIM_REGION) %>% summarise(n = n())

sub_storm_analysis = storm_analysis %>% filter(CLIM_REGION == "Southwest") 
thresh_quant = sub_storm_analysis %>% pull(DAMAGE_BOX_COX) %>%
  quantile(0.95)

thresh_quant = thresh_quant[[1]]


res_gamma_GPD = sub_storm_analysis %>%
  pull(DAMAGE_BOX_COX) %>%
  fgammagpdcon(phiu = F)
  # fgammagpdcon(useq = 7, fixedu = T)
  # fmgammagpdcon(M = 2, phiu = F)

res_gamma_GPD[-(1:2)]

par(mfrow = c(1,1))
sub_storm_analysis %>% pull(DAMAGE_BOX_COX) %>% VGAM::meplot()
abline(v = res_gamma_GPD$u, lty = "dashed", col = 'red')

evmix.diag(res_gamma_GPD, N = 250, upperfocus = F)
evmix.diag(res_gamma_GPD, N = 250, upperfocus = T)
par(mfrow = c(1,1))

## Comparison with a gamm distribution

gamma_fit = sub_storm_analysis %>% pull(DAMAGE_BOX_COX) %>% fitdistrplus::fitdist(distr = "gamma", method = "mle")

plot(gamma_fit)


### DAMAGE VISUALISATION ON THE MAP ###

## Geoplotting the residuals
boundaries_US = read.csv("./Data/US_map/boundaries_US_mesh.csv")
boundaries_US = boundaries_US[nrow(boundaries_US):1,]

custom_colors <- c("blue", "green", "red")

ggplot() +
  geom_path(data = boundaries_US, mapping = aes(x = X_CART, y = Y_CART))+
  geom_point(data = storm_analysis, mapping = aes(x = X_CART, y = Y_CART, col = DAMAGE_LOG)) +
  scale_color_gradientn(colors = custom_colors)


