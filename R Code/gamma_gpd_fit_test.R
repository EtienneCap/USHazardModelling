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
  # filter(YEAR > 1996) %>%
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
  mutate(DAMAGE_LOG = log10(DAMAGE)) %>%
  mutate(DAMAGE_LN = log(DAMAGE)) %>%
  mutate(STATE_FIPS = factor(STATE_FIPS)) %>%
  # group_by(REGION, EVENT_CAT) %>%
  mutate(DAMAGE_SCALED = scale(DAMAGE, center = F)[,1]) %>%
  mutate(n = n()) %>%
  filter(!is.nan(DAMAGE_SCALED)) %>%
  # ungroup() %>%
  mutate(DAMAGE_BOX_COX = sae::bxcx(DAMAGE, lambda_bxcx)) %>%
  mutate(min_BXCX = min(DAMAGE_BOX_COX)) %>%
  mutate(DAMAGE_BOX_COX = DAMAGE_BOX_COX - min(DAMAGE_BOX_COX)) %>%
  mutate(MEAN = mean(DAMAGE)) %>%
  mutate(MED = median(DAMAGE)) %>%
  mutate(BEGIN_DATE_TIME = as.Date(BEGIN_DATE_TIME)) %>%
  mutate(SCALER = sqrt(mean(DAMAGE^2))) %>%
  filter(!is.na(REGION)) %>%
  filter(DAMAGE > 500) %>%
  mutate(YEAR = YEAR - min(YEAR)) %>%
  mutate(PERIOD = YEAR %/% 10)

par(mfrow = c(1,1))
storm_analysis %>% pull(DAMAGE_BOX_COX) %>% hist(breaks = 30)
storm_analysis %>% pull(DAMAGE_BOX_COX) %>% VGAM::meplot()

### ASSESSING TIME IMPACT ###

quant_reg = rq(DAMAGE ~ YEAR + CUSTOM_REGION , data = storm_analysis, tau = 0.80)
summary(quant_reg)

plot(quant_reg$model$YEAR, quant_reg$fitted.values %>% log10())

storm_analysis["THRESH_QUANT"] = quant_reg$fitted.values

### DAMAGE MODELLING ###

storm_analysis %>% group_by(CUSTOM_REGION) %>% summarise(n = n())

## Threshold selection

regions = storm_analysis %>% pull(CUSTOM_REGION) %>% unique()

par(mfrow = c(3,2))

thresh_list = c("North - Ohio Valley" = 12, "South" = 13, "Northwest" = 8.5, "Southeast" = 10, "Northeast" = 10, "West" = 8)

for (reg in regions){
  storm_analysis %>% filter(CUSTOM_REGION == reg) %>% pull(DAMAGE_BOX_COX) %>% evmix::mrlplot(main = paste0("Region: ", reg), legend.loc = NULL, try.thresh = c(thresh_list[[reg]]))
}

## Bulk fitting (gamma distribution)


# gamma_fits = rep(NA, length(thresh_list))
# names(gamma_fits) = names(thresh_list)
par(mfrow = c(3,2))
for (reg in regions){
  data = storm_analysis %>%
    filter(CUSTOM_REGION == reg) %>%
    filter(DAMAGE_BOX_COX <= thresh_list[[reg]])
  
  gamma_fit = gamlss::gamlss(DAMAGE_BOX_COX ~ YEAR, sigma.formula = ~ YEAR, data = data, family = gamlss.dist::GA())
  print("######################################")
  print(reg)
  print(paste0("Number of rows: ", data %>% nrow()))
  summary(gamma_fit)
  plot(gamma_fit, main = paste0("Region: ", reg))
  print("######################################")
  saveRDS(gamma_fit, file = paste0("./R Code/Result_gamma_fit/", reg, ".rds"))
}

## Bulk fitting (truncated gamma distribution)
alpha = 0.05
for (reg in regions){
  data_filtered = storm_analysis %>%
    filter(CUSTOM_REGION == reg) %>%
    filter(DAMAGE_BOX_COX <= thresh_list[[reg]])
  
  
  # Define the truncated gamma density function
  # dtruncgamma <- function(x, shape, rate, a = 0, b = thresh_list[[reg]]) {
  #   dat = data.frame(x = x, shape = shape, rate = rate)
  #   n = nrow(dat)
  #   sapply(1:n, function(i){
  #     if (dat$x[i] < a || dat$x[i] > b) {
  #       return(-Inf)
  #     }
  #     else{
  #       return(truncdist::dtrunc(dat$x[i], "gamma", shape = dat$shape[i], rate = dat$rate[i], a = a, b = b))
  #     }
  #   })
  # }
  
  dtruncgamma_log = function(x, shape, rate, upper.bound = thresh_list[[reg]]){
    int.const = pgamma(upper.bound, shape = shape, rate = rate, log.p = T)
    lik = dgamma(x, shape = shape, rate = rate, log = T) - int.const
    lik[x > upper.bound | x< 0] = -Inf
    return(lik)
  }
  
  loglik_truncgamma_cov <- function(params, data = data_filtered$DAMAGE_BOX_COX, time = data_filtered$YEAR) {
    alpha_0 <- params[1]
    alpha_1 <- params[2]
    beta_0 <- params[3]
    beta_1 <- params[4]
    
    shape <- exp(alpha_0 + alpha_1 * time)
    rate <- exp(beta_0 + beta_1 * time)
    
    sum(dtruncgamma_log(data, shape = shape, rate = rate))
  }
  
  init_params <- c(alpha_0 = 0.5, alpha_1 = 0.02, beta_0 = 0, beta_1 = 0.02)
  
  # Fit the truncated gamma distribution using maximum likelihood estimation
  fit <- optim(init_params, loglik_truncgamma_cov, control = list(fnscale = -1, maxit = 1e5))
  
  hessian = numDeriv::hessian(function(x){loglik_truncgamma_cov(x)}, fit$par)
  cov_matrix <- solve(-hessian*data_filtered %>% nrow())
  
  std_errors <- sqrt(diag(cov_matrix))
  
  # Confidence intervals
  conf_intervals <- data.frame(
    Estimate = fit$par,
    lower.bound = fit$par + qnorm(alpha/2) * std_errors,
    upper.bound = fit$par + qnorm(1-alpha/2) * std_errors)
  
  conf_intervals["signif"] = conf_intervals$`lower.bound` * conf_intervals$`upper.bound` > 0
  
  
  
  print("######################################")
  print(reg)
  print(paste0("Number of rows: ", data_filtered %>% nrow()))
  print(fit)
  print(conf_intervals)
  print("######################################")
  saveRDS(fit, file = paste0("./R Code/Result_gamma_fit/", reg, ".rds"))
}

# Fit assessement

par(mfrow = c(3,2))
for (reg in regions){
  data_filtered = storm_analysis %>%
    filter(CUSTOM_REGION == reg) %>%
    filter(DAMAGE_BOX_COX <= thresh_list[[reg]])
  
  hist(data_filtered %>% pull(DAMAGE_BOX_COX), freq = F)
  
  fit = readRDS(paste0("./R Code/Result_gamma_fit/", reg, ".rds"))
  
  shape <- exp(fit$par[["alpha_0"]])
  rate <- exp(fit$par[["beta_0"]])
  
  r = range(data_filtered %>% pull(DAMAGE_BOX_COX))
  x = seq(r[1], r[2], length.out = 250)
  lines(x,dgamma(x, shape = shape, rate = rate), col = 'red')
}



## Tail fitting (Generalised Pareto distribution)
par(mfrow = c(3,2))
for (reg in regions){
  data = storm_analysis %>%
    filter(CUSTOM_REGION == reg)
  
  GPD_fit = extRemes::fevd(x = data$DAMAGE_BOX_COX, 
                               data = data.frame(data),
                               threshold = thresh_list[[reg]],
                               type = "GP", 
                               scale.fun = ~ YEAR,
                               shape.fun = ~ 1,
                               method = "MLE")
  print("######################################")
  print(reg)
  print(paste0("Number of rows: ", data %>% nrow()))
  summary(GPD_fit)
  print(ci.fevd(GPD_fit, type = "parameter", alpha = 0.05))
  
  # plot(GPD_fit, type = 'probprob')
  plot(GPD_fit, type = 'qq', main = paste0("Region: ", reg))
  print("######################################")
  saveRDS(GPD_fit, file = paste0("./R Code/Result_GPD_fit/", reg, ".rds"))
  
}










sub_storm_analysis = storm_analysis %>%
  # filter(PERIOD == 2) %>%
  filter(CUSTOM_REGION == "North/Ohio Valley")

# par(mfrow = c(3, 1))
# storm_analysis %>% filter(PERIOD == 0) %>% pull(DAMAGE_LN)%>% hist(breaks = 30)
# storm_analysis %>% filter(PERIOD == 1) %>% pull(DAMAGE_LN) %>% hist(breaks = 30)
# storm_analysis %>% filter(PERIOD == 2) %>% pull(DAMAGE_LN) %>% hist(breaks = 30)


res_gamma_GPD = sub_storm_analysis %>%
  pull(DAMAGE_BOX_COX) %>%
  fgammagpdcon(phiu = F)
  # fgammagpdcon(useq = 7, fixedu = T)
  # fmgammagpdcon(M = 2, phiu = F) 
  # fgammagpdcon(useq = storm_analysis$THRESH_QUANT, fixedu = TRUE)

res_gamma_GPD[-(1:2)]

par(mfrow = c(1,1))
sub_storm_analysis %>% pull(DAMAGE_BOX_COX) %>% VGAM::meplot()
abline(v = res_gamma_GPD$u, lty = "dashed", col = 'red')

evmix.diag(res_gamma_GPD, N = 250, upperfocus = F)
evmix.diag(res_gamma_GPD, N = 250, upperfocus = T)
par(mfrow = c(1,1))


test = sub_storm_analysis %>% filter(DAMAGE_BOX_COX > res_gamma_GPD$u)

lin_model = lm(DAMAGE_BOX_COX ~ YEAR, data = sub_storm_analysis)
summary(lin_model)




densplot(res_gamma_GPD)

## Comparison with a gamma fit

gamma_fit = gamlss::gamlss(DAMAGE_BOX_COX ~ YEAR + CLIM_REGION, sigma.formula = ~ YEAR + CLIM_REGION, data = storm_analysis, family = gamlss.dist::GA())

summary(gamma_fit)

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
  scale_color_gradientn(colors = custom_colors)+
  theme_minimal()


