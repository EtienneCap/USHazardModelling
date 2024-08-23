rm(list = ls())
library(cowplot)
library(tidyverse)
library(extRemes)
library(quantreg)
library(evmix)

storm_data = read.csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

peril = "Hurricane"

if (peril == 'Hail'){
  nrow = 2
  ncol = 3
} else if (peril == "Hurricane"){
  nrow = 1
  ncol = 2
} else {
  stop("Please select either Hail or Hurricane for the variable peril")
}

storm_analysis_raw = storm_data %>%
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

storm_analysis_raw = storm_analysis_raw[complete.cases(storm_analysis_raw %>% dplyr::select(c("EPISODE_ID", "X_CART", "Y_CART"))),]

min_year = storm_analysis_raw %>% pull(YEAR) %>% min()

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
lambda_bxcx = 1

storm_analysis = storm_analysis_raw %>% 
  mutate(quart = quarter(BEGIN_DATE_TIME)) %>%
  mutate(MONTH = month(BEGIN_DATE_TIME)) %>%
  left_join(US_regions %>% dplyr::select(c("STATE", "REGION", "CLIM_REGION", "CUSTOM_REGION")), by = "STATE") %>%
  mutate(REGION = factor(REGION)) %>%
  mutate(DAMAGE_LOG = log10(DAMAGE)) %>%
  mutate(DAMAGE_LN = log(DAMAGE)) %>%
  mutate(STATE_FIPS = factor(STATE_FIPS)) %>%
  mutate(n = n()) %>%
  mutate(BEGIN_DATE_TIME = as.Date(BEGIN_DATE_TIME)) %>%
  filter(!is.na(REGION)) %>%
  mutate(YEAR = YEAR - min_year) %>%
  mutate(DAMAGE_BOX_COX = DAMAGE)



## Pre-visualising the data
par(mfrow = c(1,1))
storm_analysis %>% pull(DAMAGE) %>% hist(breaks = 30)
storm_analysis %>% pull(DAMAGE_LOG) %>% hist(breaks = 30)
storm_analysis %>% pull(DAMAGE) %>% VGAM::meplot()
storm_analysis %>% pull(DAMAGE) %>% VGAM::meplot()


### DAMAGE MODELLING ###

storm_analysis %>% group_by(CUSTOM_REGION) %>% summarise(n = n())

## Threshold selection

regions = storm_analysis %>% group_by(CUSTOM_REGION) %>% filter(n() > 15) %>% pull(CUSTOM_REGION) %>% unique()

mean_excess <- function(u, data) {
  excess <- data[data > u] - u  # Calculate the excesses over the threshold u
  mean(excess)  # Return the mean of these excesses
}


storm_analysis %>% group_by(CUSTOM_REGION) %>% summarise(quantile = quantile(DAMAGE, 0.8))

if (peril == "Hurricane"){
  thresh_list = c("South" = 1e7, "Southeast" = 1e8)
} else if (peril == "Hail"){
  thresh_list = c("North - Ohio Valley" = 7.6e5, "South" =7.3e5, "Northwest" = 2.9e5, "Southeast" = 2.1e5, "Northeast" = 1.4e5, "West" = 1.8e5)
} 



plot_list = vector("list", length = length(regions))
  
i = 1

par(mfrow = c(1,3))
for (reg in regions){

  u = storm_analysis %>% filter(CUSTOM_REGION == reg) %>% pull(DAMAGE)
  me_u = sapply(u, function(x){mean_excess(x, u)})
  x = u[u>= thresh_list[[reg]]]
  y = me_u[u>= thresh_list[[reg]]]
  regression = lm(y ~ x)
  slope = regression$coefficients[2]
  intercept = regression$coefficients[1]

  p = ggplot() +
    geom_line(data = data.frame(u = u, me_u = me_u), mapping = aes(x = u, y = me_u), col = "steelblue", linewidth = 0.75) +
    geom_point(data = data.frame(u = u, me_u = me_u), mapping = aes(x = u, y = me_u), col = "steelblue", size = 1) +
    geom_vline(xintercept = thresh_list[[reg]], linetype = "dashed", color = "red", linewidth = 0.75) +
    geom_line(data = data.frame(x = x, y = intercept + slope * x), mapping = aes(x = x, y = y), color = "black", linewidth = 0.6, linetype = "dotted")+
    labs(title = paste0("Mean excess plot: ",reg), x = "Threshold u", y = "Mean excess", subtitle = paste0("(n = ", length(u), ")"))+
    scale_x_log10(limits = c(NA, NA))+
    theme_minimal()

  plot_list[[i]] = p
  i = i+1
  # storm_analysis %>% filter(CUSTOM_REGION == reg) %>% pull(DAMAGE) %>% evmix::mrlplot(main = paste0("Region: ", reg))#, try.thresh = c(thresh_list[[reg]]))
}

p = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, scale = 1) 
p



## Bulk fitting (truncated gamma distribution)
alpha = 0.05
for (reg in regions){
  data_filtered = storm_analysis %>%
    filter(CUSTOM_REGION == reg) %>%
    filter(DAMAGE <= thresh_list[[reg]]) %>%
    filter(DAMAGE > 0)

  
  dtruncgamma_log = function(x, shape, rate, upper.bound = thresh_list[[reg]]){
    int.const = pgamma(upper.bound, shape = shape, rate = rate, log.p = T)
    lik = dgamma(x, shape = shape, rate = rate, log = T) - int.const
    lik[x > upper.bound | x< 0] = -Inf
    return(lik)
  }
  
  loglik_truncgamma_cov <- function(params, data = data_filtered$DAMAGE, time = data_filtered$YEAR) {
    alpha_0 <- params[1]
    alpha_1 <- params[2]
    beta_0 <- params[3]
    # beta_1 <- params[4]
    
    shape <- exp(alpha_0 + alpha_1 * time)
    rate <- exp(beta_0) # + beta_1 * time)
    
    sum(dtruncgamma_log(data, shape = shape, rate = rate))
  }
  
  init_params <- c(alpha_0 = 0.5, alpha_1 = 0.02, beta_0 = 0)#, beta_1 = 0.02)
  
  # Fit the truncated gamma distribution using maximum likelihood estimation
  fit <- optim(init_params, loglik_truncgamma_cov, control = list(fnscale = -1, maxit = 1e5))
  
  hessian = numDeriv::hessian(function(x){loglik_truncgamma_cov(x)}, fit$par)
  cov_matrix <- solve(-hessian)
  
  std_errors <- sqrt(diag(cov_matrix))
  
  # Confidence intervals
  conf_intervals <- data.frame(
    Estimate = fit$par,
    lower.bound = fit$par + qnorm(alpha/2) * std_errors,
    upper.bound = fit$par + qnorm(1-alpha/2) * std_errors,
    pm = qnorm(1-alpha/2) * std_errors,
    std_errors = std_errors)
  
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

plot_list = vector("list", length = length(regions))
i=1
for (reg in regions){
  data_filtered = storm_analysis %>%
    filter(CUSTOM_REGION == reg) %>%
    filter(DAMAGE <= thresh_list[[reg]])
  
  fit = readRDS(paste0("./R Code/Result_gamma_fit/", reg, ".rds"))
  
  shape <- exp(fit$par[["alpha_0"]])
  rate <- exp(fit$par[["beta_0"]])
  r = range(data_filtered %>% pull(DAMAGE))
  x = seq(r[1], r[2], length.out = 250)
  # lines(x,dgamma(x, shape = shape, rate = rate), col = 'red')
  
  p = ggplot()+
    geom_histogram(data = data_filtered,
                   mapping = aes(x = DAMAGE, y = ..density..),
                   fill = "lightgray",
                   color = "black",
                   alpha = 0.5)+
    geom_line(data = data.frame(x = x, y = dgamma(x, shape = shape, rate = rate)),
              mapping = aes(x = x, y = y),
              col = 'red')+
    labs(title = paste0("Region: ", reg), x = "log(Damage)", y = "Density", subtitle = paste0("(n = ", data_filtered %>% nrow(), ")"))+
    theme_minimal()
  
  if (peril == "Hurricane"){
    p = p + scale_y_continuous(limits = c(NA, max(8e-7)))
  } else if (peril == "Hail"){
    p = p + scale_y_continuous(limits = c(NA, max(8e-5)))
  }
    
  plot_list[[i]] = p
  i = i+1
}

p = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, scale = 1) 
p

## Tail fitting (Generalised Pareto distribution)
par(mfrow = c(nrow,ncol))
for (reg in regions){
  data = storm_analysis %>%
    filter(CUSTOM_REGION == reg)
  
  GPD_fit = extRemes::fevd(x = data$DAMAGE, 
                               data = data.frame(data),
                               threshold = thresh_list[[reg]],
                               type = "GP", 
                               scale.fun = ~ YEAR,
                               shape.fun = ~ 1,
                               use.phi = T,
                               method = "MLE")
  print("######################################")
  print(reg)
  print(paste0("Number of rows: ", data %>% nrow()))
  summary(GPD_fit)
  print(ci.fevd(GPD_fit, type = "parameter", alpha = 0.05))
  
  # plot(GPD_fit, type = 'probprob')
  plot(GPD_fit, type = 'qq', main = paste0("Q-Q plot (", reg,")"), col = 'red', pch = 16) ## Exponential scaling to handle non-stationarity (minus thresh and divided by sigma_t)
  mtext(paste0("(n = ", data %>% filter(DAMAGE >= thresh_list[[reg]]) %>% nrow(),")"), line = 0.5, side = 3, adj = 0.9, cex = 0.75)
  print("######################################")
  saveRDS(GPD_fit, file = paste0("./R Code/Result_GPD_fit/", reg, ".rds"))
  
}


## Plotting the fitted densities/CDF throughout the years

if (peril == "Hail"){
  titles = c("North - Ohio Valley" = "North - Ohio Valley (n.s.)", "South" = "South (n.s.)", "Northwest" = "Northwest (n.s.)", "Southeast" = "Southeast", "Northeast" = "Northeast", "West" = "West (n.s.)")
  x = seq(range(storm_analysis %>% pull(DAMAGE))[1], range(storm_analysis %>% pull(DAMAGE))[2] *1.2, length.out = 3000)
} else if (peril == "Hurricane"){
  titles =  c("South" = "South (n.s.)", "Southeast" = "Southeast")
  x = 10^(seq(0, 13, length.out = 1000))
}



gamma_GPD_dist = function(x, alpha, beta, thresh, sigma, xi, cdf = FALSE){
  y = rep(-5, length(x))
  gamma_domain = x <= thresh
  
  if (cdf == FALSE){
    y[gamma_domain] = dgamma(x[gamma_domain], shape = alpha, rate = beta)
    y[!gamma_domain] = (1-pgamma(thresh, shape = alpha, rate = beta)) * devd(x[!gamma_domain]-thresh, scale = sigma, shape = xi, threshold = thresh, type = 'GP')
  } else {
    y[gamma_domain] = pgamma(x[gamma_domain], shape = alpha, rate = beta)
    y[!gamma_domain] = (1-last(y[gamma_domain])) * pevd(x[!gamma_domain], scale = sigma, shape = xi, threshold = thresh, type = 'GP') + last(y[gamma_domain])
  }
  
  return(y)
}
years = seq(0, 45)
n_years = storm_analysis %>% pull(YEAR) %>% unique() %>% length()

# Densities
plot_list = vector('list', length(regions))
i=1
for (reg in regions){
  data = storm_analysis %>% filter(CUSTOM_REGION == reg)
  n = nrow(data)
  gamma_fit = readRDS(paste0("./R Code/Result_gamma_fit/", reg, ".rds"))
  GPD_fit = readRDS(paste0("./R Code/Result_GPD_fit/", reg, ".rds"))
  
  densities = data.frame()
  for (year in years){
    shape <- exp(gamma_fit$par[["alpha_0"]] + year * gamma_fit$par[["alpha_1"]])
    rate <- exp(gamma_fit$par[["beta_0"]])
    
    thresh = thresh_list[[reg]]
    sigma = exp(GPD_fit$results$par[["phi0"]] + year * GPD_fit$results$par[["phi1"]])
    xi = GPD_fit$results$par[["shape"]]
    
    y = gamma_GPD_dist(x, shape, rate, thresh, sigma, xi, cdf = FALSE)
    
    densities = rbind(densities, data.frame(x = x, y = y, year = year+min_year))
  }
  
  p = ggplot()+
    geom_histogram(data = data, mapping = aes(x = DAMAGE, y = ..density..), fill = "lightgray", color = "black", alpha = 0.5) +
    geom_line(data = densities, mapping= aes(x = x, y = y, color = year, group = year))+
    # scale_color_viridis_c(option = "plasma")+
    labs(x = "Damage", y = "Density")+
    # scale_x_log10()+
    guides(color = guide_legend(title = "Year"))+
    theme_minimal()
  
  plot_list[[i]] = p
  i = i+1
}

p = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, scale = 1) 
p

# CDFs
plot_list = vector('list', length(regions))

i=1
for (reg in regions){
  data = storm_analysis %>% filter(CUSTOM_REGION == reg)
  n = nrow(data)
  

  gamma_fit = readRDS(paste0("./R Code/Result_gamma_fit/", reg, ".rds"))
  GPD_fit = readRDS(paste0("./R Code/Result_GPD_fit/", reg, ".rds"))
  
  cdfs = data.frame()
  for (year in years){
    shape <- exp(gamma_fit$par[["alpha_0"]] + year * gamma_fit$par[["alpha_1"]])
    rate <- exp(gamma_fit$par[["beta_0"]])
    
    thresh = thresh_list[[reg]]
    sigma = exp(GPD_fit$results$par[["phi0"]] + year * GPD_fit$results$par[["phi1"]])
    xi = GPD_fit$results$par[["shape"]]
    
    y = gamma_GPD_dist(x, shape, rate, thresh, sigma, xi, cdf = TRUE)
    
    cdfs = rbind(cdfs, data.frame(x = x, y = y, year = year+min_year))
  }
  
  p = ggplot()+
    geom_line(data = cdfs, mapping= aes(x = x, y = y, color = year, group = year))+
    # scale_color_viridis_c(option = "plasma")+
    labs(x = "Damage", y = "CDF", title = titles[[reg]])+
    guides(color = guide_legend(title = "Year"))+
    scale_x_log10() +
    theme_minimal()
  
  plot_list[[i]] = p
  i = i+1
}

p = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, scale = 1) 
p

## Return levels

m = seq(5,25, by = 1)

plot_list = vector("list", length = length(regions))
i=1
for (reg in regions){
  data = storm_analysis %>% filter(CUSTOM_REGION == reg)
  n = nrow(data)
  gamma_fit = readRDS(paste0("./R Code/Result_gamma_fit/", reg, ".rds"))
  GPD_fit = readRDS(paste0("./R Code/Result_GPD_fit/", reg, ".rds"))
  return_levels = data.frame()
  n_y = n / n_years
  for (year in years){
    # n_y = data %>% filter(YEAR == year) %>% nrow()
    shape <- exp(gamma_fit$par[["alpha_0"]] + year * gamma_fit$par[["alpha_1"]])
    rate <- exp(gamma_fit$par[["beta_0"]])
    
    thresh = thresh_list[[reg]]
    sigma = exp(GPD_fit$results$par[["phi0"]] + year * GPD_fit$results$par[["phi1"]])
    xi = GPD_fit$results$par[["shape"]]
    
    pu = (data %>% filter(DAMAGE > thresh) %>% nrow()) / data %>% nrow()
    
    if (!(pu > 1/min(m*n_y))){
      warning("return level too lwo, please revise it")
      n_y = n / n_years
    }
    
    # rl = exp(thresh + (sigma / xi) * ((m*n_y * pu)^xi -1) + min_BXCX)
    rl = thresh + (sigma / xi) * ((m*n_y * pu)^xi -1)
  return_levels = rbind(return_levels, data.frame(rl = rl, m = m, year = year+min_year))
    
  }
  
  p = ggplot()+
    geom_line(data = return_levels, mapping= aes(x = m, y = rl, color = year, group = year))+
    # scale_color_viridis_c(option = "plasma")+
    labs(x = "Return period (years)", y = "Damage", title = titles[[reg]])+
    scale_y_log10()+
    guides(color = guide_legend(title = "Year"))+
    theme_minimal()
  
  plot_list[[i]] = p
  i = i+1
}

p = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, scale = 1) 
p


### DAMAGE VISUALISATION ON THE MAP ###

## Geoplotting the damages
boundaries_US = read.csv("./Data/US_map/boundaries_US_mesh.csv")
boundaries_US = boundaries_US[nrow(boundaries_US):1,]

custom_colors <- c("blue", "green", "red")

ggplot() +
  geom_path(data = boundaries_US, mapping = aes(x = X_CART, y = Y_CART))+
  geom_point(data = storm_analysis, mapping = aes(x = X_CART, y = Y_CART, col = DAMAGE_LOG)) +
  scale_color_gradientn(colors = custom_colors)+
  theme_minimal()

### Example gamma/GPD mixture
CDF = TRUE

x = seq(0, 15, length.out = 1000)

y = gamma_GPD_dist(x, alpha = 2, beta = 0.5, thresh = 10, sigma = 1, xi = 3, cdf = CDF)


p = ggplot()+
  geom_line(data = data.frame(x = x, y = y), mapping = aes(x = x ,y = y), col = "black", linewidth = 0.75)+
  labs(title = "Gamma/GPD mixture", x = "x", y = "CDF")+
  geom_vline(xintercept = 10, linetype = "dashed", linewidth = 0.5, col = 'red')+
  theme_minimal()

if (!CDF){
  p = p + 
    ggpubr::geom_bracket(xmin = 0, xmax = 9.8, y.position = 0.21, label = "Bulk (gamma distribution)")+
    ggpubr::geom_bracket(xmin = 10.2, xmax = 15, y.position = 0.21, label = "Tail (GPD)")+
    ylim(0, 0.23) 
} else{
  p = p + 
    ggpubr::geom_bracket(xmin = 0, xmax = 9.8, y.position = 1.03, label = "Bulk (gamma distribution)")+
    ggpubr::geom_bracket(xmin = 10.2, xmax = 15, y.position = 1.03, label = "Tail (GPD)")+
    ylim(0, 1.05) 
}

p



