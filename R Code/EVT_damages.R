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

storm_analysis_raw = storm_analysis_raw[complete.cases(storm_analysis_raw %>% select(c("EPISODE_ID", "X_CART", "Y_CART"))),]


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

US_regions = read.csv("./Data/mapping_tables/US_STATE.csv") %>%
  mutate(State = toupper(State)) %>%
  rename(c("REGION" = "Region", "STATE" = "State"))


lambdas = seq(0, 0.3, length.out = 10)

taus = seq(0.7, 0.99, length.out = 15)


AICs = matrix(nrow = length(lambdas), ncol = length(taus))
log_lik = AICs
BICs = AICs
AICs_corrected = AICs

lambdas = c(1)
taus = c(0.90)

for (lambda_bxcx in lambdas){
  i = which(lambdas == lambda_bxcx)
  for (tau in taus){
    j = which(taus == tau)
    
    storm_analysis = storm_analysis_raw %>% 
      mutate(quart = quarter(BEGIN_DATE_TIME)) %>%
      mutate(MONTH = month(BEGIN_DATE_TIME)) %>%
      # left_join(GDP_state %>% select(c("GeoFips", "GDP", "YEAR")), by = c("YEAR" = "YEAR", "STATE_FIPS" = "GeoFips")) %>%
      left_join(SOI_data, by = c("YEAR", "MONTH")) %>%
      left_join(HPI_state %>% select(c("yr", "period", "place_name", "HPI")), by = c("YEAR" = "yr", "quart" = "period", "STATE" = "place_name")) %>%
      left_join(US_regions %>% select(c("STATE", "REGION")), by = "STATE") %>%
      mutate(REGION = factor(REGION)) %>%
      mutate(DAMAGE_LOG = log10(DAMAGE_PROPERTY)) %>%
      mutate(DAMAGE_LN = log(DAMAGE_PROPERTY)) %>%
      mutate(STATE_FIPS = factor(STATE_FIPS)) %>%
      group_by(REGION, EVENT_CAT) %>%
      mutate(DAMAGE_SCALED = scale(DAMAGE_PROPERTY, center = FALSE)[,1]) %>%
      mutate(thresh_log = quantile(DAMAGE_LOG, 0.9),
             thresh = quantile(DAMAGE_PROPERTY, 0.9),
             thresh_ln = quantile(DAMAGE_LN, 0.9),
             n = n()) %>%
      # mutate(thresh_scaled = 5 * sd(DAMAGE_SCALED)) %>%
      filter(!is.nan(DAMAGE_SCALED)) %>%
      # ungroup() %>%
      mutate(DAMAGE_BOX_COX = sae::bxcx(DAMAGE_PROPERTY/sd(DAMAGE_PROPERTY), lambda_bxcx)) %>%
      # mutate(SD = sd(DAMAGE_PROPERTY)) %>%
      mutate(MEAN = mean(DAMAGE_PROPERTY)) %>%
      mutate(MED = median(DAMAGE_PROPERTY)) %>%
      mutate(BEGIN_DATE_TIME = as.Date(BEGIN_DATE_TIME)) %>%
      mutate(SCALER = sqrt(mean(DAMAGE_PROPERTY^2))) %>%
      filter(!is.na(REGION)) 
      # ungroup() %>%
      # mutate(test_split = runif(n()) > 0.8)
    
    # storm_analysis = storm_analysis %>% filter(!test_split)
    # storm_analysis_test = storm_analysis %>% filter(test_split)
    
    ## Quantile regression
    
    quant_reg = rq(DAMAGE_BOX_COX ~ REGION, data = storm_analysis, tau = tau)
    summary(quant_reg)  
    quant_reg
    
    storm_analysis = storm_analysis %>%
      ungroup() %>%
      mutate(thresh_quant = quant_reg$fitted.values)
    
    print(paste0("Lambda = ", lambda_bxcx, ", Tau = ", tau))
    print(storm_analysis %>% mutate(test = DAMAGE_BOX_COX > thresh_quant) %>% pull() %>% sum())
    
    
    res_GPD_cov = fevd(x = storm_analysis$DAMAGE_BOX_COX, 
                       data = data.frame(storm_analysis),
                       threshold = storm_analysis$thresh_quant,
                       type = "GP", 
                       scale.fun = ~ REGION,
                       shape.fun = ~ REGION,
                       method = "MLE")
    
    k = (res_GPD_cov$results$par %>% length())
    n = storm_analysis %>% mutate(test = DAMAGE_BOX_COX > thresh_quant) %>% pull() %>% sum()
    
    AICs[i,j] = res_GPD_cov$results$value + 2*k
    log_lik[i,j] = res_GPD_cov$results$value 
    BICs[i,j] = res_GPD_cov$results$value + log(nrow(storm_analysis))*(res_GPD_cov$par.models$scale %>% length() -1)
    AICs_corrected[i,j] = AICs[i,j] + 2*k*(k+1)/(n-k-1)
    
  }
}

df <- reshape2::melt(AICs_corrected) %>%
  mutate(lambda = lambdas[Var1],
         tau = taus[Var2]) 
custom_colors <- c("blue", "green", "red")
ggplot(df, aes(lambda, tau, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = custom_colors)+
  # scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Color Plot of Matrix", x = "lambda", y = "tau") 
  theme_minimal()

storm_analysis %>% pull(DAMAGE_BOX_COX) %>% hist(breaks = 25)

  

## Linear regression
# simple_OLS = lm(LN_DAMAGE_PROPERTY ~ log(GDP), data = storm_analysis)
# summary(simple_OLS)

## GLM
# glm_fit = glm(LOG_DAMAGE_PROPERTY ~  STATE_FIPS, data = storm_analysis, family = inverse.gaussian())
# summary(glm_fit)

## GLM with EVD

# x = storm_analysis %>%
#   group_by(YEAR, STATE_FIPS) %>%
#   mutate(BOOL = DAMAGE_PROPERTY == max(DAMAGE_PROPERTY)) %>%
#   filter(BOOL)
  # summarise(DAMAGE_PROPERTY = max(DAMAGE_PROPERTY), LOG_DAMAGE_PROPERTY = max(log(DAMAGE_PROPERTY)), GDP = first(GDP), HPI = mean(HPI)) #%>%
  # mutate(GDP = scale(GDP)) %>%
  # mutate(DAMAGE_PROPERTY = scale(DAMAGE_PROPERTY)) %>%
  # mutate(LOG_DAMAGE_PROPERTY= scale(LOG_DAMAGE_PROPERTY))

# res_GEV_cov = fevd(x = LOG_DAMAGE_PROPERTY, data.frame(x), location.fun = ~HPI + log(GDP) + SOI , type = 'GEV', method = "MLE")
# res_GEV_cov
# ci.fevd(res_GEV_cov, type = "parameter")
# plot(res_GEV_cov)
# 
# x$test = (x$LOG_DAMAGE_PROPERTY - res_GEV_cov$results$par[1] - res_GEV_cov$results$par[2]*log10(x$GDP)) / res_GEV_cov$results$par[3]
# 
# res_GEV_nocov = fevd(x = LOG_DAMAGE_PROPERTY, data.frame(x %>% group_by(YEAR) %>% summarise(LOG_DAMAGE_PROPERTY = max(LOG_DAMAGE_PROPERTY))), type = "GEV")
# res_GEV_nocov
# ci.fevd(res_GEV_nocov, type = "parameter")
# plot(res_GEV_nocov)

## GLM with POT

storm_analysis %>% pull(DAMAGE_SCALED) %>% plot()
storm_analysis %>% pull(thresh_scaled) %>% lines(type = 'l', col = 'red')

VGAM::meplot(storm_analysis %>% filter(STATE == "KENTUCKY") %>% pull(DAMAGE_BOX_COX))

res_GPD_cov = fevd(x = storm_analysis$DAMAGE_BOX_COX, data = data.frame(storm_analysis), threshold = storm_analysis$thresh_quant, type = "GP", scale.fun = ~ 1, method = "MLE")
res_GPD_cov
ci.fevd(res_GPD_cov, type = "parameter", alpha = 0.1)
plot(res_GPD_cov, type = 'probprob')
plot(res_GPD_cov, type = 'qq')
plot(res_GPD_cov, type = "qq2")
plot(res_GPD_cov, type = 'rl', rperiods = c(1.5))
plot(res_GPD_cov, type = 'hist', hist.args = c("breaks" = 25, "freq" = F), main = "", xlab= "x")
plot(res_GPD_cov, type = "density")
# plot(res_GPD_cov, rperiods = c(10))
# par(mfrow = c(1,1))

regions = storm_analysis %>% pull(REGION) %>% unique() %>% sort()

### Find good coefficients for the good region


reg = "South"
i = which(regions == reg) - 1
q = storm_analysis %>% filter(REGION == reg) %>% pull(DAMAGE_BOX_COX) %>% quantile(tau)
data = storm_analysis %>% filter(REGION == reg, DAMAGE_BOX_COX >= thresh_quant) %>%
  pull(DAMAGE_BOX_COX)

print(length(data))


scale = res_GPD_cov$results$par[1]
if ((i !=0) & (res_GPD_cov$par.models$term.names$scale %>% length() == 1)){
  scale = scale + res_GPD_cov$results$par[[paste0("sigma", i)]]
}

if (res_GPD_cov$par.models$term.names$shape%>% length() != 1){
  shape = res_GPD_cov$results$par[["shape"]]
}else{
  shape = res_GPD_cov$results$par[["xi0"]]
}

if ((i !=0) & (res_GPD_cov$par.models$term.names$shape%>% length() == 1)){
  shape = shape + res_GPD_cov$results$par[[paste0("xi", i)]]
}

x = seq(0, 6, length.out = 250)
plot(x+q, ecdf(data)(x+q), main = paste0(reg, ", n = ", length(data)), ylim = c(0,1))
lines(x+q,
      pevd(x+q,
           threshold = q,
           scale = scale,
           shape = shape,
           type = "GP"),
        col = "red",
      lwd = 2)


m = 10 * 365.25
exced_rate = storm_analysis %>% filter(REGION == reg, DAMAGE_BOX_COX > thresh_quant) %>% nrow() / storm_analysis %>% nrow() 
return_level = predict(quant_reg, newdata = data.frame(REGION = regions[i+1])) + (scale / shape) * ((m * exced_rate)^shape-1)



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
  geom_point(data = storm_analysis, mapping = aes(x = X_CART, y = Y_CART, col = DAMAGE_LOG)) +
  scale_color_gradientn(colors = custom_colors)


