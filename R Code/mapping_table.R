library(tidyverse)
NWS_forecast_zone <- read.csv("Data/mapping_tables/NWS_forecast_zone.csv") 
US_State = read.csv("Data/mapping_tables/US_STATE.csv")
CZ_FIPS_coordinates <- read.csv("Data/mapping_tables/CZ_FIPS_coordinates.csv")

NWS_forecast_zone = NWS_forecast_zone %>%
  left_join(US_State %>% select(c("Abbreviation", "FIPS")), by = c("STATE" = "Abbreviation")) %>%
  mutate(US_FIPS = FIPS * 1000 + as.numeric(ZONE)) %>%
  mutate(CZ_TYPE = "Z") %>%
  rename(c("lat" = "LAT", "long" = "LON")) %>%
  filter(!is.na(US_FIPS))

FIPS_coordinates = CZ_FIPS_coordinates %>% rbind(NWS_forecast_zone %>% select(c("CZ_TYPE", "US_FIPS", "lat", "long")))

write.csv(FIPS_coordinates, "Data/mapping_tables/FIPS_coordinates_mapping_2.csv", row.names = F)                                                 
