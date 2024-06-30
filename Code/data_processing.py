


if __name__ == "__main__":
    import pandas as pd

    print("Loading datasets and mapping tables")

    # Loading tables
    FIPS_mapping = pd.read_csv("./Data/mapping_tables/FIPS_coordinate_mapping.csv")
    storm_data = pd.read_csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv")

    print("Starting the data processing routine.")

    storm_data["US_FIPS"] = (storm_data["STATE_FIPS"] * 1000 + storm_data["CZ_FIPS"]).astype("Int64").astype('str')
    FIPS_mapping["US_FIPS"] = FIPS_mapping["US_FIPS"].astype("Int64").astype('str')

    # Mapping the FIPS code to their coordinates
    storm_data = storm_data.join(FIPS_mapping[["CZ_TYPE", "US_FIPS", "lat", "long"]].set_index(["US_FIPS", "CZ_TYPE"]), on = ["US_FIPS", "CZ_TYPE"])
    storm_data["average_lat"] = (storm_data["BEGIN_LAT"] + storm_data["END_LAT"])/2
    storm_data["average_long"] = (storm_data["BEGIN_LON"] + storm_data["END_LON"])/2
    storm_data["lat_fin"] = storm_data["average_lat"].where(~storm_data["average_lat"].isna(), storm_data["BEGIN_LAT"])
    storm_data["long_fin"] = storm_data["average_long"].where(~storm_data["average_long"].isna(), storm_data["BEGIN_LON"])
    storm_data["lat_fin"] = storm_data["lat_fin"].where(~storm_data["lat_fin"].isna(), storm_data["lat"])
    storm_data["long_fin"] = storm_data["long_fin"].where(~storm_data["long_fin"].isna(), storm_data["long"])
    storm_data.drop(columns=["average_lat", "average_long", "lat", "long"], inplace = True)
    storm_data.rename(columns = {"lat_fin":'lat', "long_fin":"long"}, inplace = True)

    # Discounting amounts to take inflation into account
    storm_data["DISCOUNT_DAMAGE_PROPERTY"] = storm_data["INFLATION_INDEX"] * storm_data["DAMAGE_PROPERTY"]
    storm_data["DISCOUNT_DAMAGE_CROPS"] = storm_data["INFLATION_INDEX"] * storm_data["DAMAGE_CROPS"]

    # Narrowing the categories of events
    storm_events = storm_data["EVENT_TYPE"].unique()
    storm_events.sort()
    mapping_table_events = {"EVENT_TYPE": list(storm_events), "EVENT_CAT": ["Astronomical Low Tide", "Avalanche", "Blizzard", "Submersion", "Wind Chill", "Debris Flow", "Fog", "Fog", "Drought", "Dust Devil", "Storm", "Heat", "Wind Chill", "Flood", "Flood", "Fog", "Frost", "High Wind", "Flood", "Hail", "Hail", "Heat", "Heavy Rain", "Heavy Snow", "High Surf", "High Wind", "Hurricane", "Hurricane", "Blizzard", 'Heavy Snow', "Flood", "Lightning", "Fog", "Hail", "High Wind", "Hurricane", "Lightning", "High Wind", "Thunderstorm", "Tropical Storm", "Tropical Storm", "Northern Lights", "Rip Current", "Submersion", "Heavy Snow", "Submersion", "Submersion", "High Wind", "Thunderstorm", "Thunderstorm", "Thunderstorm", "Thunderstorm", "Thunderstorm", "Thunderstorm", "Thunderstorm", "Thunderstorm", "Thunderstorm", "Tornado", "Tornado", "Thunderstorm", "Tornado", "Tropical Depression", "Tropical Storm", "Submersion", "Volcanic Ash", "Volcanic Ash", "Waterspout", "Wildfire", "Blizzard", "Winter Weather"]}
    mapping_table_events = pd.DataFrame(mapping_table_events)

    storm_data = storm_data.join(mapping_table_events.set_index("EVENT_TYPE"), on  = "EVENT_TYPE")

    print("Data processing done. Saving data in CSV file.")

    if any(storm_data["EVENT_CAT"].isna()):
        print("/!\ Warning: some of the meteorological event type were not succesfully mapped to a broader category.")

    storm_data.to_csv("./Data/Prod_datasets/Storm_events_details_full_clean.csv", index = False)

    print("Data successfully processed and saved in CSV file in ./Data/Prod_datasets/")


