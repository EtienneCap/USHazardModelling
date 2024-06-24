import pandas as pd
import geopandas as gpd



if __name__ == "__main__":
    ## Mapping the location using the FIPS encoding

    FIPS_coor = pd.read_json("./Data/US_map/fips_map.json").T.reset_index()
    FIPS_coor.rename(columns = {"index": "US_FIPS"}, inplace = True)
    FIPS_coor["US_FIPS"] = FIPS_coor["US_FIPS"].astype("str")
    FIPS_coor["CZ_TYPE"] = "C"

    loc = gpd.read_file("./Data/US_map/NWS Forecast zone/z_05mr24.shp").drop(columns = "geometry").rename(columns = {"LON": "long", "LAT":"lat"}).drop_duplicates(subset = "STATE_ZONE")
    loc["CZ_TYPE"] = "Z"

    State_FIPS = pd.read_csv("./Data/US_map/StatesFIPSCodes.csv").drop(columns = "STATENS").rename(columns = {"STUSAB":"STATE"})

    loc = loc.join(State_FIPS.set_index("STATE"), on = "STATE").dropna(axis = 0, subset = "STATE_FIPS")

    loc = loc[["STATE", "STATE_FIPS", "CZ_TYPE","ZONE", "NAME", "lat", "long"]]
    loc["US_FIPS"] = (loc["STATE_FIPS"] * 1000 + loc["ZONE"].astype(float)).astype("Int64").astype('str')


    FIPS_mapping = pd.concat([loc[["CZ_TYPE","US_FIPS", "lat", "long"]], FIPS_coor[["CZ_TYPE","US_FIPS", "lat", "long"]]])

    FIPS_mapping.to_csv("./Data/mapping_tables/FIPS_coordinate_mapping.csv", index = False)

    print("FIPS mapping table created and stored in ./Data/Prod_datasets/FIPS_coordinate_mapping.csv")