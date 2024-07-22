
def latlong_to_cartesian(storm_dataset, lat = "lat", long = "long", R = 6371):
    # No projection
    storm_dataset['X_CART'] = R * np.cos(np.radians(storm_dataset[lat])) * np.cos(np.radians(storm_dataset[long]))
    storm_dataset["Y_CART"] = R * np.cos(storm_dataset[lat]) * np.sin(storm_dataset[long])
    storm_dataset["Z_CART"] = R * np.sin(storm_dataset[lat])

    return storm_dataset



def latlon_to_cartesian_albers(storm_dataset, lat = "lat", long = "long"):
    # Albers Equal-Area projection
    albers_equal_area = pyproj.Proj(proj='aea', lat_1=29.5, lat_2=45.5, lat_0=24.44, lon_0=-124.615, datum='WGS84')
    _, max_north = albers_equal_area(-95.1533, 49.39)
    max_east, _ = albers_equal_area(-66.94, 44.815)
    x, y = albers_equal_area(storm_dataset[long], storm_dataset[lat])
    storm_dataset["X_CART"] = x/max_east
    storm_dataset["Y_CART"] = y/max_north
    return storm_dataset


def extract_points(geometry):
    if isinstance(geometry, Polygon):
        return list(geometry.exterior.coords)
    elif isinstance(geometry, MultiPolygon):
        points = []
        for polygon in geometry.geoms:
            points.extend(list(polygon.exterior.coords[:-1]))
        return points
    else:
        return []

if __name__ == "__main__":
    import pandas as pd
    import numpy as np
    import pyproj
    import geoplot as gplt
    import geopandas as gpd
    from shapely import Point, Polygon, MultiPolygon

    print("Loading datasets and mapping tables")

    # Loading tables
    FIPS_mapping = pd.read_csv("./Data/mapping_tables/FIPS_coordinate_mapping_2.csv")
    storm_data = pd.read_csv("./Data/Prod_datasets/Storm_events_details_full_raw.csv")

    print("Starting the data processing routine.")

    storm_data["US_FIPS"] = (storm_data["STATE_FIPS"] * 1000 + storm_data["CZ_FIPS"]).astype("Int64").astype('str')
    FIPS_mapping["US_FIPS"] = FIPS_mapping["US_FIPS"].astype("Int64").astype('str')

    # # Mapping the FIPS code to their coordinates
    storm_data = storm_data.join(FIPS_mapping[["CZ_TYPE", "US_FIPS", "lat", "long"]].set_index(["US_FIPS", "CZ_TYPE"]), on = ["US_FIPS", "CZ_TYPE"])
    storm_data["average_lat"] = (storm_data["BEGIN_LAT"] + storm_data["END_LAT"])/2
    storm_data["average_long"] = (storm_data["BEGIN_LON"] + storm_data["END_LON"])/2
    storm_data["lat_fin"] = storm_data["average_lat"].where(~storm_data["average_lat"].isna(), storm_data["BEGIN_LAT"])
    storm_data["long_fin"] = storm_data["average_long"].where(~storm_data["average_long"].isna(), storm_data["BEGIN_LON"])
    storm_data["lat_fin"] = storm_data["lat_fin"].where(~storm_data["lat_fin"].isna(), storm_data["lat"])
    storm_data["long_fin"] = storm_data["long_fin"].where(~storm_data["long_fin"].isna(), storm_data["long"])
    storm_data.drop(columns=["average_lat", "average_long", "lat", "long"], inplace = True)
    storm_data.rename(columns = {"lat_fin":'lat', "long_fin":"long"}, inplace = True)

    ## Adding the Cartesian coordinates
    storm_data = latlon_to_cartesian_albers(storm_data, lat = "lat", long = "long")

    # Creating the USA border maps for the US mesh
    path = gplt.datasets.get_path("contiguous_usa")
    contiguous_usa = gpd.read_file(path)

    contiguous_usa["COUNTRY"] = "USA"
    contiguous_usa = contiguous_usa.dissolve(by = "COUNTRY")

    all_points = []

    for geometry in contiguous_usa.geometry:
        points = extract_points(geometry)
        all_points.extend(points)

    longitude = [point[0] for point in all_points]
    latitude = [point[1] for point in all_points]

    boundaries_mesh = pd.DataFrame({"long": longitude, "lat": latitude}).iloc[:-1, :]
    boundaries_mesh = latlon_to_cartesian_albers(boundaries_mesh, lat = "lat", long = "long")

    # Creating the evaluation mesh for the fdaPDE package
    usa_polygon = contiguous_usa.geometry.iloc[0]
    min_x, min_y, max_x, max_y = usa_polygon.bounds

    spacing = .25 # degree spacing for grid points
    x_coords = np.arange(np.floor(min_x), np.ceil(max_x), spacing)
    y_coords = np.arange(np.floor(min_y), np.ceil(max_y), spacing)

    grid_points = []
    for x in x_coords:
        for y in y_coords:
            point = Point(x, y)
            if usa_polygon.contains(point):
                grid_points.append((point.x, point.y))

    # Create a DataFrame from grid points
    eval_mesh = pd.DataFrame({"long": [point[0] for point in grid_points], "lat": [point[1] for point in grid_points]})
    eval_mesh = latlon_to_cartesian_albers(eval_mesh, lat = "lat", long = "long")
    
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
    boundaries_mesh.to_csv("./Data/US_map/boundaries_US_mesh.csv", index = False)
    eval_mesh.to_csv("./Data/US_map/evaluation_US_mesh.csv", index = False)

    print("Data successfully processed and saved in CSV file in ./Data/Prod_datasets/")


