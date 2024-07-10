import pandas as pd
import os
import numpy as np
from tqdm import tqdm
import re

import ftplib
import gzip
import shutil

import pyproj
import geoplot as gplt
import geopandas as gpd
from shapely import Point, Polygon, MultiPolygon


def run_file(filepath):
    with open(filepath) as file:
        exec(file.read())


if __name__ == "__main__":
    # Downloading the NOAA files and loading them to clean them
    run_file("./Python Code/data_loading_pre_processing.py")

    # Processing the data
    run_file("./Python Code/data_processing.py")
