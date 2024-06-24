import os


def run_file(filepath):
    with open(filepath) as file:
        exec(file.read())


if __name__ == "__main__":
    # Downloading the NOAA files and loading them to clean them
    run_file("./Code/data_loading_pre_processing.py")

    # Processing the data
    run_file("./Code/data_processing.py")
