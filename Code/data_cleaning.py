import pandas as pd
import os
import numpy as np
from tqdm import tqdm
import re

import ftplib
import gzip
import shutil




def clean_detailed_dataset(data):
    data.drop(columns = ["BEGIN_YEARMONTH", "BEGIN_DAY", "BEGIN_TIME", "END_YEARMONTH", "END_DAY", "END_TIME"], inplace = True)
    data["BEGIN_DATE_TIME"] = pd.to_datetime(data["BEGIN_DATE_TIME"], format = "mixed")
    data["END_DATE_TIME"] = pd.to_datetime(data["END_DATE_TIME"], format = "mixed")

    def convert_amounts(col_name, df = data):
        df[col_name] = df[col_name].astype(str)
        df[col_name] = df[col_name].apply(lambda x: re.sub(r"[^\w\s\.]", "", x))
        amount = df[col_name].str.split("(K|M|B|H|k|m|b|h)", regex = True, expand = True)
        print(amount)
        amount.replace(to_replace = [None, "", np.nan, "nan"], value = 0, inplace = True)

        def abreviation(val):
            if type(val) != type("str") or val != val: #Nan values
                return 0
            if val.lower() == "h":
                return 1e2
            if val.lower() == "k":
                return 1e3
            if val.lower() == "m":
                return 1e6
            if val.lower() == "b":
                return 1e9
        
        def force_convert_float(val):
            try:
                return float(val)
            except:
                return 0
            
        print(amount[0])

        if amount.shape[1] > 1:
            amount[1] = amount[1].apply(abreviation)

            amount[0] = amount[0].apply(force_convert_float)
            amount[1] = pd.to_numeric(amount[1])

            return(amount[0] * amount[1])
        
        else:
            return(pd.to_numeric(amount[0]))
        
    data['DAMAGE_PROPERTY'] = convert_amounts("DAMAGE_PROPERTY")
    data['DAMAGE_CROPS'] = convert_amounts("DAMAGE_CROPS")

    data = data.join(inflation_traj.set_index("YEAR"), on = "YEAR")

    return(data)

   


def clean_fatalities_dataset(data):
    data.drop(columns = ["FAT_YEARMONTH", "FAT_DAY", "FAT_TIME", "EVENT_YEARMONTH"], inplace = True)

    data["FATALITY_DATE"] = pd.to_datetime(data["FATALITY_DATE"], format = "mixed")
    
    return(data)



if __name__ == "__main__":

    ### Downloading data

    ftp = "ftp://ftp.ncei.noaa.gov/pub/data/swdi/stormevents/csvfiles/"
    hostname = "ftp.ncei.noaa.gov"
    remote_dir = "pub/data/swdi/stormevents/csvfiles/"

    ftp_server = ftplib.FTP(hostname, "", "")
    ftp_server.encoding = "utf-8"

    ftp_server.login()

    ftp_server.cwd("/pub/data/swdi/stormevents/csvfiles/")

    files = ftp_server.nlst()

    for file in tqdm(files):
        if file != "legacy":
            ftp_server.retrbinary('RETR {}'.format(file), open("./NOAA_Storm_events_dataset_zip/{}".format(file), "wb").write)
            if ".gz" in file:
                with gzip.open("./NOAA_Storm_events_dataset_zip/{}".format(file), 'rb') as file_in:
                    with open(re.sub(".gz", "", "./NOAA_Storm_events_dataset/{}".format(file)), 'wb') as file_out:
                        shutil.copyfileobj(file_in, file_out)
            else:
                with open("./NOAA_Storm_events_dataset_zip/{}".format(file), 'rb') as file_in:
                    with open(re.sub(".gz", "", "./NOAA_Storm_events_dataset/{}".format(file)), 'wb') as file_out:
                        shutil.copyfileobj(file_in, file_out)
            os.remove("./NOAA_Storm_events_dataset_zip/{}".format(file))

    ftp_server.quit()

    ### Cleaning data

    ## Inlfation adjustement
    inflation_mapping = pd.read_csv("Prod_datasets/inflation_rates_US.csv")
    inflation_mapping["INFLATION_INDEX"] = np.append(np.array([1]), np.cumprod(inflation_mapping[" value"]/100 +1)[1:])[::-1]
    inflation_mapping["YEAR"] = inflation_mapping["date"].apply(lambda x: int(x[0:4]))

    inflation_traj = inflation_mapping[["YEAR", "INFLATION_INDEX"]]

    inflation_traj.to_csv("Prod_datasets/inflation_traj_US.csv", index = False)

    ## Cleaning datasets

    path = "NOAA_Storm_events_dataset"

    for file_name in tqdm(os.listdir(path)):
        print(file_name)
        _, ext = os.path.splitext(path+"/"+file_name)
        if ext == ".csv":
            data = pd.read_csv(path+"/"+file_name)
            if "details" in file_name:
                data = clean_detailed_dataset(data)
                name = "StormEvents_details_" + file_name[29:44]
            elif "fatalities" in file_name:
                data = clean_fatalities_dataset(data)
                name = "StormEvents_fatalities_" + file_name[32:47]
            elif "locations" in file_name:
                data = data
                name = "StormEvents_locations_" + file_name[31:46]
            else:
                print("Unknown dataset type: {}".format(file_name))

            data.to_csv("NOAA_Storm_events_clean/"+name+".csv") 
            os.remove(path+"/"+file_name)

    ### Fusion of the yearly datasets
    
    path = "NOAA_Storm_events_clean"
    datasets = [i for i in os.listdir(path) if "details" in i]

    type_dict = {'CATEGORY': 'str',
                    'TOR_OTHER_CZ_STATE': 'str',
                    'TOR_OTHER_CZ_FIPS': 'str',
                    'BEGIN_RANGE': 'str',
                    'MAGNITUDE' : 'str',
                    'EVENT_NARRATIVE': 'str',
                    'FLOOD_CAUSE': 'str'}
    
    storm_data = pd.DataFrame()

    for i, dataset  in enumerate(datasets):
        data = pd.read_csv(path + "/" + dataset, dtype=type_dict)
        storm_data = pd.concat([storm_data, data], axis = 0)
        os.remove(path+"/"+dataset)
    
    storm_data.drop(columns="Unnamed: 0", inplace = True)
    storm_data.to_csv("Prod_datasets/Storm_events_details_full_clean.csv", index = False)