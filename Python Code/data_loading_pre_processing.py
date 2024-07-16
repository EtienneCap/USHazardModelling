def clean_detailed_dataset(data, inflation_traj):
    data["BEGIN_DATE"] = data["BEGIN_YEARMONTH"].astype(str).str.slice(0,4) + "-" +  data["BEGIN_YEARMONTH"].astype(str).str.slice(4,6) + "-" + data["BEGIN_DAY"].astype(str) + " " + data["BEGIN_DATE_TIME"].str.slice(9, 18)
    data["END_DATE"] = data["END_YEARMONTH"].astype(str).str.slice(0,4) + "-" +  data["END_YEARMONTH"].astype(str).str.slice(4,6) + "-" + data["END_DAY"].astype(str) + " " + data["END_DATE_TIME"].str.slice(9, 18)
    data.drop(columns = ["BEGIN_YEARMONTH", "BEGIN_DAY", "BEGIN_TIME", "END_YEARMONTH", "END_DAY", "END_TIME"], inplace = True)
    data["BEGIN_DATE_TIME"] = pd.to_datetime(data["BEGIN_DATE"], format = "%Y-%m-%d %H:%M:%S")
    data["END_DATE_TIME"] = pd.to_datetime(data["END_DATE"], format = "%Y-%m-%d %H:%M:%S")

    # Copy of the original data
    data["DAMAGE_PROPERTY_ORIGINAL"] = data["DAMAGE_PROPERTY"]
    data["DAMAGE_CROPS_ORIGINAL"] = data["DAMAGE_CROPS"]


    def convert_amounts(col_name, df = data):
        df[col_name] = df[col_name].astype(str)
        df[col_name] = df[col_name].apply(lambda x: re.sub(r"[^\w\s\.]", "", x))
        amount = df[col_name].str.split("(K|M|B|H|k|m|b|h)", regex = True, expand = True)
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

    import pandas as pd
    import os
    import numpy as np
    from tqdm import tqdm
    import re
    import sys

    import ftplib
    import gzip
    import shutil

    ### Downloading data

    print("Starting downlading and unzipping data.")

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
            try:
                ftp_server.retrbinary('RETR {}'.format(file), open("./Data/NOAA_Storm_events_dataset_zip/{}".format(file), "wb").write)
            except error as e:
                print(e)
                print(file)
                sys.exit("Error in file downloading")
            if ".gz" in file:
                with gzip.open("./Data/NOAA_Storm_events_dataset_zip/{}".format(file), 'rb') as file_in:
                    with open(re.sub(".gz", "", "./Data/NOAA_Storm_events_dataset/{}".format(file)), 'wb') as file_out:
                        shutil.copyfileobj(file_in, file_out)
            else:
                with open("./Data/NOAA_Storm_events_dataset_zip/{}".format(file), 'rb') as file_in:
                    with open(re.sub(".gz", "", "./Data/NOAA_Storm_events_dataset/{}".format(file)), 'wb') as file_out:
                        shutil.copyfileobj(file_in, file_out)
            os.remove("./Data/NOAA_Storm_events_dataset_zip/{}".format(file))
    print("Downlading and unzipping data finished.")

    ftp_server.quit()

    ### Cleaning/augmenting data

    ## Inlfation adjustement

    print("Importation of the inflation information")

    inflation_mapping = pd.read_csv("./Data/Prod_datasets/inflation_rates_US.csv")
    inflation_mapping["INFLATION_INDEX"] = np.append(np.array([1]), np.cumprod(inflation_mapping[" value"]/100 +1)[1:])[::-1]
    inflation_mapping["YEAR"] = inflation_mapping["date"].apply(lambda x: int(x[0:4]))

    inflation_traj = inflation_mapping[["YEAR", "INFLATION_INDEX"]]

    inflation_traj.to_csv("./Data/mapping_tables/inflation_traj_US.csv", index = False)

    ## Cleaning datasets

    print("Pre-processing the datasets.")

    path = "./Data/NOAA_Storm_events_dataset"

    for file_name in (pbar := tqdm(os.listdir(path), desc = "File:")):
        pbar.set_description(f" {file_name}")
        pbar.refresh()
        _, ext = os.path.splitext(path+"/"+file_name)
        if ext == ".csv":
            data = pd.read_csv(path+"/"+file_name)
            
            if "details" in file_name:
                data = clean_detailed_dataset(data, inflation_traj)
                name = "StormEvents_details_" + file_name[29:44]
                if(any(data.duplicated())):
                    raise "Stop duplicates"
            elif "fatalities" in file_name:
                data = clean_fatalities_dataset(data)
                name = "StormEvents_fatalities_" + file_name[32:47]
            elif "locations" in file_name:
                data = data
                name = "StormEvents_locations_" + file_name[31:46]
            else:
                print("Unknown dataset type: {}".format(file_name))
            os.remove(path+"/"+file_name)
            data.to_csv("./Data/NOAA_Storm_events_clean/"+name+".csv") 
            
    

    ### Fusion of the yearly datasets

    print("Merging the datasets together.")
    
    path = "./Data/NOAA_Storm_events_clean"
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
        os.remove(path+"/"+dataset)
        storm_data = pd.concat([storm_data, data], axis = 0)
        

    print("Pre-processing the datasets done.")
    try:
        storm_data.drop(columns="Unnamed: 0", inplace = True)
    except:
        a = 0
    print("Saving the dataset into a CSV file.")
    storm_data.to_csv("./Data/Prod_datasets/Storm_events_details_full_raw.csv", index = False)

    print("Data cleaning and pre-processing finished.")