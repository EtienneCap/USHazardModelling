# USHazardModelling
Master thesis of my MSc Statisics @ Imperial College London

## Project description

This project focuses on the modelisation of the damages caused by natural disaster in the US using the [Storm Event Dataset](https://www.ncdc.noaa.gov/stormevents/) from the NOAA. It aim to include:
- Modelisation of the space-time properties of some types of frequent natural disasters in the US.
- Estimation of the distribution of the damages caused by these natural catastrophes. 

## Package versions

### R (version at least 4.3.1)
R packages requirements:
  - tidyverse (2.0.0)
  - ggplot2 (3.4.3)
  - cowplot (1.1.3)
  - spatstat.geom (3.2-9)
  - sparr (2.3-10)
  - fdaPDE (1.1-19)
### Python (version at least 3.11.5)
Python libraries requirements:
  - numpy (1.26.3)
  - pandas (2.1.4)
  - pyproj (3.6.1)
  - geopandas (0.14.4)
  - geoplot (0.5.1)
  - shapely (2.0.4)
  - tqdm (4.65.0)
    
## Getting started

This project was coded with Python and R. To get started, first set your working directory at the root of this repository by using `cd [your path to this repository]` in bash terminal.

### Step 1: Data processing

The first step uses Python and involves downloading the dataset, processing it and saving it. With bash, excute `Python 'Python Code/main.py'`.

This script downloads the Storm Event Dataset from the NCEI FTP servers and processes it according to the pre-processing steps described in the dissertation.

### Step 2: Exploratory data analysis

This step is performed with Python. To replicate it, just run the notebook in `Notebooks/exploratory_data_analysis.ipynb`.

### Step 3: Statistical modelling

This part is coded in R.

#### Creating the mapping tables (optionnal)

You can re-create the mapping tables (they are already uploaded in the GitHub repository) that will be used in the R scripts by executing the file `R Code/mapping_table.R`. You can execute this file with the bash command `Rscript 'R Code/mapping_table.R'`.

#### Intensity estimation

To perform the intensity estimation step, you can execute the files `R Code/US_Hazard_STKDE.R` for hail or `R Code/US_Hazard_STDE-PDE.R` for hurricanes. These script will genereate intensity plots from the contiguous US.

#### Damages' distribution

To estimate the conditional distribution of the damages, execute the file `R Code/gamma_gpd_fit.R`. Don't forget to select the appropriate type of event you want to process (the variable `peril` can be set to either `"Hail"` or `"Hurricane"`). This script will generate the CDF plots, return levels plots and diagnostic plots displayed in the dissertation.

## Notes

- The file `R Code/Playground_STDE.R` is just an example of intensity estimation using the STDE-PDE method. It can be used as an example to play with to see the impact of the different parameters like the sample size, the smoothing coefficients or the shape of the true intensity function.
- Do not change the folder structure of the repository, as otherwise some paths to data are hardcoded in the scripts and would need to be changed if the folder are moved to another place. In particular, the folder `/Data` can be left like this and doesn't need any modification other than being downloaded from this repository.