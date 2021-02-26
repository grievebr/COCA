# COCA
This contains the code for the project "Climate Induced Habitat Changes of Fish Stocks"
It will contain a repository and link to view in a Jupyter notebook.
Project is near completion, but code will be added within the next week or so

Code Descriptions:
HSI_Loop.m - Calculates Habitat Suitability Index and Availability for species ensemble
ROMS_Downscaling.m - Downscales ROMS models to same resolution as the benthic data. Also regrids habitat models and ROMS temperature biases
habitat_model_cv.R - Creates, cross validates, and projects habitat models for every combination of available variables
Benthic_raster.R - Analyzes GeoTiffs and Shapefiles of benthic features to transform them into a model-usable dataframe
Trawl_data_clean.R - Cleans and combines NOAA/State Bottom Trawl Survey data. Also supplements variables required for swept area calculation
trawl_roms_comparison.R - Calculates ROMS bottom-water temperature bias against NEFSC observations

Data samples: 
alltrawl - NEFSC and state bottom trawl survey and all associated variables required to calculate swept area and habitat models. Only a sample because I do not have permission to share, especially the state data
NWA-SZ.HCob05T_avg_2010-04-21.nc - 7km ROMS cobalt bottom-water temperature (C) for April 21 2010. 
NWA_grid.nc - Gridfile necessary to interpret curvelinear ROMS model. 
habvar_sample.csv - sample of benthic classification values used to predict completed habitat models  
