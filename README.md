# COCA
This contains sample code for the project "Climate Induced Habitat Changes of Fish Stocks" published in ICES Journal of Marine Science October 2022.



Code Descriptions:

	HSI_Loop.m - Calculates Habitat Suitability Index and Availability for species ensemble

	ROMS_Downscaling.m - Downscales ROMS models to same resolution as the benthic data. Also regrids habitat models and ROMS temperature biases

	habitat_model_cv.R - Creates, cross validates, and projects habitat models for every combination of available variables	

	Benthic_raster.R - Analyzes GeoTiffs and Shapefiles of benthic features to transform them into a model-usable dataframe

	Trawl_data_clean.R - Cleans and combines NOAA/State Bottom Trawl Survey data. Also supplements variables required for swept area calculation

	trawl_roms_comparison.R - Calculates ROMS bottom-water temperature bias against NEFSC observations


Data samples: 

	alltrawl_sample.csv - Sample of NEFSC and state bottom trawl survey and all associated variables required to calculate swept area/night correction and create/project habitat models. Any column names in all caps conforms to standard NEFSC naming conventions.

	bwt_20100421.mat - Downscaled .mat file of ROMS ocean model for April 21 2010. 

	habvar_sample.csv - sample of benthic classification values used to predict completed habitat models  
