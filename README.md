# Mres-Project
Data and code archives for Mres Project of Mres programme: Living Planet with Computationel Methods in Ecology and Evolution  
_code_ of the project are listed in their corresponding folders<br> 
_data_ of the project can be downloaded from the following link:<br>
data/Land Cover: https://maps.elie.ucl.ac.be/CCI/viewer/download.php<br>
data/Macroclimate/chelsa02: https://envicloud.wsl.ch/#/?bucket=https%3A%2F%2Fos.unil.cloud.switch.ch%2Fchelsa02%2F&prefix=chelsa%2Fglobal%2Fbioclim%2F<br>
data/Macroclimate/ERA5: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview<br>
data/Microclimate(weather data): https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview<br>
data/Microclimate/MOD15A2H_061-20260520_184007: https://www.earthdata.nasa.gov/data/catalog/lpcloud-mod15a2h-006<br>
data/Microclimate/MCD43A3_061-20260607_163715: https://www.earthdata.nasa.gov/data/catalog/lpcloud-mcd43a3-006<br>
data/occurrence: https://www.gbif.org/occurrence/search?occurrenceStatus=present

# Data collection, analysis and modelling (Workflow)
1. Rscript Code/_Model_Input_Prepare_except_Micro.R_
   -> Acquire environmental predictors including Macroclimate (Macro_Bio05, Macro_Bio06, Macro_SMmax, Macro_SMmin), Macroclimate-based Climatic Position (Macro_CP) and Land cover
   -> Acquire thinned presence and absence records of one of the focal species. Ready for extracting Microclimate and Microclimate-based CP
2. Rscript HPC/_MicroBioclim_Final.R_
   -> Collect dtm, vegetation, soil information and local weather data, parallelly run on HPC (bash: qsub _run_chunk.pbs_). These input variables can be used to run the _runbioclim()_ model, and then microclimate data of each species is acquired (Micro_Bio05, Micro_Bio06, Micro_SMmax, Micro_SMmin)
3. Rscript Code/_microcp_finalmerge.R_
   -> Integrate the Micro and Macro datasets for the complete input variables for the Species Distribution Modelling
4. Rscript Code/_Spatial_Correlation.R_
   -> Prepare the training and test indices (5 folds) of each species for SDM fit via 3 different algorithms (Generalized Linear Model, Maxent and Random Forest)
5. Rscript Code/_ThreeModels.R_
   -> Fit SDMs with (GLM, Maxent and RF) * (Macroclimate, Microclimate and Land cover variables), and store the AUC score of each combination and species
6. Rscript Code/_Correlation.R_
   -> Calculate the Pearson corrlations among environmental variables across 10 species, to check the explained variation

# Result Visualization
1. Rscript Code/_MapVisualization.R_
   -> Visualize the Fig. 1-4
