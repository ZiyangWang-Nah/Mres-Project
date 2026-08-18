library(rgbif)
library(terra)
library(predicts)
library(sf)
library(dplyr)
library(tidyr)
library(dismo)
library(pROC)
library(randomForest)


### Download GBIF species occurrence records in Year 1981-2010 in whole neotropical region
# 1. Get the Taxon Key for Birds (Names replaceable here)
taxon_key <- name_backbone(name='Pitangus sulphuratus')$speciesKey

# 2. Define your Polygon (WKT format)
my_polygon <- 'POLYGON((-123.3374 -60.34388, -21.58016 -60.34388, -21.58016 31.11561, -73.50268 31.11561, -123.3374 31.11561, -123.3374 -60.34388))'

# 3. Trigger the download request
download_request <- occ_download(
  pred("taxonKey", taxon_key),
  pred_gte("year", 1981),         
  pred_lte("year", 2010),        
  pred_within(my_polygon),
  pred_in("basisOfRecord", c("MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "LIVING_SPECIMEN")),
  format = "SIMPLE_CSV", 
  user = "cienverdad", 
  pwd = "Andyness10", 
  email = "ziyangw480@gmail.com"
)

# 4. Check the status
occ_download_wait(download_request)

# 5. Download and Import the data
species_1981_2010 <- occ_download_get(download_request) %>%
  occ_download_import()

write.csv(species_1981_2010,"../data/occurrence/Pitangus sulphuratus.csv")
species_1981_2010 <- read.csv("../data/occurrence/Pitangus sulphuratus.csv",row.names = NULL)

presence <- species_1981_2010[, c("decimalLongitude", "decimalLatitude")]
presence_sf <- st_as_sf(presence, 
                        coords = c("decimalLongitude", "decimalLatitude"), 
                        crs = 4326) # For absence point creation by backgroundPOINTS

### Get Macroclimate data and MacroCP, along with tolerance limits
# 0. Prepare the bioclimatic (including climatic position) and land cover layer. Scale them in the same spatial resolution. Scale the presence pts 
# SETUP EXTENT
neotrop_ext <- ext(-123.3374, -21.58016, -60.34388, 31.11561)

# PREPARE CLIMATE LAYERS (Tmin, Tmax, Pmin, Pmax) 
temp_files <- c(
  "../data/Macroclimate/chelsa02/chelsa/global/bioclim/bio05/1981-2010/CHELSA_bio05_1981-2010_V.2.1.tif",
  "../data/Macroclimate/chelsa02/chelsa/global/bioclim/bio06/1981-2010/CHELSA_bio06_1981-2010_V.2.1.tif"
)

# Extract the 5 bioclimatic vars within the neotropical region
temp_stack <- rast(temp_files) %>% crop(neotrop_ext)
names(temp_stack) <- c("Tmax", "Tmin")

# 1. Prepare the soil moisture files
nc_file <- "../data/Macroclimate/ERA5/data_stream-moda.nc"
sds(nc_file)

# Check the data structure and the time 
sm_raw <- rast(nc_file)
layer_names <- names(sm_raw)
# Extract seconds prior to 'valid_time=' 
time_seconds <- as.numeric(
  sub(".*valid_time=", "", layer_names)
)

# Convert to POSIXct
dates <- as.POSIXct(
  time_seconds,
  origin = "1970-01-01",
  tz = "UTC"
)
# Extract year and month
years  <- as.integer(format(dates, "%Y"))
months <- as.integer(format(dates, "%m"))

# Calculate the mean of 12 months over the 30 years
sm_monthly_clim <- rast(
  lapply(1:12, function(m) {
    
    month_idx <- which(
      months == m
    )
    
    mean(
      sm_raw[[month_idx]],
      na.rm = TRUE
    )
  })
)
names(sm_monthly_clim) <- month.abb

SMmax <- app(
  sm_monthly_clim,
  max,
  na.rm = TRUE
)

SMmin <- app(
  sm_monthly_clim,
  min,
  na.rm = TRUE
)

SMmax_1km <- resample(
  SMmax,
  temp_stack[[1]],
  method = "bilinear"
)

SMmin_1km <- resample(
  SMmin,
  temp_stack[[1]],
  method = "bilinear"
)

names(SMmax_1km) <- "SMmax"
names(SMmin_1km) <- "SMmin"

# Construct the new bioclim stack
bioclim_stack <- c(
  temp_stack,
  SMmax_1km,
  SMmin_1km
)
writeRaster(bioclim_stack, "../result/Model Input/bioclim_stack_1981-2010.tif", overwrite = TRUE)
bioclim_stack <- rast("../result/Model Input/bioclim_stack_1981-2010.tif")

# 2. RESAMPLE LAND COVER (Modal)
landcover <- rast("../data/Land Cover/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7.tif") %>% crop(neotrop_ext)
unique_values <- unique(landcover)
print(unique_values)

# Reclassify the original 30 classes into 5 big classes for smaller df in the model
reclass_matrix_fixed <- matrix(c(
  # 1). Forest
  50,1, 60,1, 61,1, 62,1, 70,1, 80,1, 90,1, 100,1, 160,1, 170,1,
  
  # 2). Semi-natural 
  40,2, 110,2, 120,2, 121,2, 122,2, 130,2, 150,2, 152,2, 153,2, 180,2,
  
  # 3). Cropland 
  10,3, 11,3, 12,3, 20,3, 30,3,
  
  # 4). Urban 
  190,4,
  
  # 5). Barren/Ice 
  200,5, 201,5, 220,5,
  
  # Water body is NA
  210, NA
), ncol=2, byrow=TRUE)
lc_grouped <- classify(landcover,reclass_matrix_fixed)

# Use "modal" for categorical data when downscaling
landcover_1km <- resample(lc_grouped, bioclim_stack[[1]], method = "mode")
names(landcover_1km) <- "landcover"
writeRaster(landcover_1km, "../result/landcover_1km_1981-2010.tif", overwrite = TRUE)
landcover_1km <- rast("../result/Model Input/landcover_1km_1981-2010.tif")

# 3. Reduce the presence records into '1pt in one 1km * 1km cell'
cells <- cellFromXY(bioclim_stack[[1]], st_coordinates(presence_sf))
unique_cells <- unique(cells)
presence_thinned <- presence_sf[which(!duplicated(cells)), ]

presence_thinned_pts <- presence_thinned %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  rename(lon = X, lat = Y)

bioclim_stack_slim <- mask(
  bioclim_stack,
  landcover_1km
)
writeRaster(
  bioclim_stack_slim,
  "../result/Model Input/bioclim_stack_slim_1981-2010.tif",
  overwrite = TRUE
)
bioclim_stack_slim <- rast("../result/Model Input/bioclim_stack_slim_1981-2010.tif")

# 4.1 Convert the presenced_thinned_points into sf objects 
pts_sf_range <- st_as_sf(presence_thinned_pts, coords = c("lon", "lat"), crs = 4326)

# 4.2 Generate the geographical MCP via thinned presence records to represent the species range
species_range_polygon <- st_convex_hull(st_union(pts_sf_range))

# 4.3 Use MCP-polygon to crop and mask the environmental predictors
bioclim_in_species_range <- crop(bioclim_stack_slim, species_range_polygon)
bioclim_in_species_range <- mask(bioclim_in_species_range, vect(species_range_polygon))

# 4.4 Vectorize all values of environmental predictors across the species range
all_pixels_in_range <- as.data.frame(bioclim_in_species_range, na.rm = TRUE)

# 4.5 Extract the maximum and minimum value
limits <- all_pixels_in_range %>%
  summarise(across(everything(), list(
    min = ~quantile(.x, 0),   
    max = ~quantile(.x, 1)    
  )))
rm(all_pixels_in_range)
gc()

# 5. CALCULATE CLIMATIC POSITION (CP)
# Calculate the denominators (Range of tolerance)
thermal_range <- limits$Tmax_max - limits$Tmin_min
soilmoist_range <-
  limits$SMmax_max -
  limits$SMmin_min

cp_tmin <- (bioclim_stack_slim[["Tmin"]] - limits$Tmin_min) / thermal_range
cp_tmax <- (bioclim_stack_slim[["Tmax"]] - limits$Tmin_min) / thermal_range
cp_smmin <- (bioclim_stack_slim[["SMmin"]] - limits$SMmin_min) / soilmoist_range
cp_smmax <- (bioclim_stack_slim[["SMmax"]] - limits$SMmin_min) / soilmoist_range

# Combine into a stack and name them
cp_stack_slim <- c(cp_tmin, cp_tmax, cp_smmin, cp_smmax)
names(cp_stack_slim) <- c("CP_Tmin", "CP_Tmax", "CP_SMmin", "CP_SMmax")
writeRaster(cp_stack_slim, "../result/Model Input/cp_stack_slim_1981-2010.tif", overwrite = TRUE)
cp_stack_slim <- rast("../result/Model Input/cp_stack_slim_1981-2010.tif")

# Combine into final environmental stack, and remove areas in the water body for a save of RAM
final_stack_slim <- c(cp_stack_slim, landcover_1km)
names(final_stack_slim) <- c("CP_Tmin", "CP_Tmax", "CP_SMmin", "CP_SMmax", 
                             "landcover")
writeRaster(final_stack_slim, "../result/Model Input/final_stack_slim_1981-2010.tif", overwrite = TRUE)
final_stack_slim <- rast("../result/Model Input/final_stack_slim_1981-2010.tif")

### Create the background points as pueso-absent records, using random points generated from Automaxent() function. Then, combine the presence and absence records
# Load the stable versios for the specific functions
lapply(installed.packages(lib.loc="./Libraries/AutoMaxent")[,"Package"],require,lib.loc="./Libraries/AutoMaxent",character.only=TRUE)

# Load the BackgroundPOINTS function
MaxEnt_Functions <- "../data/Functions/AutoMaxent" %>% list.files(pattern="BackgroundPOINTS.R",full.names = T)
lapply(MaxEnt_Functions, function(s) source(s))

# Create the pseudo-absence data
study_area_mask <- final_stack_slim[[1]] %>% 
  as.polygons() %>% 
  st_as_sf() %>% 
  st_union()

set.seed(0309)
r_bk <- backgroundPOINTS(presence=presence_thinned, # Distribution points of the species
                         background_n = nrow(presence_thinned_pts), # Number of background points
                         range_samp = study_area_mask, # Range information of species
                         TrainTest = 1,
                         buffer.dist = 0.5,
                         weights.p="Random")
bg_pts <- r_bk$Train %>% 
  st_coordinates() %>%      
  as.data.frame() %>%       
  rename(lon = X, lat = Y)  

### Combine the presence and absence records, and corresponding bioclimatic vars and land cover data into one df (for running GLM)
presence_thinned_pts$pa <- 1
bg_pts$pa <- 0
all_points <- rbind(presence_thinned_pts,bg_pts)

# Extract env variables out of raster corresponding to the all_points but with coordinate and absence-presence missing
env_extract <- terra::extract(x = final_stack_slim, 
                              y = all_points[, c("lon", "lat")],
                              ID = FALSE) 

# Extract env variables out of raster corresponding to the all_points but with coordinate and absence-presence missing
env_extract_raw <- terra::extract(x = bioclim_stack_slim, 
                            y = all_points[, c("lon", "lat")],
                            ID = FALSE)

# Add the point and climate datasets together
points_climate <- cbind(all_points, env_extract,env_extract_raw) 
points_climate <- points_climate %>% 
  drop_na() 
points_climate$landcover <- as.factor(points_climate$landcover)

points_climate$pa_factor <- as.factor(
  ifelse(points_climate$pa == 1, "presence", "background")
)
write.csv(points_climate,"../result/Model Input/points_climate_Pitangus_sulphuratus.csv")
points_climate <- read.csv("Occurrence/points_climate.csv", row.names = NULL)
