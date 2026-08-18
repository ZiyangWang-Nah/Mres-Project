# Conduct the spatial autocorrelation after merging macro micro macrocp and microcp.
install.packages("blockCV")
library(blockCV)
library(sf)

point_finale = fread("D:/Paper9/result/Model Input/final/Pitangus_sulphuratus_FULL_Macro_Micro_CP_merged.csv")
pts_sf <- st_as_sf(point_finale, 
                   coords = c("lon", "lat"), 
                   crs = 4326)

# Calculate the range of spatial autocorrelation
sac <- cv_spatial_autocor(
  x = pts_sf,
  column = "pa",       
  plot = TRUE
)
cat("空间自相关范围：", sac$range, "米\n")

spatial_folds <- cv_spatial(
  x        = pts_sf,
  column   = "pa",
  k        = 5,
  size     = sac$range,  
  seed     = 42
)

# Use these folds instead of a random 80/20 split to train and evaluate all models
# Extract the fold index, with [[1]] as the training dataset while [[2]] as the test dataset
folds_list <- spatial_folds$folds_list

saveRDS(folds_list,"foldslist_Pitangus_sulphuratus.rds")