library(data.table)
library(dplyr)

# 1. fread the Microclimate and Macroclimate data
micro_dt <- fread("../result/Model Input/Pitangus_sulphuratus_MICRO_Tmax_Tmin.csv")
macro_dt <- fread("../result/Model Input/points_climate_Pitangus_sulphuratus.csv")

# 2. Filter the microclimate data
micro_clean <- micro_dt[status == "ok" & !is.na(Micro_Bio5) & !is.na(Micro_Bio6)]

# Ensure the correct data format 
micro_clean[, `:=`(
  Micro_Bio5  = as.numeric(Micro_Bio5),
  Micro_Bio6  = as.numeric(Micro_Bio6),
  Micro_Bio13 = as.numeric(Micro_Bio13),
  Micro_Bio14 = as.numeric(Micro_Bio14),
  lon         = as.numeric(lon),
  lat         = as.numeric(lat)
)]

# 3. Calculate the Thermal and soil moisure limit (Global Min / Max). "limits" Extracted from Model_Input_Prepare...R
T_min_min <- limits$Tmin_min
T_max_max <- limits$Tmax_max
SM_min_min <- limits$SMmin_min
SM_max_max <- limits$SMmax_max

# 4. Calculate Micro Climatic Position (Micro_CP)
micro_clean[, `:=`(
  Micro_CP_Tmax = (Micro_Bio5  - T_min_min) / (T_max_max - T_min_min),
  Micro_CP_Tmin = (Micro_Bio6  - T_min_min) / (T_max_max - T_min_min),
  Micro_CP_SMmax = (Micro_Bio13  - SM_min_min) / (SM_max_max - SM_min_min),
  Micro_CP_SMmin = (Micro_Bio14  - SM_min_min) / (SM_max_max - SM_min_min)
)]

# fwrite(micro_clean, "Ara_macao_MICRO_Tmax_Tmin_with_CP.csv")
# cat("-> 带有 Micro CP 的微气候中间文件已保存为：Ara_macao_MICRO_Tmax_Tmin_with_CP.csv\n")

# ==============================================================================
# 5. Get the intercept part of Macro and Micro data
# ==============================================================================

# Ensure the coordinate precision consistency
macro_dt[, `:=`(
  lon = round(as.numeric(lon), 6),
  lat = round(as.numeric(lat), 6)
)]

micro_clean[, `:=`(
  lon = round(as.numeric(lon), 6),
  lat = round(as.numeric(lat), 6)
)]

macro_unique <- unique(macro_dt, by = c("lon", "lat"))
micro_unique <- unique(micro_clean, by = c("lon", "lat"))

# Select the needed cols
micro_sub <- micro_unique[, .(lon, lat, Micro_Bio5, Micro_Bio6, Micro_Bio13, Micro_Bio14,
                             Micro_CP_Tmax, Micro_CP_Tmin, Micro_CP_SMmax, Micro_CP_SMmin)]

final_combined <- inner_join(macro_unique, micro_sub, by = c("lon", "lat"))

# 6. Rename the colnames of macroclimate for consistency with microcliamte data 
setnames(final_combined, 
         old = c("Tmax", "Tmin", "SMmax", "SMmin", "CP_Tmax", "CP_Tmin", "CP_SMmax", "CP_SMmin"),
         new = c("Macro_Bio5", "Macro_Bio6", "Macro_Bio13", "Macro_Bio14", 
                 "Macro_CP_Tmax", "Macro_CP_Tmin", "Macro_CP_SMmax", "Macro_CP_SMmin"),
         skip_absent = TRUE)

if ("V1" %in% names(final_combined)) final_combined[, V1 := NULL]
if ("" %in% names(final_combined)) final_combined[[1]] <- NULL

# 7. Export the input variables for SDMs
output_final_path <- "../result/Model Input/final/Pitangus_sulphuratus_FULL_Macro_Micro_CP_merged.csv"
fwrite(final_combined, output_final_path)
