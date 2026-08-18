# PACKAGES
library(data.table)
library(dplyr)
library(lubridate)
library(parallel)
library(microclimf)
# library(elevatr)
library(terra)
library(ncdf4)
library(mcera5)
library(microclimdata)

# INPUT POINTS
pts <- fread("../data/points_climate_Pitangus_sulphuratus.csv")   # Species are replaceable here

stopifnot(all(c("lon","lat","pa") %in% names(pts)))
colnames(pts)[1] <- "id"

if (!("id" %in% names(pts))) {
  pts$id <- seq_len(nrow(pts))
}
 
pts <- pts %>%
  mutate(
    id = as.integer(id),
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    pa  = as.integer(pa)
  )

start_idx <- as.integer(Sys.getenv("START_IDX", "1"))
end_idx   <- as.integer(Sys.getenv("END_IDX", as.character(nrow(pts))))

end_idx <- min(end_idx, nrow(pts))

pts <- pts[start_idx:end_idx, ]

# FUNCTION: your single-point microclimate pipeline and its prerequisites 
# data Preparation for run_one_point_pipeline()
landcover <- rast("../data/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7.tif") 
soil_type_master <- rast("../data/soil_master/neotrop_soil_type_500m.tif")

target_crs <- "EPSG:32717"

from_vals <- c(
  10, 11, 12, 20, 30, 40,
  50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90,
  100, 110, 120, 121, 122, 130, 140,
  150, 152, 153, 160, 170, 180, 190,
  200, 201, 202, 210, 220
)

to_vals <- c(
  13, 13, 13, 13, 15, 15,
  2,  4,  4,  4,  1,  1,  1,  3,  3,  3,  5,
  8,  9,  6,  6,  7, 10, 16,
  16, 16, 16, 12, 12, 12, 14,
  16, 16, 16, NA, NA
)

credentials <- data.frame(
  username = "cienverdad",
  password = "Andynesssteam10!"
)

tme <- seq(
  as.POSIXct("2000-02-18 00:00:00", tz = "UTC"),
  as.POSIXct("2000-12-31 23:00:00", tz = "UTC"),
  by = "1 day"
)

mu <- NULL    
cnt <- NULL  
n_files_ok <- 0

out_path <- "../data/Microclimate"

clim_list <- vector("list", 12)

# prerequisites functions
process_one_lai <- function(hdf_file, template_r) {
  lai_path <- paste0(
    'HDF4_EOS:EOS_GRID:"', hdf_file, '":MOD_Grid_MOD15A2H:Lai_500m'
  )
  r <- rast(lai_path)
  # Clean the obvious deviation
  r[r > 10] <- NA
  # Project the template to the LAI coordinate system 
  template_lai <- project(template_r, crs(r))
  # Inspect the intercept
  if (is.null(intersect(ext(r), ext(template_lai)))) {
    warning("No overlap: ", basename(hdf_file))
    return(NULL)
  }
  r_crop <- crop(r, template_lai)
  r_proj <- project(r_crop, template_r, method = "bilinear")
  r_out <- resample(r_proj, template_r, method = "bilinear")
  return(r_out)
}

parse_modis_date <- function(filename) {
  x <- regmatches(basename(filename), regexpr("A[0-9]{7}", basename(filename)))
  year <- as.integer(substr(x, 2, 5))
  doy  <- as.integer(substr(x, 6, 8))
  as.Date(doy - 1, origin = paste0(year, "-01-01"))
}

read_wsa_one <- function(f, template = NULL) {
  p <- paste0('HDF4_EOS:EOS_GRID:"', f, '":MOD_Grid_BRDF:Albedo_WSA_shortwave')
  r <- try(rast(p), silent = TRUE)
  if (inherits(r, "try-error")) return(NULL)
  
  # Quality control
  r[r >= 32766] <- NA
  r[r < 0] <- NA
  r[r > 1] <- NA
  
  if (!is.null(template)) {
    r <- project(r, template, method = "bilinear")
    r <- crop(r, template)
    r <- resample(r, template, method = "bilinear")
  }
  
  # Skip when all NAs
  n_ok <- global(!is.na(r), "sum", na.rm = TRUE)[1,1]
  if (is.na(n_ok) || n_ok == 0) return(NULL)
  
  r
}

lonlat_to_modis_tile <- function(lon, lat) {
  R <- 6371007.181
  tile_size <- 1111950.5196666666
  xmin <- -20015109.354
  ymax <-  10007554.677
  
  lon_rad <- lon * pi / 180
  lat_rad <- lat * pi / 180
  
  x <- R * lon_rad * cos(lat_rad)
  y <- R * lat_rad
  
  h <- floor((x - xmin) / tile_size)
  v <- floor((ymax - y) / tile_size)
  
  h <- pmax(0, pmin(35, h))
  v <- pmax(0, pmin(17, v))
  
  sprintf("h%02dv%02d", h, v)
}

build_modis_index <- function(files, product_prefix = NULL) {
  bn <- basename(files)
  
  # e.g. MOD15A2H.A2000049.h05v10.061.xxxxx.hdf
  tile <- sub(paste0("^", product_prefix, "\\.A\\d{7}\\.(h\\d{2}v\\d{2})\\..*$"), "\\1", bn)
  doy  <- sub(paste0("^", product_prefix, "\\.A(\\d{7})\\..*$"), "\\1", bn)
  
  ok <- grepl("^h\\d{2}v\\d{2}$", tile) & grepl("^\\d{7}$", doy)
  
  data.frame(
    file = files[ok],
    tile = tile[ok],
    yday = doy[ok],
    stringsAsFactors = FALSE
  )
}

get_point_files_by_tile <- function(lon, lat, modis_index_df) {
  tile <- lonlat_to_modis_tile(lon, lat)
  out <- modis_index_df[modis_index_df$tile == tile, , drop = FALSE]
  out[order(out$yday), , drop = FALSE]
}

lai_dir <- "../data/Microclimate/MOD15A2H_061-20260520_184007"
lai_files <- list.files(lai_dir, pattern="\\.hdf$", full.names=TRUE)
lai_index <- build_modis_index(lai_files, product_prefix = "MOD15A2H")

gr_dir <- "../data/Microclimate/MCD43A3_061-20260607_163715"
gr_files <- list.files(gr_dir, pattern="\\.hdf$", full.names=TRUE)
gr_index <- build_modis_index(gr_files, product_prefix = "MCD43A3")

run_one_point_pipeline <- function(lon, lat, custom_temp) {
  # make a 2km *2km range centered on the presence point
  pt <- vect(
    data.frame(lon = lon, lat = lat),
    geom = c("lon", "lat"),
    crs = "EPSG:4326"
  )
  
  # Unify the coordinate system with the ESA CCI Land cover
  pt_proj <- project(pt, "EPSG:32717")
  xy <- crds(pt_proj)
  
  e <- ext(
    xy[1] - 1000, xy[1] + 1000,
    xy[2] - 1000, xy[2] + 1000
  )
  
  # Construct 2km * 2km ref SpatRaster
  r_base <- rast(e, resolution = 500, crs = "EPSG:32717")
  values(r_base) <- 1
  r_base_degree <- project(r_base, "EPSG:4326")
  
  ### runpointmodel input1：dtm_base
  dtm_raw <- microclimdata::dem_download(r_base, msk = FALSE, zeroasna = TRUE)
  dtm_res <- resample(dtm_raw, r_base, method = "bilinear")
  if (!compareGeom(dtm_res, r_base, stopOnError = FALSE)) {
    ext(dtm_res) <- ext(r_base)
  }
  dtm_base <- wrap(dtm_res)
  
  ### runpointmodel input2：Vegp_fixed
  landcover_proj <- project(landcover, r_base, method = "near")
  landcover_test <-  crop(landcover_proj, ext(r_base) + 2000)
  
  # filter or convert all non-land rasters
  vals_lc <- values(landcover_test)
  legal_lc_mask <- (vals_lc >= 10 & vals_lc <= 190) & !is.na(vals_lc)
  legal_vals <- vals_lc[legal_lc_mask]
  # Calculate the most advantaged land cover type in the local scale 
  if (length(legal_vals) > 0) {
    fallback_lc <- as.integer(names(sort(table(legal_vals), decreasing = TRUE))[1])
  } else {
    # If the full range falls into non-land type, then manually assign the classic shrubland type
    fallback_lc <- 130
  }
  invalid_lc_mask <- is.na(vals_lc) | vals_lc == 210 | vals_lc == 211 | vals_lc == 220 | vals_lc == 0 | vals_lc > 220
  
  if (any(invalid_lc_mask)) {
    vals_lc[invalid_lc_mask] <- fallback_lc
    values(landcover_test) <- vals_lc
  }
  
  
  landcover_test <- resample(landcover_test, r_base, method = "near")
  landcover_test <- mask(landcover_test, r_base)
  if (!compareGeom(landcover_test, r_base, stopOnError = FALSE)) {
    ext(landcover_test) <- ext(r_base)
  }
  habitat_cci <- subst(landcover_test, from = from_vals, to = to_vals)
  values(habitat_cci) <- as.integer(values(habitat_cci))
  habitat_base_wrapped <- wrap(habitat_cci)
  
  vegp_test <- vegpfromhab(
    habitats = habitat_base_wrapped,
    lat = lat,
    long = lon,
    tme =  as.POSIXlt( as.POSIXct(sprintf("2000-%02d-15 12:00:00", 1:12), tz = "UTC")),
    clump0 = TRUE
  )
  
  # Replace the PAI value with the LAI from MOD15A2H_061
  lai_point_df <- get_point_files_by_tile(lon, lat, lai_index)
  lai_hdf_files <- lai_point_df$file
  # Proceed the points in cluster across years
  lai_list <- lapply(lai_hdf_files, process_one_lai, template_r = r_base)
  lai_list <- lai_list[!sapply(lai_list, is.null)]
  lai_dates <- as.Date(sapply(lai_hdf_files, parse_modis_date))
  months_i <- lubridate::month(lai_dates) # data in Jan is unavailable 
  
  pai_monthly <- rast(lapply(2:12, function(m) {
    idx <- which(months_i == m)
    if (length(idx) == 0) {
      r <- r_base
      values(r) <- NA
      return(r)
    }
    rr <- rast(lai_list[idx])
    app(rr, mean, na.rm = TRUE)
  }))
  names(pai_monthly) <- paste0("lyr.", 2:12)
  
  # Supplement the Jan data
  pai_jan <- pai_monthly[[1]] 
  pai_monthly_12 <- c(pai_jan, pai_monthly) 
  names(pai_monthly_12) <- paste0("lyr.", 1:12) # 检查 
  
  # Any cell with the valid value will be considered in the global mean value 
  global_backup_val <- global(mean(pai_monthly_12, na.rm = TRUE), "mean", na.rm = TRUE)[1,1]
  # If all cells are NA, then manually assign 0.5 in prevention of error (but this is extremely rare)
  if (is.na(global_backup_val) || global_backup_val < 0.01) global_backup_val <- 0.5
  
  pai_list_fixed <- lapply(1:12, function(m) {
    lyr <- pai_monthly_12[[m]]
    na_count <- global(is.na(lyr), "sum")[1,1]
    if (na_count == ncell(lyr)) {
      values(lyr) <- global_backup_val
      message(sprintf("⚠️ 警告：检测到第 %d 月份为全 NA 空壳（可能2月全盲），已自动使用背景均值 %.3f 熨平防崩！", m, global_backup_val))
    }
    return(lyr)
  })
  
  pai_monthly_12 <- rast(pai_list_fixed)
  names(pai_monthly_12) <- paste0("lyr.", 1:12)
  
  # Supplement the clump data (the default is 0)
  clump_r <- rast(r_base, nlyrs = 12)
  values(clump_r) <- 0
  names(clump_r) <- paste0("lyr.", 1:12)
  
  vegp_fixed <- vegp_test
  vegp_fixed$pai <- wrap(pai_monthly_12)
  vegp_fixed$clump <- wrap(clump_r)
  # str(vegp_fixed)
  
  ### runpointmodel Input3: soilc_base
  # Construct soiltype
  soil_type_crop <- crop(soil_type_master, ext(r_base_degree) + 0.01)
  soil_type_proj <- project(soil_type_crop, target_crs, method = "near")
  
  soil_type_final <- resample(soil_type_proj, r_base, method = "near")
  if (!compareGeom(soil_type_final, r_base, stopOnError = FALSE)) {
    ext(soil_type_final) <- ext(r_base)
  }

  vals_soil <- values(soil_type_final)
  
  if (any(is.na(vals_soil)) || any(vals_soil < 1) || any(vals_soil > 11)) {
    legal_vals <- vals_soil[vals_soil >= 1 & vals_soil <= 11 & !is.na(vals_soil)]
    
    if (length(legal_vals) > 0) {
      # If found valid cells, fill the rest with the most frequent one
      fallback_soil <- as.integer(names(sort(table(legal_vals), decreasing = TRUE))[1])
    } else {
      # If the cells of the whole range fall into non-land type, manually assign 4 (Loam) to them
      fallback_soil <- 4
    }
    
    vals_soil[is.na(vals_soil) | vals_soil < 1 | vals_soil > 11] <- fallback_soil
    values(soil_type_final) <- vals_soil
  }
  
  names(soil_type_final) <- "soiltype"
  
  # Calculate the ground reflectance 
    gr_df <- get_point_files_by_tile(lon, lat, gr_index)
    gr_files <- gr_df$file

  for (i in seq_along(gr_files)) {
    f <- gr_files[i]
    # cat(sprintf("[%d/%d] %s\n", i, length(gr_files), basename(f)))

    r <- read_wsa_one(f, template = r_base)
    if (is.null(r)) next

    if (is.null(mu)) {
      mu  <- r
      cnt <- ifel(!is.na(r), 1, 0)
      n_files_ok <- 1
    } else {
      valid <- !is.na(r)
      mu  <- ifel(valid, (mu * cnt + r) / (cnt + 1), mu)
      cnt <- ifel(valid, cnt + 1, cnt)
      n_files_ok <- n_files_ok + 1
    }
    if (i %% 20 == 0) {
      mu  <- writeRaster(mu,  paste0(custom_temp, "/mu_tmp.tif"),  overwrite = TRUE)
      cnt <- writeRaster(cnt, paste0(custom_temp,"/cnt_tmp.tif"), overwrite = TRUE)
      gc()
    }
  }

  n_total <- ncell(mu)
  n_na <- global(is.na(mu), "sum", na.rm = TRUE)[1,1]

  if (n_na == 0) {
    groundr_final <- mu
    message("All cells valid: keep spatially varying groundr.")
  } else {
    # If there is missing value, fill the whole map with the non-NA values
    gmean <- global(mu, "mean", na.rm = TRUE)[1,1]
    groundr_final <- r_base
    values(groundr_final) <- gmean
    message(sprintf("NA cells detected (%d/%d): use uniform mean groundr = %.4f",
                    n_na, n_total, gmean))
  }

  print(global(groundr_final, c("min","mean","max"), na.rm=TRUE))
  
  # Combine the soilc_base                    
  soilc_base <- list(
    soiltype = wrap(soil_type_final),
    groundr  = wrap(groundr_final)
  )
  class(soilc_base) <- "soilcharac"
  
  ### runpointmodel Input 4: weather
  for (m in 1:12) {
    nc_file <- sprintf("%s/month%02d/era5_2000_%02d_v2.nc",
                       out_path, m, m)
    if (!file.exists(nc_file)) {
      stop(sprintf("文件不存在: %s", nc_file))
    }
    
    st_time <- as.POSIXct(sprintf("2000-%02d-01 00:00:00", m), tz = "UTC")
    en_time <- as.POSIXct(sprintf("2000-%02d-%02d 23:00:00",
                                  m, lubridate::days_in_month(as.Date(sprintf("2000-%02d-01", m)))),
                          tz = "UTC")
    
    clim_i <- mcera5::extract_clim(
      nc         = nc_file,
      long       = lon,
      lat        = lat,
      start_time = st_time,
      end_time   = en_time,
      format     = "microclimf"
    )
    
    # Keep all cols that microclimf need
    target_cols <- c("obs_time", "temp", "relhum", "pres", "swdown",
                     "difrad", "lwdown", "windspeed", "winddir", "precip")
    
    clim_i <- clim_i[, intersect(target_cols, names(clim_i))]
    clim_i <- clim_i[, target_cols]
    
    if ("swdown" %in% names(clim_i)) {
      clim_i$swdown <- pmax(0, clim_i$swdown)
    }
    if ("difrad" %in% names(clim_i)) {
      clim_i$difrad <- pmax(0, clim_i$difrad)
    }
    
    clim_list[[m]] <- clim_i
  }
  
  climdata_2000 <- bind_rows(clim_list) %>%
    arrange(obs_time)
  pure_climdata <- as.data.frame(climdata_2000)
  
  ### Run the single-point microclimate Bioclimatic variables 
  micropoint <- runbioclim(
    climdata = pure_climdata,
    reqhgt   = 0.05,
    dtm      = dtm_base,
    vegp     = vegp_fixed,
    soilc    = soilc_base
  )
  extracted_vals <- terra::extract(micropoint, pt_proj, ID = FALSE)
  
  return(list(
    Bio5  = as.numeric(extracted_vals$bio5),
    Bio6  = as.numeric(extracted_vals$bio6),
    Bio13 = as.numeric(extracted_vals$bio13),
    Bio14 = as.numeric(extracted_vals$bio14)
  ))
}

run_one_point <- function(id, lon, lat, pa, custom_temp) {
  out <- tryCatch({
    microp <- run_one_point_pipeline(lon = lon, lat = lat, custom_temp = custom_temp)
    
    data.table(
      id = id, lon = lon, lat = lat, pa = pa,
      Micro_Bio5  = microp$Bio5,   # Max Temperature
      Micro_Bio6  = microp$Bio6,   # Min Temperature
      Micro_Bio13 = microp$Bio13,  # Wettest Soil Moisture
      Micro_Bio14 = microp$Bio14,  # Driest Soil Moisture
      status = "ok",
      error_msg = ""
    )
  }, error = function(e) {
    cat(sprintf("ERROR at point id=%d (lon=%.4f, lat=%.4f): %s\n", 
                id, lon, lat, e$message), file = stderr())
    print(e)
    cat(sprintf("\n[DEBUG] Error at ID %d. checking geometry:\n", id), file = stderr())
    cat(sprintf("Message: %s\n", e$message), file = stderr())
    
    data.table(
      id = id, lon = lon, lat = lat, pa = pa,
      Micro_Bio5  = NA_real_,
      Micro_Bio6  = NA_real_,
      Micro_Bio13 = NA_real_,
      Micro_Bio14 = NA_real_,
      status = "fail",
      error_msg = as.character(e$message) 
    )
  })
  return(out)
}

# CHUNKING (200 points per task)
chunk_size <- 200L
pts$chunk_id <- ((seq_len(nrow(pts)) - 1L) %/% chunk_size) + 1L
chunks <- split(pts, pts$chunk_id)

run_chunk <- function(df_chunk, out_file) {
  current_chunk_id <- df_chunk$chunk_id[1]
  my_temp_dir <- sprintf(
    "terra_tmp_batch_%06d_%06d_chunk_%06d",
    start_idx,
    end_idx,
    current_chunk_id
  )
  dir.create(my_temp_dir, showWarnings = FALSE, recursive = TRUE)
  terraOptions(memfrac = 0.4, tempdir = my_temp_dir, progress = 0)
  
  total_pts <- nrow(df_chunk)
  
  for (i in seq_len(total_pts)) {
    row <- df_chunk[i, ]
    
    cat(sprintf("[%s] Chunk %06d: Starting point %d/%d (ID: %d) ...\n", 
                format(Sys.time(), "%H:%M:%S"), current_chunk_id, i, total_pts, row$id))
    
    # Run the runbioclim model for the current single point
    res_point <- run_one_point(
      id = row$id, lon = row$lon, lat = row$lat, pa = row$pa,
      custom_temp = my_temp_dir
    )
    
    # Examine and fwrite
    if (!file.exists(out_file)) {
      # If the current point is the first of the chunk, then create a file and fwrite in the 1st row
      fwrite(res_point, out_file, append = FALSE, col.names = TRUE)
    } else {
      # If not, then append the data to the existed files
      fwrite(res_point, out_file, append = TRUE, col.names = FALSE)
    }
    
    flush.console()
  }
  
  unlink(my_temp_dir, recursive = TRUE)
  return(out_file)
}

# PARELLEL RUN ()
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  k <- 1                    # For the local test
} else {
  k <- as.integer(args[1])  # PBS INdex = chunk id
}

cat("=== Running Chunk", k, "===\n")
out_file <- sprintf("../result/chunk_outputs/batch_%06d_%06d_chunk_%06d.csv",start_idx,
                    end_idx, k)

if (!file.exists(out_file)) {
  run_chunk(chunks[[k]], out_file)
}

message("Done! Output saved to: ", out_file)

# # MERGE ALL POINT RESULTS (Excecute the section after points of the species are all finished)
# files <- list.files("../result/chunk_outputs", pattern = "^batch_\\d+_\\d+_chunk_\\d+\\.csv$", full.names = TRUE)
# 
# all_res <- rbindlist(lapply(files, fread), fill = TRUE)
# all_res_clean <- all_res[!is.na(Micro_Bio5)]
# all_res_clean <- unique(all_res_clean, by = "id")
# all_res_clean <- all_res_clean[order(id)]
# 
# # # Save Tmax/Tmin of the species
# fwrite(all_res_clean, "../result/chunk_outputs/Pitangus_sulphuratus_Final.csv")
# qsub -v START_IDX=3001,END_IDX=6000 run_chunk.pbs
