library(data.table)
library(dplyr)
library(purrr)
library(tidyr)

# 1. Read input data of 10 species
files <- list.files(
  "D:/Paper9/result/Model Input/final/",
  pattern = "_FULL_Macro_Micro_CP_merged\\.csv$",
  full.names = TRUE
)

all_species <- lapply(files, function(f) {
  
  dat <- fread(f)
  
  dat$species <- basename(f) |>
    sub("_FULL_Macro_Micro_CP_merged\\.csv$", "", x = _)
  
  dat
}) |>
  rbindlist(fill = TRUE)

# 2. Define the four environmental variables of predictor sets
var_sets <- list(
  
  Macro = c(
    Tmax  = "Macro_Bio5",
    Tmin  = "Macro_Bio6",
    SMmax = "Macro_Bio13",
    SMmin = "Macro_Bio14"
  ),
  
  MacroCP = c(
    Tmax  = "Macro_CP_Tmax",
    Tmin  = "Macro_CP_Tmin",
    SMmax = "Macro_CP_SMmax",
    SMmin = "Macro_CP_SMmin"
  ),
  
  Micro = c(
    Tmax  = "Micro_Bio5",
    Tmin  = "Micro_Bio6",
    SMmax = "Micro_Bio13",
    SMmin = "Micro_Bio14"
  ),
  
  MicroCP = c(
    Tmax  = "Micro_CP_Tmax",
    Tmin  = "Micro_CP_Tmin",
    SMmax = "Micro_CP_SMmax",
    SMmin = "Micro_CP_SMmin"
  )
)

# 3. Generate 6 pairs among the 4 sets 
set_pairs <- combn(
  names(var_sets),
  2,
  simplify = FALSE
)

# 4. Calculate 24 correlations for each species 
cor_results <- map_dfr(
  unique(all_species$species),
  
  function(sp) {
    
    dat_sp <- all_species %>%
      filter(species == sp)
    
    map_dfr(set_pairs, function(pair) {
      
      set1 <- pair[1]
      set2 <- pair[2]
      
      map_dfr(names(var_sets[[set1]]), function(env_var) {
        
        xvar <- var_sets[[set1]][env_var]
        yvar <- var_sets[[set2]][env_var]
        
        tmp <- dat_sp %>%
          dplyr::select(
            x = all_of(xvar),
            y = all_of(yvar)
          ) %>%
          drop_na()
        
        test <- cor.test(
          tmp$x,
          tmp$y,
          method = "pearson"
        )
        
        tibble(
          species = sp,
          predictor_set_1 = set1,
          predictor_set_2 = set2,
          environmental_variable = env_var,
          n = nrow(tmp),
          r = unname(test$estimate),
          p_value = test$p.value
        )
      })
    })
  }
)

cor_summary <- cor_results %>%
  group_by(
    predictor_set_1,
    predictor_set_2,
    environmental_variable
  ) %>%
  summarise(
    mean_r   = mean(r, na.rm = TRUE),
    median_r = median(r, na.rm = TRUE),
    min_r    = min(r, na.rm = TRUE),
    max_r    = max(r, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(cor_summary,"../result/correlation/cor_summary.csv")
