library(sf)
library(ggplot2)
library(dplyr)
library(maptiles)
library(tidyterra)
library(ggspatial)
library(shadowtext)

library(purrr)
library(lme4)
library(lmerTest)
library(emmeans)

library(data.table)
library(randomForest)
library(viridis)

neotropics <- st_read(
  "../result/Visualization/Fig1/NeotropicMap_Geo.shp"
)

neotropics <- st_transform(
  neotropics,
  crs = 4326
)

################# 绘制Fig1
species_order <- c(
  "Ara macao",
  "Tapera naevia",
  "Ardea cocoi",
  "Aramus guarauna",
  "Rupornis magnirostris",
  "Zonotrichia capensis",
  "Furnarius rufus",
  "Athene cunicularia",
  "Vanellus chilensis",
  "Pitangus sulphuratus"
)

species_first6 <- species_order[1:6]
species_last4  <- species_order[7:10]

points <- fread(
  "../result/Model Input/final/Ara_macao_FULL_Macro_Micro_CP_merged.csv"
)

Ara_presence <- points %>%
  filter(pa == 1) %>%
  dplyr::select(lon, lat) %>%
  mutate(
    species = "Ara macao",
    habitat = "Forest"
  )

all_presence <- bind_rows(
  Ara_presence,
  Aramus_presence,
  Ardea_presence,
  tapera_presence,
  Rupornis_presence,
  Zonotrichia_presence,
  Furnarius_presence,
  Athene_presence,
  Vanellus_presence,
  Pitangus_presence
)

presence_first6 <- all_presence %>% 
  filter(species %in% species_first6)

presence_first6$species <- factor(
  presence_first6$species,
  levels = species_first6
)

presence_last4 <- all_presence %>% 
  filter(species %in% species_last4)

presence_last4$species <- factor(
  presence_last4$species,
  levels = species_last4
)

habitat_colors <- c(
  "Forest"         = "#009E73",
  "Shrubland"      = "#D55E00",
  "Grassland"      = "#E69F00",
  "Wetland"        = "#0072B2",
  "Human Modified" = "#CC79A7"
)

p_species <- ggplot() +
  
  geom_sf(
    data = neotropics,
    fill = "grey92",
    color = "grey35",
    linewidth = 0.25
  ) +
  
  geom_point(
    data = presence_last4,
    aes(
      x = lon,
      y = lat,
      color = habitat
    ),
    size = 0.6,
    alpha = 0.65
  ) +
  
  scale_color_manual(
    name = "Major Habitat Category",
    values = habitat_colors,
    breaks = c(
      #"Forest",
      "Shrubland",
      "Grassland",
      #"Wetland",
      "Human Modified"
    )
  ) +
  
  coord_sf(
    xlim = c(-125, -20),
    ylim = c(-61, 33),
    expand = FALSE
  ) +
  
  facet_wrap(
    ~ species,
    nrow = 2
  ) +
  
  theme_void() +
  
  theme(
    strip.text = element_text(
      face = "italic",
      size = 14
    ),
    
    panel.border = element_rect(
      color = "grey30",
      fill = NA,
      linewidth = 0.5
    ),
    
    panel.spacing = unit(
      4,
      "mm")
  #   ,
  #   
  #   legend.position = "right",
  #   
  #   legend.title = element_text(
  #     size = 13,
  #     face = "bold"
  #   ),
  #   
  #   legend.text = element_text(
  #     size = 12
  #   )
  # ) +
  # guides(
  #   color = guide_legend(
  #     override.aes = list(
  #       size = 4,
  #       alpha = 1
  #     )
  #   )
  )
p_species
ggsave(
  "../result/Visualization/Fig1/Figure1_last4.png",
  p_species,
  width = 8,
  height = 11,
  dpi = 600
)

################# Fig2
raw_results_Pitangus_sulphuratus <- read.csv("../result/Model Output/Pitangus_sulphuratus_RawResults.csv")

# Convert to long format
results_long_Pitangus_sulphuratus <- raw_results_Pitangus_sulphuratus %>%
  dplyr::select(-1) %>%   
  pivot_longer(
    cols = starts_with("auc_"),
    names_to = c("algorithm", "config"),
    names_pattern = "auc_([^_]+)_(.*)",
    values_to = "AUC"
  ) %>%
  mutate(
    species = "Pitangus sulphuratus",
    
    algorithm = recode(
      algorithm,
      "glm" = "glm",
      "maxent" = "maxent",
      "rf" = "random forest"
    )
  ) %>%
  dplyr::select(species, fold, config, algorithm, AUC)

results_long_Pitangus_sulphuratus <- results_long_Pitangus_sulphuratus%>%
  mutate(
    species   = factor(species),
    fold      = factor(fold),
    algorithm = factor(algorithm),
    config    = as.character(config)
  )

calculate_delta <- function(data, pair_table) {
  
  purrr::pmap_dfr(
    pair_table,
    function(Comparison, model_A, model_B) {
      
      tmp <- data %>%
        filter(config %in% c(model_A, model_B)) %>%
        select(
          species,
          fold,
          algorithm,
          config,
          AUC
        ) %>%
        pivot_wider(
          names_from = config,
          values_from = AUC
        )
      if (!all(c(model_A, model_B) %in% names(tmp))) {
        return(NULL)
      }
      
      tmp %>%
        transmute(
          species,
          fold,
          algorithm,
          Comparison = Comparison,
          Delta_AUC =
            .data[[model_A]] -
            .data[[model_B]]
        ) %>%
        filter(!is.na(Delta_AUC))
    }
  )
}

# Colors of algorithms kept consistent across the 3 panels
algorithm_colors <- c(
  "glm"           = "#D55E00",
  "maxent"        = "#009E73",
  "random forest" = "#0072B2"
)

# Each species = different shapes of point
species_shapes <- c(
  "Ara macao"       = 16,
  "Ardea cocoi"     = 17,
  "Aramus guarauna" = 15,
  "Tapera naevia"   = 18,
  "Rupornis magnirostris" = 10,
  "Zonotrichia capensis" = 11,
  "Furnarius rufus" = 12,
  "Athene cunicularia" = 13,
  "Vanellus chilensis" = 14,
  "Pitangus sulphuratus" = 19
)

# Species names shown in the figure legend 
species_labels <- c(
  "Ara macao"       = "Ara macao",
  "Ardea cocoi"     = "Ardea cocoi",
  "Aramus guarauna" = "Aramus guarauna",
  "Tapera naevia"   = "Tapera naevia",
  "Rupornis magnirostris" = "Rupornis magnirostris",
  "Zonotrichia capensis" = "Zonotrichia capensis",
  "Furnarius rufus" = "Furnarius rufus",
  "Athene cunicularia" = "Athene cunicularia",
  "Vanellus chilensis" = "Vanellus chilensis",
  "Pitangus sulphuratus" = "Pitangus sulphuratus"
)

# FIg2A
pairs_A <- tibble(
  Comparison = c(
    "Macro",
    "Macro CP",
    "Micro",
    "Micro CP"
  ),
  
  model_A = c(
    "Lc_add_Macro",
    "Lc_add_MacroCp",
    "Lc_add_Micro",
    "Lc_add_MicroCp"
  ),
  
  model_B = c(
    "Macro",
    "MacroCp",
    "Micro",
    "MicroCp"
  )
)
panelA_Pitangus_sulphuratus <- calculate_delta(
  results_long_Pitangus_sulphuratus,
  pairs_A
)
panelA_final <- rbind(panelA_Ara_macao,panelA_Aramus_guarauna,panelA_Ardea_cocoi,panelA_Tapera_naevia,panelA_Rupornis_magnirostris,panelA_Zonotrichia_capensis,panelA_Furnarius_rufus,panelA_Athene_cunicularia, panelA_Vanellus_chilensis,panelA_Pitangus_sulphuratus)

model_A <- lmer(
  Delta_AUC ~
    Comparison * algorithm +
    (1 | species) +
    (1 | species:fold),
  data = panelA_final
)
anova(model_A)
summary(model_A)

emm_A <- emmeans(
  model_A,
  ~ Comparison | algorithm
)
emm_A_df <- as.data.frame(
  confint(emm_A)
)

raw_A <- panelA_final %>%
  group_by(
    species,
    algorithm,
    Comparison
  ) %>%
  summarise(
    Delta_AUC =
      mean(
        Delta_AUC,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

pA <- ggplot() +
  
  # Species-level raw points
  geom_point(
    data = raw_A,
    aes(
      x = Delta_AUC,
      y = Comparison,
      color = algorithm,
      shape = species,
      group = algorithm
    ),
    alpha = 0.45,
    size = 2.4,
    position = position_jitterdodge(
      jitter.width = 0,
      jitter.height = 0.07,
      dodge.width = 0.55
    )
  ) +
  
  # 95% CI
  geom_errorbar(
    data = emm_A_df,
    aes(
      y = Comparison,
      xmin = lower.CL,
      xmax = upper.CL,
      color = algorithm,
      group = algorithm
    ),
    width = 0,
    linewidth = 1,
    position = position_dodge(
      width = 0.5
    )
  ) +
  
  # Estimated mean
  geom_point(
    data = emm_A_df,
    aes(
      x = emmean,
      y = Comparison,
      color = algorithm,
      group = algorithm
    ),
    size = 3.5,
    shape = 19,
    position = position_dodge(
      width = 0.5
    )
  ) +
  
  # Zero reference
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  scale_colour_manual(
    name = "Modelling algorithm",
    values = algorithm_colors,
    limits = names(algorithm_colors),
    drop = FALSE
  ) +
  
  scale_shape_manual(
    name = "Species",
    values = species_shapes,
    labels = species_labels,
    drop = FALSE
  ) +
  
  labs(
    tag = "A",
    x = expression(
      Delta * "AUC (Climate + LC - Climate)"
    ),
    y = NULL,
    color = "algorithm"
  ) +
  
  theme_classic(base_size = 13)
pA
ggsave(
  "../result/Visualization/Fig2/Figure2A.png",
  pA,
  width = 11,
  height = 7,
  dpi = 600
)

# FIg2B
pairs_B <- tibble(
  Comparison = c(
    "Macro",
    "Macro CP",
    "Micro",
    "Micro CP"
  ),
  
  model_A = c(
    "Lc_int_Macro",
    "Lc_int_MacroCp",
    "Lc_int_Micro",
    "Lc_int_MicroCp"
  ),
  
  model_B = c(
    "Lc_add_Macro",
    "Lc_add_MacroCp",
    "Lc_add_Micro",
    "Lc_add_MicroCp"
  )
)
panelB_Pitangus_sulphuratus <- calculate_delta(
  results_long_Pitangus_sulphuratus,
  pairs_B
)
panelB_final <- rbind(panelB_Ara_macao,panelB_Aramus_guarauna,panelB_Ardea_cocoi,panelB_Tapera_naevia,panelB_Rupornis_magnirostris,panelB_Zonotrichia_capensis,panelB_Furnarius_rufus,panelB_Athene_cunicularia, panelB_Vanellus_chilensis,panelB_Pitangus_sulphuratus)
model_B <- lmer(
  Delta_AUC ~
    Comparison * algorithm +
    (1 | species) +
    (1 | species:fold),
  data = panelB_final
)
anova(model_B)
summary(model_B)
emm_B <- emmeans(
  model_B,
  ~ Comparison | algorithm
)
emm_B_df <- as.data.frame(
  confint(emm_B)
)
raw_B <- panelB_final %>%
  group_by(
    species,
    algorithm,
    Comparison
  ) %>%
  summarise(
    Delta_AUC =
      mean(
        Delta_AUC,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

pB <- ggplot() +
  
  # Species-level raw points
  geom_point(
    data = raw_B,
    aes(
      x = Delta_AUC,
      y = Comparison,
      color = algorithm,
      shape = species,
      group = algorithm
    ),
    alpha = 0.45,
    size = 2.4,
    position = position_jitterdodge(
      jitter.width = 0,
      jitter.height = 0.07,
      dodge.width = 0.55
    )
  ) +
  
  # 95% CI
  geom_errorbar(
    data = emm_B_df,
    aes(
      y = Comparison,
      xmin = lower.CL,
      xmax = upper.CL,
      color = algorithm,
      group = algorithm
    ),
    width = 0,
    linewidth = 1,
    position = position_dodge(
      width = 0.5
    )
  ) +
  
  # Estimated mean
  geom_point(
    data = emm_B_df,
    aes(
      x = emmean,
      y = Comparison,
      color = algorithm,
      group = algorithm
    ),
    size = 3.5,
    shape = 19,
    position = position_dodge(
      width = 0.5
    )
  ) +
  
  # Zero reference
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  scale_colour_manual(
    name = "Modelling algorithm",
    values = algorithm_colors,
    limits = names(algorithm_colors),
    drop = FALSE
  ) +
  
  scale_shape_manual(
    name = "Species",
    values = species_shapes,
    labels = species_labels,
    drop = FALSE
  ) +
  
  labs(
    tag = "B",
    x = expression(
      Delta * "AUC (Climate × LC − Climate + LC)"
    ),
    y = NULL,
    color = "algorithm"
  ) +
  
  theme_classic(base_size = 13)
pB
ggsave(
  "../result/Visualization/Fig2/Figure2B.png",
  pB,
  width = 11,
  height = 7,
  dpi = 600
)

# FIg2C
pairs_C <- tibble(
  Comparison = c(
    "Macro",
    "Macro + LC",
    "Macro × LC",
    "Micro",
    "Micro + LC",
    "Micro × LC"
  ),
  
  model_A = c(
    "MacroCp",
    "Lc_add_MacroCp",
    "Lc_int_MacroCp",
    "MicroCp",
    "Lc_add_MicroCp",
    "Lc_int_MicroCp"
  ),
  
  model_B = c(
    "Macro",
    "Lc_add_Macro",
    "Lc_int_Macro",
    "Micro",
    "Lc_add_Micro",
    "Lc_int_Micro"
  )
)
panelC_Pitangus_sulphuratus <- calculate_delta(
  results_long_Pitangus_sulphuratus,
  pairs_C
)
panelC_final <- rbind(panelC_Ara_macao,panelC_Aramus_guarauna,panelC_Ardea_cocoi,panelC_Tapera_naevia,panelC_Rupornis_magnirostris,panelC_Zonotrichia_capensis,panelC_Furnarius_rufus,panelC_Athene_cunicularia, panelC_Vanellus_chilensis,panelC_Pitangus_sulphuratus)
model_C <- lmer(
  Delta_AUC ~
    Comparison * algorithm +
    (1 | species) +
    (1 | species:fold),
  data = panelC_final
)
anova(model_C)
summary(model_C)
emm_C <- emmeans(
  model_C,
  ~ Comparison | algorithm
)
emm_C_df <- as.data.frame(
  confint(emm_C)
)
raw_C <- panelC_final %>%
  group_by(
    species,
    algorithm,
    Comparison
  ) %>%
  summarise(
    Delta_AUC =
      mean(
        Delta_AUC,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

pC <- ggplot() +
  
  # Species-level raw points
  geom_point(
    data = raw_C,
    aes(
      x = Delta_AUC,
      y = Comparison,
      color = algorithm,
      shape = species,
      group = algorithm
    ),
    alpha = 0.45,
    size = 2.4,
    position = position_jitterdodge(
      jitter.width = 0,
      jitter.height = 0.07,
      dodge.width = 0.55
    )
  ) +
  
  # 95% CI
  geom_errorbar(
    data = emm_C_df,
    aes(
      y = Comparison,
      xmin = lower.CL,
      xmax = upper.CL,
      color = algorithm,
      group = algorithm
    ),
    width = 0,
    linewidth = 1,
    position = position_dodge(
      width = 0.5
    )
  ) +
  
  # Estimated mean
  geom_point(
    data = emm_C_df,
    aes(
      x = emmean,
      y = Comparison,
      color = algorithm,
      group = algorithm
    ),
    size = 3.5,
    shape = 19,
    position = position_dodge(
      width = 0.5
    )
  ) +
  
  # Zero reference
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  scale_colour_manual(
    name = "Modelling algorithm",
    values = algorithm_colors,
    limits = names(algorithm_colors),
    drop = FALSE
  ) +
  
  scale_shape_manual(
    name = "Species",
    values = species_shapes,
    labels = species_labels,
    drop = FALSE
  ) +
  
  labs(
    tag = "C",
    x = expression(
      Delta * "AUC (CP − absolute climate)"
    ),
    y = NULL,
    color = "algorithm"
  ) +
  
  theme_classic(base_size = 13)
pC
ggsave(
  "../result/Visualization/Fig2/Figure2C.png",
  pC,
  width = 11,
  height = 7,
  dpi = 600
)

############## Fig3 
all_results <- rbind(results_long_Ara_macao,results_long_Aramus_guarauna, results_long_Ardea_cocoi, results_long_Tapera_naevia, results_long_Athene_cunicularia, results_long_Furnarius_rufus, results_long_Pitangus_sulphuratus, results_long_Rupornis_magnirostris, results_long_Vanellus_chilensis, results_long_Zonotrichia_capensis)
model_summary <- all_results %>%
  group_by(species, algorithm, config) %>%
  summarise(
    n_folds = n(),
    mean_AUC = mean(AUC, na.rm = TRUE),
    SD = sd(AUC, na.rm = TRUE),
    SE = SD / sqrt(n_folds),
    t_crit = qt(0.975, df = n_folds - 1),
    CI_lower = mean_AUC - t_crit * SE,
    CI_upper = mean_AUC + t_crit * SE,
    .groups = "drop"
  )
best_models <- model_summary %>%
  group_by(species) %>%
  slice_max(
    order_by = mean_AUC,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# Ordered by mean AUC
best_models <- best_models %>%
  mutate(
    species = forcats::fct_reorder(species, mean_AUC)
  )

FIg3 <- ggplot(
  best_models,
  aes(
    x = mean_AUC,
    y = species,
    colour = algorithm
  )
) +
  geom_errorbar(
    aes(
      xmin = CI_lower,
      xmax = CI_upper
    ),
    height = 0.2,
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  # geom_text(
  #   aes(label = config),
  #   hjust = 0.5,
  #   vjust = -1,
  #   nudge_x = 0.015,
  #   colour = "black",
  #   size = 5
  # ) + 
  scale_x_continuous(
    limits = c(0.5, 1),
    breaks = seq(0.5, 1, 0.1)
  ) +
  scale_colour_manual(
    name = "Modelling algorithm",
    values = algorithm_colors,
    limits = names(algorithm_colors),
    drop = FALSE
  ) + 
  labs(
    x = "Mean AUC (95% CI)",
    y = NULL,
    colour = "algorithm"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )
ggsave(
  "../result/Visualization/Fig3/Figure3.png",
  FIg3,
  width = 11,
  height = 7,
  dpi = 600
)

############### FIg4 
points_climate_Zonotrichia_capensis= fread("D:/Paper9/result/Model Input/final/Zonotrichia_capensis_FULL_Macro_Micro_CP_merged.csv")
# For RF
n_balance <- min(sum(points_climate_Zonotrichia_capensis$pa == 1), 
                 sum(points_climate_Zonotrichia_capensis$pa == 0))  # Take the smaller one

bg_idx_down   <- sample(which(points_climate_Zonotrichia_capensis$pa == 0), size = n_balance, replace = FALSE)
pres_idx_down <- sample(which(points_climate_Zonotrichia_capensis$pa == 1), size = n_balance, replace = FALSE)

train_down_Zonotrichia_capensis <- points_climate_Zonotrichia_capensis[c(pres_idx_down, bg_idx_down), ]
train_down_Zonotrichia_capensis$pa_factor <- as.factor(
  ifelse(train_down_Zonotrichia_capensis$pa == 1, "presence", "background")
)

# For Maxent
train_Zonotrichia_capensis_p <- points_climate_Zonotrichia_capensis[points_climate_Zonotrichia_capensis$pa == 1, 
                                   c("Macro_Bio5","Macro_Bio6","Macro_Bio13","Macro_Bio14","landcover")]
train_Zonotrichia_capensis_b  <- points_climate_Zonotrichia_capensis[points_climate_Zonotrichia_capensis$pa == 0, 
                                    c("Macro_Bio5","Macro_Bio6","Macro_Bio13","Macro_Bio14","landcover")]

# ── 32. RF Macroclimate * landcover ──────────────────────────────────────
rf_Zonotrichia_capensis <- randomForest(
  pa_factor ~ Macro_Bio5 + Macro_Bio6 + Macro_Bio13 + Macro_Bio14 + landcover,
  data     = train_down_Zonotrichia_capensis,
  ntree    = 500,
  mtry     = 2,
  importance = TRUE,
  classwt    = c("background" = 1, "presence" = 1)
)
# ── 33. RF MacroCp * landcover ──────────────────────────────────────
rf_Aramus_guarauna <- randomForest(
  pa_factor ~ Macro_CP_Tmin + Macro_CP_Tmax + Macro_CP_SMmax + Macro_CP_SMmin + landcover,
  data     = train_down_Aramus_guarauna,
  ntree    = 500,
  mtry     = 2,
  importance = TRUE,
  classwt    = c("background" = 1, "presence" = 1)
)
# ── 23. maxent Macroclimate * Landcover ──────────────────────────────────────
tryCatch({
  maxent_Zonotrichia_capensis <- maxent(
    x       = rbind(train_Zonotrichia_capensis_p, train_Zonotrichia_capensis_b),
    p       = c(rep(1, nrow(train_Zonotrichia_capensis_p)), 
                rep(0, nrow(train_Zonotrichia_capensis_b))),
    factors =  "landcover",
    args    = c("randomtestpoints=0", "betamultiplier=1",
                "autofeature=false",
                "linear=true",
                "quadratic=true",
                "hinge=true",
                "product=true",        # turn on Product (Cross-product)
                "threshold=false")
  )
})
# ── 7. GLM MacroCp * Landcover ──────────────────────────────────────
glm_Zonotrichia_capensis <- glm(
  pa ~ (Macro_CP_Tmin + Macro_CP_Tmax + Macro_CP_SMmin + Macro_CP_SMmax) * landcover,
  data   = points_climate_Zonotrichia_capensis,
  family = binomial()
)

# Prepare the absolute climatic stack
bioclim_stack_slim <- rast("../result/Model Input/bioclim_stack_slim_1981-2010.tif")
landcover_1km <- rast("../result/Model Input/landcover_1km_1981-2010.tif") 

# If the best model does not involve with CP stack, then directly jump to 'prediction_stack'
species_1981_2010 <- read.csv("../data/occurrence/Zonotrichia capensis.csv",row.names = NULL)

presence <- species_1981_2010[, c("decimalLongitude", "decimalLatitude")]
presence_sf <- st_as_sf(presence, 
                        coords = c("decimalLongitude", "decimalLatitude"), 
                        crs = 4326) # For absence point creation by backgroundPOINTS
bioclim_stack <- rast("../result/Model Input/bioclim_stack_1981-2010.tif")
cells <- cellFromXY(bioclim_stack[[1]], st_coordinates(presence_sf))
unique_cells <- unique(cells)
presence_thinned <- presence_sf[which(!duplicated(cells)), ]

presence_thinned_pts <- presence_thinned %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  rename(lon = X, lat = Y)
pts_sf_range <- st_as_sf(presence_thinned_pts, coords = c("lon", "lat"), crs = 4326)
species_range_polygon <- st_convex_hull(st_union(pts_sf_range))
bioclim_in_species_range <- crop(bioclim_stack_slim, species_range_polygon)
bioclim_in_species_range <- mask(bioclim_in_species_range, vect(species_range_polygon))
all_pixels_in_range <- as.data.frame(bioclim_in_species_range, na.rm = TRUE)
limits <- all_pixels_in_range %>%
  summarise(across(everything(), list(
    min = ~quantile(.x, 0),   # 整个分布区里的绝对最低温
    max = ~quantile(.x, 1)    # 整个分布区里的绝对最高温
  )))
rm(all_pixels_in_range)
gc()
thermal_range <- limits$Tmax_max - limits$Tmin_min
soilmoist_range <-
  limits$SMmax_max -
  limits$SMmin_min

cp_tmin <- (bioclim_stack_slim[["Tmin"]] - limits$Tmin_min) / thermal_range
cp_tmax <- (bioclim_stack_slim[["Tmax"]] - limits$Tmin_min) / thermal_range
cp_smmin <- (bioclim_stack_slim[["SMmin"]] - limits$SMmin_min) / soilmoist_range
cp_smmax <- (bioclim_stack_slim[["SMmax"]] - limits$SMmin_min) / soilmoist_range
cp_stack_slim <- c(cp_tmin, cp_tmax, cp_smmin, cp_smmax)
names(cp_stack_slim) <- c("CP_Tmin", "CP_Tmax", "CP_SMmin", "CP_SMmax")
final_stack_slim <- c(cp_stack_slim, landcover_1km)
names(final_stack_slim) <- c("CP_Tmin", "CP_Tmax", "CP_SMmin", "CP_SMmax", 
                             "landcover")



prediction_stack <- c(
  final_stack_slim,
  landcover_1km
)
prediction_df <- as.data.frame(
  prediction_stack,
  xy = TRUE,
  cells = TRUE,
  na.rm = TRUE
)
names(prediction_df) <- c(
  "cell",
  "lon",
  "lat",
  "Macro_CP_Tmin",
  "Macro_CP_Tmax",
  "Macro_CP_SMmin",
  "Macro_CP_SMmax",
  "landcover"
)
maxent_pred_df <- prediction_df %>%
  dplyr::select(
    Macro_Bio5,
    Macro_Bio6,
    Macro_Bio13,
    Macro_Bio14,
    landcover
  )
prediction_df$probability <- predict(
  glm_Zonotrichia_capensis,
  newdata = prediction_df,
  type = "response"
)
# [, "presence"]
prob_raster <- prediction_stack[[1]]
values(prob_raster) <- NA
values(prob_raster)[
  prediction_df$cell
] <- prediction_df$probability

p_prob <- ggplot() +
  
  # Prediction raster
  geom_spatraster(
    data = prob_raster
  ) +
  
  # Probability colour scale
  scale_fill_viridis_c(
    name = "Occurrence\nprobability",
    limits = c(0, 1),
    option = "viridis",
    na.value = "white"
  ) +
  
  # Neotropical boundary
  geom_sf(
    data = neotropics,
    fill = NA,
    color = "black",
    linewidth = 0.3
  ) +
  
  coord_sf(
    xlim = c(-123.3374, -21.58016),
    ylim = c(-60.34388, 31.11561),
    expand = FALSE
  ) +
  
  labs(
    title = "Zonotrichia capensis",
    subtitle = "Generalised Linear Models | Macro Climatic Position * Landcover",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "italic",
      size = 15
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    axis.title = element_text(
      size = 13
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    )
  )
p_prob
ggsave(
  "../result/Visualization/Fig4/Zonotrichia_capensis_GLM.png",
  p_prob,
  width = 7,
  height = 6,
  dpi = 600
)

