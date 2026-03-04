


# N-mixture Model Analysis for Coyote, Bobcat, and Raccoon
# Date: 10/30/2025

library(AICcmodavg)
library(unmarked)
library(tidyverse)


# DATA LOADING ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set data path
data_path <- ".../data"  

# Load detection matrices
load(file.path(data_path, "detection_matrices.RData"))

# Load site covariates
site_covariates <- read_csv(file.path(data_path, "site_covariates.csv"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# CREATE UNMARKED FRAME OBJECTS FOR MODELLING ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define target species and seasons
target_species <- c("Coyote", "Bobcat", "Raccoon")
seasons <- c("summer", "fall", "winter", "spring")

# Create unmarkedFramePCount objects for each species-season combination
for(species in target_species) {
  for(season in seasons) {
    # Get detection matrix
    detection_matrix <- all_species_matrices[[species]][[season]]
    
    # Create UMF name
    umf_name <- paste0(species, "_", season, "_umf")
    
    # Create unmarkedFramePCount object
    umf_obj <- unmarkedFramePCount(
      y = detection_matrix,
      siteCovs = site_covariates
    )
    
    # Assign to global environment
    assign(umf_name, umf_obj)
  }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# MODEL SELECTION FUNCTION ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to run final model set for a species-season combination
run_final_models <- function(species, season, env_vars, species_vars, K = 50) {
  # Get UMF object
  umf_name <- paste0(species, "_", season, "_umf")
  umf <- get(umf_name)
  
  # Generate all combinations of environmental variables
  env_combos <- list()
  for(i in 0:length(env_vars)) {
    if(i == 0) {
      env_combos[[length(env_combos) + 1]] <- "1"  
    } else {
      combs <- combn(env_vars, i, simplify = FALSE)
      for(combo in combs) {
        env_combos[[length(env_combos) + 1]] <- paste(combo, collapse = " + ")
      }
    }
  }
  
  # Generate all combinations of prey species variables
  species_combos <- list()
  for(i in 0:length(species_vars)) {
    if(i == 0) {
      species_combos[[length(species_combos) + 1]] <- "1" 
    } else {
      combs <- combn(species_vars, i, simplify = FALSE)
      for(combo in combs) {
        species_combos[[length(species_combos) + 1]] <- paste(combo, collapse = " + ")
      }
    }
  }
  
  # Detection formulas
  detection_formulas <- c("~1", "~feat")
  
  # Create final formulas
  model_formulas <- list()
  for(det_formula in detection_formulas) {
    for(env_combo in env_combos) {
      for(sp_combo in species_combos) {
        if(env_combo == "1" && sp_combo == "1") {
          model_formulas[[length(model_formulas) + 1]] <- paste0(det_formula, " ~1")
        } else if(env_combo == "1") {
          model_formulas[[length(model_formulas) + 1]] <- paste0(det_formula, " ~", sp_combo)
        } else if(sp_combo == "1") {
          model_formulas[[length(model_formulas) + 1]] <- paste0(det_formula, " ~", env_combo)
        } else {
          model_formulas[[length(model_formulas) + 1]] <- paste0(det_formula, " ~", env_combo, " + ", sp_combo)
        }
      }
    }
  }
  
  # Remove duplicates
  model_formulas <- unique(model_formulas)
  
  # Run models
  models <- list()
  model_names <- character()
  
  # Initial updates
  cat(paste0("\nRunning final models for ", species, " in ", season, " season\n"))
  cat(paste0("Total models to run: ", length(model_formulas), "\n"))
  
  for(i in 1:length(model_formulas)) {
    form <- as.formula(model_formulas[[i]])
    
    # Progress update every 50 models
    if(i %% 50 == 0) {
      cat(paste0("  Progress: ", i, "/", length(model_formulas), " models completed\n"))
    }
    
    tryCatch({
      models[[i]] <- pcount(form, data = umf, K = K)
      model_names[i] <- model_formulas[[i]]
    }, error = function(e) {
      cat(paste0("  Error in model ", i, ": ", conditionMessage(e), "\n"))
    })
    
    # Force garbage collection
    gc()
  }
  
  # Remove NULL models
  valid_indices <- !sapply(models, is.null)
  models <- models[valid_indices]
  model_names <- model_names[valid_indices]
  
  # Calculate AICc table
  aic_table <- aictab(cand.set = models, modnames = model_names, second.ord = TRUE)
  
  # Return results
  return(list(
    models = models,
    aic_table = aic_table
  ))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# RUN ALL MODELS ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These are the variables for each species-season that had substantial contribution in the initial selection process and 
# progressed to this final model set, prey species variables are separated only for organizational purposes
# Note: Raccoon models do not include prey species activity variables

# Coyote summer
cat("\n=== COYOTE SUMMER ===\n")
coyote_summer_env_vars <- c("sm", "open", "pp_nothin", "pp_har", "tsf_summer_500m_scaled")
coyote_summer_species_vars <- c("Smammal_summer_rate_scaled", "Turkey_summer_rate_scaled", "Wild_Pig_summer_rate_scaled")
coyote_summer_final <- run_final_models("Coyote", "summer", coyote_summer_env_vars, coyote_summer_species_vars, K = 50)

# Coyote fall
cat("\n=== COYOTE FALL ===\n")
coyote_fall_env_vars <- c("sm", "open", "pp_har", "tsf_fall_500m_scaled", "unq_fall_500m_scaled")
coyote_fall_species_vars <- c("Smammal_fall_rate_scaled", "Wild_Pig_fall_rate_scaled", "White_tailed_Deer_fall_rate_scaled")
coyote_fall_final <- run_final_models("Coyote", "fall", coyote_fall_env_vars, coyote_fall_species_vars, K = 50)

# Coyote winter
cat("\n=== COYOTE WINTER ===\n")
coyote_winter_env_vars <- c("bd", "sm", "open", "pp_nothin")
coyote_winter_species_vars <- c("Smammal_winter_rate_scaled")
coyote_winter_final <- run_final_models("Coyote", "winter", coyote_winter_env_vars, coyote_winter_species_vars, K = 50)

# Coyote spring
cat("\n=== COYOTE SPRING ===\n")
coyote_spring_env_vars <- c("sm", "open", "pp_nothin")
coyote_spring_species_vars <- c("Wild_Pig_spring_rate_scaled")
coyote_spring_final <- run_final_models("Coyote", "spring", coyote_spring_env_vars, coyote_spring_species_vars, K = 50)

# Bobcat summer
cat("\n=== BOBCAT SUMMER ===\n")
bobcat_summer_env_vars <- c("sm")
bobcat_summer_species_vars <- c()
bobcat_summer_final <- run_final_models("Bobcat", "summer", bobcat_summer_env_vars, bobcat_summer_species_vars, K = 50)

# Bobcat fall
cat("\n=== BOBCAT FALL ===\n")
bobcat_fall_env_vars <- c("bd", "open", "unq_fall_500m_scaled", "pp_nothin")
bobcat_fall_species_vars <- c("Smammal_fall_rate_scaled", "Wild_Pig_fall_rate_scaled")
bobcat_fall_final <- run_final_models("Bobcat", "fall", bobcat_fall_env_vars, bobcat_fall_species_vars, K = 50)

# Bobcat winter
cat("\n=== BOBCAT WINTER ===\n")
bobcat_winter_env_vars <- c("pp_har", "tsf_winter_500m_scaled", "unq_winter_500m_scaled")
bobcat_winter_species_vars <- c()
bobcat_winter_final <- run_final_models("Bobcat", "winter", bobcat_winter_env_vars, bobcat_winter_species_vars, K = 50)

# Bobcat spring
cat("\n=== BOBCAT SPRING ===\n")
bobcat_spring_env_vars <- c("open", "pp_nothin", "pp_har", "tsf_spring_500m_scaled")
bobcat_spring_species_vars <- c("Smammal_spring_rate_scaled", "White_tailed_Deer_spring_rate_scaled", "Deer_fawn_spring_rate_scaled")
bobcat_spring_final <- run_final_models("Bobcat", "spring", bobcat_spring_env_vars, bobcat_spring_species_vars, K = 50)

# Raccoon summer
cat("\n=== RACCOON SUMMER ===\n")
raccoon_summer_env_vars <- c("sm", "open", "pp_har", "tsf_summer_500m_scaled", "unq_summer_500m_scaled", "bd")
raccoon_summer_species_vars <- c()
raccoon_summer_final <- run_final_models("Raccoon", "summer", raccoon_summer_env_vars, raccoon_summer_species_vars, K = 50)

# Raccoon fall
cat("\n=== RACCOON FALL ===\n")
raccoon_fall_env_vars <- c("sm", "open", "pp_har", "tsf_fall_500m_scaled", "unq_fall_500m_scaled")
raccoon_fall_species_vars <- c()
raccoon_fall_final <- run_final_models("Raccoon", "fall", raccoon_fall_env_vars, raccoon_fall_species_vars, K = 50)

# Raccoon winter
cat("\n=== RACCOON WINTER ===\n")
raccoon_winter_env_vars <- c("bd", "sm", "open", "pp_nothin", "pp_har", "tsf_winter_500m_scaled", "unq_winter_500m_scaled")
raccoon_winter_species_vars <- c()
raccoon_winter_final <- run_final_models("Raccoon", "winter", raccoon_winter_env_vars, raccoon_winter_species_vars, K = 50)

# Raccoon spring
cat("\n=== RACCOON SPRING ===\n")
raccoon_spring_env_vars <- c("bd", "sm", "open", "pp_nothin", "unq_spring_500m_scaled")
raccoon_spring_species_vars <- c()
raccoon_spring_final <- run_final_models("Raccoon", "spring", raccoon_spring_env_vars, raccoon_spring_species_vars, K = 100)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# VIEW TOP MODELS ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n=== TOP MODELS BY SPECIES AND SEASON ===\n")

# Coyote
cat("\nCOYOTE - Summer:\n"); print(head(coyote_summer_final$aic_table, 1))
cat("\nCOYOTE - Fall:\n"); print(head(coyote_fall_final$aic_table, 1))
cat("\nCOYOTE - Winter:\n"); print(head(coyote_winter_final$aic_table, 1))
cat("\nCOYOTE - Spring:\n"); print(head(coyote_spring_final$aic_table, 1))

# Bobcat
cat("\nBOBCAT - Summer:\n"); print(head(bobcat_summer_final$aic_table, 1))
cat("\nBOBCAT - Fall:\n"); print(head(bobcat_fall_final$aic_table, 1))
cat("\nBOBCAT - Winter:\n"); print(head(bobcat_winter_final$aic_table, 1))
cat("\nBOBCAT - Spring:\n"); print(head(bobcat_spring_final$aic_table, 1))

# Raccoon
cat("\nRACCOON - Summer:\n"); print(head(raccoon_summer_final$aic_table, 1))
cat("\nRACCOON - Fall:\n"); print(head(raccoon_fall_final$aic_table, 1))
cat("\nRACCOON - Winter:\n"); print(head(raccoon_winter_final$aic_table, 1))
cat("\nRACCOON - Spring:\n"); print(head(raccoon_spring_final$aic_table, 1))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# VISUALIZE MODEL RESULTS - PLOTTING BETAS ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This creates one plot per species, showing beta estimates for each variable in the top model

# Set up plotting parameters
season_colors <- c("Summer" = "#d62728", "Fall" = "#ff7f0e", 
                   "Winter" = "#1f77b4", "Spring" = "#2ca02c")
season_shapes <- c("Summer" = 16, "Fall" = 17, "Winter" = 15, "Spring" = 18)
season_order <- c("Summer", "Fall", "Winter", "Spring")

# Define covariate order for consistent x-axis
covariate_order <- c("Distance to\nOpen", 
                     "Distance to\nUpland Hardwood", 
                     "Distance to\nBottomland Hardwood",
                     "Distance to Recently\nThinned Pine Plantation", 
                     "Distance to Not Recently\nThinned Pine Plantation",
                     "Time Since Fire", 
                     "Pyrodiversity", 
                     "Small Prey Activity", 
                     "Turkey Activity", 
                     "Pig Activity", 
                     "Adult Deer Activity", 
                     "Fawn Activity")


## Coyote beta estimates ----
coyote_summer <- data.frame(
  feature = c("Distance to\nUpland Hardwood", "Distance to\nOpen", "Time Since Fire", 
              "Small Prey Activity", "Turkey Activity"),
  estimate = c(-0.472, 0.152, -0.258, 0.418, 0.442),
  se = c(0.0729, 0.0695, 0.0963, 0.0490, 0.0376),
  season = "Summer"
)

coyote_fall <- data.frame(
  feature = c("Distance to\nUpland Hardwood", "Distance to\nOpen", "Time Since Fire", 
              "Pyrodiversity", "Small Prey Activity", "Pig Activity"),
  estimate = c(-0.278, 0.327, -0.298, -0.143, 0.213, 0.131),
  se = c(0.0727, 0.0642, 0.0915, 0.076, 0.0472, 0.0498),
  season = "Fall"
)

coyote_winter <- data.frame(
  feature = c("Distance to\nBottomland Hardwood", "Distance to\nUpland Hardwood", 
              "Distance to\nOpen", "Distance to Not Recently\nThinned Pine Plantation", 
              "Small Prey Activity"),
  estimate = c(0.235, -0.500, 0.309, -0.347, 0.168),
  se = c(0.0692, 0.0882, 0.0624, 0.0800, 0.0446),
  season = "Winter"
)

coyote_spring <- data.frame(
  feature = c("Distance to\nUpland Hardwood", "Distance to\nOpen", 
              "Distance to Not Recently\nThinned Pine Plantation", "Pig Activity"),
  estimate = c(-0.472, 0.524, -0.475, 0.401),
  se = c(0.1071, 0.0763, 0.0997, 0.0483),
  season = "Spring"
)

# Combine and prepare data
coyote_all <- rbind(coyote_summer, coyote_fall, coyote_winter, coyote_spring)
coyote_all$lower <- coyote_all$estimate - qnorm(0.975) * coyote_all$se
coyote_all$upper <- coyote_all$estimate + qnorm(0.975) * coyote_all$se
coyote_all$season <- factor(coyote_all$season, levels = season_order)
coyote_all$feature <- factor(coyote_all$feature, levels = covariate_order)

# Create plot
coyote_plot <- ggplot(coyote_all, aes(x = feature, y = estimate, 
                                      color = season, shape = season)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.3) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 1)) + 
  scale_color_manual(values = season_colors, breaks = season_order) +
  scale_shape_manual(values = season_shapes, breaks = season_order) +
  labs(title = "Coyote Beta Estimates by Season",
       x = "", 
       y = "Beta Estimate",
       color = "Season",
       shape = "Season") +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "top")

print(coyote_plot)


## Bobcat beta estimates ----
bobcat_spring <- data.frame(
  feature = c("Distance to\nOpen", "Distance to Not Recently\nThinned Pine Plantation", 
              "Distance to Recently\nThinned Pine Plantation", "Small Prey Activity", 
              "Adult Deer Activity", "Fawn Activity"),
  estimate = c(0.6167, 0.5221, -0.8539, 0.3073, 0.6370, -1.5412),
  se = c(0.195, 0.207, 0.303, 0.147, 0.204, 0.519),
  season = "Spring"
)

bobcat_summer <- data.frame(
  feature = c("Distance to\nUpland Hardwood"),
  estimate = c(-0.637),
  se = c(0.303),
  season = "Summer"
)

bobcat_fall <- data.frame(
  feature = c("Distance to\nOpen", "Small Prey Activity", "Pig Activity"),
  estimate = c(0.331, 0.523, -0.469),
  se = c(0.1257, 0.0884, 0.1756),
  season = "Fall"
)

bobcat_winter <- data.frame(
  feature = c("Distance to Recently\nThinned Pine Plantation", 
              "Time Since Fire", "Pyrodiversity"),
  estimate = c(0.616, -0.442, -0.299),
  se = c(0.113, 0.188, 0.143),
  season = "Winter"
)

# Combine and prepare data
bobcat_all <- rbind(bobcat_spring, bobcat_summer, bobcat_fall, bobcat_winter)
bobcat_all$lower <- bobcat_all$estimate - qnorm(0.975) * bobcat_all$se
bobcat_all$upper <- bobcat_all$estimate + qnorm(0.975) * bobcat_all$se
bobcat_all$season <- factor(bobcat_all$season, levels = season_order)
bobcat_all$feature <- factor(bobcat_all$feature, levels = covariate_order)

# Create plot
bobcat_plot <- ggplot(bobcat_all, aes(x = feature, y = estimate, 
                                      color = season, shape = season)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.3) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(limits = c(-3, 1.5), breaks = seq(-3, 1, by = 1)) + 
  scale_color_manual(values = season_colors, breaks = season_order) +
  scale_shape_manual(values = season_shapes, breaks = season_order) +
  labs(title = "Bobcat Beta Estimates by Season",
       x = "", 
       y = "Beta Estimate",
       color = "Season",
       shape = "Season") +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "top")

print(bobcat_plot)


## Raccoon beta estimates ----
raccoon_spring <- data.frame(
  feature = c("Distance to\nBottomland Hardwood", "Distance to\nUpland Hardwood", 
              "Distance to\nOpen", "Distance to Not Recently\nThinned Pine Plantation", 
              "Pyrodiversity"),
  estimate = c(-0.738, -0.411, -0.515, 1.044, -0.779),
  se = c(0.1071, 0.0565, 0.0766, 0.0461, 0.0676),
  season = "Spring"
)

raccoon_summer <- data.frame(
  feature = c("Distance to\nUpland Hardwood", "Distance to Recently\nThinned Pine Plantation", 
              "Time Since Fire", "Pyrodiversity"),
  estimate = c(-0.508, 0.227, -0.683, -0.603),
  se = c(0.0588, 0.0511, 0.1028, 0.0819),
  season = "Summer"
)

raccoon_fall <- data.frame(
  feature = c("Distance to\nUpland Hardwood", "Distance to\nOpen", 
              "Distance to Recently\nThinned Pine Plantation", 
              "Time Since Fire", "Pyrodiversity"),
  estimate = c(-0.667, -0.257, 0.480, -0.721, -1.067),
  se = c(0.0678, 0.0840, 0.0513, 0.0873, 0.0976),
  season = "Fall"
)

raccoon_winter <- data.frame(
  feature = c("Distance to\nOpen", "Distance to Not Recently\nThinned Pine Plantation", 
              "Distance to Recently\nThinned Pine Plantation"),
  estimate = c(-0.526, 0.341, -0.192),
  se = c(0.1091, 0.0518, 0.0789),
  season = "Winter"
)

# Combine and prepare data
raccoon_all <- rbind(raccoon_spring, raccoon_summer, raccoon_fall, raccoon_winter)
raccoon_all$lower <- raccoon_all$estimate - qnorm(0.975) * raccoon_all$se
raccoon_all$upper <- raccoon_all$estimate + qnorm(0.975) * raccoon_all$se
raccoon_all$season <- factor(raccoon_all$season, levels = season_order)
raccoon_all$feature <- factor(raccoon_all$feature, levels = covariate_order)

# Create plot
raccoon_plot <- ggplot(raccoon_all, aes(x = feature, y = estimate, 
                                        color = season, shape = season)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5), width = 0.3) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-1, 1, by = 1)) + 
  scale_color_manual(values = season_colors, breaks = season_order) +
  scale_shape_manual(values = season_shapes, breaks = season_order) +
  labs(title = "Raccoon Beta Estimates by Season",
       x = "", 
       y = "Beta Estimate",
       color = "Season",
       shape = "Season") +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "top")

print(raccoon_plot)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




