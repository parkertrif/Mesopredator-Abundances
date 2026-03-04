# Mesopredator Abundances

## Overview
This repository contains the data and code for analyzing mesopredator abundance using N-mixture models. The analysis examines seasonal patterns in relative abundance of coyotes (*Canis latrans*), bobcats (*Lynx rufus*), and raccoons (*Procyon lotor*) in relation to landcover, forest management, and prey activity covariates using camera trap data.


## Repository Structure
```
├── data/
│   ├── detection_matrices.RData
│   └── site_covariates.csv
├── scripts/
│   └── MA_Code.R
└── README.md
```

## Data Description

### detection_matrices.RData
Contains detection matrices for each target species across four seasons (summer, fall, winter, spring). Each matrix contains weekly counts of the number of individuals detected at each camera location over a 9-week sampling season.

**Structure:**
- Nested list: `all_species_matrices[[species]][[season]]`
- Species: Coyote, Bobcat, Raccoon
- Seasons: summer (Jul 27 2023 - Sep 28 2023), fall (Oct 15 2023 - Dec 17 2023), winter (Dec 29 2023 - Mar 1 2024), spring (May 1 2024 - Jul 3 2024)
- Each matrix: rows = camera locations, columns = weekly sampling occasions
- NA values indicate camera non-operational periods

### site_covariates.csv
Covariates for each camera location.

**Columns:**

*Site identifiers:*
- `placename`: Unique camera location identifier
- `feat`: Trail feature type (categorical: primary road, secondary road, field)

*Landcover covariates (all centered scaled):*
- `pine`: Distance to pine 
- `bd`: Distance to bottomland hardwood
- `sm`: Distance to upland hardwood
- `open`: Distance to open 

*Forest management covariates (all centered scaled):*
- `pp_nothin`: Distance to pine plantation not recently thinned
- `pp_har`: Distance to recently thinned pine plantation
- `tsf_[season]_500m_scaled`: Time since fire by season in 500m buffer
- `unq_[season]_500m_scaled`: Number of unique fires events (pyrodiversity) by season in 500m buffer

*Prey species activity covariates (detections per trap night: all centered and scaled):*
- `Smammal_[season]_rate_scaled`: Small prey detection rate (squirrel, rabbit, armadillo)
- `Turkey_[season]_rate_scaled`: Wild turkey detection rate
- `Wild_Pig_[season]_rate_scaled`: Wild pig detection rate
- `White_tailed_Deer_[season]_rate_scaled`: Adult white-tailed deer detection rate
- `Deer_fawn_[season]_rate_scaled`: Fawn white-tailed deer detection rate (spring and summer only)
  

## Analysis 

The analysis script (`MA_Code.R`) performs the following steps:

1. **Data Loading**: Imports detection matrices and site covariates

2. **Model Setup**: Creates unmarkedFramePCount objects for each species-season combination

3. **Model Selection**: Runs model selection using N-mixture models with:
   - Detection covariates: intercept-only or feature type
   - Abundance covariates: combinations of landcover, management, and prey activity variables
   - Upper bound on abundance (K): 50 for most models, 100 for raccoon spring

4. **Model Ranking**: Identifies top-ranked model using AICc for each species-season combination

5. **Visualization**: Creates coefficient plots showing beta estimates and 95% confidence intervals for abundance covariates across seasons, one plot per species

## Requirements

### R Version
R version >= 4.0.0

### R Packages
```r
library(AICcmodavg)   # Model selection and comparison
library(unmarked)     # N-mixture modeling
library(tidyverse)    # Data manipulation and visualization
```

## Reproducibility & Publishing
This repository is prepared to be published on Dryad alongside the associated manuscript. It contains all data and code needed to reproduce key results, following best practices for reproducible open science.
