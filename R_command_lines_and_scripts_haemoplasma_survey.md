# **Haemoplasma epidemiological survey : R command lines and script**

EN TRAVAUX - A MODIFIER
We analyzed data from 175 wild sloths captured between 1994 and 1995 during the flooding of the Petit Saut Dam (5°03′43″ N, 53°03′00″ O) on the Sinnamary River (French Guiana, South America). The clinical data include the following variables for each examined sloth: 
- `species` : Sloth species (Bt: *Bradypus tridactylus*; Cd: *Choloepus didactylus*)
- `sex` : Sex of the sloth (F: Female; M: Male)
- `age_class` : Age category (A: Adult; J: Juvenile)
- `season` : Season of capture (W: Wet; D: Dry)
- `weight` : Body weight (quantitative variable, in kg)
- `total_length` : Total body length (quantitative variable, in cm)
- `wither_height` : Height at the withers (quantitative variable, in cm)
- `neck_size` : Neck circumference (quantitative variable, in cm)
- `temperature` : Body temperature (quantitative variable, in °C)
- `hematocrit` : Hematocrit level (quantitative variable, in %)
- `health_condition` : Overall health status (G: Good; D: Deteriorated)
- `haemoplasma` : Infection status with haemotropic mycoplasmas (0: Uninfected; 1: Infected)
- `anaplasma` : Infection status with *Anaplasma* (0: Uninfected; 1: Infected)
- `tick` : Presence of ticks in the fur (0: Absent; 1: Present)
- `microfilaria` : Infection status with microfilariae (0: Uninfected; 1: Infected)
- `trypanosome` : Infection status with trypanosomes (0: Uninfected; 1: Infected)
- `babesia` : Infection status with _Babesia_ (0: Uninfected; 1: Infected)
- `bloodparasite` : Combined infection status for blood parasites/pathogens (_Anaplasma_ + microfilariae + trypanosome + _Babesia_, but excluding haemotropic mycoplasmas; 0: Uninfected; 1: Infected)
  
Details about all the experimental methods and measures are available in the related manuscript.


## Table of contents 
- [Step 1. Retrieving the data](#step-1-retrieving-the-data)
- [Step 2. Prepare the data for analysis](#step-2-prepare-the-data-for-analysis)
- [Step 3. Calculate *Anaplasma* infection prevalence](#step-3-calculate-anaplasma-infection-prevalence)
- [Step 4. Test whether _Anaplasma_ infection prevalence in _Bradypus tridactylus_ (Bt) is influenced by sex, age, season, ticks and blood parasites (GLM model 1)](#step-4-test-whether-anaplasma-infection-prevalence-in-bradypus-tridactylus-bt-is-influenced-by-sex-age-season-ticks-and-blood-parasites-glm-model-1)
- [Step 5. Test whether _Anaplasma_ infection prevalence in _Choloepus didactylus_ (Cd) is influenced by sex, age, season, ticks and blood parasites (GLM model 2)](#step-5-test-whether-anaplasma-infection-prevalence-in-choloepus-didactylus-cd-is-influenced-by-sex-age-season-ticks-and-blood-parasites-glm-model-2)
- [Step 6. Test whether the proportion of sloths carrying ticks and blood parasites vary between seasons](#step-6-test-whether-the-proportion-of-sloths-carrying-ticks-and-blood-parasites-vary-between-seasons)
- [Step 7. Impact of _Anaplasma_ infections on Scale Mass Index (SMI) (GLM models 3 and 4)](#step-7-impact-of-anaplasma-infections-on-scale-mass-index-smi-glm-models-3-and-4)
- [Step 8. Impact of _Anaplasma_ infections on neck circumference (GLM models 5 and 6)](#step-8-impact-of-anaplasma-infections-on-neck-circumference-glm-models-5-and-6)
- [Step 9. Impact of _Anaplasma_ infections on hematocrit levels (GLM models 7, 8 and 9)](#step-9-impact-of-anaplasma-infections-on-hematocrit-levels-glm-models-7-8-and-9)
- [Step 10. Impact of _Anaplasma_ infections on body temperature (CLRM models 10 and 11)](#step-10-impact-of-anaplasma-infections-on-body-temperature-clrm-models-10-and-11)
- [Step 11. Impact of _Anaplasma_ infections on general health condition](#step-11-impact-of-anaplasma-infections-on-general-health-condition)
- [Step 12. Impact of _Anaplasma_ infections on female reproductive status](#step-12-impact-of-anaplasma-infections-on-female-reproductive-status)

## Step 1. Retrieving the data

All veterinary clinical data for the two sloth species are available here: (https://github.com/olivierduron/Haemoplasma_sloth_infections/blob/main/data_haemoplasma_sloth.csv)

This database will be referred to as `data_haemoplasma` throughout the R command lines and scripts provided below. It corresponds to the dataset provided in Table S1 of the related manuscript.

Load the dataset directly from the GitHub repository to R:
```
data_haemoplasma <- read.csv ("https://raw.githubusercontent.com/olivierduron/Haemoplasma_sloth_infections/main/data_haemoplasma_sloth.csv", sep = "\t")
```


## Step 2. Prepare the data for analysis

Convert categorical variables into factors:
```
data_haemoplasma$haemoplasma      <- as.factor(data_haemoplasma$haemoplasma)
data_haemoplasma$anaplasma      <- as.factor(data_haemoplasma$anaplasma)
data_haemoplasma$species        <- as.factor(data_haemoplasma$species)
data_haemoplasma$season         <- as.factor(data_haemoplasma$season)
data_haemoplasma$sex            <- as.factor(data_haemoplasma$sex)
data_haemoplasma$age            <- as.factor(data_haemoplasma$age)
data_haemoplasma$tick           <- as.factor(data_haemoplasma$tick)
data_haemoplasma$microfilaria   <- as.factor(data_haemoplasma$microfilaria)
data_haemoplasma$trypanosome    <- as.factor(data_haemoplasma$trypanosome)
data_haemoplasma$babesia        <- as.factor(data_haemoplasma$babesia)
data_haemoplasma$bloodparasite  <- as.factor(data_haemoplasma$bloodparasite)
```

Load libraries for analysis: 
```
library(binom)
library(dplyr)
library(MASS)
library(ggplot2)
library(patchwork)
library(smatr)
library(lmtest)
library(akima)
library(pwr)
library(survival)
library(RColorBrewer)
```

## Step 3. Calculate haemoplasma infection prevalence
Calculate haemoplasma infection prevalence and 95% confidence interval for _Bradypus tridactylus_ (Bt) and _Choloepus didactylus_ (Cd):

```
prevalence_results <- data_haemoplasma %>% group_by(species) %>% summarise(n = n(), positives = sum(haemoplasma == 1), prevalence = positives / n, conf_low = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$lower, conf_high = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$upper)
print(prevalence_results)
```

Results are:
```
A tibble: 2 × 6
species     n positives prevalence conf_low conf_high
  <fct>   <int>     <int>      <dbl>    <dbl>     <dbl>
1 Bt         92         4     0.0435   0.0120     0.108
2 Cd         83        68     0.819    0.720      0.895
```

