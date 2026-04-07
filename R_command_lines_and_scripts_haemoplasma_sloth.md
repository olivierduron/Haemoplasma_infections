# **Haemoplasma veterinary evaluations : R command lines and script**

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

Test if `haemoplasma` is influenced by sloth `species`:
```
chisq.test(table(data_haemoplasma$haemoplasma, data_haemoplasma$species))
```

Results are:
```
Pearson's Chi-squared test with Yates' continuity correction
data:  table(data_haemoplasma$haemoplasma, data_haemoplasma$species)
X-squared = 105.27, df = 1, p-value < 2.2e-16
```

## Step 4. Test whether haemoplasma infection prevalence in _Bradypus tridactylus_ (Bt) is influenced by sex, age, season, ticks and other blood parasites (GLM model 1)
Create a subset `data_Bt` containing only records for _Bradypus tridactylus_ (Bt):

```
data_Bt <- subset(data_haemoplasma, species == "Bt")
```

Fit a GLM to test whether `haemoplasma` is influenced by interactions among `sex`, `age`, `season`, `tick`, and `bloodparasite` in Bt:
```
model_1 <- glm(haemoplasma ~ sex * age * season * tick * bloodparasite, data = data_Bt, family = binomial)
```

Fit a GLM to test whether `haemoplasma` infection prevalence is influenced by additive effects of `sex`, `age`, `season`, `tick`, and `bloodparasite` in Bt:
```
model_1a <- glm(haemoplasma ~ sex + age + season + tick + bloodparasite, data = data_Bt, family = binomial)
```

Compare the additive model (model_1a) to the interaction model (model_1) using a likelihood ratio test:
```
anova(model_1a, model_1, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: haemoplasma ~ sex + age + season + tick + bloodparasite
Model 2: haemoplasma ~ sex * age * season * tick * bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        86     25.402                     
2        74     19.733 12   5.6681   0.9319
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_1, model_1a)
```

Results are:
```
         df      AIC
model_1  18 55.73344
model_1a  6 37.40157
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_1a, test = "Chisq")
```

Results are:
```
Single term deletions
Model: haemoplasma ~ sex + age + season + tick + bloodparasite
              Df Deviance    AIC    LRT Pr(>Chi)  
<none>             25.402 37.402                  
sex            1   25.416 35.416 0.0144  0.90461  
age            1   25.509 35.508 0.1069  0.74370  
season         1   29.313 39.313 3.9112  0.04797 *
tick           1   27.412 37.412 2.0103  0.15623  
bloodparasite  1   27.707 37.707 2.3050  0.12896  
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_1a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        37.402   0.00000
sex           35.416  -1.98564
age           35.508  -1.89310
season        39.313   1.91118
tick          37.412   0.01031
bloodparasite 37.707   0.30503
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model_null <- glm(haemoplasma ~ 1, data = data_Bt, family = binomial)
model_sex <- glm(haemoplasma ~ sex, data = data_Bt, family = binomial)
model_age <- glm(haemoplasma ~ age, data = data_Bt, family = binomial)
model_season <- glm(haemoplasma ~ season, data = data_Bt, family = binomial)
model_tick <- glm(haemoplasma ~ tick, data = data_Bt, family = binomial)
model_bloodparasite <- glm(haemoplasma ~ bloodparasite, data = data_Bt, family = binomial)
anova(model_null, model_sex, test="Chisq")
anova(model_null, model_age, test="Chisq")
anova(model_null, model_season, test="Chisq")
anova(model_null, model_tick, test="Chisq")
anova(model_null, model_bloodparasite, test="Chisq")
aics <- AIC(model_null, model_sex, model_age, model_season, model_tick, model_bloodparasite)
aic_null <- aics["model_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        91     32.907                     
2        90     32.890  1 0.017826   0.8938
> anova(model_null, model_age, test="Chisq")
---
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ age
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        91     32.907                     
2        90     32.450  1  0.45735   0.4989
> anova(model_null, model_season, test="Chisq")
---
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ season
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        91     32.907                       
2        90     29.111  1   3.7967  0.05135 .
---
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ tick
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        91     32.907                     
2        90     31.659  1   1.2483   0.2639
---
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        91     32.907                     
2        90     31.430  1   1.4775   0.2242
---
                         AIC delta_AIC_vs_null
model_null          34.90746         0.0000000
model_sex           36.88964         1.9821737
model_age           36.45012         1.5426530
model_season        33.11076        -1.7967066
model_tick          35.65919         0.7517249
model_bloodparasite 35.42993         0.5224647
```

Tests for associations between `haemoplasma` and the presence of blood parasites (`anaplasma`, `microfilaria`, `trypanosome`, `babesia`) considered separately in Bt:
```
fisher.test(table(data_Bt$haemoplasma, data_Bt$anaplasma))
fisher.test(table(data_Bt$haemoplasma, data_Bt$microfilaria))  
fisher.test(table(data_Bt$haemoplasma, data_Bt$trypanosome))  
fisher.test(table(data_Bt$haemoplasma, data_Bt$babesia))
```

Results are:
```
Fisher's Exact Test for Count Data
data:  table(data_Bt$haemoplasma, data_Bt$anaplasma)
p-value = 0.6245
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.03989675 8.28720044
sample estimates:
odds ratio 
  0.575116 
---
data:  table(data_Bt$haemoplasma, data_Bt$microfilaria)
p-value = 0.2962
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.000000 2.964093
sample estimates:
odds ratio 
         0 
---
data:  table(data_Bt$haemoplasma, data_Bt$trypanosome)
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.00000 64.56697
sample estimates:
odds ratio 
         0 
---
data:  table(data_Bt$haemoplasma, data_Bt$babesia)
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
   0 Inf
sample estimates:
odds ratio 
         0 
```

Display the proportion of Bt sloths infected by `haemoplasma` in wet and dry `season`:
```
table_haemoplasma_season_Bt <- table(data_Bt$haemoplasma, data_Bt$season)
table_haemoplasma_season_Bt
```

Results are:
```
     D  W
  0 34 54
  1  0  4
```

Tests for associations between `haemoplasma` and `season` in Bt:
```
fisher.test(table(data_Bt$haemoplasma, data_Bt$season))  
```

Results are:
```
Fisher's Exact Test for Count Data
data:  table(data_Bt$haemoplasma, data_Bt$season)
p-value = 0.2927
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.3909555       Inf
sample estimates:
odds ratio 
       Inf 
```

## Step 5. Test whether haemoplasma infection prevalence in _Choloepus didactylus_ (Cd) is influenced by sex, age, season, ticks and other blood parasites (GLM model 2)

Create a subset `data_Cd` containing only records for _Choloepus didactylus_ (Cd):
```
data_Cd <- subset(data_haemoplasma, species == "Cd")
```

Fit a GLM to test whether `haemoplasma` is influenced by interactions among `sex`, `age`, `season`, `tick`, and `bloodparasite` in Cd:
```
model_2 <- glm(haemoplasma ~ sex * age * season * tick * bloodparasite, data = data_Cd, family = binomial)
```

Fit a GLM to test whether `haemoplasma` infection prevalence is influenced by additive effects of `sex`, `age`, `season`, `tick`, and `bloodparasite` in Cd:
```
model_2a <- glm(haemoplasma ~ sex + age + season + tick + bloodparasite, data = data_Cd, family = binomial)
```

Compare the additive model (model_2a) to the interaction model (model_2) using a likelihood ratio test:
```
anova(model_2a, model_2, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: haemoplasma ~ sex + age + season + tick + bloodparasite
Model 2: haemoplasma ~ sex * age * season * tick * bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        77     67.088                     
2        57     46.644 20   20.445   0.4304
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_2, model_2a)
```

Results are:
```
         df      AIC
model_2  26 98.64370
model_2a  6 79.08835
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_2a, test = "Chisq")
```

Results are:
```
Single term deletions
Model: haemoplasma ~ sex + age + season + tick + bloodparasite
              Df Deviance    AIC    LRT Pr(>Chi)  
<none>             67.088 79.088                  
sex            1   67.094 77.094 0.0058  0.93954  
age            1   67.806 77.806 0.7181  0.39677  
season         1   70.172 80.172 3.0832  0.07910 .
tick           1   68.851 78.851 1.7625  0.18431  
bloodparasite  1   71.140 81.140 4.0514  0.04413 *
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_2a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        79.088   0.00000
sex           77.094  -1.99425
age           77.806  -1.28192
season        80.172   1.08325
tick          78.851  -0.23745
bloodparasite 81.140   2.05144
```

Compare the null model (model2_null) to univariate models using likelihood ratio tests and AIC:
```
model2_null <- glm(haemoplasma ~ 1, data = data_Cd, family = binomial)
model2_sex <- glm(haemoplasma ~ sex, data = data_Cd, family = binomial)
model2_age <- glm(haemoplasma ~ age, data = data_Cd, family = binomial)
model2_season <- glm(haemoplasma ~ season, data = data_Cd, family = binomial)
model2_tick <- glm(haemoplasma ~ tick, data = data_Cd, family = binomial)
model2_bloodparasite <- glm(haemoplasma ~ bloodparasite, data = data_Cd, family = binomial)
anova(model2_null, model2_sex, test="Chisq")
anova(model2_null, model2_age, test="Chisq")
anova(model2_null, model2_season, test="Chisq")
anova(model2_null, model2_tick, test="Chisq")
anova(model2_null, model2_bloodparasite, test="Chisq")
aics <- AIC(model2_null, model2_sex, model2_age, model2_season, model2_tick, model2_bloodparasite)
aic_null <- aics["model2_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        82     78.433                     
2        81     77.984  1  0.44913   0.5027
---
Analysis of Deviance Table
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ age
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        82     78.433                     
2        81     78.149  1   0.2835   0.5944
---
Analysis of Deviance Table
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ season
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        82     78.433                       
2        81     73.855  1   4.5779  0.03239 *
---
Analysis of Deviance Table
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ tick
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        82     78.433                       
2        81     75.076  1   3.3566  0.06694 .
---
Analysis of Deviance Table
Model 1: haemoplasma ~ 1
Model 2: haemoplasma ~ bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        82     78.433                       
2        81     73.490  1   4.9434  0.02619 *
---
                          AIC delta_AIC_vs_null
model2_null          80.43299          0.000000
model2_sex           81.98386          1.550870
model2_age           82.14949          1.716502
model2_season        77.85506         -2.577930
model2_tick          79.07643         -1.356560
model2_bloodparasite 77.48955         -2.943435
```

Tests for associations between `haemoplasma` and the presence of blood parasites (`anaplasma`, `microfilaria`, `trypanosome`, `babesia`) considered separately in Cd:
```
fisher.test(table(data_Cd$haemoplasma, data_Cd$anaplasma))
fisher.test(table(data_Cd$haemoplasma, data_Cd$microfilaria))  
fisher.test(table(data_Cd$haemoplasma, data_Cd$trypanosome))  
fisher.test(table(data_Cd$haemoplasma, data_Cd$babesia))
```

Results are:
```
Fisher's Exact Test for Count Data
data:  table(data_Cd$haemoplasma, data_Cd$anaplasma)
p-value = 0.02172
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  1.126489 28.191553
sample estimates:
odds ratio 
  4.689964 
---
data:  table(data_Cd$haemoplasma, data_Cd$microfilaria)
p-value = 0.4456
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
   0.3757521 137.1929723
sample estimates:
odds ratio 
  2.969683 
---
data:  table(data_Cd$haemoplasma, data_Cd$trypanosome)
p-value = 0.3306
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.002632914 17.565042872
sample estimates:
odds ratio 
 0.2146914 
---
data:  table(data_Cd$haemoplasma, data_Cd$babesia)
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.1447197 66.6934646
sample estimates:
odds ratio 
  1.350276 
```

Display the proportion of Cd sloths infected by `haemoplasma` and `bloodparasite`:
```
table_haemoplasma_bloodparasite_Cd <- table(data_Cd$haemoplasma, data_Cd$bloodparasite)
table_haemoplasma_bloodparasite_Cd
```

Results are:
```
bloodparasite   0  1
haemoplasma  0 10  5
             1 24 44
```

Display the proportion of Cd sloths infected by `haemoplasma` in `anaplasma`:
```
table_haemoplasma_anaplasma_Cd <- table(data_Cd$haemoplasma, data_Cd$anaplasma)
table_haemoplasma_anaplasma_Cd
```

Results are:
```
anaplasma        0  1
haemoplasma   0 12  3
              1 31 37
```

Display the proportion of Cd sloths infected by `haemoplasma` in wet and dry `season`:
```
table_haemoplasma_season_Cd <- table(data_Cd$haemoplasma, data_Cd$season)
table_haemoplasma_season_Cd
```

Results are:
```
     D  W
  0 14  1
  1 47 21
```

Tests for associations between `haemoplasma` and `season` in Cd:
```
fisher.test(table(data_Cd$haemoplasma, data_Cd$season))  
```

Results are:
```
Fisher's Exact Test for Count Data
data:  table(data_Cd$haemoplasma, data_Cd$season)
p-value = 0.06059
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
   0.8295181 276.5385457
sample estimates:
odds ratio 
  6.158366  
```

ATTENTION ESSAI A GARDER OU PAS

Fit a GLM to test whether `haemoplasma` infection prevalence is influenced by additive effects of `anaplasma` and `season` in Cd:
```
model_2b <- glm(haemoplasma ~ anaplasma + season, data = data_Cd, family = binomial)
```

Compare model_2b to the complete additive model (model_2a) using a likelihood ratio test:
```
anova(model_2b, model_2a, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: haemoplasma ~ anaplasma + season
Model 2: haemoplasma ~ sex + age + season + tick + bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        80     69.017                     
2        77     67.088  3    1.929   0.5873
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_2a, model_2b)
```

Results are:
```
         df      AIC
model_2a  6 79.08835
model_2b  3 75.01733
```

Perform drop-one-term analysis on the model_2b additive model:
```
res <- drop1(model_2b, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
haemoplasma ~ anaplasma + season
          Df Deviance    AIC    LRT Pr(>Chi)  
<none>         69.017 75.017                  
anaplasma  1   73.855 77.855 4.8377  0.02784 *
season     1   72.229 76.229 3.2117  0.07311 .
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_2b)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
             AIC delta_AIC
<none>    75.017    0.0000
anaplasma 77.855    2.8377
season    76.229    1.2117
```

Compare the the model_2b additive model (haemoplasma ~ anaplasma + season) to the model2_anaplasma univariate model (haemoplasma ~ anaplasma) using likelihood ratio tests and AIC:
```
model2_anaplasma <- glm(haemoplasma ~ anaplasma, data = data_Cd, family = binomial)
anova(model_2b, model2_anaplasma, test="Chisq")
aics <- AIC(model_2b, model2_anaplasma)
aic_null <- aics["model_2b", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: haemoplasma ~ anaplasma + season
Model 2: haemoplasma ~ anaplasma
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        80     69.017                       
2        81     72.229 -1  -3.2117  0.07311 .
---
                      AIC delta_AIC_vs_null
model_2b         75.01733          0.000000
model2_anaplasma 76.22900          1.211667
```

Compare the the model2_anaplasma univariate model (haemoplasma ~ anaplasma) to the null model using likelihood ratio tests and AIC:
```
anova(model2_anaplasma, model2_null, test="Chisq")
aics <- AIC(model2_anaplasma, model2_null)
aic_null <- aics["model2_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: haemoplasma ~ anaplasma
Model 2: haemoplasma ~ 1
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        81     72.229                       
2        82     78.433 -1   -6.204  0.01275 *
---
                      AIC delta_AIC_vs_null
model2_anaplasma 76.22900         -4.203986
model2_null      80.43299          0.000000
```

FIN DE L ESSAI A RETIRER OU PAS

## Step 6. Impact of haemaplasma infections on Scale Mass Index (SMI) (GLM models 3) in Bt
The Scaled Mass Index (SMI) was used as a body condition indicator that standardizes individual `weight` to `body_length`, using an allometric scaling relationship. SMI was calculated following Peig & Green (2009) (https://doi.org/10.1111/j.1600-0706.2009.17643.x).

Function to calculate SMI for adult Bt:
```
data_adult_Bt <- subset(data_Bt, age == "A")
sma_model_Bt <- sma(log(weight) ~ log(total_length), data = data_adult_Bt)
b <- coef(sma_model_Bt)[2]
L0 <- mean(data_adult_Bt$total_length, na.rm = TRUE)
data_adult_Bt$SMI <- data_adult_Bt$weight * (L0 / data_adult_Bt$total_length)^b
```

Fit a GLM to test whether SMI is influenced by interactions among `haemoplasma`, `bloodparasite`, `sex` and `season` in Bt:
```
model_3 <- glm(SMI ~ haemoplasma * bloodparasite * season * sex, data = data_adult_Bt, family = gaussian(link = "identity"))
```

Fit a GLM to test whether SMI is influenced by additive effects of `haemoplasma`, `bloodparasite`, `sex` and `season` in Bt:
```
model_3a <- glm(SMI ~ haemoplasma + bloodparasite + season + sex, data = data_adult_Bt, family = gaussian(link = "identity"))
```

Compare the additive model (model_3a) to the interaction model (model_3) using a likelihood ratio test:
```
anova(model_3a, model_3, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: SMI ~ haemoplasma + bloodparasite + season + sex
Model 2: SMI ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        78     23.641                     
2        73     22.818  5  0.82298   0.7564
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_3, model_3a)
```

Results are:
```
         df      AIC
model_3  11 150.3660
model_3a  6 143.3069
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_3a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
SMI ~ haemoplasma + bloodparasite + season + sex
              Df Deviance    AIC scaled dev.  Pr(>Chi)    
<none>             23.641 143.31                          
haemoplasma    1   26.804 151.73     10.4232 0.0012444 ** 
bloodparasite  1   23.642 141.31      0.0045 0.9463824    
season         1   23.655 141.36      0.0515 0.8205475    
sex            1   27.489 153.82     12.5176 0.0004031 ***
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_3a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        143.31    0.0000
haemoplasma   151.73    8.4232
bloodparasite 141.31   -1.9955
season        141.36   -1.9485
sex           153.82   10.5176
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model3_null <- glm(SMI ~ 1, data = data_adult_Bt, family = gaussian(link = "identity"))
model3_haemoplasma <- glm(SMI ~ haemoplasma, data = data_adult_Bt, family = gaussian(link = "identity"))
model3_bloodparasite <- glm(SMI ~ bloodparasite, data = data_adult_Bt, family = gaussian(link = "identity"))
model3_season <- glm(SMI ~ season, data = data_adult_Bt, family = gaussian(link = "identity"))
model3_sex <- glm(SMI ~ sex, data = data_adult_Bt, family = gaussian(link = "identity"))
anova(model3_null, model3_haemoplasma, test="Chisq")
anova(model3_null, model3_bloodparasite, test="Chisq")
anova(model3_null, model3_season, test="Chisq")
anova(model3_null, model3_sex, test="Chisq")
aics <- AIC(model3_null, model3_haemoplasma, model3_bloodparasite, model3_season, model3_sex)
aic_null <- aics["model3_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: SMI ~ 1
Model 2: SMI ~ haemoplasma
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)   
1        82     30.838                        
2        81     27.594  1   3.2438 0.002031 **
---
Analysis of Deviance Table
Model 1: SMI ~ 1
Model 2: SMI ~ bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        82     30.838                     
2        81     30.564  1  0.27455   0.3937
---
Analysis of Deviance Table
Model 1: SMI ~ 1
Model 2: SMI ~ season
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        82     30.838                      
2        81     30.829  1 0.0087678   0.8794
---
Analysis of Deviance Table
Model 1: SMI ~ 1
Model 2: SMI ~ sex
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        82     30.838                          
2        81     26.881  1   3.9571 0.0005542 ***
---
                          AIC delta_AIC_vs_null
model3_null          157.3666          0.000000
model3_haemoplasma   150.1419         -7.224666
model3_bloodparasite 158.6243          1.257743
model3_season        159.3430          1.976398
model3_sex           147.9681         -9.398467
```

Fit a GLM to test whether SMI is influenced by additive effects of `haemoplasma` and `sex` in Bt:
```
model_3b <- glm(SMI ~ haemoplasma + sex, data = data_adult_Bt, family = gaussian(link = "identity"))
```

Compare the additive model (model_3b) to the full additive model (model_3a) using a likelihood ratio test:
```
anova(model_3b, model_3a, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: SMI ~ haemoplasma + sex
Model 2: SMI ~ haemoplasma + bloodparasite + season + sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        80     23.657                     
2        78     23.641  2 0.015892   0.9741
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_3a, model_3b)
```

Results are:
```
         df      AIC
model_3a  6 143.3069
model_3b  4 139.3626
```

Compare the the model3_haemoplasma univariate model (SMI ~ haemoplasma) to the the 'haemoplasma' + 'sex' additive model (model_3b) using likelihood ratio tests and AIC:
```
anova(model_3b, model3_haemoplasma, test="Chisq")
aics <- AIC(model_3b, model3_haemoplasma)
aic_haemoplasma <- aics["model3_haemoplasma", "AIC"]
aics$delta_AIC_vs_haemoplasma <- aics$AIC - aic_haemoplasma
print(aics[, c("AIC", "delta_AIC_vs_haemoplasma")])
```

Results are:
```
Analysis of Deviance Table
Model 1: SMI ~ haemoplasma + sex
Model 2: SMI ~ haemoplasma
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        80     23.657                          
2        81     27.594 -1  -3.9377 0.0002631 ***
---
                        AIC delta_AIC_vs_haemoplasma
model_3b           139.3626                -10.77931
model3_haemoplasma 150.1419                  0.00000
```

Compare the the model3_sex univariate model (SMI ~ sex) to the 'haemoplasma' + 'sex' additive model (model_3b) using likelihood ratio tests and AIC:
```
anova(model_3b, model3_sex, test="Chisq")
aics <- AIC(model_3b, model3_sex)
aic_sex <- aics["model3_sex", "AIC"]
aics$delta_AIC_vs_sex <- aics$AIC - aic_sex
print(aics[, c("AIC", "delta_AIC_vs_sex")])
```

Results are:
```
Analysis of Deviance Table
Model 1: SMI ~ haemoplasma + sex
Model 2: SMI ~ sex
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        80     23.657                          
2        81     26.881 -1  -3.2244 0.0009596 ***
---
                AIC delta_AIC_vs_sex
model_3b   139.3626        -8.605512
model3_sex 147.9681         0.000000
```

Fit a GLM to test whether SMI is influenced by interaction effect of `haemoplasma` and `sex` in Bt:
```
model_3c <- glm(SMI ~ haemoplasma * sex, data = data_adult_Bt, family = gaussian(link = "identity"))
```

Compare the additive model (model_3b) to the interactive model (model_3c) using a likelihood ratio test:
```
anova(model_3c, model_3b, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: SMI ~ haemoplasma * sex
Model 2: SMI ~ haemoplasma + sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        79     23.647                     
2        80     23.657 -1 -0.00961   0.8578
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_3b, model_3c)
```

Results are:
```
         df      AIC
model_3b  4 139.3626
model_3c  5 141.3289
```

Assess residual normality and heteroscedasticity of the 'haemoplasma' + 'sex' additive model (model_3b):
```
shapiro.test(model_3b$residuals)
bptest(model_3b)
```

Results are:
```
Shapiro-Wilk normality test
data:  model_3b$residuals
W = 0.9883, p-value = 0.6612
---
studentized Breusch-Pagan test
data:  model_3b
BP = 4.0982, df = 2, p-value = 0.1288
```

Calculation of mean and standard error of SMI by `haemoplasma` infection status and `sex` for Bt:
```
data_adult_Bt %>% 
  group_by(sex, haemoplasma) %>% 
  summarise(
    n = sum(!is.na(SMI)),
    mean = mean(SMI, na.rm = TRUE),
    se = sd(SMI, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) %>% 
  mutate(SMI = sprintf("%.2f ± %.2f", mean, se))
```

Results are:
```
  sex   haemoplasma     n  mean     se SMI        
  <fct> <fct>       <int> <dbl>  <dbl> <chr>      
1 F     0              39  4.43 0.0749 4.43 ± 0.07
2 F     1               2  3.56 0.0284 3.56 ± 0.03
3 M     0              40  4.88 0.0985 4.88 ± 0.10
4 M     1               2  3.90 0.298  3.90 ± 0.30
```

Generate SMI chart for Bt:
```
clean_data <- data_adult_Bt %>%
  filter(
    !is.na(weight), !is.na(total_length), !is.na(SMI),
    is.finite(weight), is.finite(total_length), is.finite(SMI)
  ) %>%
  mutate(
    sex_infect = case_when(
      sex == "M" & haemoplasma == 0 ~ "Male, uninfected",
      sex == "M" & haemoplasma == 1 ~ "Male, infected",
      sex == "F" & haemoplasma == 0 ~ "Female, uninfected",
      sex == "F" & haemoplasma == 1 ~ "Female, infected",
      TRUE ~ NA_character_
    )
  )
levels_order <- c("Male, uninfected", "Male, infected", "Female, uninfected", "Female, infected")
clean_data <- clean_data %>%
  mutate(
    sex_infect = factor(sex_infect, levels = levels_order),
    point_size = case_when(
      sex_infect %in% c("Male, uninfected", "Male, infected") ~ 3.25,
      TRUE ~ 4  # taille normale pour les cercles
    )
  )
interp_data <- with(clean_data, akima::interp(
  x = weight,
  y = total_length,
  z = SMI,
  duplicate = "mean",
  extrap = FALSE
))
interp_df <- expand.grid(
  x = interp_data$x,
  y = interp_data$y
)
interp_df$z <- as.vector(interp_data$z)
legend_point_sizes <- c(3.25, 3.25, 4, 4) / 2
ggplot() +
  geom_contour_filled(data = interp_df, aes(x = x, y = y, z = z)) +
  geom_point(
    data = clean_data,
    aes(
      x = weight,
      y = total_length,
      shape = sex_infect,
      size = point_size
    ),
    color = "black",
    stroke = 1
  ) +
  scale_fill_brewer(palette = "Green", name = "SMI level") +
  scale_shape_manual(
    name = expression(paste("Haemoplasma", " infection status")),
    values = c(
      "Male, uninfected" = 0,
      "Male, infected" = 12,
      "Female, uninfected" = 1,
      "Female, infected" = 10
    )
  ) +
  scale_size_identity(guide = "none") + 
  guides(
    shape = guide_legend(override.aes = list(size = legend_point_sizes))
  ) +
  labs(
    x = "Body mass (kg)",
    y = "Total Length (cm)",
    title = expression(paste("Scale Mass Index (SMI) of ", italic("Bradypus tridactylus")))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
```

## Step 7. Impact of haemaplasma infections on Scale Mass Index (SMI) (GLM models 4) in Cd
Function to calculate SMI for adult Cd:
```
data_adult_Cd <- subset(data_Cd, age == "A")
sma_model_Cd <- sma(log(weight) ~ log(total_length), data = data_adult_Cd)
b <- coef(sma_model_Cd)[2]
L0 <- mean(data_adult_Cd$total_length, na.rm = TRUE)
data_adult_Cd$SMI <- data_adult_Cd$weight * (L0 / data_adult_Cd$total_length)^b
```

Fit a GLM to test whether SMI is influenced by interactions among `haemoplasma`, `bloodparasite`, `sex`, and `season` in Cd:
```
model_4 <- glm(SMI ~ haemoplasma * bloodparasite * season * sex, data = data_adult_Cd, family = gaussian(link = "identity"))
```

Fit a GLM to test whether SMI is influenced by additive effects of `haemoplasma`, `bloodparasite`, `sex`, and `season` in Cd:
```
model_4a <- glm(SMI ~ haemoplasma + bloodparasite + season + sex, data = data_adult_Cd, family = gaussian(link = "identity"))
```

Compare the additive model (model_4a) to the interaction model (model_4) using a likelihood ratio test:
```
anova(model_4a, model_4, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: SMI ~ haemoplasma + bloodparasite + season + sex
Model 2: SMI ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        52     27.506                     
2        45     25.399  7    2.106   0.8102
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_4, model_4a)
```

Results are:
```
         df      AIC
model_4  13 141.6847
model_4a  6 132.2251
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_4a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
SMI ~ haemoplasma + bloodparasite + season + sex
              Df Deviance    AIC scaled dev. Pr(>Chi)  
<none>             27.506 132.22                       
haemoplasma    1   27.514 130.24      0.0178  0.89379  
bloodparasite  1   30.401 135.93      5.7058  0.01691 *
season         1   27.725 130.68      0.4539  0.50048  
sex            1   27.625 130.47      0.2465  0.61958  
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_4a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        132.22    0.0000
haemoplasma   130.24   -1.9822
bloodparasite 135.93    3.7058
season        130.68   -1.5461
sex           130.47   -1.7535
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model4_null <- glm(SMI ~ 1, data = data_adult_Cd, family = gaussian(link = "identity"))
model4_haemoplasma <- glm(SMI ~ haemoplasma, data = data_adult_Cd, family = gaussian(link = "identity"))
model4_bloodparasite <- glm(SMI ~ bloodparasite, data = data_adult_Cd, family = gaussian(link = "identity"))
model4_season <- glm(SMI ~ season, data = data_adult_Cd, family = gaussian(link = "identity"))
model4_sex <- glm(SMI ~ sex, data = data_adult_Cd, family = gaussian(link = "identity"))
anova(model4_null, model4_haemoplasma, test="Chisq")
anova(model4_null, model4_bloodparasite, test="Chisq")
anova(model4_null, model4_season, test="Chisq")
anova(model4_null, model4_sex, test="Chisq")
aics <- AIC(model4_null, model4_haemoplasma, model4_bloodparasite, model4_season, model4_sex)
aic_null <- aics["model4_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: SMI ~ 1
Model 2: SMI ~ haemoplasma
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        56     30.572                      
2        55     30.571  1 0.0004533   0.9772
---
Analysis of Deviance Table
Model 1: SMI ~ 1
Model 2: SMI ~ bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        56     30.572                       
2        55     27.933  1   2.6388  0.02264 *
---
Analysis of Deviance Table
Model 1: SMI ~ 1
Model 2: SMI ~ season
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        56     30.572                     
2        55     30.422  1  0.15022   0.6023
---
Analysis of Deviance Table
Model 1: SMI ~ 1
Model 2: SMI ~ sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        56     30.572                     
2        55     30.554  1 0.017332   0.8598
---
AIC delta_AIC_vs_null
model4_null          130.2494          0.000000
model4_haemoplasma   132.2486          1.999155
model4_bloodparasite 127.1041         -3.145367
model4_season        131.9687          1.719234
model4_sex           132.2171          1.967676
```

Assess residual normality and heteroscedasticity:
```
shapiro.test(model4_bloodparasite$residuals)
bptest(model4_bloodparasite)
```

Results are:
```
Shapiro-Wilk normality test
data:  model4_bloodparasite$residuals
W = 0.97176, p-value = 0.2024
---
studentized Breusch-Pagan test
data:  model4_bloodparasite
BP = 0.37341, df = 1, p-value = 0.5412
```

Post hoc power analyses for SMI tests in Cd (for full interaction model):
```
n <- nrow(na.omit(data_adult_Cd[, c("SMI", "haemoplasma", "bloodparasite", "season", "sex")]))
k <- 15
pwr.f2.test(u = k, v = n - k - 1, f2 = 0.30, sig.level = 0.05)
pwr.f2.test(u = k, v = n - k - 1, f2 = 0.20, sig.level = 0.05)
```

Results are:
```
Multiple regression power calculation 
u = 15
v = 41
f2 = 0.3
sig.level = 0.05
power = 0.5978527
---
Multiple regression power calculation 
u = 15
v = 41
f2 = 0.2
sig.level = 0.05
power = 0.3970745
```

Post hoc power analyses for SMI tests in Cd (for `SMI` ~ `bloodparasite` and adding `haemoplasma`):
```
n <- nrow(na.omit(data_adult_Cd[, c("SMI", "bloodparasite" "haemoplasma")]))
k <- 3
pwr.f2.test(u = k, v = n - k - 1, f2 = 0.30, sig.level = 0.05)
pwr.f2.test(u = k, v = n - k - 1, f2 = 0.20, sig.level = 0.05)
```

Results are:
```
Multiple regression power calculation 
u = 3
v = 53
f2 = 0.3
sig.level = 0.05
power = 0.9321324
---
Multiple regression power calculation 
u = 3
v = 53
f2 = 0.2
sig.level = 0.05
power = 0.787044
```

AUTRE ESSAI BLOODPARASITE EN 1 PAR 1 A VOIR SI GARDER OU PAS

Fit a GLM to test whether SMI is influenced by interactions among `haemoplasma`, `anaplasma`, `microfilaria`, `trypanosome`, `babesia`, `sex` and `season` in Cd:
```
model_4bis <- glm(SMI ~ haemoplasma * anaplasma * microfilaria * trypanosome * babesia * season * sex, data = data_adult_Cd, family = gaussian(link = "identity"))
```

Fit a GLM to test whether SMI is influenced by additive effects of `haemoplasma`, `anaplasma`, `microfilaria`, `trypanosome`, `babesia`, `sex` and `season` in Cd:
```
model_4bisa <- glm(SMI ~ haemoplasma + anaplasma + microfilaria + trypanosome + babesia + season + sex, data = data_adult_Cd, family = gaussian(link = "identity"))
```

Compare the additive model (model_4bisa) to the interaction model (model_4bis) using a likelihood ratio test:
```
anova(model_4bisa, model_4bis, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: SMI ~ haemoplasma + anaplasma + microfilaria + trypanosome + babesia + season + sex
Model 2: SMI ~ haemoplasma * anaplasma * microfilaria * trypanosome * babesia * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        49     28.333                     
2        34     20.363 15   7.9702   0.5786
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_4bis, model_4bisa)
```

Results are:
```
            df      AIC
model_4bis  24 151.0877
model_4bisa  9 139.9153
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_4bisa, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
SMI ~ haemoplasma + anaplasma + microfilaria + trypanosome + babesia + season + sex
             Df Deviance    AIC scaled dev. Pr(>Chi)
<none>            28.333 139.91                     
haemoplasma   1   28.335 137.92     0.00399   0.9496
anaplasma     1   28.885 139.01     1.09981   0.2943
microfilaria  1   29.568 140.35     2.43203   0.1189
trypanosome   1   28.452 138.15     0.23902   0.6249
babesia       1   28.431 138.11     0.19698   0.6572
season        1   28.647 138.54     0.62787   0.4281
sex           1   28.430 138.11     0.19435   0.6593
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_4bisa)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                AIC delta_AIC
<none>       139.91   0.00000
haemoplasma  137.92  -1.99601
anaplasma    139.01  -0.90019
microfilaria 140.35   0.43203
trypanosome  138.15  -1.76098
babesia      138.11  -1.80302
season       138.54  -1.37213
sex          138.11  -1.80565
```

FIN ESSAI 4bis MAIS IL FAUDRA SANS DOUTE CREER UNE VARIABLE  HAEMOPLASMA + BLOODPARASITE :(

Generate SMI chart for Cd:
```
clean_data <- data_adult_Cd %>%
  filter(
    !is.na(weight), !is.na(total_length), !is.na(SMI),
    is.finite(weight), is.finite(total_length), is.finite(SMI)
  ) %>%
  mutate(
    sex_infect = case_when(
      sex == "M" & haemoplasma == 0 ~ "Male, uninfected",
      sex == "M" & haemoplasma == 1 ~ "Male, infected",
      sex == "F" & haemoplasma == 0 ~ "Female, uninfected",
      sex == "F" & haemoplasma == 1 ~ "Female, infected",
      TRUE ~ NA_character_
    )
  )
levels_order <- c("Male, uninfected", "Male, infected", "Female, uninfected", "Female, infected")
clean_data <- clean_data %>%
  mutate(
    sex_infect = factor(sex_infect, levels = levels_order),
    point_size = case_when(
      sex_infect %in% c("Male, uninfected", "Male, infected") ~ 3.25,
      TRUE ~ 4  # taille normale pour les cercles
    )
  )
interp_data <- with(clean_data, akima::interp(
  x = weight,
  y = total_length,
  z = SMI,
  duplicate = "mean",
  extrap = FALSE
))
interp_df <- expand.grid(
  x = interp_data$x,
  y = interp_data$y
)
interp_df$z <- as.vector(interp_data$z)
legend_point_sizes <- c(3.25, 3.25, 4, 4) / 2
ggplot() +
  geom_contour_filled(data = interp_df, aes(x = x, y = y, z = z)) +
  geom_point(
    data = clean_data,
    aes(
      x = weight,
      y = total_length,
      shape = sex_infect,
      size = point_size
    ),
    color = "black",
    stroke = 1
  ) +
  scale_fill_brewer(palette = "YlOrBr", name = "SMI level") +
  scale_shape_manual(
    name = expression(paste("Haemoplasma", " infection status")),
    values = c(
      "Male, uninfected" = 0,
      "Male, infected" = 12,
      "Female, uninfected" = 1,
      "Female, infected" = 10
    )
  ) +
  scale_size_identity(guide = "none") + 
  guides(
    shape = guide_legend(override.aes = list(size = legend_point_sizes))
  ) +
  labs(
    x = "Body mass (kg)",
    y = "Total Length (cm)",
    title = expression(paste("Scale Mass Index (SMI) of ", italic("Choloepus didactylus")))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
```

## Step 8. Impact of heamoplasma infections on neck circumference (GLM models 5 and 6)

Fit a GLM to test whether neck circumference is influenced by interactions among `haemoplasma`, `bloodparasite`, `sex`, and `season` in Bt:
```
model_5 <- glm(log(neck_size) ~ haemoplasma * bloodparasite * season * sex, data = data_adult_Bt, family = gaussian(link = "identity"))
```

Fit a GLM to test whether SMI is influenced by additive effects of `anaplasma`, `sex`, and `season` in Bt:
```
model_5a <- glm(log(neck_size) ~ haemoplasma + bloodparasite + season + sex, data = data_adult_Bt, family = gaussian(link = "identity"))
```

Compare the additive model (model_5a) to the interaction model (model_5) using a likelihood ratio test:
```
anova(model_5a, model_5, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: log(neck_size) ~ haemoplasma + bloodparasite + season + sex
Model 2: log(neck_size) ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        56    0.47867                     
2        52    0.45529  4 0.023386   0.6143
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_5, model_5a)
```

Results are:
```
         df       AIC
model_5  10 -105.6494
model_5a  6 -110.5939
```


Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_5a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
log(neck_size) ~ haemoplasma + bloodparasite + season + sex
              Df Deviance     AIC scaled dev. Pr(>Chi)
<none>            0.47867 -110.59                     
haemoplasma    1  0.48450 -111.86     0.73846   0.3902
bloodparasite  1  0.48917 -111.27     1.32339   0.2500
season         1  0.48060 -112.35     0.24482   0.6207
sex            1  0.49638 -110.38     2.21631   0.1366
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_5a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                  AIC delta_AIC
<none>        -110.59   0.00000
haemoplasma   -111.86  -1.26154
bloodparasite -111.27  -0.67661
season        -112.35  -1.75518
sex           -110.38   0.21631
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model5_null <- glm(log(neck_size) ~ 1, data = data_adult_Bt, family = gaussian(link = "identity"))
model5_haemoplasma <- glm(log(neck_size) ~ haemoplasma, data = data_adult_Bt, family = gaussian(link = "identity"))
model5_bloodparasite <- glm(log(neck_size) ~ bloodparasite, data = data_adult_Bt, family = gaussian(link = "identity"))
model5_season <- glm(log(neck_size) ~ season, data = data_adult_Bt, family = gaussian(link = "identity"))
model5_sex <- glm(log(neck_size) ~ sex, data = data_adult_Bt, family = gaussian(link = "identity"))
anova(model5_null, model5_haemoplasma, test="Chisq")
anova(model5_null, model5_bloodparasite, test="Chisq")
anova(model5_null, model5_season, test="Chisq")
anova(model5_null, model5_sex, test="Chisq")
aics <- AIC(model5_null, model5_haemoplasma, model5_bloodparasite, model5_season, model5_sex)
aic_null <- aics["model5_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ haemoplasma
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        60    0.52778                     
2        59    0.51330  1 0.014479    0.197
---
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ bloodparasite
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        60    0.52778                     
2        59    0.51039  1 0.017391   0.1562
---
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ season
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        60    0.52778                      
2        59    0.52413  1 0.0036529   0.5214
---
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        60    0.52778                       
2        59    0.49900  1 0.028786  0.06505 .
---
model5_null          -112.6360        0.00000000
model5_haemoplasma   -112.3329        0.30312861
model5_bloodparasite -112.6798       -0.04384486
model5_season        -111.0596        1.57633716
model5_sex           -114.0572       -1.42123531
```

Fit a linear model to test the null hypothesis (`neck_size` ~ 1) in adult Bt, assessing model fit and checking residual normality:
```
model_5b <- glm(log(neck_size) ~ 1, data = data_adult_Bt, family = gaussian(link = "identity"))
anova(model_5b, model_5, test = "Chisq")
AIC(model_5b, model_5)
shapiro.test(model_5b$residuals)
```

Results are:
```
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        60    0.52778                     
2        52    0.45529  8 0.072498   0.4066
---
         df       AIC
model_5b  2 -112.6360
model_5  10 -105.6494
---
Shapiro-Wilk normality test
data:  model_5b$residuals
W = 0.96651, p-value = 0.09328
```

Fit a GLM to test whether neck circumference is influenced by interactions among `haemoplasma`, `bloodparasite`, `sex`, and `season` in Cd:
```
model_6 <- glm(log(neck_size) ~ haemoplasma * bloodparasite * season * sex, data = data_adult_Cd, family = gaussian(link = "identity"))
```

Fit a GLM to test whether SMI is influenced by additive effects of `haemoplasma`, `bloodparasite`, `sex`, and `season` in Cd:
```
model_6a <- glm(log(neck_size) ~ haemoplasma + bloodparasite + season + sex, data = data_adult_Cd, family = gaussian(link = "identity"))
```

Compare the additive model (model_6a) to the interaction model (model_6) using a likelihood ratio test:
```
anova(model_6a, model_6, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: log(neck_size) ~ haemoplasma + bloodparasite + season + sex
Model 2: log(neck_size) ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        43    0.55930                     
2        36    0.42106  7  0.13824   0.1067
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_6, model_6a)
```

Results are:
```
         df       AIC
model_6  13 -65.11804
model_6a  6 -65.49074
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_6a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
log(neck_size) ~ haemoplasma + bloodparasite + season + sex
              Df Deviance     AIC scaled dev. Pr(>Chi)
<none>            0.55930 -65.491                     
haemoplasma    1  0.55935 -67.486     0.00426   0.9480
bloodparasite  1  0.56108 -67.338     0.15253   0.6961
season         1  0.57932 -65.802     1.68826   0.1938
sex            1  0.57644 -66.042     1.44873   0.2287
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_6a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                  AIC delta_AIC
<none>        -65.491   0.00000
haemoplasma   -67.486  -1.99574
bloodparasite -67.338  -1.84747
season        -65.802  -0.31174
sex           -66.042  -0.55127
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model6_null <- glm(log(neck_size) ~ 1, data = data_adult_Cd, family = gaussian(link = "identity"))
model6_haemoplasma <- glm(log(neck_size) ~ haemoplasma, data = data_adult_Cd, family = gaussian(link = "identity"))
model6_bloodparasite <- glm(log(neck_size) ~ bloodparasite, data = data_adult_Cd, family = gaussian(link = "identity"))
model6_season <- glm(log(neck_size) ~ season, data = data_adult_Cd, family = gaussian(link = "identity"))
model6_sex <- glm(log(neck_size) ~ sex, data = data_adult_Cd, family = gaussian(link = "identity"))
anova(model6_null, model6_haemoplasma, test="Chisq")
anova(model6_null, model6_bloodparasite, test="Chisq")
anova(model6_null, model6_season, test="Chisq")
anova(model6_null, model6_sex, test="Chisq")
aics <- AIC(model6_null, model6_haemoplasma, model6_bloodparasite, model6_season, model6_sex)
aic_null <- aics["model6_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ haemoplasma
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        47    0.60752                      
2        46    0.60512  1 0.0024025   0.6691
---
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ bloodparasite
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        47    0.60752                      
2        46    0.60334  1 0.0041827   0.5723
---
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ season
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        47    0.60752                     
2        46    0.58020  1 0.027321   0.1411
---
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        47    0.60752                     
2        46    0.58267  1 0.024851   0.1613
---
                           AIC delta_AIC_vs_null
model6_null          -69.52144       0.000000000
model6_haemoplasma   -67.71164       1.809799276
model6_bloodparasite -67.85305       1.668385770
model6_season        -69.73015      -0.208715492
model6_sex           -69.52621      -0.004777113
```

Fit a linear model to test the null hypothesis (`neck_size` ~ 1) in adult Cd, assessing model fit and checking residual normality:
```
model_6b <- glm(log(neck_size) ~ 1, data = data_adult_Cd, family = gaussian(link = "identity"))
anova(model_6b, model_6, test = "Chisq")
AIC(model_6b, model_6)
shapiro.test(model_6b$residuals)
```

Results are:
```
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        47    0.60752                     
2        36    0.42106 11  0.18645   0.1433
---
         df       AIC
model_6b  2 -69.52144
model_6  13 -65.11804
---
Shapiro-Wilk normality test
data:  model_6b$residuals
W = 0.97242, p-value = 0.3137
```

## Step 9. Impact of _Anaplasma_ infections on hematocrit levels (GLM models 7, 8 and 9)
Fit a GLM to test whether `hematocrit` is influenced by interactions among `haemoplasma`, `bloodparasite`, `sex`, and `season` in Bt:
```
model_7 <- glm(hematocrit ~ haemoplasma * bloodparasite * season * sex, data = data_adult_Bt, family = Gamma(link = "log"))
```

Fit a GLM to test whether `hematocrit` is influenced by additive effects of haemoplasma`, `bloodparasite`, `sex`, and `season` in Bt:
```
model_7a <- glm(hematocrit ~ haemoplasma + bloodparasite + season + sex, data = data_adult_Bt, family = Gamma(link = "log"))
```

Compare the additive model (model_7a) to the interaction model (model_7) using a likelihood ratio test:
```
anova(model_7a, model_7, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ haemoplasma + bloodparasite + season + sex
Model 2: hematocrit ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        79     1.3166                     
2        74     1.1985  5  0.11806   0.2065
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_7, model_7a)
```

Results are:
```
         df      AIC
model_7  11 517.8352
model_7a  6 515.7467
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_7a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
hematocrit ~ haemoplasma + bloodparasite + season + sex
              Df Deviance    AIC scaled dev. Pr(>Chi)  
<none>             1.3166 515.75                       
haemoplasma    1   1.3209 514.00      0.2550  0.61358  
bloodparasite  1   1.3258 514.29      0.5474  0.45937  
season         1   1.3238 514.18      0.4285  0.51273  
sex            1   1.3825 517.67      3.9191  0.04774 *
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_7a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        515.75    0.0000
haemoplasma   514.00   -1.7450
bloodparasite 514.29   -1.4526
season        514.18   -1.5715
sex           517.67    1.9191
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model7_null <- glm(hematocrit ~ 1, data = data_adult_Bt, family = Gamma(link = "log"))
model7_haemoplasma <- glm(hematocrit ~ haemoplasma, data = data_adult_Bt, family = Gamma(link = "log"))
model7_bloodparasite <- glm(hematocrit ~ bloodparasite, data = data_adult_Bt, family = Gamma(link = "log"))
model7_season <- glm(hematocrit ~ season, data = data_adult_Bt, family = Gamma(link = "log"))
model7_sex <- glm(hematocrit ~ sex, data = data_adult_Bt, family = Gamma(link = "log"))
anova(model7_null, model7_haemoplasma, test="Chisq")
anova(model7_null, model7_bloodparasite, test="Chisq")
anova(model7_null, model7_season, test="Chisq")
anova(model7_null, model7_sex, test="Chisq")
aics <- AIC(model7_null, model7_haemoplasma, model7_bloodparasite, model7_season, model7_sex)
aic_null <- aics["model7_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ haemoplasma
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        83     1.4015                     
2        82     1.3936  1 0.007882   0.4994
---
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ bloodparasite
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        83     1.4015                      
2        82     1.3950  1 0.0064542   0.5418
---
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ season
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        83     1.4015                      
2        82     1.3926  1 0.0089277   0.4709
---
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        83     1.4015                       
2        82     1.3423  1 0.059141  0.05883 .
---
                          AIC delta_AIC_vs_null
model7_null          513.0105          0.000000
model7_haemoplasma   514.5354          1.524936
model7_bloodparasite 514.6217          1.611190
model7_season        514.4722          1.461711
model7_sex           511.3790         -1.631481
```

Fit a linear model to test the null hypothesis (`hematocrit` ~ 1) in adult Bt, assessing model fit:
```
model_7b <- glm(hematocrit ~ 1, data = data_adult_Bt, family = Gamma(link = "log"))
anova(model_7b, model_7, test = "Chisq")
AIC(model_7b, model_7)
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        83     1.4015                     
2        74     1.1985  9  0.20297   0.1932
---
         df      AIC
model_7b  2 513.0105
model_7  11 517.8352
```

Generate QQ plot for model_7b to assess normality of residuals:
```
res_dev <- residuals(model_7b, type = "deviance")
qqnorm(res_dev, main = "QQ Plot (Deviance residuals)")
qqline(res_dev, col = "red")
```

Calculation of mean and standard error of `hematocrit` by `haemoplasma` for Bt:
```
data_adult_Bt %>%
  group_by(haemoplasma) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results are:
```
A tibble: 2 × 3
  haemoplasma mean_hematocrit se_hematocrit
  <fct>                 <dbl>         <dbl>
1 0                      39.0         0.569
2 1                      40.8         3.04 
```

Fit a GLM to test whether `hematocrit` is influenced by interactions among `haemoplasma`, `bloodparasite`, `sex`, and `season` in Cd:
```
model_8 <- glm(hematocrit ~ haemoplasma * bloodparasite * season * sex, data = data_adult_Cd, family = Gamma(link = "log"))
```

Fit a GLM to test whether `hematocrit` is influenced by additive effects of `haemoplasma`, `bloodparasite`, `sex`, and `season` in Cd:
```
model_8a <- glm(hematocrit ~ haemoplasma + bloodparasite + season + sex, data = data_adult_Cd, family = Gamma(link = "log"))
```

Compare the additive model (model_8a) to the interaction model (model_8) using a likelihood ratio test:
```
anova(model_8a, model_8, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ haemoplasma + bloodparasite + season + sex
Model 2: hematocrit ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        55     1.1832                     
2        49     1.1262  6 0.057072   0.8368
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_8, model_8a)
```

Results are:
```
         df      AIC
model_8  12 393.2948
model_8a  6 384.2705
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_8a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
hematocrit ~ haemoplasma + bloodparasite + season + sex
              Df Deviance    AIC scaled dev. Pr(>Chi)   
<none>             1.1832 384.27                        
haemoplasma    1   1.1883 382.54      0.2645 0.607020   
bloodparasite  1   1.2285 384.64      2.3705 0.123646   
season         1   1.3174 389.30      7.0307 0.008013 **
sex            1   1.2199 384.19      1.9221 0.165622   
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_8a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        384.27    0.0000
haemoplasma   382.54   -1.7355
bloodparasite 384.64    0.3705
season        389.30    5.0307
sex           384.19   -0.0779
```

Fit a linear model to test the model_8b (`hematocrit` ~ `season`) in adult Cd, assessing model fit:
```
model_8b <- glm(hematocrit ~ season, data = data_adult_Cd, family = Gamma(link = "log"))
anova(model_8b, model_8, test = "Chisq")
AIC(model_8b, model_8)
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ season
Model 2: hematocrit ~ haemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        58     1.2872                     
2        49     1.1262  9  0.16103    0.552
---
         df      AIC
model_8b  3 383.3404
model_8  12 393.2948
```

Generate QQ plot for model_8b to assess normality of residuals:
```
res_dev <- residuals(model_8b, type = "deviance")
qqnorm(res_dev, main = "QQ Plot (Deviance residuals)")
qqline(res_dev, col = "red")
```

Calculation of mean and standard error of `hematocrit` by `haemoplasma` for Cd:
```
data_adult_Cd %>%
  group_by(haemoplasma) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results are:
```
A tibble: 2 × 3
  haemoplasma mean_hematocrit se_hematocrit
  <fct>                 <dbl>         <dbl>
1 0                      39.1         1.48 
2 1                      38.6         0.798
```

Calculation of mean and standard error of `hematocrit` by `season` for Cd:
```
data_adult_Cd %>%
  group_by(season) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results are:
```
A tibble: 2 × 3
  season mean_hematocrit se_hematocrit
  <fct>            <dbl>         <dbl>
1 D                 37.6         0.892
2 W                 41.1         0.944
```

Create violin plots for `hematocrit`
```
label_style <- element_text(size = 28, face = "bold") 

pA <- ggplot(data_adult_Bt, aes(x = factor(haemoplasma, levels = c(0, 1),
                                labels = c("Uninfected", "Infected")),
                                y = hematocrit)) +
  geom_violin(fill = "darkolivegreen3", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = NULL,  
       y = "Hematocrit (%)",   
       title = "A") +
  scale_y_continuous(limits = c(20, 60)) +
  scale_x_discrete(labels = NULL) + 
  theme_minimal() +
  theme(plot.title = label_style)

pB <- ggplot(data_adult_Bt, aes(x = factor(season, levels = c("D", "W"),
                                labels = c("Dry", "Wet")),
                                y = hematocrit)) +
  geom_violin(fill = "darkolivegreen3", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = NULL,
       y = NULL,
       title = "B") +
  scale_y_continuous(limits = c(20, 60)) +
  scale_x_discrete(labels = NULL) +  
  theme_minimal() +
  theme(plot.title = label_style)

pC <- ggplot(data_adult_Bt, aes(x = factor(sex, levels = c("M", "F"),
                                labels = c("Male", "Female")),
                                y = hematocrit)) +
  geom_violin(fill = "darkolivegreen3", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = NULL,
       y = NULL,
       title = "C") +
  scale_y_continuous(limits = c(20, 60)) +
  scale_x_discrete(labels = NULL) +  
  theme_minimal() +
  theme(plot.title = label_style)

pD <- ggplot(data_adult_Cd, aes(x = factor(haemoplasma, levels = c(0, 1),
                                labels = c("Uninfected", "Infected")),
                                y = hematocrit)) +
  geom_violin(fill = "goldenrod1", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = expression(paste("Haemoplasma", " infection status")),
       y = "Hematocrit (%)",
       title = "D") +
  scale_y_continuous(limits = c(10, 60)) +
  theme_minimal() +
  theme(plot.title = label_style)

pE <- ggplot(data_adult_Cd, aes(x = factor(season, levels = c("D", "W"),
                                labels = c("Dry", "Wet")),
                                y = hematocrit)) +
  geom_violin(fill = "goldenrod1", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = "Season",
       y = NULL,
       title = "E") +
  scale_y_continuous(limits = c(10, 60)) +
  theme_minimal() +
  theme(plot.title = label_style)

pF <- ggplot(data_adult_Cd, aes(x = factor(sex, levels = c("M", "F"),
                                labels = c("Male", "Female")),
                                y = hematocrit)) +
  geom_violin(fill = "goldenrod1", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = "Sex",
       y = NULL,
       title = "F") +
  scale_y_continuous(limits = c(10, 60)) +
  theme_minimal() +
  theme(plot.title = label_style)
final_plot <- (pA | pB | pC) / (pD | pE | pF)
print(final_plot)
```

## Step 10. Impact of Haemoplasma infections on body temperature (CLRM models 10 and 11)
Convert `temperature` to numeric, handle left-censored values (<32°C) for analysis in Bt:
```
data_adult_Bt <- data_adult_Bt %>%
  mutate(
    temperature_numeric = as.numeric(ifelse(temperature == "< 32.00", 32, temperature)),
    censored = ifelse(temperature == "< 32.00", TRUE, FALSE)
  )
```

Create a left-censored Surv object (temp) for `temperature` in Bt:
```
temp <- Surv(data_adult_Bt$temperature_numeric,
                  event = !data_adult_Bt$censored,
                  type = "left")
```

Fit Gaussian survival regression models to test the effects of `haemoplasma`, `bloodparasite`, `season`, `sex` on `temperature`, with (model_10) and without (model_10b) interactions in Bt:
```
model_10 <- survreg(temp ~ haemoplasma * bloodparasite * season * sex, data = data_adult_Bt, dist = "gaussian")
model_10a <- survreg(temp ~ haemoplasma + bloodparasite + season + sex, data = data_adult_Bt, dist = "gaussian")

```

Compare models using ANOVA and AIC to evaluate the contribution of interaction terms:
```
anova(model_10a, model_10, test = "Chisq")
AIC(model_10a, model_10)
```

Results are:
```
                                       Terms Resid. Df    -2*LL Test Df  Deviance  Pr(>Chi)
1 haemoplasma + bloodparasite + season + sex        27 79.74336      NA        NA        NA
2 haemoplasma * bloodparasite * season * sex        16 79.05558    = 11 0.6877884 0.9999927
---
          df       AIC
model_10a  6  91.74336
model_10  17 113.05558
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_10a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
temp ~ haemoplasma + bloodparasite + season + sex
              Df    AIC    LRT Pr(>Chi)  
<none>           91.743                  
haemoplasma    1 91.269 1.5261   0.2167  
bloodparasite  1 89.952 0.2084   0.6480  
season         1 90.275 0.5320   0.4658  
sex            1 94.666 4.9228   0.0265 *
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_10a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        91.743   0.00000
haemoplasma   91.269  -0.47387
bloodparasite 89.952  -1.79159
season        90.275  -1.46797
sex           94.666   2.92284
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model10_null <- survreg(temp ~ 1, data = data_adult_Bt, dist = "gaussian")
model10_haemoplasma <- survreg(temp ~ haemoplasma, data = data_adult_Bt, dist = "gaussian")
model10_bloodparasite <- survreg(temp ~ bloodparasite, data = data_adult_Bt, dist = "gaussian")
model10_season <- survreg(temp ~ season, data = data_adult_Bt, dist = "gaussian")
model10_sex <- survreg(temp ~ sex, data = data_adult_Bt, dist = "gaussian")
anova(model10_null, model10_haemoplasma, test="Chisq")
anova(model10_null, model10_bloodparasite, test="Chisq")
anova(model10_null, model10_season, test="Chisq")
anova(model10_null, model10_sex, test="Chisq")
aics <- AIC(model10_null, model10_haemoplasma, model10_bloodparasite, model10_season, model10_sex)
aic_null <- aics["model10_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
        Terms Resid. Df    -2*LL Test Df  Deviance  Pr(>Chi)
1           1        31 85.36015      NA        NA        NA
2 haemoplasma        30 84.81430    =  1 0.5458493 0.4600186
---
          Terms Resid. Df    -2*LL Test Df   Deviance  Pr(>Chi)
1             1        31 85.36015      NA         NA        NA
2 bloodparasite        30 85.27469    =  1 0.08546037 0.7700297
---
   Terms Resid. Df    -2*LL Test Df    Deviance  Pr(>Chi)
1      1        31 85.36015      NA          NA        NA
2 season        30 85.35333    =  1 0.006827965 0.9341446
---
  Terms Resid. Df    -2*LL Test Df Deviance   Pr(>Chi)
1     1        31 85.36015      NA       NA         NA
2   sex        30 81.66861    =  1 3.691544 0.05468897
---
                           AIC delta_AIC_vs_null
model10_null          89.36015          0.000000
model10_haemoplasma   90.81430          1.454151
model10_bloodparasite 91.27469          1.914540
model10_season        91.35333          1.993172
model10_sex           87.66861         -1.691544
```

Generate a QQ-plot of deviance residuals from model11_null to visually assess normality:
```
resid_temp <- residuals(model10_null, type = "deviance")
qqnorm(resid_temp)
qqline(resid_temp, col = "red", lwd = 1)
```

Convert `temperature` to numeric, handle left-censored values (<32°C) for sanalysis in Cd:
```
data_adult_Cd <- data_adult_Cd %>%
  mutate(
    temperature_numeric = as.numeric(ifelse(temperature == "< 32.00", 32, temperature)),
    censored = ifelse(temperature == "< 32.00", TRUE, FALSE)
  )
```

Create a left-censored Surv object (temp) for `temperature` in Cd:
```
temp <- Surv(data_adult_Cd$temperature_numeric,
                  event = !data_adult_Cd$censored,
                  type = "left")
```

Fit Gaussian survival regression models to test the effects of `haemoplasma`, `bloodparasite`, `season`, `sex` on `temperature`, with (model_11) and without (model_11b) interactions in Cd:
```
model_11 <- survreg(temp ~ haemoplasma * bloodparasite * season * sex, data = data_adult_Cd, dist = "gaussian")
model_11a <- survreg(temp ~ haemoplasma + bloodparasite + season + sex, data = data_adult_Cd, dist = "gaussian")
```

Compare models using ANOVA and AIC to evaluate the contribution of interaction terms:
```
anova(model_11a, model_11, test = "Chisq")
AIC(model_11a, model_11)
```

Results are:
```
                                       Terms Resid. Df    -2*LL Test Df Deviance  Pr(>Chi)
1 haemoplasma + bloodparasite + season + sex        13 58.22666      NA       NA        NA
2 haemoplasma * bloodparasite * season * sex         2 54.04089    = 11  4.18577 0.9641642
---
          df      AIC
model_11a  6 70.22666
model_11  17 88.04089
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_11a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
temp ~ haemoplasma + bloodparasite + season + sex
              Df    AIC      LRT Pr(>Chi)
<none>           70.227                  
haemoplasma    1 68.390 0.163699   0.6858
bloodparasite  1 68.235 0.008139   0.9281
season         1 68.269 0.042462   0.8367
sex            1 68.405 0.178242   0.6729
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_11a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        70.227    0.0000
haemoplasma   68.390   -1.8363
bloodparasite 68.235   -1.9919
season        68.269   -1.9575
sex           68.405   -1.8218
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model11_null <- survreg(temp ~ 1, data = data_adult_Cd, dist = "gaussian")
model11_haemoplasma <- survreg(temp ~ haemoplasma, data = data_adult_Cd, dist = "gaussian")
model11_bloodparasite <- survreg(temp ~ bloodparasite, data = data_adult_Cd, dist = "gaussian")
model11_season <- survreg(temp ~ season, data = data_adult_Cd, dist = "gaussian")
model11_sex <- survreg(temp ~ sex, data = data_adult_Cd, dist = "gaussian")
anova(model11_null, model11_haemoplasma, test="Chisq")
anova(model11_null, model11_bloodparasite, test="Chisq")
anova(model11_null, model11_season, test="Chisq")
anova(model11_null, model11_sex, test="Chisq")
aics <- AIC(model11_null, model11_haemoplasma, model11_bloodparasite, model11_season, model11_sex)
aic_null <- aics["model11_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
        Terms Resid. Df    -2*LL Test Df  Deviance  Pr(>Chi)
1           1        17 58.73649      NA        NA        NA
2 haemoplasma        16 58.42303    =  1 0.3134589 0.5755654
---
          Terms Resid. Df    -2*LL Test Df   Deviance  Pr(>Chi)
1             1        17 58.73649      NA         NA        NA
2 bloodparasite        16 58.67821    =  1 0.05828631 0.8092253
---
   Terms Resid. Df    -2*LL Test Df   Deviance  Pr(>Chi)
1      1        17 58.73649      NA         NA        NA
2 season        16 58.64685    =  1 0.08964222 0.7646325
---
  Terms Resid. Df    -2*LL Test Df  Deviance  Pr(>Chi)
1     1        17 58.73649      NA        NA        NA
2   sex        16 58.60173    =  1 0.1347589 0.7135479
---
                           AIC delta_AIC_vs_null
model11_null          62.73649          0.000000
model11_haemoplasma   64.42303          1.686541
model11_bloodparasite 64.67821          1.941714
model11_season        64.64685          1.910358
model11_sex           64.60173          1.865241
```

Generate a QQ-plot of deviance residuals from model11_null to visually assess normality:
```
resid_temp <- residuals(model11_null, type = "deviance")
qqnorm(resid_temp)
qqline(resid_temp, col = "red", lwd = 1)
```

## Step 11. Impact of Haemoplasma infections on general health condition 
Test the association between `haemoplasma` and `health_condition` in Bt:
```
table_health_condition_haemoplasma_Bt <- table(data_Bt$haemoplasma, data_Bt$health_condition)
table_health_condition_haemoplasma_Bt
fisher.test(table_health_condition_haemoplasma_Bt)
```

Results are:
```
     D  G
  0 10 78
  1  0  4
---
Fisher's Exact Test for Count Data
data:  table_health_condition_haemoplasma_Bt
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.0750502       Inf
sample estimates:
odds ratio 
       Inf 
```

Test the association between `haemoplasma` and `health_condition` in Cd:
```
table_health_condition_haemoplasma_Cd <- table(data_Cd$haemoplasma, data_Cd$health_condition)
table_health_condition_haemoplasma_Cd
fisher.test(table_health_condition_haemoplasma_Cd)
```

Results are:
```
     D  G
  0  1 14
  1  4 64
---
Fisher's Exact Test for Count Data
data:  table_health_condition_haemoplasma_Cd
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.02167699 12.75597732
sample estimates:
odds ratio 
  1.140982 
```

## Step 12. Impact of Haemoplasma_ infections on female reproductive status 
Test the association between `haemoplasma` and `female_reproductive_status` in Bt:
```
table_haemoplasma_infection_female_Bt <- table(data_Bt$haemoplasma, data_Bt$female_reproductive_status)
table_haemoplasma_infection_female_Bt
fisher.test(table_haemoplasma_infection_female_Bt)
```

Results are:
```
    Female lactating with a young Female non pregnant non lactating Pregnant female
  0                             8                                27               6
  1                             0                                 2               0
---
Fisher's Exact Test for Count Data
data:  table_haemoplasma_infection_female_Bt
p-value = 1
alternative hypothesis: two.sided
```

Test the association between `haemoplasma` and `female_reproductive_status` in Cd:
```
table_haemoplasma_infection_female_Cd <- table(data_Cd$haemoplasma, data_Cd$female_reproductive_status)
table_haemoplasma_infection_female_Cd
fisher.test(table_haemoplasma_infection_female_Cd)
```

Results are:
```
    Female lactating with a young Female non pregnant non lactating Pregnant female
  0                             2                                 7               1
  1                             5                                32               2
---
Fisher's Exact Test for Count Data
data:  table_haemoplasma_infection_female_Cd
p-value = 0.5051
alternative hypothesis: two.sided
```
