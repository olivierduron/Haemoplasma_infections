# **Hemoplasma epidemiological survey : R command lines and script**

We analyzed data on 626 individuals belonging to 44 species of wild mammals captured in French Guiana. The epidemiological dataset includes the following variables for each sampled individual:
- `species` : Species identity (one of the 44 wild mammal species included in this study)
- `order` : Species taxonomic order
- `hemoplasma` : Infection status with hemotropic mycoplasmas (0: Uninfected; 1: Infected)
- `sex` : Sex of the individual (M: Male; F: Female)
- `vertical_stratum` : Primary habitat use in the forest vertical strata (Ground, Canopy, Mixed)
- `activity` : Activity rhythm (Nocturnal, Diurnal)
- `diet` : Dietary category (Phytophage, Omnivore, Insectivore, Carnivore)
- `sociality` : Social organization (Solitary, Group)
- `anaplasmataceae` : Infection status with bacteria of the Anaplasmataceae family (*Anaplasma*, *Ehrlichia* and *Allocryptoplasma*) (0: Uninfected; 1: Infected)
- `apicomplexa` : Infection status with blood parasites, including piroplasmids (*Babesia* and *Theileria*) and haemogregarines (*Hepatozoon* and *Hemolivia*) (0: Uninfected; 1: Infected)
- `trypanosoma` : Infection status with trypanosomes (0: Uninfected; 1: Infected)
- `filaria` : Infection status with microfilariae (0: Uninfected; 1: Infected)
   
Details about all the experimental methods are available in the related manuscript.

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

The epidemiological dataset for mammals is available [here](data_hemoplasma_stat.csv)

This database will be referred to as `data_hemoplasma_stat` throughout the R command lines and scripts provided below. It corresponds to part of the dataset provided in Table S1 of the related manuscript.

Load the dataset directly from the GitHub repository to R and explore the dataset by summarizing the different modalities and their frequencies for each variable:
```
data_hemoplasma_stat <- read.csv2("https://raw.githubusercontent.com/olivierduron/Hemoplasma_infections/main/data_hemoplasma_stat.csv")
data_hemoplasma_stat
str(data_hemoplasma_stat)
get_modalities <- function(x) {sort(table(x), decreasing = TRUE)}
lapply(data_hemoplasma_stat, get_modalities)
```

## Step 2. Prepare the data for analysis

Convert categorical variables into factors:
```
data_hemoplasma_stat$species        <- as.factor(data_hemoplasma_stat$species)
data_hemoplasma_stat$order           <- as.factor(data_hemoplasma_stat$order)
data_hemoplasma_stat$hemoplasma      <- as.factor(data_hemoplasma_stat$hemoplasma)
data_hemoplasma_stat$sex   <- as.factor(data_hemoplasma_stat$sex)
data_hemoplasma_stat$vertical_stratum   <- as.factor(data_hemoplasma_stat$vertical_stratum)
data_hemoplasma_stat$activity   <- as.factor(data_hemoplasma_stat$activity)
data_hemoplasma_stat$diet   <- as.factor(data_hemoplasma_stat$diet)
data_hemoplasma_stat$sociality   <- as.factor(data_hemoplasma_stat$sociality)
data_hemoplasma_stat$anaplasmataceae      <- as.factor(data_hemoplasma_stat$anaplasmataceae)
data_hemoplasma_stat$apicomplexa         <- as.factor(data_hemoplasma_stat$apicomplexa)
data_hemoplasma_stat$trypanosoma    <- as.factor(data_hemoplasma_stat$trypanosoma)
data_hemoplasma_stat$filaria   <- as.factor(data_hemoplasma_stat$filaria)
```

Load libraries for analysis: 
```
library(dplyr)
library(ggplot2)
library(scales)
library(ggthemes)
library(lme4)
library(car)
library(emmeans)
library(brms)
library(patchwork)
```

## Step 3. Test whether hemoplasma detection depends on species sampling effort

Prepare species-level dataset and visualization (Fig. S1) : 
```
species_summary <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(as.numeric(as.character(hemoplasma)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    n_negative = n_sampled - n_positive,
    prevalence = n_positive / n_sampled
  ) %>%
  rowwise() %>%
  mutate(
    ci = list(binom.test(n_positive, n_sampled)$conf.int),
    ci_low = ci[[1]],
    ci_high = ci[[2]]
  ) %>%
  ungroup()
species_summary

#Fig S1
model <- lm(prevalence ~ n_sampled, data = species_summary)
newdata <- data.frame(
  n_sampled = seq(min(species_summary$n_sampled),
                  max(species_summary$n_sampled),
                  length.out = 200)
)
pred <- predict(model, newdata, interval = "confidence")
newdata$fit <- pred[, "fit"]
newdata$lwr <- pred[, "lwr"]
newdata$upr <- pred[, "upr"]
p <- ggplot(species_summary, aes(x = n_sampled, y = prevalence)) +
  geom_ribbon(
    data = newdata,
    aes(x = n_sampled, ymin = lwr, ymax = upr),
    fill = "grey70",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = newdata,
    aes(x = n_sampled, y = fit),
    color = "blue",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal() +
   labs(
    x = "Sample size per species",
    y = "Hemoplasma prevalence"
  )
print(p)
pdf("Fig_S1_hemoplasma_sampling_effect.pdf", width = 7, height = 5)
print(p)
dev.off()
```

Results are: 
```
# A tibble: 44 × 8
   species                 n_sampled n_positive n_negative prevalence ci_low ci_high
 1 Alouatta_macconnelli        22         20          2     0.909   0.708    0.989 
 2 Bradypus_tridactylus       108          4        104     0.0370  0.0102   0.0921
 3 Cabassous_unicinctus         2          0          2     0       0        0.842 
 4 Caluromys_philander          5          0          5     0       0        0.522 
 5 Cebus_apella                 1          0          1     0       0        0.975 
 6 Choloepus_didactylus        90         72         18     0.8     0.702    0.877 
 7 Coendou_melanurus            1          0          1     0       0        0.975 
 8 Coendou_sp                   3          1          2     0.333   0.00840  0.906 
 9 Cyclopes_didactylus          1          0          1     0       0        0.975 
10 Dasypus_novemcinctus        15          5         10     0.333   0.118    0.616 
# 34 more rows
```

Spearman correlation test :
```
cor.test(
  species_summary$n_sampled,
  species_summary$prevalence,
  method = "spearman"
)
```

Results are : 
```
Spearman's rank correlation rho
data:  species_summary$n_sampled and species_summary$prevalence
S = 9760.9, p-value = 0.03915
alternative hypothesis: true rho is not equal to 0
sample estimates:
     rho 
0.312126 
```

## Step 4. Test infection prevalence variation across mammalian orders (+ test sampling bias effect; GLMM with species random effect : model #1) : 

Create and visualize contingency tables order : 
```
df_species <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),                             
    n_infected = sum(hemoplasma == 1, na.rm = TRUE), 
    prevalence = n_infected / n  
  )
df_species <- data_hemoplasma_stat %>%
  group_by(species, order) %>%
  summarise(
    n = n(),
    n_infected = sum(as.numeric(as.character(hemoplasma)) == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    infected = ifelse(n_infected > 0, 1, 0),
    .groups = "drop"
  )
df_order <- df_species %>%
  group_by(order) %>%
  summarise(
    infected_species = sum(infected),
    uninfected_species = n() - sum(infected),
    .groups = "drop"
  )
contingency_table <- df_order %>%
  select(infected_species, uninfected_species) %>%
  as.matrix()
rownames(contingency_table) <- df_order$order
contingency_table
```

Synthesis tables are :
```
Order                infected_species uninfected_species
Carnivora                      3                  3
Cingulata                      1                  1
Didelphimorphia                4                  4
Pilosa                         2                  2
Primates                       2                  3
Rodentia                       8                 11
```

Data preparation :
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  mutate(
    hemoplasma = as.numeric(as.character(hemoplasma))
  ) %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup()
```

Full GLMM: 
```
mod_full <- glmer(
  hemoplasma ~ order + log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 1e5))
)
```

Model selection (LRT) :
```
mod_no_order <- update(mod_full, . ~ . - order)
mod_no_logn  <- update(mod_full, . ~ . - log_n)
anova(mod_full, mod_no_order, test = "Chisq")
anova(mod_full, mod_no_logn, test = "Chisq")
```

Results are: 
```
> anova(mod_full, mod_no_order, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod_no_order: hemoplasma ~ log_n + (1 | species)
mod_full: hemoplasma ~ order + log_n + (1 | species)
             npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
mod_no_order    3 416.94 430.20 -205.47    410.94                       
mod_full        8 414.76 450.11 -199.38    398.76 12.185  5    0.03233 *
> anova(mod_full, mod_no_logn, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod_no_logn: hemoplasma ~ order + (1 | species)
mod_full: hemoplasma ~ order + log_n + (1 | species)
            npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
mod_no_logn    7 416.19 447.13 -201.10    402.19                       
mod_full       8 414.76 450.11 -199.38    398.76 3.4387  1    0.06369 .
```

AIC comparison : 
```
AIC(mod_full, mod_no_order, mod_no_logn)
```

Results are :
```
             df      AIC
mod_full      8 414.7550
mod_no_order  3 416.9406
mod_no_logn   7 416.1937
```

Summarize full model :
```
summary(mod_full)
```

Results are :
```
Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: hemoplasma ~ order + log_n + (1 | species)
   Data: data_hemoplasma_stat
Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))

      AIC       BIC    logLik -2*log(L)  df.resid 
    414.8     450.1    -199.4     398.8       606 

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-3.0444 -0.2310 -0.1769  0.1224  4.7470 

Random effects:
 Groups  Name        Variance Std.Dev.
 species (Intercept) 3.347    1.829   
Number of obs: 614, groups:  species, 44

Fixed effects:
                     Estimate Std. Error z value Pr(>|z|)   
(Intercept)           -1.3932     1.1807  -1.180  0.23800   
orderCingulata        -1.8148     2.0681  -0.878  0.38021   
orderDidelphimorphia  -2.5358     1.5494  -1.637  0.10170   
orderPilosa           -2.9356     2.0342  -1.443  0.14899   
orderPrimates          1.0400     1.7047   0.610  0.54182   
orderRodentia         -3.8349     1.4124  -2.715  0.00662 **
log_n                  0.6804     0.3915   1.738  0.08225 . 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) ordrCn ordrDd ordrPl ordrPr ordrRd
orderCinglt -0.419                                   
ordrDdlphmr -0.527  0.501                            
orderPilosa -0.316  0.445  0.651                     
orderPrimts -0.515  0.431  0.598  0.526              
orderRodent -0.599  0.534  0.751  0.676  0.640       
log_n       -0.333 -0.249 -0.434 -0.587 -0.277 -0.409
```

Predicted infection probabilities with confidence interval per order :
```
emm_prob <- emmeans(mod_full, ~ order, type = "response")
prob_df <- as.data.frame(emm_prob)
prob_df <- prob_df %>%
  mutate(
    percent = prob * 100,
    lower_percent = asymp.LCL * 100,
    upper_percent = asymp.UCL * 100
  )
emm_prob
```

Results are :
```
            order       prob         SE  df   asymp.LCL  asymp.UCL   percent lower_percent upper_percent
1       Carnivora 0.32124636 0.25153810 Inf 0.047018484 0.81949870 32.124636     4.7018484     81.949870
2       Cingulata 0.14399309 0.21028873 Inf 0.005903137 0.82654332 14.399309     0.5903137     82.654332
3 Didelphimorphia 0.10004408 0.07972107 Inf 0.019224258 0.38667771 10.004408     1.9224258     38.667771
4          Pilosa 0.13682287 0.14390491 Inf 0.014342216 0.63326195 13.682287     1.4342216     63.326195
5        Primates 0.72790324 0.24377576 Inf 0.193356424 0.96759081 72.790324    19.3356424     96.759081
6        Rodentia 0.02417335 0.01635629 Inf 0.006324077 0.08794235  2.417335     0.6324077      8.794235
```

Tukey-adjusted post-hoc comparisons of estimated marginal means :
```
pairs(emm_prob, adjust = "tukey")
```

Results are :
```
 contrast                    estimate   SE  df z.ratio p.value
 Carnivora - Cingulata         1.0345 2.05 Inf   0.505  0.9960
 Carnivora - Didelphimorphia   1.4487 1.44 Inf   1.007  0.9159
 Carnivora - Pilosa            1.0939 1.66 Inf   0.659  0.9863
 Carnivora - Primates         -1.7321 1.67 Inf  -1.039  0.9050
 Carnivora - Rodentia          2.9500 1.33 Inf   2.222  0.2273
 Cingulata - Didelphimorphia   0.4142 1.89 Inf   0.220  0.9999
 Cingulata - Pilosa            0.0594 2.05 Inf   0.029  1.0000
 Cingulata - Primates         -2.7665 2.05 Inf  -1.347  0.7586
 Cingulata - Rodentia          1.9155 1.80 Inf   1.065  0.8952
 Didelphimorphia - Pilosa     -0.3548 1.46 Inf  -0.244  0.9999
 Didelphimorphia - Primates   -3.1807 1.46 Inf  -2.180  0.2473
 Didelphimorphia - Rodentia    1.5013 1.07 Inf   1.409  0.7218
 Pilosa - Primates            -2.8259 1.67 Inf  -1.694  0.5358
 Pilosa - Rodentia             1.8561 1.34 Inf   1.383  0.7372
 Primates - Rodentia           4.6820 1.34 Inf   3.485  0.0065
Results are given on the log odds ratio (not the response) scale. 
P value adjustment: tukey method for comparing a family of 6 estimates 
```

Create a plot of hemoplasma prevalence by species and mammalian order :
```
df_species <- data_hemoplasma_stat %>%
  group_by(species, order) %>%
  summarise(
    n = n(),
    n_infected = sum(as.numeric(as.character(hemoplasma)) == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  )
mod_glmm <- glmer(
  hemoplasma ~ order + log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 1e5))
)
emm_prob <- emmeans(mod_glmm, ~ order, type = "response")
prob_df <- as.data.frame(emm_prob)
prob_df <- prob_df %>%
  arrange(prob) %>%
  mutate(order = factor(order, levels = order))
df_species <- df_species %>%
  mutate(order = factor(order, levels = levels(prob_df$order)))
order_colors <- c(
  "Primates" = "#0072B2",
  "Pilosa" = "#E69F00",
  "Cingulata" = "#009E73",
  "Rodentia" = "#D55E00",
  "Carnivora" = "#CC79A7",
  "Didelphimorphia" = "#F0E442"
)
set.seed(1)
p <- ggplot() +
  geom_jitter(
    data = df_species,
    aes(x = order, y = prevalence, color = order, size = n),
    width = 0.35,
    height = 0,
    alpha = 0.5
  ) +
  geom_segment(
    data = prob_df,
    aes(
      x = as.numeric(order) - 0.25,
      xend = as.numeric(order) + 0.25,
      y = asymp.LCL,
      yend = asymp.LCL,
      color = order
    ),
    linewidth = 1
  ) +
  geom_segment(
    data = prob_df,
    aes(
      x = as.numeric(order) - 0.25,
      xend = as.numeric(order) + 0.25,
      y = prob,
      yend = prob,
      color = order
    ),
    linewidth = 1.5
  ) +
  geom_segment(
    data = prob_df,
    aes(
      x = as.numeric(order) - 0.25,
      xend = as.numeric(order) + 0.25,
      y = asymp.UCL,
      yend = asymp.UCL,
      color = order
    ),
    linewidth = 1
  ) +
  scale_color_manual(values = order_colors) +
  scale_size_continuous(range = c(2, 10), name = "Sample size") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Mammalian order",
    y = "Hemoplasma prevalence"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )
print(p)
ggsave(
  filename = "Fig_2_Hemoplasma_prevalence_by_order.pdf",
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)
```

## Step 5. Test infection prevalence variation across sex in infected mammalian species (GLMM with species random effect : model #2) : 

Prepare data :
```
species_infected <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(infected = any(hemoplasma == 1, na.rm = TRUE)) %>%
  filter(infected) %>%
  pull(species)
data_inf <- data_hemoplasma_stat %>%
  filter(species %in% species_infected)
data_inf <- data_inf %>%
  mutate(
    hemoplasma = as.numeric(as.character(hemoplasma))
  ) %>%
  filter(!is.na(sex)) %>%   # IMPORTANT pour LRT
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup()
```

GLMM full model :
```
mod_sex_inf <- glmer(
  hemoplasma ~ sex + log_n + (1 | species),
  family = binomial,
  data = data_inf,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 1e5))
)
summary(mod_sex_inf)
```

XXX : 
```
# ----------------------------
# 4. Test effet sexe (LRT correct)
# ----------------------------
mod_no_sex_inf <- update(mod_sex_inf, . ~ . - sex)

anova(mod_sex_inf, mod_no_sex_inf, test = "Chisq")

# ----------------------------
# 5. Probabilités interprétables (emmeans)
# ----------------------------
emm_sex_inf <- emmeans(mod_sex_inf, ~ sex, type = "response")

prob_sex_df <- as.data.frame(emm_sex_inf)

print(prob_sex_df)
```







## Step 5. Test infection prevalence variation across species ecological traits (GLMM with species random effect : model #2) : 

Prepare data (sampling effort) :
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup()
```

Fit models : 
```
# Null model (random effect only)
mod_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat
)

# Sampling effort model
mod_sampling <- glmer(
  hemoplasma ~ log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 1e5))
)

# Full ecological model
mod_traits_full <- glmer(
  hemoplasma ~ log_n + vertical_stratum + activity + diet + sociality + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)
```

Model summaries: 
```
summary(mod_traits_full)
```

Results are: 
```
Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: hemoplasma ~ log_n + vertical_stratum + activity + diet + sociality +      (1 | species)
   Data: data_hemoplasma_stat
Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))

      AIC       BIC    logLik -2*log(L)  df.resid 
    421.8     466.0    -200.9     401.8       604 

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-2.9825 -0.2388 -0.1417  0.1252  4.8857 

Random effects:
 Groups  Name        Variance Std.Dev.
 species (Intercept) 4.676    2.162   
Number of obs: 614, groups:  species, 44

Fixed effects:
                       Estimate Std. Error z value Pr(>|z|)  
(Intercept)             -0.2649     2.4589  -0.108   0.9142  
log_n                    0.7135     0.3841   1.858   0.0632 .
vertical_stratumGround   0.2483     1.1325   0.219   0.8265  
vertical_stratumMixed   -1.9842     2.2823  -0.869   0.3846  
activityNocturnal       -3.1289     1.6083  -1.946   0.0517 .
dietInsectivore         -1.2706     2.7936  -0.455   0.6492  
dietOmnivore            -1.1872     2.4231  -0.490   0.6242  
dietPhytophage          -0.7991     2.6323  -0.304   0.7615  
socialitySolitary       -0.2983     1.4131  -0.211   0.8328  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) log_n  vrtc_G vrtc_M actvtN dtInsc dtOmnv dtPhyt
log_n       -0.137                                                 
vrtcl_strtG -0.500 -0.008                                          
vrtcl_strtM -0.284  0.007  0.324                                   
actvtyNctrn  0.362 -0.142 -0.270 -0.082                            
dietInsctvr -0.555 -0.101  0.112  0.008 -0.507                     
dietOmnivor -0.735 -0.206  0.292  0.046 -0.586  0.781              
dietPhytphg -0.812 -0.170  0.460  0.190 -0.551  0.714  0.894       
socltySltry -0.621  0.013  0.209  0.273 -0.530  0.264  0.392  0.460
```

Likelihood Ratio Tests (LRT) :
```
anova(mod_null, mod_sampling, mod_traits_full, test = "Chisq")
```

Results:
```
Data: data_hemoplasma_stat
Models:
mod_null: hemoplasma ~ 1 + (1 | species)
mod_sampling: hemoplasma ~ log_n + (1 | species)
mod_traits_full: hemoplasma ~ log_n + vertical_stratum + activity + diet + sociality + (1 | species)
                npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod_null           2 416.38 425.22 -206.19    412.38                     
mod_sampling       3 416.94 430.20 -205.47    410.94 1.4399  1     0.2302
mod_traits_full   10 421.76 465.96 -200.88    401.76 9.1816  7     0.2399
```

AIC comparison : 
```
model_set <- list(
  null = mod_null,
  sampling = mod_sampling,
  traits = mod_traits_full
)
AIC_table <- data.frame(
  model = names(model_set),
  AIC = sapply(model_set, AIC)
)
AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table
```

Results are:
```
            model      AIC delta_AIC
null         null 416.3804 0.0000000
sampling sampling 416.9406 0.5601158
traits     traits 421.7590 5.3785552
```

Single-term deletion:
```
mod_drop <- drop1(mod_traits_full, test = "Chisq")
mod_drop
```

Results are:
```
Single term deletions
Model:
hemoplasma ~ log_n + vertical_stratum + activity + diet + sociality + 
    (1 | species)
                 npar    AIC    LRT Pr(Chi)  
<none>                421.76                 
log_n               1 423.37 3.6160 0.05723 .
vertical_stratum    2 418.86 1.0977 0.57760  
activity            1 423.09 3.3279 0.06812 .
diet                3 416.11 0.3554 0.94930  
sociality           1 419.80 0.0445 0.83297  
```

Predicted probabilities and 95% CI:
```
emm_stratum <- emmeans(mod_traits_full, ~ vertical_stratum, type = "response")
prob_stratum <- as.data.frame(emm_stratum) %>%
  mutate(
    percent = prob * 100,
    lower = asymp.LCL * 100,
    upper = asymp.UCL * 100
  )
emm_activity <- emmeans(mod_traits_full, ~ activity, type = "response")
prob_activity <- as.data.frame(emm_activity) %>%
  mutate(
    percent = prob * 100,
    lower = asymp.LCL * 100,
    upper = asymp.UCL * 100
  )
emm_diet <- emmeans(mod_traits_full, ~ diet, type = "response")
prob_diet <- as.data.frame(emm_diet) %>%
  mutate(
    percent = prob * 100,
    lower = asymp.LCL * 100,
    upper = asymp.UCL * 100
  )
emm_social <- emmeans(mod_traits_full, ~ sociality, type = "response")
prob_social <- as.data.frame(emm_social) %>%
  mutate(
    percent = prob * 100,
    lower = asymp.LCL * 100,
    upper = asymp.UCL * 100
  )
prob_stratum
prob_activity
prob_diet
prob_social
```

Results are: 
```
  vertical_stratum       prob        SE  df   asymp.LCL asymp.UCL   percent      lower    upper
1           Canopy 0.41118131 0.2977796 Inf 0.058983632 0.8861024 41.118131  5.8983632 88.61024
2           Ground 0.47232960 0.2609295 Inf 0.103147778 0.8744772 47.232960 10.3147778 87.44772
3            Mixed 0.08759731 0.1765953 Inf 0.001261804 0.8794559  8.759731  0.1261804 87.94559

   activity       prob         SE  df   asymp.LCL asymp.UCL   percent      lower    upper
1   Diurnal 0.65174222 0.32496720 Inf 0.101611655 0.9687158 65.174222 10.1611655 96.87158
2 Nocturnal 0.07570553 0.09187215 Inf 0.006209163 0.5177779  7.570553  0.6209163 51.77779

         diet      prob        SE  df   asymp.LCL asymp.UCL  percent     lower    upper
1   Carnivore 0.4691611 0.5834336 Inf 0.008880171 0.9886597 46.91611 0.8880171 98.86597
2 Insectivore 0.1987451 0.3171466 Inf 0.004978985 0.9247860 19.87451 0.4978985 92.47860
3    Omnivore 0.2123682 0.1722342 Inf 0.034594109 0.6698367 21.23682 3.4594109 66.98367
4  Phytophage 0.2844348 0.2717129 Inf 0.028219860 0.8447450 28.44348 2.8219860 84.47450

  sociality      prob        SE  df  asymp.LCL asymp.UCL  percent    lower    upper
1     Group 0.3124687 0.3073413 Inf 0.02679051 0.8823988 31.24687 2.679051 88.23988
2  Solitary 0.2522081 0.2257055 Inf 0.03129654 0.7788041 25.22081 3.129654 77.88041
```

Create a plot of hemoplasma prevalence by species ecological traits :
```
df_species <- data_hemoplasma_stat %>%
  group_by(species, vertical_stratum, activity, diet, sociality) %>%
  summarise(
    n = n(),
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  ) %>%
  mutate(
    prevalence = ifelse(is.nan(prevalence), 0, prevalence)
  )

fix_emm <- function(df, var){

  df <- as.data.frame(df)

  df[[var]] <- as.character(df[[var]])

  df$prob  <- as.numeric(df$prob)
  df$lower <- as.numeric(df$asymp.LCL)
  df$upper <- as.numeric(df$asymp.UCL)

  df
}

prob_stratum  <- fix_emm(prob_stratum, "vertical_stratum")
prob_activity  <- fix_emm(prob_activity, "activity")
prob_diet      <- fix_emm(prob_diet, "diet")
prob_social    <- fix_emm(prob_social, "sociality")

make_panel <- function(df, pred, var, palette, title){

  df[[var]] <- factor(df[[var]])
  pred[[var]] <- factor(pred[[var]])

  ggplot() +

    # species points
    geom_jitter(
      data = df,
      aes(x = .data[[var]], y = prevalence,
          color = .data[[var]], size = n),
      width = 0.20,
      alpha = 0.65,
      na.rm = TRUE
    ) +

    # LOWER CI
    geom_segment(
      data = pred,
      aes(
        x = as.numeric(.data[[var]]) - 0.15,
        xend = as.numeric(.data[[var]]) + 0.15,
        y = lower,
        yend = lower,
        color = .data[[var]]
      ),
      linewidth = 0.9,
      na.rm = TRUE
    ) +

    # MEAN
    geom_segment(
      data = pred,
      aes(
        x = as.numeric(.data[[var]]) - 0.15,
        xend = as.numeric(.data[[var]]) + 0.15,
        y = prob,
        yend = prob,
        color = .data[[var]]
      ),
      linewidth = 1.3,
      na.rm = TRUE
    ) +

    # UPPER CI
    geom_segment(
      data = pred,
      aes(
        x = as.numeric(.data[[var]]) - 0.15,
        xend = as.numeric(.data[[var]]) + 0.15,
        y = upper,
        yend = upper,
        color = .data[[var]]
      ),
      linewidth = 0.9,
      na.rm = TRUE
    ) +

    scale_color_manual(values = palette, drop = FALSE) +
    scale_size_continuous(range = c(2, 9), name = "Sample size") +

    scale_y_continuous(
      labels = percent_format(accuracy = 1)
    ) +

    labs(
      x = NULL,
      y = "Hemoplasma prevalence",
      title = title
    ) +

    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}
pA <- make_panel(df_species, prob_stratum, "vertical_stratum", stratum_colors, "A")
pB <- make_panel(df_species, prob_activity, "activity", activity_colors, "B")
pC <- make_panel(df_species, prob_diet, "diet", diet_colors, "C")
pD <- make_panel(df_species, prob_social, "sociality", social_colors, "D")
figure1 <- (pA | pB) / (pC | pD)
figure1
pdf(file = file.path(getwd(), "Figure1_Hemoplasma.pdf"),
    width = 12, height = 9, useDingbats = FALSE)
print(figure1)
dev.off()
```

yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy





Class size effects on hemoplasma infection prevalence (GLMM with species random effect, model #3) : 
```
df_body <- data %>%
  group_by(species, body_size) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )
# Full model (body size only)
mod_body <- glmer(
  hemoplasma ~ body_size + (1 | species),
  data = data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

# Null model (only random effect)
mod_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  data = data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

# Summary
summary(mod_body)

# LRT (global test of body_size effect)
anova(mod_null, mod_body, test = "Chisq")

# AIC comparison
AIC_table <- data.frame(
  model = c("body_size + species", "null (species only)"),
  AIC = c(AIC(mod_body), AIC(mod_null))
)
AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table


# Post-hoc (if needed)
library(emmeans)
emmeans(mod_body, ~ body_size, type = "response")

```

Results are : 
```
Data: data
Models:
mod_null: hemoplasma ~ 1 + (1 | species)
mod_body: hemoplasma ~ body_size + (1 | species)
         npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod_null    2 416.38 425.22 -206.19    412.38                     
mod_body    4 417.54 435.22 -204.77    409.54 2.8398  2     0.2417

                model      AIC delta_AIC
1 body_size + species 417.5406  1.160212
2 null (species only) 416.3804  0.000000
```

Predicted infection probabilities  with confidence interval per mammalian class size :
```
emm_body <- emmeans(mod_body, ~ body_size, type = "response")
body_pred <- as.data.frame(emm_body)
body_pred
```

Results are : 
```
body_size      prob        SE  df   asymp.LCL asymp.UCL
 Large     0.5705841 0.5767685 Inf 0.013002934 0.9925935
 Medium    0.0857229 0.0578205 Inf 0.021605423 0.2847425
 Small     0.0258783 0.0243023 Inf 0.003999273 0.1494869
```


XXXXXXXXXXXXXXXXXXXXXXXXXX









## Step 3. Hemoplasma prevalence analysis across host species

Summarize by species:

```
df_species <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),                             
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),  # infected individuals
    prevalence = n_infected / n  
  )
print(df_species, n = Inf)
```

Results are:
```
# A tibble: 44 × 4
   species                       n n_infected prevalence
   <fct>                     <int>      <int>      <dbl>
 1 Alouatta_macconnelli         22         20     0.909 
 2 Bradypus_tridactylus        108          4     0.0370
 3 Cabassous_unicinctus          2          0     0     
 4 Caluromys_philander           5          0     0     
 5 Cebus_apella                  1          0     0     
 6 Choloepus_didactylus         90         72     0.8   
 7 Coendou_melanurus             1          0     0     
 8 Coendou_sp                    3          1     0.333 
 9 Cyclopes_didactylus           1          0     0     
10 Dasypus_novemcinctus         15          5     0.333 
11 Didelphis_marsupialis        51         22     0.431 
12 Eira_barbara                  4          0     0     
13 Felis_wiedii                  1          0     0     
14 Galictis_vittata              4          3     0.75  
15 Holochilus_sciureus           5          1     0.2   
16 Hydrochoerus_hydrochaeris     2          0     0     
17 Hylaeamys_megacephalus       15          0     0     
18 Hylaeamys_yunganus           10          0     0     
19 Lontra_longicaudis            1          1     1     
20 Makalata_didelphoides         8          0     0     
21 Marmosa_lepida                1          1     1     
22 Marmosa_murina               20          2     0.1   
23 Marmosops_parvidens           5          0     0     
24 Mesomys_hispidus             13          0     0     
25 Metachirus_nudicaudatus       5          0     0     
26 Micoureus_demerarae          16          0     0     
27 Mus_musculus                 34          0     0     
28 Neacomys_dubosti              1          0     0     
29 Neacomys_paracou              8          0     0     
30 Nectomys_rattus               4          2     0.5   
31 Oecomys_auyantepui           16          1     0.0625
32 Oecomys_bicolor              16          0     0     
33 Oligoryzomys_fulvescens       7          1     0.143 
34 Philander_opossum            20          8     0.4   
35 Pithecia_pithecia             1          0     0     
36 Potos_flavus                  2          1     0.5   
37 Proechimys_cuvieri           18          1     0.0556
38 Proechimys_guyannensis       20          1     0.05  
39 Puma_yagouaroundi             5          0     0     
40 Rattus_rattus                19          1     0.0526
41 Saguinus_midas               41         41     1     
42 Saimiri_sciureus              1          0     0     
43 Sciurus_aestuans              1          0     0     
44 Tamandua_tetradactyla         3          0     0     
```

Scatter plot
```
ggplot(df_species, aes(x = n, y = prevalence)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(
    x = "Number of sampled individuals per species",
    y = "Hemoplasma prevalence"
  )
```

Binomial GLM with all species
```
glm_model <- glm(
  cbind(n_infected, n - n_infected) ~ n,
  family = binomial,
  data = df_species
)
summary(glm_model)

# Null model
glm_model_null <- glm(
  cbind(n_infected, n - n_infected) ~ 1,
  family = binomial,
  data = df_species
)

# Likelihood ratio test
anova_all <- anova(glm_model_null, glm_model, test = "Chisq")
print(anova_all)
```

Results are:
```
Analysis of Deviance Table
Model 1: cbind(n_infected, n - n_infected) ~ 1
Model 2: cbind(n_infected, n - n_infected) ~ n
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        43     442.01                          
2        42     428.11  1   13.903 0.0001925 ***
```

Binomial GLM filtering species with n < 5
```
df_species_filtered <- df_species %>%
  filter(n >= 5)
# GLM with filtered data
glm_model_f <- glm(
  cbind(n_infected, n - n_infected) ~ n,
  family = binomial,
  data = df_species_filtered
)
summary(glm_model_f)

# Null model for filtered data
glm_model_null_f <- glm(
  cbind(n_infected, n - n_infected) ~ 1,
  family = binomial,
  data = df_species_filtered
)

# Likelihood ratio test
anova_filtered <- anova(glm_model_null_f, glm_model_f, test = "Chisq")
print(anova_filtered)
```

Results are:
```
Analysis of Deviance Table
Model 1: cbind(n_infected, n - n_infected) ~ 1
Model 2: cbind(n_infected, n - n_infected) ~ n
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        25     419.11                          
2        24     405.19  1   13.916 0.0001912 ***
```

## Step 5. Hemoplasma prevalence by species

Prepare species-level data and keep only species with at least 1 infected individual
```
df_species <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_infected = sum(as.numeric(as.character(hemoplasma)) == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  )
df_infected <- df_species %>%
  filter(n_infected > 0)
```

Binomial GLM (for all infected species)
```
glm_all <- glm(
  cbind(n_infected, n - n_infected) ~ species,
  family = binomial,
  data = df_infected
)
glm_null_all <- glm(
  cbind(n_infected, n - n_infected) ~ 1,
  family = binomial,
  data = df_infected
)
anova_all <- anova(glm_null_all, glm_all, test = "Chisq")
print("=== All infected species ===")
print(summary(glm_all))
print(anova_all)
```

Results are:
```
Call:
glm(formula = cbind(n_infected, n - n_infected) ~ species, family = binomial, 
    data = df_infected)
Coefficients:
                                 Estimate Std. Error z value Pr(>|z|)    
(Intercept)                        2.3026     0.7416   3.105  0.00190 ** 
speciesBradypus_tridactylus       -5.5607     0.8998  -6.180 6.41e-10 ***
speciesCholoepus_didactylus       -0.9163     0.7870  -1.164  0.24434    
speciesCoendou_sp                 -2.9957     1.4318  -2.092  0.03641 *  
speciesDasypus_novemcinctus       -2.9957     0.9220  -3.249  0.00116 ** 
speciesDidelphis_marsupialis      -2.5788     0.7937  -3.249  0.00116 ** 
speciesGalictis_vittata           -1.2040     1.3723  -0.877  0.38032    
speciesHolochilus_sciureus        -3.6889     1.3416  -2.750  0.00597 ** 
speciesLontra_longicaudis         21.2635 79462.0195   0.000  0.99979    
speciesMarmosa_lepida             21.2635 79461.9966   0.000  0.99979    
speciesMarmosa_murina             -4.4998     1.0515  -4.280 1.87e-05 ***
speciesNectomys_rattus            -2.3026     1.2450  -1.849  0.06439 .  
speciesOecomys_auyantepui         -5.0106     1.2715  -3.941 8.12e-05 ***
speciesOligoryzomys_fulvescens    -4.0943     1.3102  -3.125  0.00178 ** 
speciesPhilander_opossum          -2.7081     0.8708  -3.110  0.00187 ** 
speciesPotos_flavus               -2.3026     1.5969  -1.442  0.14932    
speciesProechimys_cuvieri         -5.1358     1.2684  -4.049 5.14e-05 ***
speciesProechimys_guyannensis     -5.2470     1.2660  -4.145 3.40e-05 ***
speciesRattus_rattus              -5.1930     1.2671  -4.098 4.16e-05 ***
speciesSaguinus_midas             24.1352 52162.4892   0.000  0.99963    

(Dispersion parameter for binomial family taken to be 1)
    Null deviance: 3.0552e+02  on 19  degrees of freedom
Residual deviance: 5.0351e-10  on  0  degrees of freedom
AIC: 81.765
Number of Fisher Scoring iterations: 22

Analysis of Deviance Table
Model 1: cbind(n_infected, n - n_infected) ~ 1
Model 2: cbind(n_infected, n - n_infected) ~ species
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        19     305.52                          
2         0       0.00 19   305.52 < 2.2e-16 ***
```

Binomial GLM (for infected species with n >= 5)
```
df_infected_f <- df_infected %>% filter(n >= 5)
glm_filtered <- glm(
  cbind(n_infected, n - n_infected) ~ species,
  family = binomial,
  data = df_infected_f
)
glm_null_filtered <- glm(
  cbind(n_infected, n - n_infected) ~ 1,
  family = binomial,
  data = df_infected_f
)
anova_filtered <- anova(glm_null_filtered, glm_filtered, test = "Chisq")
print("=== Species with n >= 5 ===")
print(summary(glm_filtered))
print(anova_filtered)
```
Results are:
```
Call: glm(formula = cbind(n_infected, n - n_infected) ~ species, family = binomial, 
    data = df_infected_f)
Coefficients:
                                 Estimate Std. Error z value Pr(>|z|)    
(Intercept)                        2.3026     0.7416   3.105  0.00190 ** 
speciesBradypus_tridactylus       -5.5607     0.8998  -6.180 6.41e-10 ***
speciesCholoepus_didactylus       -0.9163     0.7870  -1.164  0.24434    
speciesDasypus_novemcinctus       -2.9957     0.9220  -3.249  0.00116 ** 
speciesDidelphis_marsupialis      -2.5788     0.7937  -3.249  0.00116 ** 
speciesHolochilus_sciureus        -3.6889     1.3416  -2.750  0.00597 ** 
speciesMarmosa_murina             -4.4998     1.0515  -4.280 1.87e-05 ***
speciesOecomys_auyantepui         -5.0106     1.2715  -3.941 8.12e-05 ***
speciesOligoryzomys_fulvescens    -4.0943     1.3102  -3.125  0.00178 ** 
speciesPhilander_opossum          -2.7081     0.8708  -3.110  0.00187 ** 
speciesProechimys_cuvieri         -5.1358     1.2684  -4.049 5.14e-05 ***
speciesProechimys_guyannensis     -5.2470     1.2660  -4.145 3.40e-05 ***
speciesRattus_rattus              -5.1930     1.2671  -4.098 4.16e-05 ***
speciesSaguinus_midas             24.1352 52162.4892   0.000  0.99963    

(Dispersion parameter for binomial family taken to be 1)
    Null deviance: 2.9957e+02  on 13  degrees of freedom
Residual deviance: 2.7045e-10  on  0  degrees of freedom
AIC: 63.069
Number of Fisher Scoring iterations: 22

Analysis of Deviance Table
Model 1: cbind(n_infected, n - n_infected) ~ 1
Model 2: cbind(n_infected, n - n_infected) ~ species
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        13     299.57                          
2         0       0.00 13   299.57 < 2.2e-16 ***
```

Scatter plot for infected species
```
plot_file <- "hemoplasma_prevalence_species.pdf"
pdf(plot_file, width = 8, height = 5)  # PDF output
ggplot(df_infected, aes(x = species, y = prevalence, size = n)) +
  geom_point(alpha = 0.5, color = "#69b3a2") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Species",
    y = "Hemoplasma prevalence",
    size = "Sample size"
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_size_continuous(range = c(2, 8))
dev.off()
cat("PDF saved to:", plot_file, "\n")
```

## Step 6. Hemoplasma prevalence by species | order

GLMM GLOBAL
```
glmm_order <- glmer(
  hemoplasma ~ order + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(optimizer = "bobyqa")
)
glmm_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(optimizer = "bobyqa")
)
anova(glmm_null, glmm_order, test = "Chisq")
summary(glmm_order)

Results are:
```
Data: data_hemoplasma_stat
Models:
glmm_null: hemoplasma ~ 1 + (1 | species)
glmm_order: hemoplasma ~ order + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
glmm_null     2 416.38 425.22 -206.19    412.38                       
glmm_order    7 416.19 447.13 -201.10    402.19 10.187  5    0.07011 .

Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: hemoplasma ~ order + (1 | species)
   Data: data_hemoplasma_stat
Control: glmerControl(optimizer = "bobyqa")

      AIC       BIC    logLik -2*log(L)  df.resid 
    416.2     447.1    -201.1     402.2       607 

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-2.9140 -0.2196 -0.1671  0.1416  4.8727 

Random effects:
 Groups  Name        Variance Std.Dev.
 species (Intercept) 3.633    1.906   
Number of obs: 614, groups:  species, 44

Fixed effects:
                     Estimate Std. Error z value Pr(>|z|)  
(Intercept)            -0.748      1.154  -0.648   0.5167  
orderCingulata         -1.034      2.046  -0.505   0.6132  
orderDidelphimorphia   -1.449      1.439  -1.007   0.3141  
orderPilosa            -1.094      1.661  -0.659   0.5102  
orderPrimates           1.732      1.668   1.039   0.2990  
orderRodentia          -2.950      1.327  -2.222   0.0263 *

Correlation of Fixed Effects:
            (Intr) ordrCn ordrDd ordrPl ordrPr
orderCinglt -0.552                            
ordrDdlphmr -0.788  0.459                     
orderPilosa -0.680  0.402  0.566              
orderPrimts -0.675  0.403  0.567  0.498       
orderRodent -0.853  0.500  0.706  0.617  0.619
```


FIGURE TRAITS
```
library(dplyr)
library(ggplot2)
library(tidyr)

# =========================
# 1. AGGREGATION PAR ESPECE
# =========================
df_species <- data_hemoplasma_stat %>%
  group_by(species, body_size, vertical_stratum, locomotion,
           activity, diet, sociality) %>%
  summarise(
    n = n(),
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  )

# =========================
# 2. FORMAT LONG
# =========================
df_long <- df_species %>%
  pivot_longer(
    cols = c(body_size, vertical_stratum, locomotion,
             activity, diet, sociality),
    names_to = "trait",
    values_to = "category"
  )

df_long$trait <- factor(df_long$trait,
                        levels = c("body_size",
                                   "vertical_stratum",
                                   "locomotion",
                                   "activity",
                                   "diet",
                                   "sociality"))

# =========================
# 3. COULEURS PAR TRAIT (SEULEMENT 6)
# =========================
trait_colors <- c(
  body_size = "#A6CEE3",
  vertical_stratum = "#B2DF8A",
  locomotion = "#FDBF6F",
  activity = "#CAB2D6",
  diet = "#FFFF99",
  sociality = "#FB9A99"
)

# =========================
# 4. PLOT
# =========================
p <- ggplot(df_long, aes(x = category, y = prevalence)) +
  
  geom_point(aes(size = n, fill = trait),
             shape = 21, color = "black", alpha = 0.85) +
  
  scale_fill_manual(values = trait_colors) +
  
  scale_size(range = c(2, 8)) +
  
  facet_wrap(~ trait, scales = "free_x", ncol = 3) +
  
  theme_classic(base_size = 14) +
  
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  
  labs(
    x = NULL,
    y = "Hemoplasma prevalence (per species)"
  )

p
ggsave(
  filename = "hemoplasma_traits_prevalence.pdf",
  plot = p,
  width = 12,
  height = 7
)
```

