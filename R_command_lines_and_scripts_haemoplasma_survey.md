# **Hemoplasma epidemiological survey: R scripts and analysis pipeline**

We analyzed data from 626 individuals belonging to 44 species of wild mammals sampled in French Guiana. The dataset includes host ecological traits and infection status with multiple blood-borne parasites.

Details of sampling and laboratory procedures are provided in the associated manuscript.

## Dataset description
Each row corresponds to one individual and includes the following variables :
- `species` : Species identity (44 wild mammal species)
- `order` : Taxonomic order
- `hemoplasma` : Infection status with hemotropic mycoplasmas (0 = Uninfected; 1 = Infected)
- `sex` : Sex of the individual (M = Male; F = Female)
- `vertical_stratum` : Habitat use in the forest strata (Ground, Canopy, Mixed)
- `activity` : Activity rhythm (Nocturnal, Diurnal)
- `diet` : Dietary category (Phytophage, Omnivore, Insectivore, Carnivore)
- `sociality` : Social organization (Solitary, Group)
- `anaplasmataceae` : Infection status with bacteria of the Anaplasmataceae family (*Anaplasma*, *Ehrlichia* and *Allocryptoplasma*) (0 = Uninfected; 1 = Infected)
- `apicomplexa` : Infection status with blood parasites, including piroplasmids (*Babesia* and *Theileria*) and haemogregarines (*Hepatozoon* and *Hemolivia*) (0 = Uninfected; 1 = Infected)
- `trypanosoma` : Infection status with trypanosomes (0 = Uninfected; 1 = Infected)
- `filaria` : Infection status with microfilariae (0 = Uninfected; 1 = Infected)
   
## Table of contents A UPDATER !!!
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

## Step 1. Data retrieval

The epidemiological dataset is available in the GitHub repository [here](data_hemoplasma_stat.csv)
```
data_hemoplasma_stat <- read.csv2("https://raw.githubusercontent.com/olivierduron/Hemoplasma_infections/main/data_hemoplasma_stat.csv")
data_hemoplasma_stat
str(data_hemoplasma_stat)
get_modalities <- function(x) {sort(table(x), decreasing = TRUE)}
lapply(data_hemoplasma_stat, get_modalities)
```

## Step 2. Data preparation

Convert categorical variables :
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

Load required libraries : 
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

## Step 3. Assess sampling effort bias on hemoplasma prevalence

We first evaluated whether infection prevalence depends on species sampling effort.

### Species-level summary
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
```

Results : 
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

### Visualization (Fig. S1)
```
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

### Spearman correlation test
```
cor.test(
  species_summary$n_sampled,
  species_summary$prevalence,
  method = "spearman"
)
```

Results : 
```
Spearman's rank correlation rho
data:  species_summary$n_sampled and species_summary$prevalence
S = 9760.9, p-value = 0.03915
alternative hypothesis: true rho is not equal to 0
sample estimates:
     rho 
0.312126 
```

### Interpretation
Hemoplasma prevalence increased weakly but significantly with sample size per species, suggesting that prevalence in mammals is likely underestimated in less sampled species.

## Step 4. Variation in hemoplasma infection across mammalian orders (GLMM model 1) 

### Contingency table
```
df_species <- data_hemoplasma_stat %>%
  group_by(species, order) %>%
  summarise(
    n = n(),
    n_infected = sum(as.numeric(as.character(hemoplasma)) == 1, na.rm = TRUE),
    infected = as.integer(n_infected > 0),
    .groups = "drop"
  )
df_order <- df_species %>%
  group_by(order) %>%
  summarise(
    infected_species = sum(infected),
    uninfected_species = n() - sum(infected),
    .groups = "drop"
  )
contingency_table <- as.matrix(df_order[, c("infected_species", "uninfected_species")])
rownames(contingency_table) <- df_order$order
contingency_table
```

Results :
```
Order                infected_species uninfected_species
Carnivora                      3                  3
Cingulata                      1                  1
Didelphimorphia                4                  4
Pilosa                         2                  2
Primates                       2                  3
Rodentia                       8                 11
```

### GLMM preparation
Sampling effort is included as a covariate (log-transformed number of individuals per species) :
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup() %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma)))
```

### GLMM (Model 1)
This model tests whether `hemoplasma` infection probability varies among mammalian orders (`order`) while controlling for differences in sampling effort (`log_n`) and accounting for species-level random effects (`1 | species`).
```
mod1_full <- glmer(
  hemoplasma ~ order + log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)
```

### Model term significance testing
Model terms were evaluated using likelihood ratio tests via single-term deletions (drop1 function with Chi-square tests).
```
res1 <- drop1(mod1_full, test = "Chisq")
res1
```

Results :
```
Single term deletions
Model:
hemoplasma ~ order + log_n + (1 | species)
       npar    AIC     LRT Pr(Chi)  
<none>      414.76                  
order     5 416.94 12.1855 0.03233 *
log_n     1 416.19  3.4387 0.06369 .
```

### Interpretation
Hemoplasma infection probability varied significantly among mammalian orders (χ² test, _p_ = 0.032), indicating a non-random distribution of infection across host taxonomic groups.

A marginal effect of sampling effort (`log_n`) was also detected (_p_ = 0.064), suggesting a weak influence of species sampling intensity on observed prevalence.

### Model comparison with null and univariate models
```
mod1_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

mod1_order <- glmer(
  hemoplasma ~ order + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

mod1_log_n <- glmer(
  hemoplasma ~ log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

anova(mod1_null, mod1_order, test = "Chisq")
anova(mod1_null, mod1_log_n, test = "Chisq")

aics <- AIC(mod1_null, mod1_order, mod1_log_n)
aic_null <- aics["mod1_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
aics[, c("AIC", "delta_AIC_vs_null")]
```

Results : 
```
> anova(mod1_null, mod1_order, test="Chisq")
Data: data_hemoplasma_stat
Models:
mod1_null: hemoplasma ~ 1 + (1 | species)
mod1_order: hemoplasma ~ order + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
mod1_null     2 416.38 425.22 -206.19    412.38                       
mod1_order    7 416.19 447.13 -201.10    402.19 10.187  5    0.07011 .

> anova(mod1_null, mod1_log_n, test="Chisq")
Data: data_hemoplasma_stat
Models:
mod1_null: hemoplasma ~ 1 + (1 | species)
mod1_log_n: hemoplasma ~ log_n + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod1_null     2 416.38 425.22 -206.19    412.38                     
mod1_log_n    3 416.94 430.20 -205.47    410.94 1.4399  1     0.2302

                AIC delta_AIC_vs_null
mod1_null  416.3804         0.0000000
mod1_order 416.1937        -0.1866964
mod1_log_n 416.9406         0.5601158
```

### Interpretation (model comparison)
The inclusion of `order` slightly improved model fit compared to the null model, although this effect was not statistically significant (_p_ = 0.07), suggesting a weak signal of taxonomic structure in infection probability.

In contrast, `log_n` did not improve model fit (_p_ = 0.23), indicating no detectable effect of sampling effort.

AIC comparisons supported these results, with minimal differences between models (ΔAIC < 1), indicating no strong support for any predictor over the null model.

Overall, results suggest a weak but consistent tendency for variation in hemoplasma infection across mammalian orders.

### Post-hoc analysis of differences between mammalian orders (model-based pairwise comparisons)
We perform post-hoc comparisons to identify which orders differ in hemoplasma infection probability.
```
emm_order <- emmeans(mod1_full, pairwise ~ order, type = "response")
emm_order
```

Results:
```
$emmeans
 order             prob     SE  df asymp.LCL asymp.UCL
 Carnivora       0.7165 0.2960 Inf    0.1267     0.978
 Cingulata       0.2916 0.3570 Inf    0.0138     0.924
 Didelphimorphia 0.1668 0.1250 Inf    0.0334     0.537
 Pilosa          0.1183 0.1300 Inf    0.0115     0.607
 Primates        0.8773 0.1400 Inf    0.3576     0.989
 Rodentia        0.0518 0.0373 Inf    0.0121     0.195
Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 

$contrasts
 contrast                    odds.ratio       SE  df null z.ratio p.value
 Carnivora / Cingulata           6.1399  12.7000 Inf    1   0.878  0.9520
 Carnivora / Didelphimorphia    12.6271  19.6000 Inf    1   1.637  0.5741
 Carnivora / Pilosa             18.8321  38.3000 Inf    1   1.443  0.7005
 Carnivora / Primates            0.3535   0.6030 Inf    1  -0.610  0.9904
 Carnivora / Rodentia           46.2877  65.4000 Inf    1   2.715  0.0722
 Cingulata / Didelphimorphia     2.0566   3.8300 Inf    1   0.387  0.9989
 Cingulata / Pilosa              3.0671   6.6300 Inf    1   0.518  0.9955
 Cingulata / Primates            0.0576   0.1170 Inf    1  -1.402  0.7260
 Cingulata / Rodentia            7.5388  13.4000 Inf    1   1.138  0.8656
 Didelphimorphia / Pilosa        1.4914   2.3300 Inf    1   0.256  0.9999
 Didelphimorphia / Primates      0.0280   0.0410 Inf    1  -2.441  0.1422
 Didelphimorphia / Rodentia      3.6658   3.8600 Inf    1   1.233  0.8206
 Pilosa / Primates               0.0188   0.0346 Inf    1  -2.158  0.2577
 Pilosa / Rodentia               2.4579   3.6800 Inf    1   0.600  0.9911
 Primates / Rodentia           130.9573 177.0000 Inf    1   3.616  0.0041
P value adjustment: tukey method for comparing a family of 6 estimates 
Tests are performed on the log odds ratio scale 
```

### Interpretation
Post-hoc pairwise comparisons showed variation in hemoplasma infection probability among mammalian orders, but most contrasts were not significant after Tukey correction.

A significant difference was found between Primates and Rodentia (_p_ = 0.0041), suggesting higher infection probabilities in Primates compared to Rodents, while all other comparisons were non-significant.

### Create a plot of hemoplasma prevalence by species and mammalian order
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

## Step 5. Variation in hemoplasma infection across sex within species where infection was detected (GLMM model 2)

### Data preparation
We restricted here the analysis to mammalian species with at least one infected individual, in order to test sex effects within relevant host species.
```
data_sex <- data_hemoplasma_stat[
  complete.cases(data_hemoplasma_stat[, c("hemoplasma", "sex", "species")]),
]

species_infected_sex <- data_sex %>%
  group_by(species) %>%
  summarise(infected = any(hemoplasma == 1, na.rm = TRUE)) %>%
  filter(infected) %>%
  pull(species)

data_inf_sex <- data_sex %>%
  filter(species %in% species_infected_sex) %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma))) %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup()
```

### GLMM (Model 2)
This model tests whether hemoplasma infection probability differs between sexes (`sex`) while controlling for sampling effort (`log_n`) and accounting for species-level random effects (`1 | species`).
```
mod2_full <- glmer(
  hemoplasma ~ sex + log_n + (1 | species),
  family = binomial,
  data = data_inf_sex,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)
```

### Model term significance testing
Model terms were evaluated using likelihood ratio tests via single-term deletions (drop1 function with Chi-square tests).
```
res2 <- drop1(mod2_full, test = "Chisq")
res2
```

Results :
```
Single term deletions
Model:
hemoplasma ~ sex + log_n + (1 | species)
       npar    AIC     LRT Pr(Chi)  
<none>      279.07                
sex       1 279.59 2.52007  0.1124
log_n     1 277.51 0.44085  0.5067
```

### Interpretation
No significant effect of `sex` on `hemoplasma` infection probability was detected (χ² = 2.52, _p_ = 0.11).

Sampling effort (`log_n`) had no detectable effect (p = 0.51).

### Model comparison with null and univariate models
```
mod2_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_inf_sex,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

mod2_sex <- glmer(
  hemoplasma ~ sex + (1 | species),
  family = binomial,
  data = data_inf_sex,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

mod2_log_n <- glmer(
  hemoplasma ~ log_n + (1 | species),
  family = binomial,
  data = data_inf_sex,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

anova(mod2_null, mod2_sex, test = "Chisq")
anova(mod2_null, mod2_log_n, test = "Chisq")

aics <- AIC(mod2_null, mod2_sex, mod2_log_n)
aic_null <- aics["mod2_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
aics[, c("AIC", "delta_AIC_vs_null")]
```

Results :
```
> anova(mod2_null, mod2_sex, test = "Chisq")
Data: data_inf_sex
Models:
mod2_null: hemoplasma ~ 1 + (1 | species)
mod2_sex: hemoplasma ~ sex + (1 | species)
          npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod2_null    2 278.07 285.65 -137.03    274.07                     
mod2_sex     3 277.51 288.89 -135.75    271.51 2.5598  1     0.1096

> anova(mod2_null, mod2_log_n, test = "Chisq")
Data: data_inf_sex
Models:
mod2_null: hemoplasma ~ 1 + (1 | species)
mod2_log_n: hemoplasma ~ log_n + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod2_null     2 278.07 285.65 -137.03    274.07                     
mod2_log_n    3 279.59 290.97 -136.79    273.59 0.4806  1     0.4881

                AIC delta_AIC_vs_null
mod2_null  278.0684         0.0000000
mod2_sex   277.5086        -0.5598415
mod2_log_n 279.5878         1.5193843
```

### Interpretation (model comparison)
The inclusion of `sex` or `log_n` did not improve model fit compared to the null model (_p_ = 0.11 and 0.49, respectively).

AIC comparison supported these results, with only a minimal improvement for the model including `sex` (ΔAIC < 1), indicating weak support for this predictor. Overall, results provide no clear evidence for sex-biased hemoplasma infection.

### Post-hoc analysis of differences between sex in infected species (model-based pairwise comparisons)

We perform post-hoc comparisons to identify if sexes differ in hemoplasma infection probability.

```
emm_sex <- emmeans(mod2_full, pairwise ~ sex, type = "response")
emm_sex
```

Results:
```
$emmeans
 sex  prob    SE  df asymp.LCL asymp.UCL
 F   0.371 0.202 Inf    0.0972     0.764
 M   0.506 0.214 Inf    0.1603     0.846
Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 

$contrasts
 contrast odds.ratio    SE  df null z.ratio p.value
 F / M         0.576 0.199 Inf    1  -1.594  0.1109
```

### Interpretation
Predicted infection probabilities were slightly higher in males than in females, but this difference was not statistically significant (_p_ = 0.11).

### Species-level sex bias in hemoplasma infection (Fisher exact tests)
To test whether `hemoplasma` infection differs between males and females within host `species`, we restricted analyses to infected `species` that had at least 20 sexed individuals to ensure minimum statistical power.


### Data preparation
```
data_clean <- data_hemoplasma_stat %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma))) %>%
  filter(!is.na(sex))

species_keep <- data_clean %>%
  group_by(species) %>%
  summarise(
    n_sexed = n(),
    any_infected = any(hemoplasma == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_sexed >= 20, any_infected) %>%
  pull(species)

write.csv(species_keep,
          "species_used_fisher_sex_bias.csv",
          row.names = FALSE)

species_sex_summary <- data_clean %>%
  filter(species %in% species_keep) %>%
  group_by(species, sex) %>%
  summarise(
    n_total = n(),
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),
    n_uninfected = sum(hemoplasma == 0, na.rm = TRUE),
    prevalence = n_infected / n_total,
    .groups = "drop"
  ) %>%
  arrange(species, sex)
species_sex_summary %>%
  arrange(species, sex)
```

Results for infected `species` that had at least 20 sexed individuals :
```
# A tibble: 8 × 6
  species               sex   n_total n_infected n_uninfected prevalence
  <fct>                 <fct>   <int>      <int>        <int>      <dbl>
1 Bradypus_tridactylus  F          43          2           41     0.0465
2 Bradypus_tridactylus  M          49          2           47     0.0408
3 Choloepus_didactylus  F          49         39           10     0.796 
4 Choloepus_didactylus  M          34         29            5     0.853 
5 Didelphis_marsupialis F          19          8           10     0.421 
6 Didelphis_marsupialis M          19         13            5     0.684 
7 Saguinus_midas        F          15         15            0     1     
8 Saguinus_midas        M          23         23            0     1     
```

### Fisher exact tests per species
```
fisher_results <- data_fisher %>%
  group_by(species) %>%
  summarise(
    tab = list(table(sex, hemoplasma)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    fisher = list(
      if (all(dim(tab) == c(2, 2))) {
        fisher.test(tab)
      } else {
        NULL
      }
    ),
    p_value = if (!is.null(fisher)) fisher$p.value else NA_real_,
    odds_ratio = if (!is.null(fisher)) as.numeric(fisher$estimate) else NA_real_
  ) %>%
  ungroup() %>%
  select(species, p_value, odds_ratio)
```

### Multiple testing correction (FDR)
```
fisher_results <- fisher_results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr")
  ) %>%
  arrange(p_adj)

fisher_results
```

Results : 
```
# A tibble: 4 × 4
  species               p_value odds_ratio  p_adj
  <fct>                   <dbl>      <dbl>  <dbl>
1 Didelphis_marsupialis   0.176      3.14   0.527
2 Choloepus_didactylus    0.573      1.48   0.860
3 Bradypus_tridactylus    1          0.874  1    
4 Saguinus_midas         NA         NA     NA    
```

### Interpretation
Fisher exact tests performed separately for each species (restricted to species with ≥20 sexed individuals and at least one infected individual) revealed no significant sex differences in hemoplasma infection after FDR correction. For _Saguinus midas_, the test could not be computed due to insufficient data structure in the contingency table (ie, all individuals are infected).

### Create a plot of hemoplasma prevalence by infected species
```
species_infected <- data_hemoplasma_stat %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma))) %>%
  group_by(species) %>%
  summarise(any_infected = any(hemoplasma == 1, na.rm = TRUE)) %>%
  filter(any_infected) %>%
  pull(species)

data_clean <- data_hemoplasma_stat %>%
  filter(
    species %in% species_infected,
    !is.na(sex)
  ) %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma)))

species_keep <- data_clean %>%
  group_by(species) %>%
  summarise(n_sexed = n()) %>%
  filter(n_sexed >= 20) %>%
  pull(species)

df_plot <- data_clean %>%
  filter(species %in% species_keep) %>%
  group_by(species, sex) %>%
  summarise(
    n = n(),
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  ) %>%
  mutate(
    species_sex = paste(species, sex)
  )

species_order <- unique(df_plot$species)

levels_ordered <- unlist(lapply(species_order, function(sp) {
  c(paste(sp, "M"),
    paste(sp, "F"),
    paste0("gap_", sp))
}))

df_plot$species_sex <- factor(df_plot$species_sex,
                              levels = levels_ordered)

species_cols <- setNames(hue_pal()(length(species_order)), species_order)

fill_cols <- c()
for (sp in species_order) {
  fill_cols[paste(sp, "M")] <- species_cols[sp]
  fill_cols[paste(sp, "F")] <- alpha(species_cols[sp], 0.5)
}

df_plot_plot <- df_plot %>%
  mutate(
    species_sex_plot = ifelse(grepl("gap_", species_sex), NA, species_sex)
  )

p_sex_species <- ggplot(df_plot_plot,
                        aes(x = species_sex_plot,
                            y = prevalence,
                            fill = species_sex)) +
  
  geom_bar(stat = "identity",
           color = "black",
           na.rm = TRUE) +
  
  scale_fill_manual(values = fill_cols,
                    na.translate = FALSE) +
  
  scale_x_discrete(drop = FALSE) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank()
  ) +
  
  labs(
    x = "Species / Sex",
    y = "Hemoplasma infection prevalence",
    fill = "Species & Sex"
  )

p_sex_species

pdf("Fig_1C_hemoplasma_sex_species.pdf", width = 12, height = 5)
print(p_sex_species)
dev.off()
```

## Step 6. Variation in hemoplasma infection across host ecological traits (GLMM Model 3) 

### Data preparation
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup()
```

### GLMM (Model 3)

This model tests whether hemoplasma infection probability varies according to host ecological traits (`vertical_stratum`, `activity`, `diet`, `sociality`), while controlling for sampling effort (`log_n`) and accounting for species-level random effects (`1 | species`).
```
mod3_full <- glmer(
  hemoplasma ~ vertical_stratum + activity + diet + sociality + log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)
```

### Model term significance testing
Model terms were evaluated using likelihood ratio tests via single-term deletions (drop1 function with Chi-square tests).
```
res3 <- drop1(mod3_full, test = "Chisq")
res3
```

Results :
```
Single term deletions
Model:
hemoplasma ~ vertical_stratum + activity + diet + sociality + log_n + (1 | species)
                 npar    AIC    LRT Pr(Chi)  
<none>                421.76                 
vertical_stratum    2 418.86 1.0977 0.57760  
activity            1 423.09 3.3279 0.06812 .
diet                3 416.11 0.3554 0.94930  
sociality           1 419.80 0.0445 0.83297
log_n               1 423.37 3.6160 0.05723 .  
```

### Interpretation
No significant effect of host ecological traits was detected on hemoplasma infection probability. Activity and sampling effort showed weak non-significant trends but did not reach statistical significance.


### Model comparison with null and univariate models
```

mod3_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_vertical_stratum <- glmer(
  hemoplasma ~ vertical_stratum + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_activity <- glmer(
  hemoplasma ~ activity + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_diet <- glmer(
  hemoplasma ~ diet + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_sociality <- glmer(
  hemoplasma ~ sociality + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_log_n <- glmer(
  hemoplasma ~ log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

anova(mod3_null, mod3_vertical_stratum, test = "Chisq")
anova(mod3_null, mod3_activity, test = "Chisq")
anova(mod3_null, mod3_diet, test = "Chisq")
anova(mod3_null, mod3_sociality, test = "Chisq")
anova(mod3_null, mod3_log_n, test = "Chisq")

aics <- AIC(mod3_null,
            mod3_vertical_stratum,
            mod3_activity,
            mod3_diet,
            mod3_sociality,
            mod3_log_n)

aic_null <- aics["mod3_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
aics[, c("AIC", "delta_AIC_vs_null")]
```

Results : 
```
> anova(mod3_null, mod3_vertical_stratum, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_vertical_stratum: hemoplasma ~ vertical_stratum + (1 | species)
                      npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod3_null                2 416.38 425.22 -206.19    412.38                     
mod3_vertical_stratum    4 419.55 437.23 -205.78    411.55 0.8301  2     0.6603

> anova(mod3_null, mod3_activity, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_activity: hemoplasma ~ activity + (1 | species)
              npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
mod3_null        2 416.38 425.22 -206.19    412.38                       
mod3_activity    3 413.14 426.40 -203.57    407.14 5.2439  1    0.02202 *

> anova(mod3_null, mod3_diet, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_diet: hemoplasma ~ diet + (1 | species)
          npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod3_null    2 416.38 425.22 -206.19    412.38                     
mod3_diet    5 419.44 441.54 -204.72    409.44 2.9409  3     0.4008

> anova(mod3_null, mod3_sociality, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_sociality: hemoplasma ~ sociality + (1 | species)
               npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod3_null         2 416.38 425.22 -206.19    412.38                     
mod3_sociality    3 417.58 430.84 -205.79    411.58 0.8027  1     0.3703

> anova(mod3_null, mod3_log_n, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_log_n: hemoplasma ~ log_n + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod3_null     2 416.38 425.22 -206.19    412.38                     
mod3_log_n    3 416.94 430.20 -205.47    410.94 1.4399  1     0.2302

                           AIC delta_AIC_vs_null
mod3_null             416.3804         0.0000000
mod3_vertical_stratum 419.5503         3.1698637
mod3_activity         413.1365        -3.2439419
mod3_diet             419.4395         3.0590596
mod3_sociality        417.5777         1.1973107
mod3_log_n            416.9406         0.5601158
```

### Interpretation
Activity showed a marginal effect in the univariate model comparison (_p_ = 0.022), but this pattern was not robust in the full multivariate model (see above).

### Marginal means (estimated infection probabilities)
We estimated marginal means (back-transformed to response scale) for each level of host ecological traits.

```
emm_vertical_stratum <- emmeans(mod3_full, ~ vertical_stratum, type = "response")
emm_activity <- emmeans(mod3_full, ~ activity, type = "response")
emm_diet <- emmeans(mod3_full, ~ diet, type = "response")
emm_sociality <- emmeans(mod3_full, ~ sociality, type = "response")

emm_vertical_stratum
emm_activity
emm_diet
emm_sociality
```

Results :
```
 vertical_stratum   prob    SE  df asymp.LCL asymp.UCL
 Canopy           0.4112 0.298 Inf   0.05898     0.886
 Ground           0.4723 0.261 Inf   0.10315     0.874
 Mixed            0.0876 0.177 Inf   0.00126     0.879

 activity    prob     SE  df asymp.LCL asymp.UCL
 Diurnal   0.6517 0.3250 Inf   0.10161     0.969
 Nocturnal 0.0757 0.0919 Inf   0.00621     0.518

 diet         prob    SE  df asymp.LCL asymp.UCL
 Carnivore   0.469 0.583 Inf   0.00888     0.989
 Insectivore 0.199 0.317 Inf   0.00498     0.925
 Omnivore    0.212 0.172 Inf   0.03459     0.670
 Phytophage  0.284 0.272 Inf   0.02822     0.845

 sociality  prob    SE  df asymp.LCL asymp.UCL
 Group     0.312 0.307 Inf    0.0268     0.882
 Solitary  0.252 0.226 Inf    0.0313     0.779

Results are averaged over the levels of: vertical_stratum, activity, diet 
Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 
```

### Post-hoc analysis of differences between levelS of host ecological traits (model-based pairwise comparisons)
```
pairs(emmeans(mod3_full, ~ vertical_stratum), adjust = "tukey")
pairs(emmeans(mod3_full, ~ activity), adjust = "tukey")
pairs(emmeans(mod3_full, ~ diet), adjust = "tukey")
pairs(emmeans(mod3_full, ~ sociality), adjust = "tukey")
```

Results:
```
> pairs(emmeans(mod3_full, ~ vertical_stratum), adjust = "tukey")
 contrast        estimate   SE  df z.ratio p.value
 Canopy - Ground   -0.248 1.13 Inf  -0.219  0.9739
 Canopy - Mixed     1.984 2.28 Inf   0.869  0.6596
 Ground - Mixed     2.233 2.20 Inf   1.017  0.5661

> pairs(emmeans(mod3_full, ~ activity), adjust = "tukey")
 contrast            estimate   SE  df z.ratio p.value
 Diurnal - Nocturnal     3.13 1.61 Inf   1.946  0.0517

> pairs(emmeans(mod3_full, ~ diet), adjust = "tukey")
 contrast                 estimate   SE  df z.ratio p.value
 Carnivore - Insectivore    1.2706 2.79 Inf   0.455  0.9687
 Carnivore - Omnivore       1.1872 2.42 Inf   0.490  0.9614
 Carnivore - Phytophage     0.7991 2.63 Inf   0.304  0.9903
 Insectivore - Omnivore    -0.0834 1.76 Inf  -0.047  1.0000
 Insectivore - Phytophage  -0.4716 2.06 Inf  -0.229  0.9958
 Omnivore - Phytophage     -0.3881 1.18 Inf  -0.328  0.9878

> pairs(emmeans(mod3_full, ~ sociality), adjust = "tukey")
 contrast         estimate   SE  df z.ratio p.value
 Group - Solitary    0.298 1.41 Inf   0.211  0.8328
```

### Interpretation
Marginal means showed broadly overlapping confidence intervals across all ecological trait categories. Tukey-adjusted pairwise comparisons revealed no significant differences among levels of vertical stratum, diet, or sociality, and only a marginal diurnal–nocturnal contrast (_p_ = 0.052).

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

    geom_jitter(
      data = df,
      aes(x = .data[[var]], y = prevalence,
          color = .data[[var]], size = n),
      width = 0.20,
      alpha = 0.65,
      na.rm = TRUE
    ) +

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
pdf(file = file.path(getwd(), "Figure2_Hemoplasma_ecological_traits.pdf"),
    width = 12, height = 9, useDingbats = FALSE)
print(figure1)
dev.off()
```

## Step 7. Association between hemoplasma infection and co-infecting blood parasites (GLMM Model 4)
### Data preparation
We tested whether `hemoplasma` infection probability is associated with co-infection by other blood parasites, including `anaplasmataceae`, `apicomplexa`, `trypanosoma`, and `filaria`. Species-level sampling effort was accounted for using `log_n`, and species identity was included as a random effect.
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup()

data_coinf <- data_hemoplasma_stat %>%
  filter(
    !is.na(hemoplasma),
    !is.na(anaplasmataceae),
    !is.na(apicomplexa),
    !is.na(trypanosoma),
    !is.na(filaria)
  )
```

## GLMM (Model 4)
This model tests whether `hemoplasma` infection probability is associated with co-infection by other blood parasites.
```
mod4_full <- glmer(
  hemoplasma ~ anaplasmataceae + apicomplexa + trypanosoma + filaria +
    log_n + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)
```

## Model term significance testing

Model terms were evaluated using likelihood ratio tests via single-term deletions.
```
res4 <- drop1(mod4_full, test = "Chisq")
res4
```

Results :
```
Single term deletions
Model:
hemoplasma ~ anaplasmataceae + apicomplexa + trypanosoma + filaria + 
    log_n + (1 | species)
                npar    AIC    LRT  Pr(Chi)   
<none>               120.13                   
anaplasmataceae    1 122.84 4.7040 0.030093 * 
apicomplexa        1 118.14 0.0109 0.916691   
trypanosoma        1 120.53 2.3976 0.121525   
filaria            1 118.14 0.0042 0.948289   
log_n              1 127.77 9.6327 0.001911 **
```

## Interpretation
`anaplasmataceae` infection was significantly associated with `hemoplasma` infection probability (χ² = 4.70, _p_ = 0.030), with higher infection probability in co-infected hosts. Sampling effort also had a strong effect (χ² = 9.63, _p_ = 0.0019), indicating lower detection probability with increasing sample size. No significant effects were detected for `apicomplexa`, `trypanosoma`, or `filaria` (all _p_ > 0.10).

## Model comparison with null and univariate models
```
mod4_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod4_anaplasma <- glmer(
  hemoplasma ~ anaplasmataceae + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

mod4_apicomplexa <- glmer(
  hemoplasma ~ apicomplexa + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

mod4_trypanosoma <- glmer(
  hemoplasma ~ trypanosoma + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

mod4_filaria <- glmer(
  hemoplasma ~ filaria + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

anova(mod4_null, mod4_anaplasma, test = "Chisq")
anova(mod4_null, mod4_apicomplexa, test = "Chisq")
anova(mod4_null, mod4_trypanosoma, test = "Chisq")
anova(mod4_null, mod4_filaria, test = "Chisq")

aics <- AIC(mod4_null,
            mod4_anaplasma,
            mod4_apicomplexa,
            mod4_trypanosoma,
            mod4_filaria)

aic_null <- aics["mod4_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
aics[, c("AIC", "delta_AIC_vs_null")]
```

Results : 
```
> anova(mod4_null, mod4_anaplasma, test = "Chisq")
Data: data_coinf
Models:
mod4_null: hemoplasma ~ 1 + (1 | species)
mod4_anaplasma: hemoplasma ~ anaplasmataceae + (1 | species)
               npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
mod4_null         2 125.55 131.92 -60.777    121.55                       
mod4_anaplasma    3 124.09 133.63 -59.044    118.09 3.4647  1    0.06269 .

> anova(mod4_null, mod4_apicomplexa, test = "Chisq")
Data: data_coinf
Models:
mod4_null: hemoplasma ~ 1 + (1 | species)
mod4_apicomplexa: hemoplasma ~ apicomplexa + (1 | species)
                 npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod4_null           2 125.55 131.92 -60.777    121.55                     
mod4_apicomplexa    3 127.47 137.02 -60.736    121.47 0.0816  1     0.7751

> anova(mod4_null, mod4_trypanosoma, test = "Chisq")
Data: data_coinf
Models:
mod4_null: hemoplasma ~ 1 + (1 | species)
mod4_trypanosoma: hemoplasma ~ trypanosoma + (1 | species)
                 npar    AIC    BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)
mod4_null           2 125.55 131.92 -60.777    121.55                    
mod4_trypanosoma    3 126.22 135.77 -60.112    120.22 1.329  1      0.249

> anova(mod4_null, mod4_filaria, test = "Chisq")
Data: data_coinf
Models:
mod4_null: hemoplasma ~ 1 + (1 | species)
mod4_filaria: hemoplasma ~ filaria + (1 | species)
             npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod4_null       2 125.55 131.92 -60.777    121.55                     
mod4_filaria    3 127.55 137.10 -60.776    121.55 0.0018  1     0.9662

                      AIC delta_AIC_vs_null
mod4_null        125.5536         0.0000000
mod4_anaplasma   124.0889        -1.4646561
mod4_apicomplexa 127.4719         1.9183628
mod4_trypanosoma 126.2246         0.6710417
mod4_filaria     127.5518         1.9981994
```

## Interpretation
Model comparisons confirmed a weak but consistent support for an association between `hemoplasma` infection and `anaplasmataceae` co-infection (ΔAIC ≈ −1.46), whereas no improvement in model fit was observed for `apicomplexa`, `trypanosoma`, or `filaria`.

## Post-hoc analysis of differences between blood parasites (model-based pairwise comparisons)
```
emm_anaplasma <- emmeans(mod4_full, ~ anaplasmataceae, type = "response")
emm_apicomplexa <- emmeans(mod4_full, ~ apicomplexa, type = "response")
emm_trypanosoma <- emmeans(mod4_full, ~ trypanosoma, type = "response")
emm_filaria <- emmeans(mod4_full, ~ filaria, type = "response")

emm_anaplasma
pairs(emmeans(mod4_full, ~ anaplasmataceae), adjust = "tukey")
emm_apicomplexa
pairs(emmeans(mod4_full, ~ apicomplexa), adjust = "tukey")
emm_trypanosoma
pairs(emmeans(mod4_full, ~ trypanosoma), adjust = "tukey")
emm_filaria
pairs(emmeans(mod4_full, ~ filaria), adjust = "tukey")
```

Results : 
```
 anaplasmataceae   prob     SE  df asymp.LCL asymp.UCL
 0               0.0635 0.0688 Inf   0.00695     0.396
 1               0.1878 0.1480 Inf   0.03345     0.607

 contrast                            estimate   SE  df z.ratio p.value
 anaplasmataceae0 - anaplasmataceae1    -1.23 0.61 Inf  -2.013  0.0442

 apicomplexa  prob     SE  df asymp.LCL asymp.UCL
 0           0.117 0.0844 Inf    0.0263     0.396
 1           0.105 0.1370 Inf    0.0068     0.669

 contrast                    estimate   SE  df z.ratio p.value
 apicomplexa0 - apicomplexa1    0.123 1.16 Inf   0.106  0.9158

 trypanosoma   prob     SE  df asymp.LCL asymp.UCL
 0           0.2947 0.1350 Inf   0.10506     0.598
 1           0.0361 0.0583 Inf   0.00141     0.499

 contrast                    estimate  SE  df z.ratio p.value
 trypanosoma0 - trypanosoma1     2.41 1.5 Inf   1.608  0.1079

 filaria  prob    SE  df asymp.LCL asymp.UCL
 0       0.113 0.100 Inf    0.0178     0.475
 1       0.109 0.111 Inf    0.0128     0.537

 contrast            estimate    SE  df z.ratio p.value
 filaria0 - filaria1   0.0426 0.657 Inf   0.065  0.9482

Results are averaged over the levels of: anaplasmataceae, apicomplexa, trypanosoma, log_n 
Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 
```

## Interpretation
`hemoplasma` infection probability was higher in hosts infected by `anaplasmataceae` (0.19 vs 0.06; _p_ = 0.044), while no significant differences were detected for `apicomplexa`, `trypanosoma`, or `filaria` (all _p _> 0.10), with strongly overlapping confidence intervals across groups.

## Create a plot of hemoplasma prevalence by co-infecting blood parasites
A FAIRE !!!



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

