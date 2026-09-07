# **Hemoplasma veterinary evaluations : R command lines and script**

We analyzed data from 175 wild sloths captured between 1994 and 1995 during the flooding of the Petit Saut Dam (5°03′43″ N, 53°03′00″ O) on the Sinnamary River (French Guiana, South America). The clinical data include the following variables for each examined sloth: 
- `species` : Sloth species (Bt: *Bradypus tridactylus*; Cd: *Choloepus didactylus*)
- `sex` : Sex of the sloth (F: female; M: male)
- `age_class` : Age category (A: adult; J: juvenile)
- `season` : Season of capture (W: wet; D: dry)
- `weight` : Body weight (quantitative variable, in kg)
- `total_length` : Total body length (quantitative variable, in cm)
- `wither_height` : Height at the withers (quantitative variable, in cm)
- `neck_size` : Neck circumference (quantitative variable, in cm)
- `temperature` : Body temperature (quantitative variable, in °C)
- `hematocrit` : Hematocrit level (quantitative variable, in %)
- `health_condition` : Overall health status (G: good; D: deteriorated)
- `female_reproductive_status` : Reproductive state of female individuals, indicating whether they were reproductively active or inactive at the time of sampling (Female non pregnant non lactating / Pregnant female / Female lactating with a young)
- `hemoplasma` : Infection status with hemotropic mycoplasmas (0 = uninfected; 1 = infected)
- `anaplasmataceae` : Infection status with bacteria of the Anaplasmataceae family (here, *Anaplasma amazonensis*) (0 = uninfected; 1 = infected)
- `apicomplexa` : Infection status with piroplasmids (here, *Babesia* sp.) (0 = uninfected; 1 = infected)
  
Details about all the experimental methods and measures are available in the related manuscript.


## Table of contents
- [Step 1. Retrieving the data](#step-1-retrieving-the-data)
- [Step 2. Prepare the data for analysis](#step-2-prepare-the-data-for-analysis)
- [Step 3. Calculate hemoplasma infection prevalence](#step-3-calculate-hemoplasma-infection-prevalence)
- [Step 4. Create the pathogens variable and sloth species subset](#step-4-create-the-pathogens-variable-and-sloth-species-subset)
- [Step 5. Impact of hemoplasma infections on Scale Mass Index (SMI) in adult Bt](#step-5-impact-of-hemoplasma-infections-on-scale-mass-index-smi-in-adult-bt)
- [Step 6. Impact of hemoplasma infections on Scale Mass Index (SMI) in adult Cd](#step-6-impact-of-hemoplasma-infections-on-scale-mass-index-smi-in-adult-cd)
- [Step 7. Impact of hemoplasma infections on neck circumference](#step-7-impact-of-hemoplasma-infections-on-neck-circumference)
- [Step 8. Impact of hemoplasma infections on hematocrit levels](#step-8-impact-of-hemoplasma-infections-on-hematocrit-levels)
- [Step 9. Impact of hemoplasma infections on body temperature](#step-9-impact-of-hemoplasma-infections-on-body-temperature)
- [Step 10. Impact of hemoplasma infections on general health_condition](#step-10-impact-of-hemoplasma-infections-on-general-health_condition)
- [Step 11. Impact of hemoplasma infections on female_reproductive_status](#step-11-impact-of-hemoplasma-infections-on-female_reproductive_status)

## Step 1. Retrieving the data

All veterinary clinical data for the two sloth species are available [here](https://github.com/olivierduron/Hemoplasma_infections/blob/main/data_hemoplasma_sloth.csv).

This database will be referred to as `data_hemoplasma` throughout the R command lines and scripts provided below. It corresponds to the dataset provided in Table S1 of the related manuscript.

Load the dataset directly from the GitHub repository to R
```
data_hemoplasma <- read.csv2(
  "https://raw.githubusercontent.com/olivierduron/Hemoplasma_infections/main/data_hemoplasma_sloth.csv",
  na.strings = c("NA", "")
)
data_hemoplasma
```


## Step 2. Prepare the data for analysis

### Convert categorical variables into factors
```
data_hemoplasma$species        <- as.factor(data_hemoplasma$species)
data_hemoplasma$season         <- as.factor(data_hemoplasma$season)
data_hemoplasma$sex            <- as.factor(data_hemoplasma$sex)
data_hemoplasma$age            <- as.factor(data_hemoplasma$age)
data_hemoplasma$hemoplasma      <- as.factor(data_hemoplasma$hemoplasma)
data_hemoplasma$anaplasmataceae      <- as.factor(data_hemoplasma$anaplasmataceae)
data_hemoplasma$apicomplexa       <- as.factor(data_hemoplasma$apicomplexa)
```

### Load libraries for analysis
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
library(emmeans)
```

## Step 3. Calculate `hemoplasma` infection prevalence
### Calculate `hemoplasma` infection prevalence and 95% confidence interval for _Bradypus tridactylus_ (Bt) and _Choloepus didactylus_ (Cd)

```
prevalence_results <- data_hemoplasma %>% group_by(species) %>% summarise(n = n(), positives = sum(hemoplasma == 1), prevalence = positives / n, conf_low = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$lower, conf_high = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$upper)
print(prevalence_results)
```

Results :
| species | n | positives | prevalence | conf_low | conf_high |
|:--------|--:|----------:|-----------:|---------:|----------:|
| Bt      | 92 | 4  | 0.0435 | 0.0120 | 0.108 |
| Cd      | 83 | 68 | 0.8190 | 0.7200 | 0.895 |

### Test if `hemoplasma` is influenced by sloth `species`:
```
fisher.test(table(data_hemoplasma$hemoplasma, data_hemoplasma$species))
```

-> Results : `hemoplasma` prevalence differed strongly between the two sloth `species`, from 4.3% (4/92; 95% CI: 1.2–10.8%) in *Bradypus tridactylus* to 81.9% (68/83; 95% CI: 72.0–89.5%) in *Choloepus didactylus* (χ²₁ = 105.27, p < 2.2 × 10⁻¹⁶).

-> Interpretation : This strong interspecific difference suggests that host species and associated ecological or evolutionary traits may constrain `hemoplasma` infection.

## Step 4. Create the `pathogens` variable and sloth `species` subset

### Create the `pathogens` variable by merging `anaplasmataceae` and `apicomplexa` (0 = uninfected ; 1 = infected by `anaplasmataceae` and/or `apicomplexa`)
```
data_hemoplasma <- data_hemoplasma %>%
  mutate(
    pathogens = ifelse(
      anaplasmataceae == 1 | apicomplexa == 1,
      1, 0
    ),
    species = factor(species)
  )
data_hemoplasma$pathogens <- as.factor(data_hemoplasma$pathogens)
```

### Convert variables:
```
data_hemoplasma$weight <- as.numeric(data_hemoplasma$weight)
data_hemoplasma$total_length <- as.numeric(data_hemoplasma$total_length)
data_hemoplasma$wither_height <- as.numeric(data_hemoplasma$wither_height)
data_hemoplasma$neck_size <- as.numeric(data_hemoplasma$neck_size)
data_hemoplasma$hematocrit <- as.numeric(data_hemoplasma$hematocrit)
```

### Create a subset `data_Bt` containing only records for _Bradypus tridactylus_ (Bt)
```
data_Bt <- subset(data_hemoplasma, species == "Bt")
table(data_Bt$hemoplasma)
```

### Create a subset `data_Cd` containing only records for _Choloepus didactylus_ (Cd):
```
data_Cd <- subset(data_hemoplasma, species == "Cd")
table(data_Cd$hemoplasma)
```

## Step 5. Impact of `hemoplasma` infections on Scale Mass Index (SMI) in adult Bt
The Scaled Mass Index (SMI) was used as a body condition indicator that standardizes individual `weight` to `body_length`, using an allometric scaling relationship. SMI was calculated following Peig & Green (2009) (https://doi.org/10.1111/j.1600-0706.2009.17643.x).

### Function to calculate SMI for adult Bt
```
data_adult_Bt <- subset(data_Bt, age == "A")
sma_model_Bt <- sma(log(weight) ~ log(total_length), data = data_adult_Bt)
b <- coef(sma_model_Bt)[2]
L0 <- mean(data_adult_Bt$total_length, na.rm = TRUE)
data_adult_Bt$SMI <- data_adult_Bt$weight * (L0 / data_adult_Bt$total_length)^b
```

### Fit a GLM to test whether SMI is influenced by `hemoplasma`, `pathogens`, `sex` and `season` in Bt
```
model_SMIBt_full <- glm(
  SMI ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)
summary(model_SMIBt_full)
model_SMIBt_3way <- glm(
  SMI ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)
summary(model_SMIBt_3way)
model_SMIBt_2way <- glm(
  SMI ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)
summary(model_SMIBt_2way)
model_SMIBt_add <- glm(
  SMI ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)
summary(model_SMIBt_add)
model_SMIBt_null <- glm(
  SMI ~ 1,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)
summary(model_SMIBt_null)
anova(
  model_SMIBt_null,
  model_SMIBt_add,
  test = "Chisq"
)
anova(
  model_SMIBt_add,
  model_SMIBt_2way,
  test = "Chisq"
)
anova(
  model_SMIBt_2way,
  model_SMIBt_3way,
  test = "Chisq"
)
anova(
  model_SMIBt_3way,
  model_SMIBt_full,
  test = "Chisq"
)
AIC_table <- AIC(
  model_SMIBt_null,
  model_SMIBt_add,
  model_SMIBt_2way,
  model_SMIBt_3way,
  model_SMIBt_full
)
AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table
drop1(
  model_SMIBt_add,
  test = "Chisq"
)
model_SMIBt_hemo <- glm(
  SMI ~ hemoplasma,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)
model_SMIBt_sex <- glm(
  SMI ~ sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)
model_SMIBt_season <- glm(
  SMI ~ season,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)
model_SMIBt_pathogens <- glm(
  SMI ~ pathogens,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)
AIC_models <- AIC(
  model_SMIBt_null,
  model_SMIBt_hemo,
  model_SMIBt_sex,
  model_SMIBt_season,
  model_SMIBt_pathogens
)
AIC_models$delta_AIC <- AIC_models$AIC - AIC(model_SMIBt_null)
AIC_models
par(mfrow = c(2, 2))
plot(model_SMIBt_add)
par(mfrow = c(1, 1))
shapiro.test(residuals(model_SMIBt_add))
plot(
  fitted(model_SMIBt_add),
  residuals(model_SMIBt_add),
  xlab = "Fitted SMI",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)
with(
  data_adult_Bt,
  table(hemoplasma, pathogens, season, sex)
)
```
-> Results : A Gaussian GLM showed that SMI was significantly associated with the additive effects of `hemoplasma` infection, other `pathogens` infection, `season` and `sex` (LRT: χ²₄ = 7.80, p < 0.001). The additive model was best supported (AIC = 141.18), with no improvement from adding two-way interactions (LRT: χ²₄ = 1.12, p = 0.437) or three-way interactions (LRT: χ²₁ = 0.31, p = 0.305). The four-way interaction model was not estimable beyond the three-way model because of data sparsity. In the final additive model, SMI was lower in `hemoplasma`-positive individuals (β = −0.96 ± 0.28 SE, LRT: χ²₁ = 11.39, p < 0.001) and higher in males (β = 0.46 ± 0.12 SE, LRT: χ²₁ = 14.24, p < 0.001), but was not associated with other `pathogens` (LRT: χ²₁ = 2.14, p = 0.144) or `season` (LRT: χ²₁ = 0.07, p = 0.798). Residuals were normally distributed (Shapiro–Wilk: W = 0.993, p = 0.952).

-> Interpretation : SMI was negatively associated with `hemoplasma` infection and positively associated with male `sex`, with no evidence for effects of other `pathogens`, `season` or interactions.


### Calculation of mean and standard error of SMI by `hemoplasma` infection status and `sex` for Bt
```
emmeans_hemoplasma <- emmeans(model_SMIBt_final, ~ hemoplasma)
emmeans_sex <- emmeans(model_SMIBt_final, ~ sex)
emmeans_sex_hemoplasma <- emmeans(model_SMIBt_final, ~ hemoplasma * sex)
summary(emmeans_hemoplasma, infer = c(TRUE, TRUE))
summary(emmeans_sex, infer = c(TRUE, TRUE))
summary(emmeans_sex_hemoplasma, infer = c(TRUE, TRUE))
pairs(emmeans_hemoplasma)
pairs(emmeans_sex)
emmeans_hemoplasma_df <- as.data.frame(summary(emmeans_hemoplasma, infer = c(TRUE, TRUE)))
emmeans_sex_df <- as.data.frame(summary(emmeans_sex, infer = c(TRUE, TRUE)))
emmeans_sex_hemoplasma_df <- as.data.frame(summary(emmeans_sex_hemoplasma, infer = c(TRUE, TRUE)))
emmeans_hemoplasma_df
emmeans_sex_df
emmeans_sex_hemoplasma_df
data_adult_Bt %>% 
  group_by(sex, hemoplasma) %>% 
  summarise(
    n = sum(!is.na(SMI)),
    mean = mean(SMI, na.rm = TRUE),
    se = sd(SMI, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) %>% 
  mutate(SMI = sprintf("%.2f ± %.2f", mean, se))
```

Results :
| sex | hemoplasma | n | mean | se | SMI |
|:----|:-----------|--:|-----:|----:|:----|
| F | 0 | 39 | 4.43 | 0.0749 | 4.43 ± 0.07 |
| F | 1 | 2 | 3.56 | 0.0284 | 3.56 ± 0.03 |
| M | 0 | 40 | 4.88 | 0.0985 | 4.88 ± 0.10 |
| M | 1 | 2 | 3.90 | 0.298 | 3.90 ± 0.30 |

### Generate SMI chart for Bt
```
clean_data <- data_adult_Bt %>%
  filter(
    !is.na(weight), !is.na(total_length), !is.na(SMI),
    is.finite(weight), is.finite(total_length), is.finite(SMI)
  ) %>%
  mutate(
    sex_infect = case_when(
      sex == "M" & hemoplasma == 0 ~ "Male, uninfected",
      sex == "M" & hemoplasma == 1 ~ "Male, infected",
      sex == "F" & hemoplasma == 0 ~ "Female, uninfected",
      sex == "F" & hemoplasma == 1 ~ "Female, infected",
      TRUE ~ NA_character_
    )
  )

levels_order <- c(
  "Male, uninfected",
  "Male, infected",
  "Female, uninfected",
  "Female, infected"
)

clean_data <- clean_data %>%
  mutate(
    sex_infect = factor(sex_infect, levels = levels_order),
    point_size = case_when(
      sex_infect %in% c("Male, uninfected", "Male, infected") ~ 3.25,
      TRUE ~ 4
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

p_SMI_Bt <- ggplot() +
  geom_contour_filled(
    data = interp_df,
    aes(x = x, y = y, z = z)
  ) +
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
  scale_fill_brewer(
    palette = "YlOrBr",
    name = "SMI level"
  ) +
  scale_shape_manual(
    name = expression(paste("Hemoplasma", " infection status")),
    values = c(
      "Male, uninfected" = 0,
      "Male, infected" = 12,
      "Female, uninfected" = 1,
      "Female, infected" = 10
    )
  ) +
  scale_size_identity(guide = "none") +
  guides(
    shape = guide_legend(
      override.aes = list(size = legend_point_sizes)
    )
  ) +
  labs(
    x = "Body mass (kg)",
    y = "Total Length (cm)",
    title = expression(
      paste(
        "Scale Mass Index (SMI) of ",
        italic("Bradypus didactylus")
      )
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    panel.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

p_SMI_Bt

ggsave(
  filename = "SMI_Bradypus_didactylus.png",
  plot = p_SMI_Bt,
  width = 22.5,
  height = 7,
  units = "in",
  dpi = 300
)
```

## Step 6. Impact of `hemoplasma` infections on Scale Mass Index (SMI) in adult Cd
### Function to calculate SMI for adult Cd
```
data_adult_Cd <- subset(data_Cd, age == "A")
sma_model_Cd <- sma(log(weight) ~ log(total_length), data = data_adult_Cd)
b <- coef(sma_model_Cd)[2]
L0 <- mean(data_adult_Cd$total_length, na.rm = TRUE)
data_adult_Cd$SMI <- data_adult_Cd$weight * (L0 / data_adult_Cd$total_length)^b
```

### Fit a GLM to test whether SMI is influenced by `hemoplasma`, `pathogens`, `sex`, and `season` in Cd
```
model_SMICd_full <- glm(
  SMI ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)
summary(model_SMICd_full)
model_SMICd_3way <- glm(
  SMI ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)
summary(model_SMICd_3way)
model_SMICd_2way <- glm(
  SMI ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)
summary(model_SMICd_2way)
model_SMICd_add <- glm(
  SMI ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)
summary(model_SMICd_add)
model_SMICd_null <- glm(
  SMI ~ 1,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)
summary(model_SMICd_null)
anova(
  model_SMICd_null,
  model_SMICd_add,
  test = "Chisq"
)
anova(
  model_SMICd_add,
  model_SMICd_2way,
  test = "Chisq"
)
anova(
  model_SMICd_2way,
  model_SMICd_3way,
  test = "Chisq"
)
anova(
  model_SMICd_3way,
  model_SMICd_full,
  test = "Chisq"
)
AIC_table <- AIC(
  model_SMICd_null,
  model_SMICd_add,
  model_SMICd_2way,
  model_SMICd_3way,
  model_SMICd_full
)
AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table
drop1(
  model_SMICd_add,
  test = "Chisq"
)
model_SMICd_hemo <- glm(
  SMI ~ hemoplasma,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)
model_SMICd_sex <- glm(
  SMI ~ sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)
model_SMICd_season <- glm(
  SMI ~ season,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)
model_SMICd_pathogens <- glm(
  SMI ~ pathogens,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)
AIC_models <- AIC(
  model_SMICd_null,
  model_SMICd_hemo,
  model_SMICd_sex,
  model_SMICd_season,
  model_SMICd_pathogens
)
AIC_models$delta_AIC <- AIC_models$AIC - AIC(model_SMICd_null)
AIC_models
par(mfrow = c(2, 2))
plot(model_SMICd_add)
par(mfrow = c(1, 1))
shapiro.test(residuals(model_SMICd_add))
plot(
  fitted(model_SMICd_add),
  residuals(model_SMICd_add),
  xlab = "Fitted SMI",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)
with(
  data_adult_Cd,
  table(hemoplasma, pathogens, season, sex)
)
```

-> Results : No evidence of an association between SMI and `hemoplasma` infection was detected in adult *C. didactylus*. The null model was best supported (AIC = 130.25), with no improvement in model fit when `hemoplasma` infection, other blood-borne `pathogens`, `season`, and `sex` were included as additive predictors (LRT: χ²₄ = 1.02, p = 0.774; ΔAIC = +6.07). None of the individual predictors significantly improved the additive model (hemoplasma: χ²₁ = 0.005, p = 0.944; other pathogens: χ²₁ = 1.61, p = 0.204; season: χ²₁ = 0.61, p = 0.434; sex: χ²₁ = 0.11, p = 0.745). Residuals of the additive model were normally distributed (Shapiro–Wilk: W = 0.980, p = 0.484).

-> Interpretation : In adult *C. didactylus*, body condition was not detectably associated with `hemoplasma` infection, co-infection with other blood-borne `pathogens`, `season`, or `sex`. The higher AIC of the additive model relative to the null model indicates that adding these predictors did not improve the explanation of variation in SMI. Overall, the data provide no evidence that `hemoplasma` infection is associated with reduced body condition in this `species`.

### Generate SMI chart for Cd
```
clean_data <- data_adult_Cd %>%
  filter(
    !is.na(weight), !is.na(total_length), !is.na(SMI),
    is.finite(weight), is.finite(total_length), is.finite(SMI)
  ) %>%
  mutate(
    sex_infect = case_when(
      sex == "M" & hemoplasma == 0 ~ "Male, uninfected",
      sex == "M" & hemoplasma == 1 ~ "Male, infected",
      sex == "F" & hemoplasma == 0 ~ "Female, uninfected",
      sex == "F" & hemoplasma == 1 ~ "Female, infected",
      TRUE ~ NA_character_
    )
  )

levels_order <- c(
  "Male, uninfected",
  "Male, infected",
  "Female, uninfected",
  "Female, infected"
)

clean_data <- clean_data %>%
  mutate(
    sex_infect = factor(sex_infect, levels = levels_order),
    point_size = case_when(
      sex_infect %in% c("Male, uninfected", "Male, infected") ~ 3.25,
      TRUE ~ 4
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

p_SMI_Cd <- ggplot() +
  geom_contour_filled(
    data = interp_df,
    aes(x = x, y = y, z = z)
  ) +
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
  scale_fill_brewer(
    palette = "YlOrBr",
    name = "SMI level"
  ) +
  scale_shape_manual(
    name = expression(paste("Hemoplasma", " infection status")),
    values = c(
      "Male, uninfected" = 0,
      "Male, infected" = 12,
      "Female, uninfected" = 1,
      "Female, infected" = 10
    )
  ) +
  scale_size_identity(guide = "none") +
  guides(
    shape = guide_legend(
      override.aes = list(size = legend_point_sizes)
    )
  ) +
  labs(
    x = "Body mass (kg)",
    y = "Total Length (cm)",
    title = expression(
      paste(
        "Scale Mass Index (SMI) of ",
        italic("Choloepus didactylus")
      )
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    panel.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

p_SMI_Cd

ggsave(
  filename = "SMI_Choloepus_didactylus.png",
  plot = p_SMI_Cd,
  width = 22.5,
  height = 7,
  units = "in",
  dpi = 300
)
```

## Step 7. Impact of `hemoplasma` infections on neck circumference

### Fit a GLM to test whether neck circumference is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Bt
```
model_5 <- glm(
  log(neck_size) ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_3way <- glm(
  log(neck_size) ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_2way <- glm(
  log(neck_size) ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_add <- glm(
  log(neck_size) ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_null <- glm(
  log(neck_size) ~ 1,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

anova(model_5_null, model_5_add, test = "Chisq")
anova(model_5_add, model_5_2way, test = "Chisq")
anova(model_5_2way, model_5_3way, test = "Chisq")
anova(model_5_3way, model_5, test = "Chisq")

AIC_table <- AIC(
  model_5_null,
  model_5_add,
  model_5_2way,
  model_5_3way,
  model_5
)
AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table

drop1(
  model_5_add,
  test = "Chisq"
)

model_5_hemo <- glm(
  log(neck_size) ~ hemoplasma,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_pathogens <- glm(
  log(neck_size) ~ pathogens,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_season <- glm(
  log(neck_size) ~ season,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_sex <- glm(
  log(neck_size) ~ sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

AIC_models <- AIC(
  model_5_null,
  model_5_hemo,
  model_5_pathogens,
  model_5_season,
  model_5_sex
)
AIC_models$delta_AIC <- AIC_models$AIC - AIC(model_5_null)
AIC_models

par(mfrow = c(2, 2))
plot(model_5_add)
par(mfrow = c(1, 1))

shapiro.test(residuals(model_5_null))

plot(
  fitted(model_5_null),
  residuals(model_5_null),
  xlab = "Fitted log(neck circumference)",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)
```
-> Results : In adult *Bradypus tridactylus*, none of the interaction models improved model fit relative to simpler models (additive vs. two-way interactions: LRT, χ²₃ = 0.01, *p* = 0.765; two-way vs. three-way interactions: χ²₁ = 0.003, *p* = 0.554), and the three-way model was equivalent to the full model because no additional parameters were estimable. The additive model did not improve on the null model (χ²₄ = 0.04, *p* = 0.298), and the null model had the lowest AIC (−112.64; ΔAIC = 0.00), compared with the additive model (−109.75; ΔAIC = 2.88). In the additive model, none of the predictors was significantly associated with log-transformed neck circumference (`hemoplasma`: χ²₁ = 0.72, *p* = 0.398; other `pathogens`: χ²₁ = 0.48, *p* = 0.487; `season`: χ²₁ = 0.22, *p* = 0.638; `sex`: χ²₁ = 2.43, *p* = 0.119). Single-predictor models provided little additional support for any predictor, although the `sex`-only model had a slightly lower AIC than the null model (ΔAIC = −1.42). Residuals of the null model showed no significant departure from normality (Shapiro–Wilk: W = 0.967, *p* = 0.093).

-> Interpretation : Neck circumference was not significantly associated with `hemoplasma` infection, other blood-borne `pathogens`, or `season` in adult *B. tridactylus*. 

### Fit a GLM to test whether neck circumference is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Cd
```
model_6 <- glm(
  log(neck_size) ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_3way <- glm(
  log(neck_size) ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_2way <- glm(
  log(neck_size) ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_add <- glm(
  log(neck_size) ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_null <- glm(
  log(neck_size) ~ 1,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

anova(model_6_null, model_6_add, test = "Chisq")
anova(model_6_add, model_6_2way, test = "Chisq")
anova(model_6_2way, model_6_3way, test = "Chisq")
anova(model_6_3way, model_6, test = "Chisq")

AIC_table <- AIC(
  model_6_null,
  model_6_add,
  model_6_2way,
  model_6_3way,
  model_6
)
AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table

drop1(
  model_6_add,
  test = "Chisq"
)

model_6_hemo <- glm(
  log(neck_size) ~ hemoplasma,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_pathogens <- glm(
  log(neck_size) ~ pathogens,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_season <- glm(
  log(neck_size) ~ season,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_sex <- glm(
  log(neck_size) ~ sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

AIC_models <- AIC(
  model_6_null,
  model_6_hemo,
  model_6_pathogens,
  model_6_season,
  model_6_sex
)
AIC_models$delta_AIC <- AIC_models$AIC - AIC(model_6_null)
AIC_models

par(mfrow = c(2, 2))
plot(model_6_add)
par(mfrow = c(1, 1))

shapiro.test(residuals(model_6_null))

plot(
  fitted(model_6_null),
  residuals(model_6_null),
  xlab = "Fitted log(neck circumference)",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)
```

-> Results : In adult *Choloepus didactylus*, adding the four predictors as additive effects did not improve model fit relative to the null model (LRT: χ²₄ = 0.05, *p* = 0.370; ΔAIC = +3.45). Adding two-way interactions provided no further improvement (χ²₆ = 0.01, *p* = 0.992), whereas the addition of three-way interactions significantly reduced residual deviance (χ²₂ = 0.09, *p* = 0.026). However, the three-way model had a substantially higher AIC than the null model (ΔAIC = +9.29), and was therefore not supported overall. In the additive model, none of the predictors was significantly associated with log-transformed neck circumference (`hemoplasma`: χ²₁ = 0.001, *p* = 0.981; other `pathogens`: χ²₁ = 0.74, *p* = 0.391; `season`: χ²₁ = 1.40, *p* = 0.237; `sex`: χ²₁ = 1.40, *p* = 0.238). Single-predictor models provided no meaningful support for any predictor, with AIC values very similar to that of the null model (ΔAIC ranging from −0.21 for `season` to +1.81 for `hemoplasma`). Residuals of the null model showed no significant departure from normality (Shapiro–Wilk: W = 0.972, *p* = 0.314).

-> Interpretation : Neck circumference was not associated with `hemoplasma` infection, other blood-borne `pathogens`, `season`, or `sex` in adult *C. didactylus*. 

## Step 8. Impact of `hemoplasma` infections on `hematocrit` levels
### Fit a GLM to test whether `hematocrit` is influenced by `hemoplasma`, `pathogens`, `sex`, and `season` in Bt
```
model_7 <- glm(
  hematocrit ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_3way <- glm(
  hematocrit ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_2way <- glm(
  hematocrit ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_add <- glm(
  hematocrit ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_null <- glm(
  hematocrit ~ 1,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

anova(model_7_null, model_7_add, test = "Chisq")
anova(model_7_add, model_7_2way, test = "Chisq")
anova(model_7_2way, model_7_3way, test = "Chisq")
anova(model_7_3way, model_7, test = "Chisq")

AIC_table <- AIC(
  model_7_null,
  model_7_add,
  model_7_2way,
  model_7_3way,
  model_7
)
AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table

drop1(
  model_7_add,
  test = "Chisq"
)

model_7_hemo <- glm(
  hematocrit ~ hemoplasma,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_pathogens <- glm(
  hematocrit ~ pathogens,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_season <- glm(
  hematocrit ~ season,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_sex <- glm(
  hematocrit ~ sex,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

AIC_models <- AIC(
  model_7_null,
  model_7_hemo,
  model_7_pathogens,
  model_7_season,
  model_7_sex
)
AIC_models$delta_AIC <- AIC_models$AIC - AIC(model_7_null)
AIC_models

par(mfrow = c(2, 2))
plot(model_7_add)
par(mfrow = c(1, 1))

plot(
  fitted(model_7_add),
  residuals(model_7_add, type = "deviance"),
  xlab = "Fitted hematocrit",
  ylab = "Deviance residuals"
)
abline(h = 0, lty = 2)
```
-> -> Results : In adult *Bradypus tridactylus*, the additive model did not improve on the null model (LRT: χ²₄ = 0.08, *p* = 0.321; ΔAIC = +3.13). Neither two-way nor three-way interactions improved fit (χ²₄ = 0.10, *p* = 0.193; χ²₁ = 0.0001, *p* = 0.928). None of the predictors was significant in the additive model (`hemoplasma`: χ²₁ = 0.46, *p* = 0.497; `pathogens`: χ²₁ = 0.18, *p* = 0.669; `season`: χ²₁ = 0.40, *p* = 0.527; `sex`: χ²₁ = 3.36, *p* = 0.067). The null model had the lowest AIC (513.01). The `sex`-only model had a slightly lower AIC (ΔAIC = −1.63), but this provided weak support.

-> Interpretation : `hematocrit` showed no significant association with `hemoplasma` infection, other blood-borne `pathogens`, `season`, or `sex` in adult *Bradypus tridactylus*.

### Calculation of mean and standard error of `hematocrit` by `hemoplasma` for Bt
```
data_adult_Bt %>%
  group_by(hemoplasma) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results :
| Hemoplasma infection | Mean hematocrit (%) | SE |
|:---------------------|--------------------:|---:|
| Negative             | 39.0                | 0.57 |
| Positive             | 40.8                | 3.04 |


### Fit a GLM to test whether `hematocrit` is influenced by `hemoplasma`, `pathogens`, `sex`, and `season` in Cd
```
model_8 <- glm(
  hematocrit ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_3way <- glm(
  hematocrit ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_2way <- glm(
  hematocrit ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_add <- glm(
  hematocrit ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_null <- glm(
  hematocrit ~ 1,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

anova(model_8_null, model_8_add, test = "Chisq")
anova(model_8_add, model_8_2way, test = "Chisq")
anova(model_8_2way, model_8_3way, test = "Chisq")
anova(model_8_3way, model_8, test = "Chisq")

AIC_table <- AIC(
  model_8_null,
  model_8_add,
  model_8_2way,
  model_8_3way,
  model_8
)
AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table

drop1(
  model_8_add,
  test = "Chisq"
)

model_8_hemo <- glm(
  hematocrit ~ hemoplasma,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_pathogens <- glm(
  hematocrit ~ pathogens,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_season <- glm(
  hematocrit ~ season,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_sex <- glm(
  hematocrit ~ sex,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

AIC_models <- AIC(
  model_8_null,
  model_8_hemo,
  model_8_pathogens,
  model_8_season,
  model_8_sex
)
AIC_models$delta_AIC <- AIC_models$AIC - AIC(model_8_null)
AIC_models

par(mfrow = c(2, 2))
plot(model_8_add)
par(mfrow = c(1, 1))

plot(
  fitted(model_8_add),
  residuals(model_8_add, type = "deviance"),
  xlab = "Fitted hematocrit",
  ylab = "Deviance residuals"
)
abline(h = 0, lty = 2)
```
-> Results : In adult *Choloepus didactylus*, the additive model improved fit relative to the null model (LRT: χ²₄ = 0.21, *p* = 0.028; ΔAIC = −1.82) and was retained. Adding two-way or three-way interactions did not improve fit (χ²₆ = 0.08, *p* = 0.692; χ²₁ = 0.003, *p* = 0.703). In the additive model, `hematocrit` was significantly associated with `season` (χ²₁ = 8.17, *p* = 0.004), but not with `hemoplasma`, other `pathogens`, or `sex` (*p* > 0.10). The `season`-only model also had the lowest AIC among single-predictor models (ΔAIC = −2.47 relative to the null).

-> Interpretation : `hematocrit` showed a seasonal pattern in adult *C. didactylus*, but was not detectably associated with `hemoplasma` infection, other `pathogens`, or `sex`.

### Calculation of mean and standard error of `hematocrit` by `hemoplasma` for Cd
```
data_adult_Cd %>%
  group_by(hemoplasma) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results :
| Hemoplasma infection | Mean hematocrit (%) | SE |
|:---------------------|--------------------:|---:|
| Negative             | 39.1                | 1.48 |
| Positive             | 38.6                | 0.80 |

### Calculation of mean and standard error of `hematocrit` by `season` for Cd
```
data_adult_Cd %>%
  group_by(season) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results :
| Season | Mean hematocrit (%) | SE |
|:-------|--------------------:|---:|
| Dry    | 37.6                | 0.89 |
| Wet    | 41.1                | 0.94 |


Create violin plots for `hematocrit`
```
label_style <- element_text(size = 28, face = "bold")

panel_theme <- theme_minimal() +
  theme(
    plot.title = label_style,
    plot.title.position = "plot",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

pC <- ggplot(data_adult_Bt, aes(
  x = factor(hemoplasma, levels = c(0, 1),
             labels = c("Uninfected", "Infected")),
  y = hematocrit
)) +
  geom_violin(fill = "darkorange2", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, shape = 21, fill = NA, color = "black", stroke = 0.7) +
  labs(
    x = "Hemoplasma infection status",
    y = "Hematocrit (%)",
    title = "C"
  ) +
  scale_y_continuous(limits = c(20, 60)) +
  panel_theme

pD <- ggplot(data_adult_Cd, aes(
  x = factor(hemoplasma, levels = c(0, 1),
             labels = c("Uninfected", "Infected")),
  y = hematocrit
)) +
  geom_violin(fill = "darkorange2", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, shape = 21, fill = NA, color = "black", stroke = 0.7) +
  labs(
    x = "Hemoplasma infection status",
    y = "Hematocrit (%)",
    title = "D"
  ) +
  scale_y_continuous(limits = c(10, 60)) +
  panel_theme

final_plot <- pC | pD

print(final_plot)

ggsave(
  "hematocrit_hemoplasma_by_species.png",
  plot = final_plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)
```

## Step 9. Impact of `hemoplasma` infections on body `temperature`
### Convert `temperature` to numeric, handle left-censored values (<32°C) for analysis, and create a left-censored Surv object (temp) for `temperature` in Bt and Cd
```
data_adult_Bt <- data_adult_Bt %>%
  mutate(
    temperature_numeric = as.numeric(ifelse(temperature == "< 32.00", 32, temperature)),
    censored = temperature == "< 32.00",
    temp_surv = Surv(temperature_numeric, event = !censored, type = "left")
  )

data_adult_Cd <- data_adult_Cd %>%
  mutate(
    temperature_numeric = as.numeric(ifelse(temperature == "< 32.00", 32, temperature)),
    censored = temperature == "< 32.00",
    temp_surv = Surv(temperature_numeric, event = !censored, type = "left")
  )

table(data_adult_Bt$censored, useNA = "ifany")
table(data_adult_Cd$censored, useNA = "ifany")
```

### Fit Gaussian survival regression models to test the effects of `hemoplasma`, `pathogens`, `season`, `sex` on `temperature` in Bt
```
model_10 <- survreg(
  temp_surv ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Bt,
  dist = "gaussian"
)

model_10_3way <- survreg(
  temp_surv ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Bt,
  dist = "gaussian"
)

model_10_2way <- survreg(
  temp_surv ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Bt,
  dist = "gaussian"
)

model_10_add <- survreg(
  temp_surv ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Bt,
  dist = "gaussian"
)

model_10_null <- survreg(
  temp_surv ~ 1,
  data = data_adult_Bt,
  dist = "gaussian"
)

anova(model_10_null, model_10_add)
anova(model_10_add, model_10_2way)
anova(model_10_2way, model_10_3way)
anova(model_10_3way, model_10)

AIC_table_10 <- AIC(
  model_10_null,
  model_10_add,
  model_10_2way,
  model_10_3way,
  model_10
)

AIC_table_10$delta_AIC <- AIC_table_10$AIC - min(AIC_table_10$AIC)

AIC_table_10

drop1(model_10_add, test="Chisq")

model_10_hemo <- survreg(
  temp_surv ~ hemoplasma,
  data=data_adult_Bt,
  dist="gaussian"
)

model_10_pathogens <- survreg(
  temp_surv ~ pathogens,
  data=data_adult_Bt,
  dist="gaussian"
)

model_10_season <- survreg(
  temp_surv ~ season,
  data=data_adult_Bt,
  dist="gaussian"
)

model_10_sex <- survreg(
  temp_surv ~ sex,
  data=data_adult_Bt,
  dist="gaussian"
)

AIC_models_10 <- AIC(
  model_10_null,
  model_10_hemo,
  model_10_pathogens,
  model_10_season,
  model_10_sex
)

AIC_models_10$delta_AIC <- AIC_models_10$AIC - AIC(model_10_null)

AIC_models_10

model_10_final <- model_10_add

summary(model_10_final)

par(mfrow=c(2,2))
plot(model_10_final)
par(mfrow=c(1,1))

plot(
  fitted(model_10_final),
  residuals(model_10_final),
  xlab="Fitted body temperature",
  ylab="Residuals"
)
abline(h=0,lty=2)

shapiro.test(residuals(model_10_final))
```
-> Results : In adult *Bradypus tridactylus*, adding the four predictors as additive effects did not significantly improve model fit relative to the null model (LRT: χ²₄ = 5.91, *p* = 0.206; ΔAIC = +2.09), and the null model had the lowest AIC (89.36). Adding two-way or three-way interactions did not improve fit (χ²₆ = 1.66, *p* = 0.948; χ²₄ = 0.13, *p* = 0.998), and the three-way and four-way models were equivalent (χ²₁ = 0, *p* = 1.00). In the additive model, body `temperature` was significantly lower in males (β = −1.65 ± 0.72 SE, *p* = 0.022), whereas `hemoplasma` infection (β = 2.06 ± 1.79 SE, *p* = 0.250), other blood-borne `pathogens` (β = 0.51 ± 0.72 SE, *p* = 0.475), and `season` (β = −0.53 ± 0.80 SE, *p* = 0.507) were not significantly associated with body `temperature`. The `sex`-only model had the lowest AIC among single-predictor models (ΔAIC = −1.69 relative to the null). Residuals did not significantly deviate from normality (Shapiro–Wilk: W = 0.958, *p* = 0.222).

-> Interpretation : Body `temperature` showed a `sex`-related difference in adult *Bradypus tridactylus*, with lower `temperature` in males. However, the additive model was not clearly supported over the null model based on AIC, and there was no detectable association between body `temperature` and `hemoplasma` infection, other blood-borne `pathogens`, or `season`.

### Fit Gaussian survival regression models to test the effects of `hemoplasma`, `pathogens`, `season`, `sex` on `temperature` in Cd
```
model_11 <- survreg(
  temp_surv ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Cd,
  dist = "gaussian"
)

model_11_3way <- survreg(
  temp_surv ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Cd,
  dist = "gaussian"
)

model_11_2way <- survreg(
  temp_surv ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Cd,
  dist = "gaussian"
)

model_11_add <- survreg(
  temp_surv ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Cd,
  dist = "gaussian"
)

model_11_null <- survreg(
  temp_surv ~ 1,
  data = data_adult_Cd,
  dist = "gaussian"
)

anova(model_11_null, model_11_add)
anova(model_11_add, model_11_2way)
anova(model_11_2way, model_11_3way)
anova(model_11_3way, model_11)

AIC_table_11 <- AIC(
  model_11_null,
  model_11_add,
  model_11_2way,
  model_11_3way,
  model_11
)

AIC_table_11$delta_AIC <- AIC_table_11$AIC - min(AIC_table_11$AIC)

AIC_table_11

drop1(model_11_add, test="Chisq")

model_11_hemo <- survreg(
  temp_surv ~ hemoplasma,
  data=data_adult_Cd,
  dist="gaussian"
)

model_11_pathogens <- survreg(
  temp_surv ~ pathogens,
  data=data_adult_Cd,
  dist="gaussian"
)

model_11_season <- survreg(
  temp_surv ~ season,
  data=data_adult_Cd,
  dist="gaussian"
)

model_11_sex <- survreg(
  temp_surv ~ sex,
  data=data_adult_Cd,
  dist="gaussian"
)

AIC_models_11 <- AIC(
  model_11_null,
  model_11_hemo,
  model_11_pathogens,
  model_11_season,
  model_11_sex
)

AIC_models_11$delta_AIC <- AIC_models_11$AIC - AIC(model_11_null)

AIC_models_11

model_11_final <- model_11_add

summary(model_11_final)

plot(
  fitted(model_11_final),
  residuals(model_11_final),
  xlab="Fitted body temperature",
  ylab="Residuals"
)
abline(h=0,lty=2)

shapiro.test(residuals(model_11_final))
```
-> Results : In adult *Choloepus didactylus*, the additive model did not improve model fit relative to the null model (LRT: χ²₄ = 0.56, *p* = 0.967; ΔAIC = +7.44), and the null model had the lowest AIC (62.74). Adding two-way or three-way interactions did not improve fit (χ²₆ = 3.55, *p* = 0.737; χ²₄ = 0, *p* = 1.00), and the three-way and four-way models were equivalent (χ²₁ = 0, *p* = 1.00). None of the predictors was significantly associated with body `temperature` in the additive model (`hemoplasma`: β = 0.34 ± 0.67 SE, *p* = 0.611; other blood-borne `pathogens`: β = −0.15 ± 0.60 SE, *p* = 0.802; `season`: β = 0.18 ± 0.65 SE, *p* = 0.782; `sex`: β = −0.19 ± 0.64 SE, *p* = 0.768). Single-predictor models provided similarly little support, with ΔAIC values ranging from +1.69 to +1.97 relative to the null. Residuals of the additive model showed no significant departure from normality (Shapiro–Wilk: W = 0.977, *p* = 0.907).

-> Interpretation : Body `temperature` showed no detectable association with `hemoplasma` infection, other blood-borne `pathogens`, `season`, or `sex` in adult *Choloepus didactylus*.

## Step 10. Impact of `hemoplasma` infections on general `health_condition` 
### Test the association between `hemoplasma` and `health_condition` in Bt
```
table_health_condition_hemoplasma_Bt <- table(data_Bt$hemoplasma, data_Bt$health_condition)
table_health_condition_hemoplasma_Bt
fisher.test(table_health_condition_hemoplasma_Bt)
```
-> Results : In adult *Bradypus tridactylus*, `hemoplasma` infection was not associated with `health_condition` (Fisher’s exact test, p = 1.00). All four `hemoplasma`-positive individuals were classified as having a good `health_condition`.

-> Interpretation : There was no detectable association between `hemoplasma` infection and `health_condition` in adult *Bradypus tridactylus*. 

### Test the association between `hemoplasma` and `health_condition` in Cd
```
table_health_condition_hemoplasma_Cd <- table(data_Cd$hemoplasma, data_Cd$health_condition)
table_health_condition_hemoplasma_Cd
fisher.test(table_health_condition_hemoplasma_Cd)
```
-> Results :In adult *Choloepus didactylus*, `hemoplasma` infection was not associated with `health_condition` (Fisher’s exact test, *p* = 1.00). Among `hemoplasma`-positive individuals, 64/68 (94.1%) were classified as having a good `health_condition`, compared with 14/15 (93.3%) among `hemoplasma`-negative individuals.

-> Interpretation : There was no detectable association between `hemoplasma` infection and `health_condition` in adult *Choloepus didactylus*.


## Step 11. Impact of `hemoplasma` infections on `female_reproductive_status`
### Test the association between `hemoplasma` and `female_reproductive_status` in Bt
```
table_hemoplasma_infection_female_Bt <- table(data_Bt$hemoplasma, data_Bt$female_reproductive_status)
table_hemoplasma_infection_female_Bt
fisher.test(table_hemoplasma_infection_female_Bt)
```
-> Results : In adult *Bradypus tridactylus* females, `hemoplasma` infection was not associated with `female_reproductive_status` (Fisher’s exact test, *p* = 1.00). The two `hemoplasma`-positive females were both classified as non-pregnant and non-lactating, while no infections were detected among lactating or pregnant females.

-> Interpretation: There was no detectable association between `hemoplasma` infection and `female_reproductive_status` in adult *Bradypus tridactylus*. However, the very small number of `hemoplasma`-positive females (n = 2) limits the power of this comparison.

### Test the association between `hemoplasma` and `female_reproductive_status` in Cd
```
table_hemoplasma_infection_female_Cd <- table(data_Cd$hemoplasma, data_Cd$female_reproductive_status)
table_hemoplasma_infection_female_Cd
fisher.test(table_hemoplasma_infection_female_Cd)
```
-> Results : In adult *Choloepus didactylus* females, `hemoplasma` infection was not associated with `female_reproductive_status` (Fisher’s exact test, *p* = 0.505). Among `hemoplasma`-positive females, 5/39 (12.8%) were lactating with a young, 32/39 (82.1%) were non-pregnant and non-lactating, and 2/39 (5.1%) were pregnant. The corresponding proportions among `hemoplasma`-negative females were 2/10 (20.0%), 7/10 (70.0%), and 1/10 (10.0%), respectively.

-> Interpretation : There was no detectable association between `hemoplasma` infection and `female_reproductive_status` in adult *Choloepus didactylus*.
