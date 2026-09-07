# **Hemoplasma epidemiological survey: R scripts and analysis pipeline**

We analyzed data from 626 individuals belonging to 44 species of wild mammals sampled in French Guiana. Details of sampling, examined traits and laboratory procedures are provided in the associated manuscript.

## Dataset description
-> For the epidemiological survey dataset ([here](data_hemoplasma_stat.csv)), each row corresponds to one individual (n = 626) and includes the following variables : 
- `species` : Species identity (44 wild mammal species)
- `order` : Taxonomic order
- `hemoplasma` : Infection status with hemotropic mycoplasmas (0 = uninfected; 1 = infected)
- `sex` : Sex of the individual (M = male; F = female)
- `anaplasmataceae` : Infection status with bacteria of the Anaplasmataceae family (*Anaplasma*, *Ehrlichia* and *Allocryptoplasma*) (0 = uninfected; 1 = infected)
- `apicomplexa` : Infection status with piroplasmids (*Babesia* and *Theileria*) and haemogregarines (*Hepatozoon* and *Hemolivia*) (0 = uninfected; 1 = infected)

-> For the life trait dataset ([here](https://github.com/olivierduron/Hemoplasma_infections/blob/main/data_mammal_traits.csv)), each row corresponds to one species (n = 44) and includes the following variables :
-  `species` : Species identity (44 wild mammal species)
-  `dietinv` : Percentage of the diet consisting of invertebrates (from 0 to 100%)
-  `dietvet` : Percentage of the diet consisting of vertebrates (from 0 to 100%)
-  `dietplant` : Percentage of the diet consisting of plants (fruits, leaves, seeds, nectar, etc; from 0 to 100%)
-  `strata` : Foraging stratum category (G = Ground level, including aquatic foraging; S = Scansorial; Ar = Arboreal)
-  `activitynocturnal` : Foraging activity at night (0 = no; 1 = yes)
-  `activitycrepuscular` : Foraging activity at at twilight (0 = no; 1 = yes)
-  `activitydiurnal` : Foraging activity at day (0 = no; 1 = yes)
-  `bodymass` : Mean adult body mass (grammes)
-  `longevity` : Mean longevity (years)
-  `femalematurity` : Mean age at maturity for females (days)
-  `littersize` : Mean litter size (number of offspring per litter)
   
## Table of contents
- [Step 1. Data retrieval and preparation](#step-1-data-retrieval-and-preparation)
- [Step 2. Species-level summary for hemoplasma](#step-2-species-level-summary-for-hemoplasma)
- [Step 3. Variation in hemoplasma infection status according to the host sex](#step-3-variation-in-hemoplasma-infection-status-according-to-the-host-sex)
- [Step 4. Variation in hemoplasma infection status according to the presence of other blood-borne pathogens (anaplasmataceae and apicomplexa)](#step-4-variation-in-hemoplasma-infection-status-according-to-the-presence-of-other-blood-borne-pathogens-anaplasmataceae-and-apicomplexa)
- [Step 5. Hemoplasma prevalence by mammalian order](#step-5-hemoplasma-prevalence-by-mammalian-order)
- [Step 6. Phylogeny of the 44 mammalian species](#step-6-phylogeny-of-the-44-mammalian-species)
- [Step 7. Hemoplasma infection prevalence and mammal trait-based analyses](#step-7-hemoplasma-infection-prevalence-and-mammal-trait-based-analyses)

## Step 1. Data retrieval and preparation

### Retrieve and examine the epidemiological dataset
```
data_hemoplasma_stat <- read.csv2("https://raw.githubusercontent.com/olivierduron/Hemoplasma_infections/main/data_hemoplasma_stat.csv")
data_hemoplasma_stat
str(data_hemoplasma_stat)
get_modalities <- function(x) {sort(table(x), decreasing = TRUE)}
lapply(data_hemoplasma_stat, get_modalities)
```

### Retrieve and examine the life trait dataset
```
data_mammal_traits <- read.csv2("https://raw.githubusercontent.com/olivierduron/Hemoplasma_infections/main/data_mammal_traits.csv")
data_mammal_traits
str(data_mammal_traits)
get_modalities <- function(x) {sort(table(x), decreasing = TRUE)}
lapply(data_mammal_traits, get_modalities)
```

### Convert categorical variables
```
data_hemoplasma_stat$species        <- as.factor(data_hemoplasma_stat$species)
data_hemoplasma_stat$order           <- as.factor(data_hemoplasma_stat$order)
data_hemoplasma_stat$sex   <- as.factor(data_hemoplasma_stat$sex)
data_hemoplasma_stat$hemoplasma      <- as.factor(data_hemoplasma_stat$hemoplasma)
data_hemoplasma_stat$anaplasmataceae      <- as.factor(data_hemoplasma_stat$anaplasmataceae)
data_hemoplasma_stat$apicomplexa         <- as.factor(data_hemoplasma_stat$apicomplexa)
data_mammal_traits$species        <- as.factor(data_mammal_traits$species)
data_mammal_traits$strata        <- as.factor(data_mammal_traits$strata)
data_mammal_traits$activitynocturnal        <- as.factor(data_mammal_traits$activitynocturnal)
data_mammal_traits$activitycrepuscular        <- as.factor(data_mammal_traits$activitydiurnal)
data_mammal_traits$activitydiurnal        <- as.factor(data_mammal_traits$activitydiurnal)
```

### Load required libraries 
```
library(binom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayestestR)
library(posterior)
library(rotl)
library(stringr)
library(ape)
library(glmmTMB)
library(picante)
library(phytools)
library(MCMCglmm)
library(scales)
library(ggthemes)
library(ggtree)
library(ggtext)
library(lme4)
library(car)
library(ggrepel)
library(emmeans)
library(brms)
library(patchwork)
```

## Step 2. `species`-level summary for `hemoplasma`

### Calculate `species`-level `hemoplasma` prevalence with 95% confidence intervals (CI; Wilson method)
```
species_summary <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(as.numeric(as.character(hemoplasma)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_negative = n_sampled - n_positive,
    prevalence = n_positive / n_sampled
  ) %>%
  rowwise() %>%
  mutate(
    ci_low = binom::binom.confint(
      n_positive, n_sampled, method = "wilson"
    )$lower,
    ci_high = binom::binom.confint(
      n_positive, n_sampled, method = "wilson"
    )$upper
  ) %>%
  ungroup() %>%
  mutate(
    ci_low = ifelse(n_positive == 0, 0, ci_low)
  )
print(species_summary, n = Inf)
```

-> Results : 
| Species | N sampled | N positive | N negative | Prevalence | 95% CI |
|:---|---:|---:|---:|---:|:---|
| *Alouatta macconnelli* | 22 | 20 | 2 | 90.9% | 72.2–97.5% |
| *Bradypus tridactylus* | 108 | 4 | 104 | 3.7% | 1.5–9.1% |
| *Cabassous unicinctus* | 2 | 0 | 2 | 0% | 0–65.8% |
| *Caluromys philander* | 5 | 0 | 5 | 0% | 0–43.4% |
| *Sapajus apella* | 1 | 0 | 1 | 0% | 0–79.3% |
| *Choloepus didactylus* | 90 | 72 | 18 | 80.0% | 70.6–87.0% |
| *Coendou melanurus* | 1 | 0 | 1 | 0% | 0–79.3% |
| *Coendou prehensilis* | 3 | 1 | 2 | 33.3% | 6.2–79.2% |
| *Cyclopes didactylus* | 1 | 0 | 1 | 0% | 0–79.3% |
| *Dasypus novemcinctus* | 15 | 5 | 10 | 33.3% | 15.2–58.3% |
| *Didelphis marsupialis* | 51 | 24 | 27 | 47.1% | 34.1–60.5% |
| *Eira barbara* | 4 | 0 | 4 | 0% | 0–49.0% |
| *Leopardus wiedii* | 1 | 0 | 1 | 0% | 0–79.3% |
| *Galictis vittata* | 4 | 3 | 1 | 75.0% | 30.1–95.4% |
| *Holochilus sciureus* | 5 | 1 | 4 | 20.0% | 3.6–62.4% |
| *Hydrochoerus hydrochaeris* | 2 | 0 | 2 | 0% | 0–65.8% |
| *Hylaeamys megacephalus* | 15 | 0 | 15 | 0% | 0–20.4% |
| *Hylaeamys yunganus* | 10 | 0 | 10 | 0% | 0–27.8% |
| *Lontra longicaudis* | 1 | 1 | 0 | 100% | 20.7–100% |
| *Makalata didelphoides* | 8 | 0 | 8 | 0% | 0–32.4% |
| *Marmosa lepida* | 1 | 1 | 0 | 100% | 20.7–100% |
| *Marmosa murina* | 20 | 3 | 17 | 15.0% | 5.2–36.0% |
| *Marmosops parvidens* | 5 | 0 | 5 | 0% | 0–43.4% |
| *Mesomys hispidus* | 13 | 0 | 13 | 0% | 0–22.8% |
| *Metachirus nudicaudatus* | 5 | 0 | 5 | 0% | 0–43.4% |
| *Micoureus demerarae* | 16 | 2 | 14 | 12.5% | 3.5–36.0% |
| *Mus musculus* | 34 | 0 | 34 | 0% | 0–10.2% |
| *Neacomys dubosti* | 1 | 0 | 1 | 0% | 0–79.3% |
| *Neacomys paracou* | 8 | 0 | 8 | 0% | 0–32.4% |
| *Nectomys rattus* | 4 | 2 | 2 | 50.0% | 15.0–85.0% |
| *Oecomys auyantepui* | 16 | 1 | 15 | 6.3% | 1.1–28.3% |
| *Oecomys bicolor* | 16 | 0 | 16 | 0% | 0–19.4% |
| *Oligoryzomys fulvescens* | 7 | 1 | 6 | 14.3% | 2.6–51.3% |
| *Philander opossum* | 20 | 9 | 11 | 45.0% | 25.8–65.8% |
| *Pithecia pithecia* | 1 | 0 | 1 | 0% | 0–79.3% |
| *Potos flavus* | 2 | 1 | 1 | 50.0% | 9.5–90.5% |
| *Proechimys cuvieri* | 18 | 2 | 16 | 11.1% | 3.1–32.8% |
| *Proechimys guyannensis* | 20 | 1 | 19 | 5.0% | 0.9–23.6% |
| *Puma yagouaroundi* | 5 | 0 | 5 | 0% | 0–43.4% |
| *Rattus rattus* | 19 | 2 | 17 | 10.5% | 2.9–31.4% |
| *Saguinus midas* | 41 | 41 | 0 | 100% | 91.4–100% |
| *Saimiri sciureus* | 1 | 0 | 1 | 0% | 0–79.3% |
| *Sciurus aestuans* | 1 | 0 | 1 | 0% | 0–79.3% |
| *Tamandua tetradactyla* | 3 | 0 | 3 | 0% | 0–56.1% |

### Visualization of `species`-level `hemoplasma` prevalence with 95% confidence intervals (CI; Wilson method), restricted to `species` with n ≥ 15  (n `species` = 16) (Fig. 1) 
```
windowsFonts(Calibri = windowsFont("Calibri"))

species_summary <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(
      as.numeric(as.character(hemoplasma)),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    n_negative = n_sampled - n_positive,
    prevalence = n_positive / n_sampled
  ) %>%
  rowwise() %>%
  mutate(
    ci_low = binom::binom.confint(
      n_positive,
      n_sampled,
      method = "wilson"
    )$lower,
    ci_high = binom::binom.confint(
      n_positive,
      n_sampled,
      method = "wilson"
    )$upper
  ) %>%
  ungroup() %>%
  mutate(
    ci_low = ifelse(n_positive == 0, 0, ci_low)
  )

print(species_summary, n = Inf)

# Species with n >= 15
species_summary_n15 <- species_summary %>%
  filter(n_sampled >= 15)

species_order <- c(
  "Alouatta_macconnelli",
  "Saguinus_midas",
  "Sapajus_apella",
  "Saimiri_sciureus",
  "Pithecia_pithecia",
  "Bradypus_tridactylus",
  "Choloepus_didactylus",
  "Cyclopes_didactylus",
  "Tamandua_tetradactyla",
  "Cabassous_unicinctus",
  "Dasypus_novemcinctus",
  "Hydrochoerus_hydrochaeris",
  "Holochilus_sciureus",
  "Hylaeamys_megacephalus",
  "Hylaeamys_yunganus",
  "Neacomys_dubosti",
  "Neacomys_paracou",
  "Nectomys_rattus",
  "Oecomys_auyantepui",
  "Oecomys_bicolor",
  "Oligoryzomys_fulvescens",
  "Makalata_didelphoides",
  "Mesomys_hispidus",
  "Proechimys_cuvieri",
  "Proechimys_guyannensis",
  "Coendou_melanurus",
  "Coendou_prehensilis",
  "Mus_musculus",
  "Rattus_rattus",
  "Sciurus_aestuans",
  "Leopardus_wiedii",
  "Puma_yagouaroundi",
  "Eira_barbara",
  "Galictis_vittata",
  "Lontra_longicaudis",
  "Potos_flavus",
  "Caluromys_philander",
  "Didelphis_marsupialis",
  "Marmosa_lepida",
  "Marmosa_murina",
  "Marmosops_parvidens",
  "Metachirus_nudicaudatus",
  "Micoureus_demerarae",
  "Philander_opossum"
)

plot_data <- species_summary_n15 %>%
  mutate(
    group = case_when(
      species %in% c(
        "Alouatta_macconnelli",
        "Saguinus_midas",
        "Sapajus_apella",
        "Saimiri_sciureus",
        "Pithecia_pithecia"
      ) ~ "Primates",
      species %in% c(
        "Bradypus_tridactylus",
        "Choloepus_didactylus",
        "Cyclopes_didactylus",
        "Tamandua_tetradactyla"
      ) ~ "Xenarthrans",
      species %in% c(
        "Cabassous_unicinctus",
        "Dasypus_novemcinctus"
      ) ~ "Armadillos",
      species %in% c(
        "Hydrochoerus_hydrochaeris",
        "Holochilus_sciureus",
        "Hylaeamys_megacephalus",
        "Hylaeamys_yunganus",
        "Neacomys_dubosti",
        "Neacomys_paracou",
        "Nectomys_rattus",
        "Oecomys_auyantepui",
        "Oecomys_bicolor",
        "Oligoryzomys_fulvescens",
        "Makalata_didelphoides",
        "Mesomys_hispidus",
        "Proechimys_cuvieri",
        "Proechimys_guyannensis",
        "Coendou_melanurus",
        "Coendou_prehensilis",
        "Mus_musculus",
        "Rattus_rattus",
        "Sciurus_aestuans"
      ) ~ "Rodents",
      species %in% c(
        "Leopardus_wiedii",
        "Puma_yagouaroundi",
        "Eira_barbara",
        "Galictis_vittata",
        "Lontra_longicaudis",
        "Potos_flavus"
      ) ~ "Carnivores",
      species %in% c(
        "Caluromys_philander",
        "Didelphis_marsupialis",
        "Marmosa_lepida",
        "Marmosa_murina",
        "Marmosops_parvidens",
        "Metachirus_nudicaudatus",
        "Micoureus_demerarae",
        "Philander_opossum"
      ) ~ "Didelphids"
    ),
    species_label = gsub(
      "_",
      " ",
      as.character(species)
    ),
    species_label = factor(
      species_label,
      levels = rev(
        gsub("_", " ", species_order)
      )
    )
  )

group_colors <- c(
  "Primates" = "#264478",
  "Xenarthrans" = "#C65911",
  "Armadillos" = "#666666",
  "Rodents" = "#D6A500",
  "Carnivores" = "#375623",
  "Didelphids" = "#4472C4"
)

species_labels <- plot_data %>%
  distinct(species_label, group) %>%
  mutate(
    label = paste0(
      "<span style='color:",
      group_colors[group],
      "'><i>",
      as.character(species_label),
      "</i></span>"
    )
  ) %>%
  select(species_label, label)

label_vector <- setNames(
  species_labels$label,
  species_labels$species_label
)

p_species_prevalence <- ggplot(
  plot_data,
  aes(
    x = prevalence * 100,
    y = species_label,
    colour = group
  )
) +
  geom_segment(
    aes(
      x = ci_low * 100,
      xend = ci_high * 100,
      y = species_label,
      yend = species_label
    ),
    linewidth = 0.8
  ) +
  geom_point(
    shape = 21,
    fill = "white",
    size = 5.5,
    stroke = 0.9
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  scale_y_discrete(
    labels = label_vector
  ) +
  scale_colour_manual(
    values = group_colors
  ) +
  labs(
    x = "Hemoplasma prevalence",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "Calibri"
    ),
    
    # Larger species names on Y axis
    axis.text.y = ggtext::element_markdown(
      size = 15,
      margin = margin(r = 10),
      family = "Calibri"
    ),
    
    # Larger X-axis labels
    axis.text.x = element_text(
      size = 13,
      family = "Calibri"
    ),
    
    # Larger X-axis title
    axis.title.x = element_text(
      size = 14,
      family = "Calibri"
    ),
    
    panel.grid.major.x = element_line(
      linewidth = 0.3,
      colour = "grey85"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    legend.position = "none"
  )

print(p_species_prevalence)

ggsave(
  filename = "species_Hemoplasma_prevalence_n15.png",
  plot = p_species_prevalence,
  width = 7,
  height = 8,
  units = "in",
  dpi = 300
)
```

### Test for the association between `hemoplasma` prevalence and sample size per `species`, for all `species` (n `species` = 44, complete dataset) (Spearman correlation test) 
```
cor.test(
  species_summary$n_sampled,
  species_summary$prevalence,
  method = "spearman",
  exact = FALSE
)
```

-> Results :
Spearman's ρ = 0.334, p = 0.027

-> Interpretation : 
A weak but significant positive association was detected between sample size and `species`-level `hemoplasma` prevalence (complete dataset).


### Test for the association between `hemoplasma` prevalence and sample size per `species`, restricted `species` with n ≥ 15 (n `species` = 16, conservative dataset) (Spearman correlation test) 
```
species_summary_n15 <- species_summary %>%
  filter(n_sampled >= 15)
cor_test_n15 <- cor.test(
  species_summary_n15$n_sampled,
  species_summary_n15$prevalence,
  method = "spearman",
  exact = FALSE
)
print(cor_test_n15)
```
-> Results:
Spearman's ρ = 0.370, p = 0.159

-> Interpretation:
No significant association was detected between sample size and `species`-level `hemoplasma` prevalence among `species` with n ≥ 15 (conservative dataset).

### Visualization of the association between `hemoplasma` prevalence and sample size per `species` (Fig. S1) 
```
windowsFonts(Calibri = windowsFont("Calibri"))
order_colors <- c(
  "Primates (monkeys)" = "#264478",
  "Pilosa (sloths, anteaters)" = "#C65911",
  "Cingulata (armadillos)" = "#666666",
  "Rodentia (rodents)" = "#D6A500",
  "Carnivora (carnivorans)" = "#375623",
  "Didelphimorphia (opossums)" = "#4472C4"
)
plot_data <- species_summary %>%
  filter(
    is.finite(n_sampled),
    is.finite(prevalence),
    n_sampled > 0
  ) %>%
  mutate(
    log_n = log10(n_sampled),
    order = case_when(
      species %in% c(
        "Alouatta_macconnelli",
        "Saguinus_midas",
        "Sapajus_apella",
        "Saimiri_sciureus",
        "Pithecia_pithecia"
      ) ~ "Primates (monkeys)",
      species %in% c(
        "Bradypus_tridactylus",
        "Choloepus_didactylus",
        "Cyclopes_didactylus",
        "Tamandua_tetradactyla"
      ) ~ "Pilosa (sloths, anteaters)",
      species %in% c(
        "Cabassous_unicinctus",
        "Dasypus_novemcinctus"
      ) ~ "Cingulata (armadillos)",
      species %in% c(
        "Hydrochoerus_hydrochaeris",
        "Holochilus_sciureus",
        "Hylaeamys_megacephalus",
        "Hylaeamys_yunganus",
        "Neacomys_dubosti",
        "Neacomys_paracou",
        "Nectomys_rattus",
        "Oecomys_auyantepui",
        "Oecomys_bicolor",
        "Oligoryzomys_fulvescens",
        "Makalata_didelphoides",
        "Mesomys_hispidus",
        "Proechimys_cuvieri",
        "Proechimys_guyannensis",
        "Coendou_melanurus",
        "Coendou_prehensilis",
        "Mus_musculus",
        "Rattus_rattus",
        "Sciurus_aestuans"
      ) ~ "Rodentia (rodents)",
      species %in% c(
        "Leopardus_wiedii",
        "Puma_yagouaroundi",
        "Eira_barbara",
        "Galictis_vittata",
        "Lontra_longicaudis",
        "Potos_flavus"
      ) ~ "Carnivora (carnivorans)",
      species %in% c(
        "Caluromys_philander",
        "Didelphis_marsupialis",
        "Marmosa_lepida",
        "Marmosa_murina",
        "Marmosops_parvidens",
        "Metachirus_nudicaudatus",
        "Micoureus_demerarae",
        "Philander_opossum"
      ) ~ "Didelphimorphia (opossums)",
      TRUE ~ NA_character_
    )
  )
conservative_species <- species_summary %>%
  filter(n_sampled >= 15) %>%
  pull(species)
print(conservative_species)
print(paste("Number of species with n >= 15:", length(conservative_species)))
species_abbr <- setNames(
  paste0(
    substr(sub("_.*", "", conservative_species), 1, 1),
    substr(sub("^[^_]+_", "", conservative_species), 1, 1)
  ),
  conservative_species
)
plot_data <- plot_data %>%
  mutate(
    species_abbr = unname(
      species_abbr[as.character(species)]
    )
  )
p <- ggplot(
  plot_data,
  aes(
    x = log_n,
    y = prevalence
  )
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    fill = "grey70",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    data = plot_data %>% filter(n_sampled < 15),
    aes(
      colour = order,
      fill = order
    ),
    shape = 21,
    size = 4,
    stroke = 0.9
  ) +
  geom_point(
    data = plot_data %>% filter(n_sampled >= 15),
    aes(
      colour = order
    ),
    shape = 21,
    fill = "white",
    size = 4,
    stroke = 1.1
  ) +
  geom_text(
    data = plot_data %>%
      filter(
        n_sampled >= 15,
        !is.na(species_abbr)
      ),
    aes(
      label = species_abbr,
      colour = order
    ),
    hjust = -0.5,
    vjust = 0.5,
    size = 4,
    family = "Calibri",
    show.legend = FALSE
  ) +
  scale_colour_manual(
    name = "Host order",
    values = order_colors
  ) +
  scale_fill_manual(
    values = order_colors,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = log10(c(1, 2, 5, 10, 20, 50, 100)),
    labels = c(1, 2, 5, 10, 20, 50, 100),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Sample size per species",
    y = "Hemoplasma prevalence"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Calibri"),
    axis.text.x = element_text(size = 13, family = "Calibri"),
    axis.text.y = element_text(size = 13, family = "Calibri"),
    axis.title.x = element_text(size = 14, family = "Calibri"),
    axis.title.y = element_text(size = 14, family = "Calibri"),
    legend.title = element_text(size = 13, family = "Calibri"),
    legend.text = element_text(size = 12, family = "Calibri"),
    legend.position = "right"
  )
print(p)
ggsave(
  filename = "hemoplasma_prevalence_vs_sample_size_n15.png",
  plot = p,
  width = 9,
  height = 5.5,
  units = "in",
  dpi = 300
)
```

## Step 3. Variation in `hemoplasma` infection status according to the host `sex`
### Test whether `hemoplasma` infection probability differs between `sex` while accounting for species-level random effects (`1 | species`) (complete dataset, 44 `species`) 
Fit the full GLMM (model 1) :
```
model_sex_data <- data_hemoplasma_stat %>%
  filter(
    !is.na(sex),
  )
model1_a <- glmer(
  hemoplasma ~ sex + (1 | species),
  data = model_sex_data,
  family = binomial
)
summary(model1_a)
model1_b <- glmer(
  hemoplasma ~ 1 + (1 | species),
  data = model_sex_data,
  family = binomial,
)
anova(
  model1_a,
  model1_b,
  test = "Chisq"
)
AIC(model1_a, model1_b)
```

Calculate the odds ratio and 95% HDI for the effect of `sex` on `hemoplasma` infection (complete dataset, 44 `species`)
```
model_sex_bayes <- brm(
  hemoplasma ~ sex + (1 | species),
  data = model_sex_data,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 1234
)
posterior_sex <- as_draws_df(
  model_sex_bayes
)
or_sex <- exp(
  posterior_sex$b_sexM
)
or_sex_results <- data.frame(
    variable = "Sex (M vs F)",
    OR = median(
    or_sex
  ),
    HDI_low = hdi(
    or_sex,
    ci = 0.95
  )$CI_low,
    HDI_high = hdi(
    or_sex,
    ci = 0.95
  )$CI_high
)
or_sex_results
```

### Test whether `hemoplasma` infection probability differs between `sex` while accounting for `species`-level random effects (1 | `species`) (conservative dataset, 16 `species`)
Fit the full GLMM (model 1_n15) :
```
model_sex_data_n15 <- data_hemoplasma_stat %>%
  filter(
    !is.na(sex),
    species %in% (
      species_summary %>%
        filter(n_sampled >= 15) %>%
        pull(species)
    )
  )

model1_a_n15 <- glmer(
  hemoplasma ~ sex + (1 | species),
  data = model_sex_data_n15,
  family = binomial
)

summary(model1_a_n15)

model1_b_n15 <- glmer(
  hemoplasma ~ 1 + (1 | species),
  data = model_sex_data_n15,
  family = binomial
)

anova(
  model1_a_n15,
  model1_b_n15,
  test = "Chisq"
)

AIC(model1_a_n15, model1_b_n15)
```
Calculate the odds ratio and 95% HDI for the effect of `sex` on `hemoplasma` infection (conservative species-level dataset, 16 `species`)
```
model_sex_data_n15 <- data_hemoplasma_stat %>%
  filter(
    !is.na(sex),
    species %in% (
      species_summary %>%
        filter(n_sampled >= 15) %>%
        pull(species)
    )
  )

model_sex_bayes_n15 <- brm(
  hemoplasma ~ sex + (1 | species),
  data = model_sex_data_n15,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 1234
)

posterior_sex_n15 <- as_draws_df(
  model_sex_bayes_n15
)

or_sex_n15 <- exp(
  posterior_sex_n15$b_sexM
)

or_sex_results_n15 <- data.frame(
  variable = "Sex (M vs F)",
  OR = median(or_sex_n15),
  HDI_low = hdi(
    or_sex_n15,
    ci = 0.95
  )$CI_low,
  HDI_high = hdi(
    or_sex_n15,
    ci = 0.95
  )$CI_high
)

or_sex_results_n15
```

-> Results : For the complete mammal dataset (44 `species`), `sex` was not significantly associated with `hemoplasma` infection status (LRT: χ²₁ = 2.33, p = 0.127). The model including `sex` had a slightly lower AIC than the null model (314.73 vs. 315.06; ΔAIC = 0.33). Males showed higher estimated odds of infection than females (OR = 1.69, 95% HDI: 0.77–3.13), but the HDI included 1. For the conservative species-level dataset (16 `species`), the association remained non-significant (LRT: χ²₁ = 2.90, p = 0.088; ΔAIC = 0.90), with males again showing higher estimated odds of infection (OR = 1.85, 95% HDI: 0.73–3.49).

-> Interpretation : There was no significant evidence for a `sex` effect on `hemoplasma` infection probability in either the complete dataset or the conservative species-level dataset. The consistent tendency towards higher infection odds in males suggests a possible weak `sex` effect, but uncertainty remained substantial.

## Step 4. Variation in `hemoplasma` infection status according to the presence of other blood-borne pathogens (`anaplasmataceae` and `apicomplexa`)

### Data preparation: Create the `pathogens` variable by merging `anaplasmataceae` and `apicomplexa` (0 = uninfected ; 1 = infected by `anaplasmataceae` and/or `apicomplexa`)
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  mutate(
    pathogens = ifelse(
      anaplasmataceae == 1 | apicomplexa == 1,
      1, 0
    ),
    species = factor(species)
  )
data_hemoplasma_stat$pathogens <- as.factor(data_hemoplasma_stat$pathogens)
```

### Test whether `hemoplasma` infection probability differs with infections by other blood-borne pathogens (`pathogens`) while accounting for species-level random effects (`1 | species`) (complete dataset, 44 `species`)
Fit the full GLMM (model 2)
```
model2_a <- glmer(
  hemoplasma ~ pathogens + (1 | species),
  data = data_hemoplasma_stat,
  family = binomial
)
summary(model2_a)
model2_b <- glmer(
  hemoplasma ~ 1 + (1 | species),
  data = data_hemoplasma_stat,
  family = binomial,
)
anova(
  model2_a,
  model2_b,
  test = "Chisq"
)
AIC(model2_a, model2_b)
```

Calculate the odds ratio and 95% HDI for the effect of `pathogens` on `hemoplasma` infection (complete dataset, 44 `species`)
```
model_pathogens_bayes <- brm(
  hemoplasma ~ pathogens + (1 | species),
  data = data_hemoplasma_stat,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 1,
  seed = 1234
)
posterior_pathogens <- as_draws_df(
  model_pathogens_bayes
)
or_pathogens <- exp(
  posterior_pathogens$b_pathogens1
)
or_pathogens_results <- data.frame(
    variable = "Pathogens (1 vs 0)",
    OR = median(
    or_pathogens
  ),
    HDI_low = hdi(
    or_pathogens,
    ci = 0.95
  )$CI_low,
  
  HDI_high = hdi(
    or_pathogens,
    ci = 0.95
  )$CI_high
)
or_pathogens_results
```

### Test whether `hemoplasma` infection probability differs with infections by other blood-borne pathogens (`pathogens`) while accounting for species-level random effects (`1 | species`) (conservative dataset, 16 `species`)
Fit the full GLMM (model 2_n15)
```
species_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    .groups = "drop"
  ) %>%
  filter(n_sampled >= 15)

model2_data_n15 <- data_hemoplasma_stat %>%
  filter(species %in% species_n15$species)

model2_a_n15 <- glmer(
  hemoplasma ~ pathogens + (1 | species),
  data = model2_data_n15,
  family = binomial
)

summary(model2_a_n15)

model2_b_n15 <- glmer(
  hemoplasma ~ 1 + (1 | species),
  data = model2_data_n15,
  family = binomial
)

anova(
  model2_a_n15,
  model2_b_n15,
  test = "Chisq"
)

AIC(
  model2_a_n15,
  model2_b_n15
)
```

Calculate the odds ratio and 95% HDI for the effect of `pathogens` on `hemoplasma` infection (conservative dataset, 16 `species`)
```
species_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    .groups = "drop"
  ) %>%
  filter(n_sampled >= 15)

model_pathogens_bayes_n15 <- brm(
  hemoplasma ~ pathogens + (1 | species),
  data = data_hemoplasma_stat %>%
    filter(species %in% species_n15$species),
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 1,
  seed = 1234
)

posterior_pathogens_n15 <- as_draws_df(
  model_pathogens_bayes_n15
)

or_pathogens_n15 <- exp(
  posterior_pathogens_n15$b_pathogens1
)

or_pathogens_results_n15 <- data.frame(
  variable = "Pathogens (1 vs 0)",
  OR = median(
    or_pathogens_n15
  ),
  HDI_low = hdi(
    or_pathogens_n15,
    ci = 0.95
  )$CI_low,
  HDI_high = hdi(
    or_pathogens_n15,
    ci = 0.95
  )$CI_high
)

or_pathogens_results_n15
```
-> Results: `pathogen` was significantly associated with `hemoplasma` infection in both the complete dataset (GLMM, LRT: χ²₁ = 7.56, p = 0.006, ΔAIC = 5.56; OR = 2.85, 95% HDI: 1.07–5.71) and the conservative dataset (χ²₁ = 7.94, p = 0.005, ΔAIC = 5.94; OR = 3.03, 95% HDI: 1.10–6.09).

-> Interpretation: Individuals positive for other blood-borne `pathogen` had approximately three-fold higher odds of `hemoplasma` infection, with consistent evidence for this association in both datasets.

### Test whether `hemoplasma` infection probability differs with infections between `anaplasmataceae` and `apicomplexa` while accounting for species-level random effects (1 | `species`) (complete dataset, 44 `species`)
Fit the full GLMM (model 3)
```
model3_a <- glmer(
  hemoplasma ~ anaplasmataceae * apicomplexa + (1 | species),
  data = data_hemoplasma_stat,
  family = binomial
)
summary(model3_a)
model3_b <- glmer(
  hemoplasma ~ anaplasmataceae + apicomplexa + (1 | species),
  data = data_hemoplasma_stat,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
anova(
  model3_b,
  model3_a,
  test = "Chisq"
)
AIC(model3_a, model3_b)
model3_c <- glmer(
  hemoplasma ~ anaplasmataceae + (1 | species),
  data = data_hemoplasma_stat,
  family = binomial,
)
anova(
  model3_b,
  model3_c,
  test = "Chisq"
)
AIC(model3_b, model3_c)
model3_d <- glmer(
  hemoplasma ~ (1 | species),
  data = data_hemoplasma_stat,
  family = binomial,
)
anova(
  model3_c,
  model3_d,
  test = "Chisq"
)
AIC(model3_c, model3_d)
```

Calculate the odds ratio and 95% HDI for the effect of `anaplasmataceae` and `apicomplexa` on `hemoplasma` infection (complete dataset, 44 `species`)
```
model3_bayes <- brm(
  hemoplasma ~ anaplasmataceae + apicomplexa + (1 | species),
  data = data_hemoplasma_stat,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 1,
  seed = 1234
)
posterior_model3 <- as_draws_df(
  model3_bayes
)
or_anaplasmataceae <- exp(
  posterior_model3$b_anaplasmataceae1
)

# Odds ratio for Apicomplexa
or_apicomplexa <- exp(
  posterior_model3$b_apicomplexa1
)
or_model3_results <- data.frame(
  variable = c(
    "Anaplasmataceae (1 vs 0)",
    "Apicomplexa (1 vs 0)"
  ),
  OR = c(
    median(or_anaplasmataceae),
    median(or_apicomplexa)
  ),
  HDI_low = c(
    hdi(
      or_anaplasmataceae,
      ci = 0.95
    )$CI_low,
    hdi(
      or_apicomplexa,
      ci = 0.95
    )$CI_low
  ),
  HDI_high = c(
    hdi(
      or_anaplasmataceae,
      ci = 0.95
    )$CI_high,
    hdi(
      or_apicomplexa,
      ci = 0.95
    )$CI_high
  )
)
or_model3_results
```

### Test whether `hemoplasma` infection probability differs with infections between `anaplasmataceae` and `apicomplexa` while accounting for species-level random effects (1 | `species`) (conservative dataset, 16 `species`)
Fit the full GLMM (model 3_n15)
```
model3_data_n15 <- data_hemoplasma_stat %>%
  filter(
    species %in% (
      species_summary %>%
        filter(n_sampled >= 15) %>%
        pull(species)
    )
  )

model3_a_n15 <- glmer(
  hemoplasma ~ anaplasmataceae * apicomplexa + (1 | species),
  data = model3_data_n15,
  family = binomial
)

summary(model3_a_n15)

model3_b_n15 <- glmer(
  hemoplasma ~ anaplasmataceae + apicomplexa + (1 | species),
  data = model3_data_n15,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

anova(
  model3_b_n15,
  model3_a_n15,
  test = "Chisq"
)

AIC(model3_a_n15, model3_b_n15)

model3_c_n15 <- glmer(
  hemoplasma ~ anaplasmataceae + (1 | species),
  data = model3_data_n15,
  family = binomial
)

anova(
  model3_b_n15,
  model3_c_n15,
  test = "Chisq"
)

AIC(model3_b_n15, model3_c_n15)

model3_d_n15 <- glmer(
  hemoplasma ~ (1 | species),
  data = model3_data_n15,
  family = binomial
)

anova(
  model3_c_n15,
  model3_d_n15,
  test = "Chisq"
)

AIC(model3_c_n15, model3_d_n15)
```

Calculate the odds ratio and 95% HDI for the effect of `anaplasmataceae` and `apicomplexa` on `hemoplasma` infection (conservative dataset, 16 `species`)
```
model3_bayes_n15 <- brm(
  hemoplasma ~ anaplasmataceae + apicomplexa + (1 | species),
  data = model3_data_n15,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 1,
  seed = 1234
)

posterior_model3_n15 <- as_draws_df(
  model3_bayes_n15
)

or_anaplasmataceae_n15 <- exp(
  posterior_model3_n15$b_anaplasmataceae1
)

or_apicomplexa_n15 <- exp(
  posterior_model3_n15$b_apicomplexa1
)

or_model3_results_n15 <- data.frame(
  variable = c(
    "Anaplasmataceae (1 vs 0)",
    "Apicomplexa (1 vs 0)"
  ),
  OR = c(
    median(or_anaplasmataceae_n15),
    median(or_apicomplexa_n15)
  ),
  HDI_low = c(
    hdi(
      or_anaplasmataceae_n15,
      ci = 0.95
    )$CI_low,
    hdi(
      or_apicomplexa_n15,
      ci = 0.95
    )$CI_low
  ),
  HDI_high = c(
    hdi(
      or_anaplasmataceae_n15,
      ci = 0.95
    )$CI_high,
    hdi(
      or_apicomplexa_n15,
      ci = 0.95
    )$CI_high
  )
)

or_model3_results_n15
```
-> Results : In the complete dataset (44 `species`), the `anaplasmataceae` × `apicomplexa` interaction did not significantly improve model fit (LRT: χ²₁ = 2.71, p = 0.100; ΔAIC = 0.71) and was therefore not retained. Adding `apicomplexa` to the model containing `anaplasmataceae` provided only weak evidence for improved fit (LRT: χ²₁ = 3.20, p = 0.074; ΔAIC = 1.20), whereas adding `anaplasmataceae` to the null model significantly improved model fit (LRT: χ²₁ = 7.01, p = 0.008; ΔAIC = 5.01). `Anaplasmataceae`-positive individuals had higher estimated odds of `hemoplasma` infection (OR = 3.09, 95% HDI: 0.85–7.24), although the HDI included 1; the corresponding effect of `apicomplexa` was also positive but uncertain (OR = 2.93, 95% HDI: 0.50–8.72). In the conservative dataset (16 `species`), the interaction was likewise not supported (LRT: χ²₁ = 2.39, p = 0.122; ΔAIC = 0.39). Adding `apicomplexa` provided weak evidence for improved fit (LRT: χ²₁ = 3.47, p = 0.063; ΔAIC = 1.47), while `anaplasmataceae` remained significantly associated with `hemoplasma` infection (LRT: χ²₁ = 7.09, p = 0.008; ΔAIC = 5.09). Estimated odds were similarly higher for `anaplasmataceae`-positive individuals (OR = 3.14, 95% HDI: 0.95–7.54) and for `Apicomplexa`-positive individuals (OR = 3.33, 95% HDI: 0.48–11.33), although both HDIs included 1.

-> Interpretation : Across both datasets, `anaplasmataceae` consistently contributed to explaining variation in `hemoplasma` infection probability, whereas evidence for an additional effect of `apicomplexa` remained weak. The similar effect estimates in the conservative dataset support the robustness of this pattern.

### Sensitivity analysis 
A leave-one-species-out analysis was further performed to assess whether the association between `pathogens` and `hemoplasma` infection was driven by any single `species`
```
species_list <- unique(
  data_hemoplasma_stat$species
)
leave_one_species <- lapply(
  species_list,
  function(sp) {    
    data_tmp <- data_hemoplasma_stat %>%
      filter(species != sp)    
    model_full <- glmer(
      hemoplasma ~ pathogens + (1 | species),
      data = data_tmp,
      family = binomial,
      control = glmerControl(
        optimizer = "bobyqa"
      )
    )    
    model_null <- glmer(
      hemoplasma ~ 1 + (1 | species),
      data = data_tmp,
      family = binomial,
      control = glmerControl(
        optimizer = "bobyqa"
      )
    )    
    coef_model <- summary(
      model_full
    )$coefficients[
      "pathogens1",
    ]    
    lrt <- anova(
      model_full,
      model_null,
      test = "Chisq"
    )    
    data.frame(
      excluded_species = sp,
      estimate = coef_model["Estimate"],
      SE = coef_model["Std. Error"],
      z = coef_model["z value"],
      p_value = coef_model["Pr(>|z|)"],
      OR = exp(
        coef_model["Estimate"]
      ),
      LRT_chisq = lrt$Chisq[2],
      LRT_p = lrt$`Pr(>Chisq)`[2],
      AIC_full = AIC(
        model_full
      ),
      AIC_null = AIC(
        model_null
      ),
      delta_AIC = AIC(
        model_null
      ) - AIC(
        model_full
      )
    )
  }
)
leave_one_species_results <- bind_rows(
  leave_one_species
)
leave_one_species_results
```

-> Table of the leave-one-species-out sensitivity analysis :
| Excluded species | Estimate | SE | z | p-value | OR | LRT χ² | LRT p | AIC full | AIC null | ΔAIC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *Caluromys philander* | 1.028 | 0.386 | 2.660 | 0.008 | 2.795 | 7.511 | 0.006 | 447.314 | 452.825 | 5.511 |
| *Didelphis marsupialis* | 0.735 | 0.427 | 1.722 | 0.085 | 2.084 | 3.056 | 0.080 | 377.008 | 378.064 | 1.056 |
| *Marmosa lepida* | 1.039 | 0.387 | 2.685 | 0.007 | 2.828 | 7.662 | 0.006 | 444.964 | 450.627 | 5.662 |
| *Marmosa murina* | 1.036 | 0.388 | 2.670 | 0.008 | 2.817 | 7.579 | 0.006 | 428.465 | 434.045 | 5.579 |
| *Marmosops parvidens* | 1.028 | 0.386 | 2.660 | 0.008 | 2.795 | 7.511 | 0.006 | 447.314 | 452.825 | 5.511 |
| *Metachirus nudicaudatus* | 1.028 | 0.386 | 2.660 | 0.008 | 2.795 | 7.511 | 0.006 | 447.314 | 452.825 | 5.511 |
| *Micoureus demerarae* | 1.035 | 0.388 | 2.669 | 0.008 | 2.815 | 7.571 | 0.006 | 433.727 | 439.298 | 5.571 |
| *Philander opossum* | 1.143 | 0.399 | 2.865 | 0.004 | 3.135 | 8.840 | 0.003 | 415.227 | 422.067 | 6.840 |
| *Bradypus tridactylus* | 1.283 | 0.420 | 3.056 | 0.002 | 3.609 | 10.123 | 0.001 | 408.205 | 416.327 | 8.123 |
| *Choloepus didactylus* | 1.033 | 0.543 | 1.905 | 0.057 | 2.811 | 3.861 | 0.049 | 355.046 | 356.906 | 1.861 |
| *Cyclopes didactylus* | 1.030 | 0.387 | 2.664 | 0.008 | 2.801 | 7.536 | 0.006 | 447.993 | 453.529 | 5.536 |
| *Tamandua tetradactyla* | 1.029 | 0.387 | 2.661 | 0.008 | 2.797 | 7.519 | 0.006 | 447.576 | 453.095 | 5.519 |
| *Cabassous unicinctus* | 1.029 | 0.387 | 2.662 | 0.008 | 2.799 | 7.525 | 0.006 | 447.753 | 453.278 | 5.525 |
| *Dasypus novemcinctus* | 0.861 | 0.399 | 2.157 | 0.031 | 2.366 | 4.872 | 0.027 | 428.330 | 431.202 | 2.872 |
| *Hydrochoerus hydrochaeris* | 1.052 | 0.389 | 2.704 | 0.007 | 2.864 | 7.787 | 0.005 | 447.491 | 453.278 | 5.787 |
| *Holochilus sciureus* | 1.036 | 0.388 | 2.670 | 0.008 | 2.817 | 7.580 | 0.006 | 441.417 | 446.997 | 5.580 |
| *Hylaeamys megacephalus* | 1.025 | 0.386 | 2.657 | 0.008 | 2.788 | 7.492 | 0.006 | 446.598 | 452.090 | 5.492 |
| *Hylaeamys yunganus* | 1.026 | 0.386 | 2.658 | 0.008 | 2.791 | 7.499 | 0.006 | 446.888 | 452.387 | 5.499 |
| *Neacomys dubosti* | 1.030 | 0.387 | 2.664 | 0.008 | 2.801 | 7.536 | 0.006 | 447.993 | 453.529 | 5.536 |
| *Neacomys paracou* | 1.027 | 0.386 | 2.659 | 0.008 | 2.792 | 7.503 | 0.006 | 447.034 | 452.537 | 5.503 |
| *Nectomys rattus* | 1.039 | 0.388 | 2.679 | 0.007 | 2.825 | 7.629 | 0.006 | 439.885 | 445.515 | 5.629 |
| *Oecomys auyantepui* | 1.033 | 0.388 | 2.665 | 0.008 | 2.809 | 7.546 | 0.006 | 438.896 | 444.442 | 5.546 |
| *Oecomys bicolor* | 1.025 | 0.386 | 2.657 | 0.008 | 2.788 | 7.491 | 0.006 | 446.549 | 452.040 | 5.490 |
| *Oligoryzomys fulvescens* | 1.035 | 0.388 | 2.669 | 0.008 | 2.815 | 7.570 | 0.006 | 440.700 | 446.270 | 5.569 |
| *Makalata didelphoides* | 1.027 | 0.386 | 2.659 | 0.008 | 2.792 | 7.503 | 0.006 | 447.034 | 452.537 | 5.503 |
| *Mesomys hispidus* | 1.026 | 0.386 | 2.657 | 0.008 | 2.789 | 7.495 | 0.006 | 446.704 | 452.198 | 5.495 |
| *Proechimys cuvieri* | 1.035 | 0.388 | 2.668 | 0.008 | 2.815 | 7.568 | 0.006 | 433.228 | 438.796 | 5.568 |
| *Proechimys guyannensis* | 1.032 | 0.387 | 2.664 | 0.008 | 2.807 | 7.540 | 0.006 | 438.392 | 443.932 | 5.540 |
| *Coendou melanurus* | 1.030 | 0.387 | 2.664 | 0.008 | 2.801 | 7.536 | 0.006 | 447.993 | 453.529 | 5.536 |
| *Coendou prehensilis* | 1.096 | 0.394 | 2.782 | 0.005 | 2.991 | 8.298 | 0.004 | 441.807 | 448.105 | 6.298 |
| *Mus musculus* | 1.023 | 0.385 | 2.655 | 0.008 | 2.782 | 7.476 | 0.006 | 445.911 | 451.387 | 5.476 |
| *Rattus rattus* | 1.035 | 0.388 | 2.668 | 0.008 | 2.814 | 7.566 | 0.006 | 433.000 | 438.566 | 5.566 |
| *Sciurus aestuans* | 1.030 | 0.387 | 2.664 | 0.008 | 2.801 | 7.536 | 0.006 | 447.993 | 453.529 | 5.536 |
| *Leopardus wiedii* | 1.058 | 0.389 | 2.719 | 0.007 | 2.881 | 7.878 | 0.005 | 447.651 | 453.529 | 5.878 |
| *Puma yagouaroundi* | 1.028 | 0.386 | 2.660 | 0.008 | 2.795 | 7.511 | 0.006 | 447.314 | 452.825 | 5.511 |
| *Eira barbara* | 1.046 | 0.389 | 2.689 | 0.007 | 2.845 | 7.699 | 0.006 | 447.249 | 452.948 | 5.699 |
| *Galictis vittata* | 1.004 | 0.389 | 2.580 | 0.010 | 2.730 | 7.048 | 0.008 | 440.708 | 445.756 | 5.048 |
| *Lontra longicaudis* | 1.039 | 0.387 | 2.685 | 0.007 | 2.828 | 7.662 | 0.006 | 444.965 | 450.627 | 5.662 |
| *Potos flavus* | 1.038 | 0.388 | 2.677 | 0.007 | 2.823 | 7.617 | 0.006 | 443.389 | 449.005 | 5.617 |
| *Alouatta macconnelli* | 1.042 | 0.386 | 2.702 | 0.007 | 2.836 | 7.752 | 0.005 | 428.508 | 434.260 | 5.752 |
| *Saguinus midas* | 1.045 | 0.382 | 2.737 | 0.006 | 2.842 | 7.930 | 0.005 | 436.904 | 442.834 | 5.930 |
| *Sapajus apella* | 1.030 | 0.387 | 2.664 | 0.008 | 2.801 | 7.536 | 0.006 | 447.993 | 453.529 | 5.536 |
| *Saimiri sciureus* | 1.030 | 0.387 | 2.664 | 0.008 | 2.801 | 7.536 | 0.006 | 447.993 | 453.529 | 5.536 |
| *Pithecia pithecia* | 1.030 | 0.387 | 2.664 | 0.008 | 2.801 | 7.536 | 0.006 | 447.993 | 453.529 | 5.536 |

-> Results: Leave-one-species-out analyses showed a consistently positive association between `hemoplasma` and `pathogens`, with odds ratios ranging from 2.08 to 3.61. The association remained supported by likelihood-ratio tests in 42/44 species-exclusion models (LRT *p* < 0.05), with highly similar estimates in most cases (OR ≈ 2.8–2.9). The strongest association was observed when *Bradypus tridactylus* was excluded (OR = 3.61, LRT *p* = 0.0015), whereas excluding *Didelphis marsupialis* reduced statistical support (OR = 2.08, LRT *p* = 0.080). 

-> Interpretation: The positive association between `hemoplasma` and `pathogens` was generally robust to the exclusion of individual host `species` and was not driven by a single mammal `species`. Excluding *Bradypus tridactylus* strengthened the estimated association, indicating that this `species` tends to attenuate the overall effect. 

### Visualization of odds ratios and 95% HDIs for `sex`, `pathogens`, `apicomplexa` and `anaplasmataceae` for complete (44 `species`) and conservative dataset (16 `species`)
```
or_results <- data.frame(
  variable = c(
    "Sex (M vs F)",
    "Pathogens (1 vs 0)",
    "Anaplasmataceae (1 vs 0)",
    "Apicomplexa (1 vs 0)"
  ),
  OR = c(
    1.69,
    2.85,
    3.09,
    2.93
  ),
  HDI_low = c(
    0.77,
    1.07,
    0.85,
    0.50
  ),
  HDI_high = c(
    3.13,
    5.71,
    7.24,
    8.72
  ),
  dataset = "Complete dataset (44 species)"
)

or_results_n15 <- data.frame(
  variable = c(
    "Sex (M vs F)",
    "Pathogens (1 vs 0)",
    "Anaplasmataceae (1 vs 0)",
    "Apicomplexa (1 vs 0)"
  ),
  OR = c(
    1.85,
    3.03,
    3.14,
    3.33
  ),
  HDI_low = c(
    0.73,
    1.10,
    0.95,
    0.48
  ),
  HDI_high = c(
    3.49,
    6.09,
    7.54,
    11.33
  ),
  dataset = "Conservative dataset (16 species)"
)

plot_or <- bind_rows(
  or_results,
  or_results_n15
) %>%
  mutate(
    variable = factor(
      variable,
      levels = c(
        "Apicomplexa (1 vs 0)",
        "Anaplasmataceae (1 vs 0)",
        "Pathogens (1 vs 0)",
        "Sex (M vs F)"
      )
    ),
    y_base = as.numeric(variable),
    y = ifelse(
      dataset == "Complete dataset (44 species)",
      y_base + 0.08,
      y_base - 0.08
    )
  )

p <- ggplot(
  plot_or,
  aes(
    x = OR,
    y = y
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_segment(
    aes(
      x = HDI_low,
      xend = HDI_high,
      y = y,
      yend = y
    ),
    colour = "black",
    linewidth = 1.2
  ) +
  geom_point(
    data = subset(
      plot_or,
      dataset == "Complete dataset (44 species)"
    ),
    aes(
      x = OR,
      y = y
    ),
    shape = 21,
    colour = "black",
    fill = "black",
    size = 5,
    stroke = 1.2
  ) +
  geom_point(
    data = subset(
      plot_or,
      dataset == "Conservative dataset (16 species)"
    ),
    aes(
      x = OR,
      y = y
    ),
    shape = 21,
    colour = "black",
    fill = "white",
    size = 5,
    stroke = 1.2
  ) +
  scale_y_continuous(
    breaks = 1:4,
    labels = c(
      "Apicomplexa (1 vs 0)",
      "Anaplasmataceae (1 vs 0)",
      "Pathogens (1 vs 0)",
      "Sex (M vs F)"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 12),
    breaks = seq(0, 12, by = 2)
  ) +
  labs(
    x = "Odds ratio (95% HDI)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "Calibri"
    ),
    axis.text.y = element_text(
      size = 14,
      family = "Calibri"
    ),
    axis.text.x = element_text(
      size = 13,
      family = "Calibri"
    ),
    axis.title.x = element_text(
      size = 15,
      family = "Calibri"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.9
    )
  )

print(p)

ggsave(
  filename = "OR_hemoplasma_complete_vs_n15.png",
  plot = p,
  width = 7,
  height = 4.8,
  units = "in",
  dpi = 300
)
```

### Visualization of `hemoplasma` prevalence with 95% confidence intervals (CI; Wilson method) in *Bradypus tridactylus* and *Didelphis marsupialis* according to `pathogens` (`anaplasmataceae` + `apicomplexa`) infection status
```
species_select <- c(
  "Didelphis_marsupialis",
  "Bradypus_tridactylus"
)
plot_data <- data_hemoplasma_stat %>%
  filter(
    species %in% species_select,
    !is.na(hemoplasma),
    !is.na(pathogens)
  ) %>%
  group_by(species, pathogens) %>%
  summarise(
    n_positive = sum(hemoplasma == 1),
    n_total = n(),
    prevalence = n_positive / n_total,
    .groups = "drop"
  ) %>%
  mutate(
    z = qnorm(0.975),
    denominator = 1 + z^2 / n_total,
    center = (
      prevalence +
        z^2 / (2 * n_total)
    ) / denominator,
    half_width = (
      z *
        sqrt(
          prevalence * (1 - prevalence) / n_total +
            z^2 / (4 * n_total^2)
        )
    ) / denominator,
    CI_low = center - half_width,
    CI_high = center + half_width,
    prevalence_percent = prevalence * 100,
    CI_low_percent = CI_low * 100,
    CI_high_percent = CI_high * 100,
    species_label = case_when(
      species == "Didelphis_marsupialis" ~ "Didelphis marsupialis",
      species == "Bradypus_tridactylus" ~ "Bradypus tridactylus"
    ),
    pathogens_label = case_when(
      pathogens == 0 ~ "infected by Anaplasmataceae/Apicomplexa",
      pathogens == 1 ~ "uninfected by Anaplasmataceae/Apicomplexa"
    ),
    group = case_when(
      species == "Bradypus_tridactylus" ~ "Xenarthrans",
      species == "Didelphis_marsupialis" ~ "Didelphids"
    )
  )
plot_data <- plot_data %>%
  mutate(
    group_y = factor(
      paste(species_label, pathogens_label),
      levels = c(
        "Bradypus tridactylus infected by Anaplasmataceae/Apicomplexa",
        "Bradypus tridactylus uninfected by Anaplasmataceae/Apicomplexa",
        "Didelphis marsupialis infected by Anaplasmataceae/Apicomplexa",
        "Didelphis marsupialis uninfected by Anaplasmataceae/Apicomplexa"
      )
    )
  )
print(
  plot_data %>%
    select(
      species_label,
      pathogens_label,
      group,
      n_positive,
      n_total,
      prevalence_percent,
      CI_low_percent,
      CI_high_percent
    )
)
sig_data <- data.frame(
  species = c(
    "Bradypus tridactylus",
    "Didelphis marsupialis"
  ),
  y = c(
    1.5,
    3.5
  ),
  label = c(
    "NS",
    "*"
  )
)
p <- ggplot(
  plot_data,
  aes(
    x = prevalence_percent,
    y = group_y,
    colour = group
  )
) +
  geom_errorbar(
    aes(
      xmin = CI_low_percent,
      xmax = CI_high_percent
    ),
    orientation = "y",
    width = 0,
    linewidth = 1.2
  ) +
  geom_point(
    shape = 21,
    fill = "white",
    size = 8,
    stroke = 1.3
  ) +
  geom_text(
    data = sig_data,
    aes(
      x = 103,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    colour = "black",
    family = "Calibri",
    size = 6,
    hjust = 0.5
  ) +
  scale_colour_manual(
    values = c(
      "Primates" = "#264478",
      "Xenarthrans" = "#C65911",
      "Armadillos" = "#666666",
      "Rodents" = "#D6A500",
      "Carnivores" = "#375623",
      "Didelphids" = "#4472C4"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 110),
    breaks = seq(0, 100, by = 20),
    expand = expansion(
      mult = c(0.02, 0.03)
    )
  ) +
  scale_y_discrete(
    labels = c(
      "<i>Bradypus tridactylus</i><br>infected by Anaplasmataceae/Apicomplexa",
      "<i>Bradypus tridactylus</i><br>uninfected by Anaplasmataceae/Apicomplexa",
      "<i>Didelphis marsupialis</i><br>infected by Anaplasmataceae/Apicomplexa",
      "<i>Didelphis marsupialis</i><br>uninfected by Anaplasmataceae/Apicomplexa"
    )
  ) +
  labs(
    x = "Hemoplasma prevalence (%)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "Calibri"
    ),
    axis.text.y = ggtext::element_markdown(
      size = 15,
      family = "Calibri",
      lineheight = 0.95
    ),
    axis.text.x = element_text(
      size = 14,
      family = "Calibri"
    ),
    axis.title.x = element_text(
      size = 16,
      family = "Calibri"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.9
    ),
    legend.position = "none"
  )
print(p)
ggsave(
  filename = "Hemoplasma_prevalence_Pathogens_Wilson_coloured.png",
  plot = p,
  width = 9,
  height = 5.5,
  units = "in",
  dpi = 300
)
```

### Fisher's exact test for association between `hemoplasma` and `pathogens` in *Bradypus tridactylus* and *Didelphis marsupialis*
```
species_select <- c(
  "Bradypus_tridactylus",
  "Didelphis_marsupialis"
)
for (sp in species_select) {  
  cat("\n========================================\n")
  cat(sp, "\n")
  cat("========================================\n")  
  data_sp <- data_hemoplasma_stat %>%
    filter(
      species == sp,
      !is.na(hemoplasma),
      !is.na(pathogens)
    )
    tab <- table(
    Hemoplasma = data_sp$hemoplasma,
    Pathogens = data_sp$pathogens
  )  
  print(tab)  
   fisher_result <- fisher.test(tab)  
  print(fisher_result)  
  cat(
    "\nOdds ratio = ",
    round(fisher_result$estimate, 3),
    "\nP-value = ",
    signif(fisher_result$p.value, 4),
    "\n"
  )
}
```

-> Results and interpretation : 
*Bradypus tridactylus* : No significant association between `pathogens` status and `haemoplasma` infection (Fisher’s exact test, OR = 0.56, 95% CI: 0.04–7.95, p = 0.619).
*Didelphis marsupialis* : `haemoplasma` infection was significantly associated with `pathogens` positivity (OR = 10.26, 95% CI: 1.15–499.19, p = 0.019), with `haemoplasma`-positive individuals ~10-fold more likely to be `pathogens`-positive.

## Step 5. `hemoplasma` prevalence by mammalian `order`
### Observed `hemoplasma` prevalence and 95% Wilson CI by `order` (complete dataset, 44 `species`) 
```
order_prevalence <- data_hemoplasma_stat %>%
  group_by(order) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    prevalence = n_positive / n_sampled,
    CI_low = binom.confint(
      n_positive,
      n_sampled,
      method = "wilson"
    )$lower,
    CI_high = binom.confint(
      n_positive,
      n_sampled,
      method = "wilson"
    )$upper
  ) %>%
  ungroup() %>%
  mutate(
    prevalence_percent = prevalence * 100,
    CI_low_percent = CI_low * 100,
    CI_high_percent = CI_high * 100
  )
order_prevalence
```

-> Results (complete dataset, 44 `species`) :
| Order | N sampled | N positive | Prevalence | 95% CI (Wilson) |
|---|---:|---:|---:|---:|
| Carnivora | 17 | 5 | 29.4% | 13.3–53.1% |
| Cingulata | 17 | 5 | 29.4% | 13.3–53.1% |
| Didelphimorphia | 123 | 39 | 31.7% | 24.1–40.4% |
| Pilosa | 202 | 76 | 37.6% | 31.2–44.5% |
| Primates | 66 | 61 | 92.4% | 83.5–96.7% |
| Rodentia | 201 | 11 | 5.47% | 3.08–9.53% |

### Observed `hemoplasma` prevalence and 95% Wilson CI by `order` (conservative dataset, 16 `species`) 
```
species_n15 <- species_summary %>%
  filter(n_sampled >= 15) %>%
  pull(species)

order_prevalence_n15 <- data_hemoplasma_stat %>%
  filter(species %in% species_n15) %>%
  group_by(order) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    prevalence = n_positive / n_sampled,
    CI_low = binom.confint(
      n_positive,
      n_sampled,
      method = "wilson"
    )$lower,
    CI_high = binom.confint(
      n_positive,
      n_sampled,
      method = "wilson"
    )$upper
  ) %>%
  ungroup() %>%
  mutate(
    prevalence_percent = prevalence * 100,
    CI_low_percent = CI_low * 100,
    CI_high_percent = CI_high * 100
  )

order_prevalence_n15
```

-> Results (conservative species-level dataset, 16 `species`) :
| Mammalian order | n sampled | n positive | Prevalence (%) | 95% CI (%) |
|---|---:|---:|---:|---:|
| Cingulata | 15 | 5 | 33.3 | 15.2–58.3 |
| Didelphimorphia | 107 | 38 | 35.5 | 27.1–44.9 |
| Pilosa | 198 | 76 | 38.4 | 31.9–45.3 |
| Primates | 63 | 61 | 96.8 | 89.1–99.1 |
| Rodentia | 138 | 6 | 4.35 | 2.01–9.16 |

### Visualization of `hemoplasma` prevalence and 95% Wilson CI by mammalian `order` for complete (44 `species`) and conservative (16 `species`) dataset
```
order_prevalence_complete <- data.frame(
  order = c(
    "Carnivora",
    "Cingulata",
    "Didelphimorphia",
    "Pilosa",
    "Primates",
    "Rodentia"
  ),
  n_sampled = c(
    17,
    17,
    123,
    202,
    66,
    201
  ),
  n_positive = c(
    5,
    5,
    39,
    76,
    61,
    11
  ),
  prevalence = c(
    29.4,
    29.4,
    31.7,
    37.6,
    92.4,
    5.47
  ),
  CI_low = c(
    13.3,
    13.3,
    24.1,
    31.2,
    83.5,
    3.08
  ),
  CI_high = c(
    53.1,
    53.1,
    40.4,
    44.5,
    96.7,
    9.53
  ),
  dataset = "Complete"
)

order_prevalence_n15 <- data.frame(
  order = c(
    "Cingulata",
    "Didelphimorphia",
    "Pilosa",
    "Primates",
    "Rodentia"
  ),
  n_sampled = c(
    15,
    107,
    198,
    63,
    138
  ),
  n_positive = c(
    5,
    38,
    76,
    61,
    6
  ),
  prevalence = c(
    33.3,
    35.5,
    38.4,
    96.8,
    4.35
  ),
  CI_low = c(
    15.2,
    27.1,
    31.9,
    89.1,
    2.01
  ),
  CI_high = c(
    58.3,
    44.9,
    45.3,
    99.1,
    9.16
  ),
  dataset = "Conservative"
)

order_levels <- c(
  "Rodentia",
  "Cingulata",
  "Carnivora",
  "Didelphimorphia",
  "Pilosa",
  "Primates"
)

order_prevalence <- bind_rows(
  order_prevalence_complete,
  order_prevalence_n15
) %>%
  mutate(
    order = factor(
      order,
      levels = order_levels
    ),
    dataset = factor(
      dataset,
      levels = c(
        "Conservative",
        "Complete"
      )
    ),
    y = as.numeric(order) +
      ifelse(
        dataset == "Complete",
        0.23,
        -0.23
      )
  )

order_colors <- c(
  "Primates" = "#264478",
  "Pilosa" = "#C65911",
  "Didelphimorphia" = "#4472C4",
  "Carnivora" = "#375623",
  "Cingulata" = "#666666",
  "Rodentia" = "#D6A500"
)

p_order_prevalence <- ggplot(
  order_prevalence,
  aes(
    x = prevalence,
    y = y
  )
) +
  geom_errorbar(
    aes(
      xmin = CI_low,
      xmax = CI_high,
      colour = order
    ),
    orientation = "y",
    width = 0,
    linewidth = 1.2
  ) +
  geom_point(
    data = subset(
      order_prevalence,
      dataset == "Complete"
    ),
    aes(
      colour = order,
      fill = order
    ),
    shape = 21,
    size = 8,
    stroke = 1.3
  ) +
  geom_point(
    data = subset(
      order_prevalence,
      dataset == "Conservative"
    ),
    aes(
      colour = order
    ),
    shape = 21,
    fill = "white",
    size = 8,
    stroke = 1.3
  ) +
  annotate(
    "point",
    x = 10,
    y = 6.65,
    shape = 21,
    fill = "black",
    colour = "black",
    size = 5,
    stroke = 1.2
  ) +
  annotate(
    "text",
    x = 14,
    y = 6.65,
    label = "Complete dataset (44 species)",
    hjust = 0,
    vjust = 0.5,
    size = 4,
    family = "Calibri"
  ) +
  annotate(
    "point",
    x = 10,
    y = 6.25,
    shape = 21,
    fill = "white",
    colour = "black",
    size = 5,
    stroke = 1.2
  ) +
  annotate(
    "text",
    x = 14,
    y = 6.25,
    label = "Conservative dataset (16 species)",
    hjust = 0,
    vjust = 0.5,
    size = 4,
    family = "Calibri"
  ) +
  scale_colour_manual(
    values = order_colors
  ) +
  scale_fill_manual(
    values = order_colors
  ) +
  scale_y_continuous(
    breaks = seq_along(order_levels),
    labels = order_levels,
    expand = expansion(
      mult = c(0.03, 0.15)
    )
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(
      mult = c(0.02, 0.03)
    )
  ) +
  labs(
    x = "Hemoplasma prevalence",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "Calibri"
    ),
    axis.text.y = element_text(
      size = 14,
      family = "Calibri"
    ),
    axis.text.x = element_text(
      size = 13,
      family = "Calibri"
    ),
    axis.title.x = element_text(
      size = 15,
      family = "Calibri"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.9
    ),
    legend.position = "none"
  )

print(p_order_prevalence)

ggsave(
  filename = "order_hemoplasma_prevalence_complete_vs_n15.png",
  plot = p_order_prevalence,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300
)
```

### Test whether `hemoplasma` infection probability differs between mammalian `order` (complete dataset, 44 `species`)
Fit the full GLMM (model 4) :
```
model4_a <- glmer(
  hemoplasma ~ order + (1 | species),
  data = data_hemoplasma_stat,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa"
  )
)
summary(model4_a)
model4_b <- glmer(
  hemoplasma ~ 1 + (1 | species),
  data = data_hemoplasma_stat,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa"
  )
)
order_test <- anova(
  model4_b,
  model4_a,
  test = "Chisq"
)
order_test
AIC(
model4_b,
  model4_a
)
```

### Post-hoc pairwise comparisons (odds ratios)
```
order_emmeans <- emmeans(
  model4_a,
  ~ order
)
order_OR <- pairs(
  order_emmeans,
  adjust = "tukey"
)
order_OR
order_OR_results <- summary(
  order_OR,
  infer = TRUE
) %>%
  as.data.frame() %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(asymp.LCL),
    CI_high = exp(asymp.UCL)
  ) %>%
  select(
    contrast,
    OR,
    CI_low,
    CI_high,
    p.value
  )
order_OR_results
```

-> Tukey-adjusted pairwise comparisons of `hemoplasma` infection odds among mammalian `order`
| Contrast | OR | 95% CI low | 95% CI high | Tukey-adjusted *p* |
|---|---:|---:|---:|---:|
| Carnivora – Cingulata | 1.646 | 0.010 | 277.081 | 0.9998 |
| Carnivora – Didelphimorphia | 1.685 | 0.050 | 57.025 | 0.9983 |
| Carnivora – Pilosa | 1.730 | 0.028 | 108.148 | 0.9990 |
| Carnivora – Primates | 0.099 | 0.002 | 6.570 | 0.6187 |
| Carnivora – Rodentia | 9.640 | 0.358 | 259.806 | 0.3654 |
| Cingulata – Didelphimorphia | 1.024 | 0.009 | 114.724 | 1.0000 |
| Cingulata – Pilosa | 1.051 | 0.006 | 184.983 | 1.0000 |
| Cingulata – Primates | 0.060 | 0.0003 | 10.997 | 0.6404 |
| Cingulata – Rodentia | 5.858 | 0.063 | 544.967 | 0.8769 |
| Didelphimorphia – Pilosa | 1.027 | 0.028 | 38.247 | 1.0000 |
| Didelphimorphia – Primates | 0.059 | 0.001 | 2.335 | 0.2411 |
| Didelphimorphia – Rodentia | 5.722 | 0.418 | 78.330 | 0.4022 |
| Pilosa – Primates | 0.057 | 0.001 | 3.927 | 0.3855 |
| Pilosa – Rodentia | 5.571 | 0.192 | 161.553 | 0.6939 |
| Primates – Rodentia | **96.957** | **3.168** | **2967.547** | **0.0019** |

### Test whether `hemoplasma` infection probability differs between mammalian `order` (conservative dataset, 16 `species`)
Fit the full GLMM (model 4_n15) :
```
model4_data_n15 <- data_hemoplasma_stat %>%
  filter(
    species %in% (
      species_summary %>%
        filter(n_sampled >= 15) %>%
        pull(species)
    )
  )

model4_a_n15 <- glmer(
  hemoplasma ~ order + (1 | species),
  data = model4_data_n15,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa"
  )
)

summary(model4_a_n15)

model4_b_n15 <- glmer(
  hemoplasma ~ 1 + (1 | species),
  data = model4_data_n15,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa"
  )
)

order_test_n15 <- anova(
  model4_b_n15,
  model4_a_n15,
  test = "Chisq"
)

order_test_n15

AIC(
  model4_b_n15,
  model4_a_n15
)
```

### Post-hoc pairwise comparisons (odds ratios) (conservative dataset, 16 `species`)
```
order_emmeans_n15 <- emmeans(
  model4_a_n15,
  ~ order
)

order_OR_n15 <- pairs(
  order_emmeans_n15,
  adjust = "tukey"
)

order_OR_n15

order_OR_results_n15 <- summary(
  order_OR_n15,
  infer = TRUE
) %>%
  as.data.frame() %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(asymp.LCL),
    CI_high = exp(asymp.UCL)
  ) %>%
  select(
    contrast,
    OR,
    CI_low,
    CI_high,
    p.value
  )

order_OR_results_n15
```

-> Tukey-adjusted pairwise comparisons of `hemoplasma` infection odds among mammalian `order` (conservative dataset, 16 `species`)
| Contrast | OR | 95% CI | Adjusted p-value |
|---|---:|---:|---:|
| Cingulata – Didelphimorphia | 1.28 | 0.021–78.86 | 0.9998 |
| Cingulata – Pilosa | 1.13 | 0.013–96.20 | >0.9999 |
| Cingulata – Primates | 0.011 | 0.00008–1.55 | 0.093 |
| Cingulata – Rodentia | 16.18 | 0.258–1014.72 | 0.353 |
| Didelphimorphia – Pilosa | 0.88 | 0.041–19.11 | >0.9999 |
| Didelphimorphia – Primates | 0.0086 | 0.00020–0.375 | 0.005 |
| Didelphimorphia – Rodentia | 12.66 | 0.940–170.40 | 0.060 |
| Pilosa – Primates | 0.0097 | 0.00016–0.607 | 0.019 |
| Pilosa – Rodentia | 14.37 | 0.656–314.63 | 0.128 |
| Primates – Rodentia | 1476.24 | 31.96–68189.05 | 0.000002 |

-> Results : Mammalian `order` significantly improved model fit in both the complete dataset (44 species; LRT: χ²₅ = 11.43, p = 0.044, ΔAIC = −1.43) and the conservative dataset (16 species; LRT: χ²₄ = 21.09, p < 0.001, ΔAIC = −13.09). In the complete dataset, only Primates and Rodentia differed significantly after Tukey correction, with substantially higher odds of `hemoplasma` infection in Primates (OR = 96.96, 95% CI: 3.17–2967.55, p = 0.0019). In the conservative dataset, Primates also showed higher odds than Didelphimorphia (OR = 116.62, 95% CI: 2.70–5036.73, p = 0.005), Pilosa (OR = 102.75, 95% CI: 1.64–6451.12, p = 0.019), and Rodentia (OR = 1476.24, 95% CI: 31.96–68189.05, p < 0.0001). All other pairwise comparisons were non-significant after Tukey correction.

-> Interpretation : `hemoplasma` infection probability varied among mammalian orders in both datasets, with the strongest and most consistent contrast involving the substantially higher prevalence in Primates, particularly compared with Rodentia.

### Visualization of odds ratios for mammalian order for complete (44 `species`) and conservative dataset (16 `species`)
```
order_OR_results <- summary(
  order_OR,
  infer = TRUE
) %>%
  as.data.frame() %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(asymp.LCL),
    CI_high = exp(asymp.UCL),
    dataset = "Complete dataset"
  ) %>%
  select(
    contrast,
    OR,
    CI_low,
    CI_high,
    p.value,
    dataset
  )

order_OR_results_n15 <- summary(
  order_OR_n15,
  infer = TRUE
) %>%
  as.data.frame() %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(asymp.LCL),
    CI_high = exp(asymp.UCL),
    dataset = "Conservative dataset"
  ) %>%
  select(
    contrast,
    OR,
    CI_low,
    CI_high,
    p.value,
    dataset
  )

contrast_levels <- c(
  "Carnivora - Cingulata",
  "Carnivora - Didelphimorphia",
  "Carnivora - Pilosa",
  "Carnivora - Primates",
  "Carnivora - Rodentia",
  "Cingulata - Didelphimorphia",
  "Cingulata - Pilosa",
  "Cingulata - Primates",
  "Cingulata - Rodentia",
  "Didelphimorphia - Pilosa",
  "Didelphimorphia - Primates",
  "Didelphimorphia - Rodentia",
  "Pilosa - Primates",
  "Pilosa - Rodentia",
  "Primates - Rodentia"
)

plot_OR <- bind_rows(
  order_OR_results,
  order_OR_results_n15
) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = rev(contrast_levels)
    ),
    y_base = as.numeric(contrast),
    y = ifelse(
      dataset == "Complete dataset",
      y_base + 0.22,
      y_base - 0.22
    )
  )

p_order_OR <- ggplot(
  plot_OR,
  aes(
    x = OR,
    y = y
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_errorbar(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    orientation = "y",
    width = 0,
    linewidth = 0.9,
    colour = "black"
  ) +
  geom_point(
    data = subset(
      plot_OR,
      dataset == "Complete dataset"
    ),
    shape = 21,
    fill = "black",
    colour = "black",
    size = 5,
    stroke = 1.1
  ) +
  geom_point(
    data = subset(
      plot_OR,
      dataset == "Conservative dataset"
    ),
    shape = 21,
    fill = "white",
    colour = "black",
    size = 5,
    stroke = 1.1
  ) +
  annotate(
    "point",
    x = 500,
    y = 15.55,
    shape = 21,
    fill = "black",
    colour = "black",
    size = 4,
    stroke = 1
  ) +
  annotate(
    "text",
    x = 700,
    y = 15.55,
    label = "Complete dataset",
    hjust = 0,
    vjust = 0.5,
    size = 3.2,
    family = "Calibri"
  ) +
  annotate(
    "point",
    x = 500,
    y = 15.15,
    shape = 21,
    fill = "white",
    colour = "black",
    size = 4,
    stroke = 1
  ) +
  annotate(
    "text",
    x = 700,
    y = 15.15,
    label = "Conservative dataset",
    hjust = 0,
    vjust = 0.5,
    size = 3.2,
    family = "Calibri"
  ) +
  scale_y_continuous(
    breaks = 1:15,
    labels = rev(contrast_levels)
  ) +
  scale_x_log10(
    breaks = c(
      0.001,
      0.01,
      0.1,
      1,
      10,
      100,
      1000,
      10000,
      100000
    ),
    labels = c(
      "0.001",
      "0.01",
      "0.1",
      "1",
      "10",
      "100",
      "1,000",
      "10,000",
      "100,000"
    ),
    expand = expansion(
      mult = c(0.02, 0.05)
    )
  ) +
  labs(
    x = "Odds ratio (log scale)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "Calibri"
    ),
    axis.text.y = element_text(
      size = 10,
      family = "Calibri"
    ),
    axis.text.x = element_text(
      size = 10,
      family = "Calibri"
    ),
    axis.title.x = element_text(
      size = 11,
      family = "Calibri"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

print(p_order_OR)

ggsave(
  filename = "order_haemoplasma_OR_complete_vs_n15.png",
  plot = p_order_OR,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300
)
```
## Step 6. Phylogeny of the 44 mammalian species
### List of mammalian species
```
mammal_species <- c(
  "Alouatta macconnelli",
  "Saguinus midas",
  "Sapajus apella",
  "Saimiri sciureus",
  "Pithecia pithecia",
  "Bradypus tridactylus",
  "Choloepus didactylus",
  "Cyclopes didactylus",
  "Tamandua tetradactyla",
  "Cabassous unicinctus",
  "Dasypus novemcinctus",
  "Hydrochoerus hydrochaeris",
  "Holochilus sciureus",
  "Hylaeamys megacephalus",
  "Hylaeamys yunganus",
  "Neacomys dubosti",
  "Neacomys paracou",
  "Nectomys rattus",
  "Oecomys auyantepui",
  "Oecomys bicolor",
  "Oligoryzomys fulvescens",
  "Makalata didelphoides",
  "Mesomys hispidus",
  "Proechimys cuvieri",
  "Proechimys guyannensis",
  "Coendou melanurus",
  "Coendou prehensilis",
  "Mus musculus",
  "Rattus rattus",
  "Sciurus aestuans",
  "Leopardus wiedii",
  "Puma yagouaroundi",
  "Eira barbara",
  "Galictis vittata",
  "Lontra longicaudis",
  "Potos flavus",
  "Caluromys philander",
  "Didelphis marsupialis",
  "Marmosa lepida",
  "Marmosa murina",
  "Marmosops parvidens",
  "Metachirus nudicaudatus",
  "Micoureus demerarae",
  "Philander opossum"
)
```
### Extract the subtree from the Open Tree Taxonomy and add branch lengths using Grafen's method
```
taxon_matches <- tnrs_match_names(
  names = mammal_species,
  context_name = "Mammals",
  do_approximate_matching = TRUE
)
taxon_matches[, c(
  "search_string",
  "unique_name",
  "score",
  "ott_id",
  "is_synonym",
  "number_matches"
)]
ott_ids <- ott_id(taxon_matches)
in_tree <- is_in_tree(ott_ids)
taxa_not_in_tree <- taxon_matches[
  !in_tree,
  c(
    "search_string",
    "unique_name",
    "score",
    "ott_id"
  )
]
taxa_not_in_tree
mammal_tree <- tol_induced_subtree(
  ott_ids = ott_ids[in_tree],
  label_format = "name"
)
mammal_tree$tip.label <- gsub(
  " ",
  "_",
  mammal_tree$tip.label
)
mammal_tree_grafen <- compute.brlen(
  mammal_tree,
  method = "Grafen"
)
mammal_tree_grafen$tip.label[
  mammal_tree_grafen$tip.label ==
    "Alouatta_seniculus_macconnelli"
] <- "Alouatta_macconnelli"
expected_species <- gsub(
  " ",
  "_",
  mammal_species
)
missing_species <- setdiff(
  expected_species,
  mammal_tree_grafen$tip.label
)
missing_species
Ntip(mammal_tree_grafen)
all(expected_species %in% mammal_tree_grafen$tip.label)
setdiff(
  mammal_tree_grafen$tip.label,
  expected_species
)
plot(
  mammal_tree_grafen,
  cex = 0.7,
  no.margin = TRUE
)
write.tree(
  mammal_tree_grafen,
  file = "mammal_phylogeny_Grafen.tre"
)
write.csv(
  taxon_matches,
  file = "OpenTree_taxon_matching.csv",
  row.names = FALSE
)
```

### Pairwise phylogenetic distances and heatmap
```
phylo_dist <- cophenetic.phylo(mammal_tree_grafen)
dim(phylo_dist)
write.csv(
  phylo_dist,
  file = "pairwise_phylogenetic_distances_Grafen.csv",
  row.names = TRUE
)
phylo_dist_long <- as.data.frame(phylo_dist) %>%
  mutate(
    species_1 = rownames(.)
  ) %>%
  pivot_longer(
    cols = -species_1,
    names_to = "species_2",
    values_to = "phylogenetic_distance"
  )
tree_order <- mammal_tree_grafen$tip.label
phylo_dist_long <- phylo_dist_long %>%
  mutate(
    species_1 = factor(
      species_1,
      levels = tree_order
    ),
    species_2 = factor(
      species_2,
      levels = rev(tree_order)
    )
  )
ggplot(
  phylo_dist_long,
  aes(
    x = species_1,
    y = species_2,
    fill = phylogenetic_distance
  )
) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "darkgreen",
    name = "Phylogenetic\ndistance"
  ) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7
    ),
    axis.text.y = element_text(
      size = 7
    ),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )
ggsave(
  "pairwise_phylogenetic_distances_heatmap.png",
  width = 10,
  height = 10,
  dpi = 300
)
```

### Test phylogenetic signal of `hemoplasma` prevalence (Pagel's lambda) in complete (44 `species`) and conservative dataset (16 `species`)
```
species_prev_complete <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n_sampled,
    .groups = "drop"
  ) %>%
  mutate(
    species = as.character(species)
  )

cat(
  "Number of species in complete dataset:",
  nrow(species_prev_complete),
  "\n"
)

cat(
  "Number of species in phylogeny:",
  length(mammal_tree_grafen$tip.label),
  "\n"
)

cat("\nSpecies in complete dataset but not in phylogeny:\n")

print(
  setdiff(
    species_prev_complete$species,
    mammal_tree_grafen$tip.label
  )
)

cat("\nSpecies in phylogeny but not in complete dataset:\n")

print(
  setdiff(
    mammal_tree_grafen$tip.label,
    species_prev_complete$species
  )
)

stopifnot(
  nrow(species_prev_complete) == 44
)

stopifnot(
  length(mammal_tree_grafen$tip.label) == 44
)

stopifnot(
  setequal(
    species_prev_complete$species,
    mammal_tree_grafen$tip.label
  )
)

prevalence_complete <- species_prev_complete$prevalence

names(prevalence_complete) <- species_prev_complete$species

prevalence_complete <- prevalence_complete[
  mammal_tree_grafen$tip.label
]

stopifnot(
  identical(
    names(prevalence_complete),
    mammal_tree_grafen$tip.label
  )
)

pagel_lambda_complete <- phylosig(
  tree = mammal_tree_grafen,
  x = prevalence_complete,
  method = "lambda",
  test = TRUE,
  nsim = 1000
)

pagel_lambda_complete


species_prev_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n_sampled,
    .groups = "drop"
  ) %>%
  filter(
    n_sampled >= 15
  ) %>%
  mutate(
    species = as.character(species)
  )

cat(
  "\nNumber of species in conservative dataset:",
  nrow(species_prev_n15),
  "\n"
)

cat(
  "Number of species in phylogeny:",
  length(mammal_tree_grafen$tip.label),
  "\n"
)

cat("\nSpecies in conservative dataset but not in phylogeny:\n")

print(
  setdiff(
    species_prev_n15$species,
    mammal_tree_grafen$tip.label
  )
)

cat("\nSpecies in phylogeny but not in conservative dataset:\n")

print(
  setdiff(
    mammal_tree_grafen$tip.label,
    species_prev_n15$species
  )
)

stopifnot(
  nrow(species_prev_n15) == 16
)

tree_n15 <- drop.tip(
  mammal_tree_grafen,
  setdiff(
    mammal_tree_grafen$tip.label,
    species_prev_n15$species
  )
)

stopifnot(
  length(tree_n15$tip.label) == 16
)

stopifnot(
  setequal(
    species_prev_n15$species,
    tree_n15$tip.label
  )
)

prevalence_n15 <- species_prev_n15$prevalence

names(prevalence_n15) <- species_prev_n15$species

prevalence_n15 <- prevalence_n15[
  tree_n15$tip.label
]

stopifnot(
  identical(
    names(prevalence_n15),
    tree_n15$tip.label
  )
)

pagel_lambda_n15 <- phylosig(
  tree = tree_n15,
  x = prevalence_n15,
  method = "lambda",
  test = TRUE,
  nsim = 1000
)

pagel_lambda_n15
```
-> Results : `hemoplasma` prevalence showed no detectable phylogenetic signal across the complete dataset of 44 mammalian `species` (Pagel’s λ = 0.00008, p = 1.00). The conservative dataset showed a strong phylogenetic signal (λ = 0.775), but with a marginal statistical support (p = 0.066).

-> Interpretation : While no phylogenetic structure was detected across all 44 species, the conservative dataset revealed a stronger, but marginal, phylogenetic signal for greater similarity in `hemoplasma` prevalence among closely related `species`.

### Visualization of the phylogenetic distribution of `hemoplasma` prevalence across 44 mammalian `species`
```
species_prevalence <- data_hemoplasma_stat %>%
  group_by(species, order) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n_sampled,
    .groups = "drop"
  ) %>%
  mutate(
    species_label = gsub("_", " ", species),
    group = case_when(
      order == "Didelphimorphia" ~ "Didelphimorphia",
      order == "Pilosa" ~ "Pilosa",
      order == "Primates" ~ "Primates",
      order == "Rodentia" ~ "Rodentia",
      order == "Cingulata" ~ "Cingulata",
      order == "Carnivora" ~ "Carnivora"
    ),
    n_class = if_else(n_sampled >= 15, "n ≥ 15", "n < 15")
  )
group_colors <- c(
  "Didelphimorphia" = "#4472C4",
  "Pilosa" = "#C65911",
  "Primates" = "#264478",
  "Rodentia" = "#D6A500",
  "Cingulata" = "#666666",
  "Carnivora" = "#375623"
)
tree_44 <- drop.tip(
  mammal_tree_grafen,
  setdiff(
    mammal_tree_grafen$tip.label,
    species_prevalence$species
  )
)
Ntip(tree_44)
p_tree <- ggtree(
  tree_44,
  layout = "rectangular"
)
p_tree <- p_tree %<+% species_prevalence
x_max <- max(
  p_tree$data$x[p_tree$data$isTip],
  na.rm = TRUE
)
p <- p_tree +
  geom_point(
    data = p_tree$data %>% filter(isTip, n_sampled < 15),
    aes(
      x = x_max + 0.25,
      y = y,
      size = prevalence,
      colour = group,
      fill = group
    ),
    shape = 21,
    stroke = 1
  ) +
  geom_point(
    data = p_tree$data %>% filter(isTip, n_sampled >= 15),
    aes(
      x = x_max + 0.25,
      y = y,
      size = prevalence,
      colour = group
    ),
    shape = 21,
    fill = "white",
    stroke = 1
  ) +
  geom_tiplab(
    aes(
      label = species_label,
      colour = group
    ),
    size = 4,
    hjust = 0,
    offset = 0.45,
    fontface = "italic",
    show.legend = FALSE
  ) +
  scale_colour_manual(
    name = "Order",
    values = group_colors,
    guide = guide_legend(
      override.aes = list(
        shape = 16,
        size = 3,
        fill = group_colors
      )
    )
  ) +
  scale_fill_manual(
    values = group_colors,
    guide = "none"
  ) +
  scale_size_continuous(
    name = "Hemoplasma\nprevalence",
    range = c(2.5, 6.5),
    limits = c(0, 1),
    labels = percent_format(accuracy = 1)
  ) +
  xlim(
    0,
    x_max + 5
  ) +
  theme_tree2() +
  theme(
    text = element_text(family = "sans"),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
p
ggsave(
  filename = "mammal_tree_Hemoplasma_prevalence.png",
  plot = p,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)
```

### Exhaustive phylogenetic clade screening `hemoplasma` prevalence (complete data set, 44 `species')
```
species_data <- data_hemoplasma_stat %>%
  mutate(
    hemoplasma = as.numeric(as.character(hemoplasma)),
    species = as.character(species)
  ) %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n_sampled,
    .groups = "drop"
  )
stopifnot(
  all(species_data$species %in% mammal_tree_grafen$tip.label)
)
cat(
  "Number of species:",
  nrow(species_data),
  "\n"
)
test_clade <- function(node, tree, species_data) {
  clade_species <- extract.clade(
    tree,
    node = node
  )$tip.label
  test_data <- species_data %>%
    mutate(
      clade = if_else(
        species %in% clade_species,
        "Clade",
        "Outside"
      ),
      clade = factor(
        clade,
        levels = c("Outside", "Clade")
      )
    )
  n_clade <- sum(test_data$clade == "Clade")
  n_outside <- sum(test_data$clade == "Outside")
  if (
    n_clade < 2 ||
    n_outside < 2
  ) {
    return(NULL)
  }
  model_null <- tryCatch(
    glm(
      cbind(
        n_positive,
        n_sampled - n_positive
      ) ~ 1,
      data = test_data,
      family = binomial
    ),
    error = function(e) NULL
  )
  model_full <- tryCatch(
    glm(
      cbind(
        n_positive,
        n_sampled - n_positive
      ) ~ clade,
      data = test_data,
      family = binomial
    ),
    error = function(e) NULL
  )
  if (
    is.null(model_null) ||
    is.null(model_full)
  ) {
    return(NULL)
  }
  LRT <- tryCatch(
    anova(
      model_null,
      model_full,
      test = "Chisq"
    ),
    error = function(e) NULL
  )
  if (is.null(LRT)) {
    return(NULL)
  }
  chi2 <- LRT$Deviance[2]
  p_value <- LRT$`Pr(>Chi)`[2]
  coef_table <- summary(model_full)$coefficients
  clade_row <- grep(
    "^clade",
    rownames(coef_table)
  )
  if (length(clade_row) != 1) {
    return(NULL)
  }
  estimate <- coef_table[
    clade_row,
    "Estimate"
  ]
  OR <- exp(estimate)
  prevalence_clade <- test_data %>%
    filter(clade == "Clade") %>%
    summarise(
      prev = sum(n_positive) /
        sum(n_sampled)
    ) %>%
    pull(prev)
  prevalence_outside <- test_data %>%
    filter(clade == "Outside") %>%
    summarise(
      prev = sum(n_positive) /
        sum(n_sampled)
    ) %>%
    pull(prev)
  tibble(
    node = node,
    n_species_clade = n_clade,
    n_species_outside = n_outside,
    n_individuals_clade = sum(
      test_data$n_sampled[
        test_data$clade == "Clade"
      ]
    ),
    n_individuals_outside = sum(
      test_data$n_sampled[
        test_data$clade == "Outside"
      ]
    ),
    prevalence_clade = prevalence_clade,
    prevalence_outside = prevalence_outside,
    OR = OR,
    chi2 = chi2,
    p_value = p_value,
    species_clade = paste(
      clade_species,
      collapse = "; "
    )
  )
}
internal_nodes <- (
  Ntip(mammal_tree_grafen) + 1
):(
  Ntip(mammal_tree_grafen) +
    Nnode(mammal_tree_grafen)
)
clade_results_list <- lapply(
  internal_nodes,
  test_clade,
  tree = mammal_tree_grafen,
  species_data = species_data
)
clade_results <- bind_rows(
  clade_results_list
)
clade_results <- clade_results %>%
  mutate(
    p_adjusted = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(
    p_adjusted,
    p_value
  )
clade_results %>%
  select(
    node,
    n_species_clade,
    n_species_outside,
    n_individuals_clade,
    n_individuals_outside,
    prevalence_clade,
    prevalence_outside,
    OR,
    chi2,
    p_value,
    p_adjusted,
    species_clade
  ) %>%
  print(n = 30)
significant_clades <- clade_results %>%
  filter(
    p_adjusted < 0.05
  )
significant_clades
suggestive_clades <- clade_results %>%
  filter(
    p_adjusted < 0.10
  )
suggestive_clades
clade_results %>%
  arrange(p_value) %>%
  select(
    node,
    n_species_clade,
    n_species_outside,
    n_individuals_clade,
    n_individuals_outside,
    prevalence_clade,
    prevalence_outside,
    OR,
    chi2,
    p_value,
    p_adjusted,
    species_clade
  ) %>%
  print(n = 20)
```
### Exhaustive phylogenetic clade screening `hemoplasma` prevalence (conservative dataset, 16 `species')
```
species_data <- data_hemoplasma_stat %>%
  mutate(
    hemoplasma = as.numeric(as.character(hemoplasma)),
    species = as.character(species)
  ) %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n_sampled,
    .groups = "drop"
  ) %>%
  filter(n_sampled >= 15)
stopifnot(
  all(species_data$species %in% mammal_tree_grafen$tip.label)
)
cat(
  "Number of species (n >= 15):",
  nrow(species_data),
  "\n"
)
test_clade <- function(node, tree, species_data) {
  clade_species <- extract.clade(
    tree,
    node = node
  )$tip.label
  test_data <- species_data %>%
    mutate(
      clade = if_else(
        species %in% clade_species,
        "Clade",
        "Outside"
      ),
      clade = factor(
        clade,
        levels = c("Outside", "Clade")
      )
    )
  n_clade <- sum(test_data$clade == "Clade")
  n_outside <- sum(test_data$clade == "Outside")
  if (
    n_clade < 2 ||
    n_outside < 2
  ) {
    return(NULL)
  }
  model_null <- tryCatch(
    glm(
      cbind(
        n_positive,
        n_sampled - n_positive
      ) ~ 1,
      data = test_data,
      family = binomial
    ),
    error = function(e) NULL
  )
  model_full <- tryCatch(
    glm(
      cbind(
        n_positive,
        n_sampled - n_positive
      ) ~ clade,
      data = test_data,
      family = binomial
    ),
    error = function(e) NULL
  )
  if (
    is.null(model_null) ||
    is.null(model_full)
  ) {
    return(NULL)
  }
  LRT <- tryCatch(
    anova(
      model_null,
      model_full,
      test = "Chisq"
    ),
    error = function(e) NULL
  )
  if (is.null(LRT)) {
    return(NULL)
  }
  chi2 <- LRT$Deviance[2]
  p_value <- LRT$`Pr(>Chi)`[2]
  coef_table <- summary(model_full)$coefficients
  clade_row <- grep(
    "^clade",
    rownames(coef_table)
  )
  if (length(clade_row) != 1) {
    return(NULL)
  }
  estimate <- coef_table[
    clade_row,
    "Estimate"
  ]
  OR <- exp(estimate)
  prevalence_clade <- test_data %>%
    filter(clade == "Clade") %>%
    summarise(
      prev = sum(n_positive) /
        sum(n_sampled)
    ) %>%
    pull(prev)
  prevalence_outside <- test_data %>%
    filter(clade == "Outside") %>%
    summarise(
      prev = sum(n_positive) /
        sum(n_sampled)
    ) %>%
    pull(prev)
  tibble(
    node = node,
    n_species_clade = n_clade,
    n_species_outside = n_outside,
    n_individuals_clade = sum(
      test_data$n_sampled[
        test_data$clade == "Clade"
      ]
    ),
    n_individuals_outside = sum(
      test_data$n_sampled[
        test_data$clade == "Outside"
      ]
    ),
    prevalence_clade = prevalence_clade,
    prevalence_outside = prevalence_outside,
    OR = OR,
    chi2 = chi2,
    p_value = p_value,
    species_clade = paste(
      clade_species[clade_species %in% species_data$species],
      collapse = "; "
    )
  )
}
internal_nodes <- (
  Ntip(mammal_tree_grafen) + 1
):(
  Ntip(mammal_tree_grafen) +
    Nnode(mammal_tree_grafen)
)
clade_results_list <- lapply(
  internal_nodes,
  test_clade,
  tree = mammal_tree_grafen,
  species_data = species_data
)
clade_results <- bind_rows(
  clade_results_list
)
clade_results <- clade_results %>%
  mutate(
    p_adjusted = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(
    p_adjusted,
    p_value
  )
clade_results %>%
  select(
    node,
    n_species_clade,
    n_species_outside,
    n_individuals_clade,
    n_individuals_outside,
    prevalence_clade,
    prevalence_outside,
    OR,
    chi2,
    p_value,
    p_adjusted,
    species_clade
  ) %>%
  print(n = Inf)
significant_clades <- clade_results %>%
  filter(
    p_adjusted < 0.05
  )
significant_clades
suggestive_clades <- clade_results %>%
  filter(
    p_adjusted < 0.10
  )
suggestive_clades
clade_results %>%
  arrange(p_value) %>%
  select(
    node,
    n_species_clade,
    n_species_outside,
    prevalence_clade,
    prevalence_outside,
    OR,
    chi2,
    p_value,
    p_adjusted,
    species_clade
  ) %>%
  print(n = 20)
```

-> Results : For the complete dataset, phylogenetically defined clades showed strong heterogeneity in Hemoplasma prevalence. After Benjamini–Hochberg correction, 22 of 42 tested clades were significantly associated with prevalence differences (padj < 0.05). The strongest association involved a primate clade of four species, with 93.8% prevalence versus 24.2% outside the clade (OR = 47.7, χ² = 128, p < 0.001, padj = 4.1 × 10⁻²⁸). Conversely, several rodent clades showed markedly lower prevalence than species outside the clade (e.g. 5.5% vs. 43.8%, OR = 0.075, χ² = 111, p < 0.001, padj < 10⁻²⁴). For the conservative dataset, the same overall pattern was observed. Thirteen of 24 tested clades remained significant after correction. The strongest association involved *Saguinus midas* + *Alouatta macconnelli*, with 96.8% prevalence versus 27.3% outside the clade (OR = 81.3, χ² = 124, p < 0.001, padj = 8.3 × 10⁻²⁸). Several rodent clades also showed substantially lower prevalence, including a seven-species clade with 4.4% prevalence versus 47.0% outside (OR = 0.051, χ² = 100, p < 0.001, padj = 8.6 × 10⁻²³).

-> Interpretation : These analyses provide strong evidence for phylogenetically structured variation in Hemoplasma prevalence, and the pattern is not driven solely by species with small sample sizes, as most major associations persisted in the conservative dataset. The recurrent contrast between high-prevalence primate lineages and low-prevalence rodent lineages suggests that host evolutionary history may contribute substantially to variation in infection probability. However, because many clade tests are overlapping and therefore non-independent, these results should be interpreted as evidence of phylogenetic clustering rather than as independent effects of each individual clade.

## Step 7. `hemoplasma` infection prevalence and mammal trait-based analyses
### Data retrieval and convert categorical variables 
The life trait dataset is available in the GitHub repository [here](https://github.com/olivierduron/Hemoplasma_infections/blob/main/data_mammal_traits.csv).
```
data_mammal_traits <- read.csv2("https://raw.githubusercontent.com/olivierduron/Hemoplasma_infections/main/data_mammal_traits.csv")
data_mammal_traits
data_mammal_traits$species        <- as.factor(data_mammal_traits$species)
data_mammal_traits$strata        <- as.factor(data_mammal_traits$strata)
data_mammal_traits$activitynocturnal        <- as.factor(data_mammal_traits$activitynocturnal)
data_mammal_traits$activitycrepuscular        <- as.factor(data_mammal_traits$activitydiurnal)
data_mammal_traits$activitydiurnal        <- as.factor(data_mammal_traits$activitydiurnal)
```

### Variation in `hemoplasma` infection prevalence according to the host diet (complete dataset, 44 `species`)
```
species_data <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, dietinv, dietvet, dietplant),
    by = "species"
  )
data_dietinv <- species_data %>%
  filter(!is.na(dietinv))
data_dietvet <- species_data %>%
  filter(!is.na(dietvet))
data_dietplant <- species_data %>%
  filter(!is.na(dietplant))
model_dietinv_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_dietinv,
  family = betabinomial(link = "logit")
)
model_dietinv <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ dietinv,
  data = data_dietinv,
  family = betabinomial(link = "logit")
)
model_dietvet_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_dietvet,
  family = betabinomial(link = "logit")
)
model_dietvet <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ dietvet,
  data = data_dietvet,
  family = betabinomial(link = "logit")
)
model_dietplant_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_dietplant,
  family = betabinomial(link = "logit")
)
model_dietplant <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ dietplant,
  data = data_dietplant,
  family = betabinomial(link = "logit")
)
lrt_dietinv <- anova(
  model_dietinv_null,
  model_dietinv
)
lrt_dietvet <- anova(
  model_dietvet_null,
  model_dietvet
)
lrt_dietplant <- anova(
  model_dietplant_null,
  model_dietplant
)
extract_beta_results <- function(model, null_model, variable, lrt) {
  coef_table <- summary(model)$coefficients$cond
  estimate <- coef_table[variable, "Estimate"]
  SE <- coef_table[variable, "Std. Error"]
  data.frame(
    variable = variable,
    n_species = nobs(model),
    estimate = estimate,
    SE = SE,
    z = coef_table[variable, "z value"],
    p_coefficient = coef_table[variable, "Pr(>|z|)"],
    LRT_chisq = lrt$Chisq[2],
    LRT_p = lrt$`Pr(>Chisq)`[2],
    AIC_null = AIC(null_model),
    AIC_model = AIC(model),
    delta_AIC = AIC(null_model) - AIC(model),
    CI_low = estimate - 1.96 * SE,
    CI_high = estimate + 1.96 * SE,
    OR = exp(estimate),
    OR_low = exp(estimate - 1.96 * SE),
    OR_high = exp(estimate + 1.96 * SE)
  )
}
results_diet <- bind_rows(
  extract_beta_results(
    model_dietinv,
    model_dietinv_null,
    "dietinv",
    lrt_dietinv
  ),
  extract_beta_results(
    model_dietvet,
    model_dietvet_null,
    "dietvet",
    lrt_dietvet
  ),
  extract_beta_results(
    model_dietplant,
    model_dietplant_null,
    "dietplant",
    lrt_dietplant
  )
)
results_diet <- results_diet %>%
  mutate(
    LRT_p_FDR = p.adjust(LRT_p, method = "BH"),
    p_coefficient_FDR = p.adjust(p_coefficient, method = "BH")
  )
print(results_diet)
```
### Variation in `hemoplasma` infection prevalence according to the host diet (conservative dataset, 16 `species`)
```
species_data_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  filter(
    n >= 15
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, dietinv, dietvet, dietplant),
    by = "species"
  )

data_dietinv_n15 <- species_data_n15 %>%
  filter(!is.na(dietinv))

data_dietvet_n15 <- species_data_n15 %>%
  filter(!is.na(dietvet))

data_dietplant_n15 <- species_data_n15 %>%
  filter(!is.na(dietplant))

model_dietinv_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_dietinv_n15,
  family = betabinomial(link = "logit")
)

model_dietinv_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ dietinv,
  data = data_dietinv_n15,
  family = betabinomial(link = "logit")
)

model_dietvet_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_dietvet_n15,
  family = betabinomial(link = "logit")
)

model_dietvet_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ dietvet,
  data = data_dietvet_n15,
  family = betabinomial(link = "logit")
)

model_dietplant_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_dietplant_n15,
  family = betabinomial(link = "logit")
)

model_dietplant_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ dietplant,
  data = data_dietplant_n15,
  family = betabinomial(link = "logit")
)

lrt_dietinv_n15 <- anova(
  model_dietinv_null_n15,
  model_dietinv_n15
)

lrt_dietvet_n15 <- anova(
  model_dietvet_null_n15,
  model_dietvet_n15
)

lrt_dietplant_n15 <- anova(
  model_dietplant_null_n15,
  model_dietplant_n15
)

extract_beta_results <- function(
  model,
  null_model,
  variable,
  lrt
) {
  coef_table <- summary(model)$coefficients$cond
  
  estimate <- coef_table[variable, "Estimate"]
  SE <- coef_table[variable, "Std. Error"]
  
  data.frame(
    variable = variable,
    n_species = nobs(model),
    estimate = estimate,
    SE = SE,
    z = coef_table[variable, "z value"],
    p_coefficient = coef_table[
      variable,
      "Pr(>|z|)"
    ],
    LRT_chisq = lrt$Chisq[2],
    LRT_p = lrt$`Pr(>Chisq)`[2],
    AIC_null = AIC(null_model),
    AIC_model = AIC(model),
    delta_AIC = AIC(null_model) - AIC(model),
    CI_low = estimate - 1.96 * SE,
    CI_high = estimate + 1.96 * SE,
    OR = exp(estimate),
    OR_low = exp(estimate - 1.96 * SE),
    OR_high = exp(estimate + 1.96 * SE)
  )
}

results_diet_n15 <- bind_rows(
  extract_beta_results(
    model_dietinv_n15,
    model_dietinv_null_n15,
    "dietinv",
    lrt_dietinv_n15
  ),
  extract_beta_results(
    model_dietvet_n15,
    model_dietvet_null_n15,
    "dietvet",
    lrt_dietvet_n15
  ),
  extract_beta_results(
    model_dietplant_n15,
    model_dietplant_null_n15,
    "dietplant",
    lrt_dietplant_n15
  )
)

results_diet_n15 <- results_diet_n15 %>%
  mutate(
    LRT_p_FDR = p.adjust(
      LRT_p,
      method = "BH"
    ),
    p_coefficient_FDR = p.adjust(
      p_coefficient,
      method = "BH"
    )
  )

print(results_diet_n15)
```
-> Results: None of the dietary composition variables was significantly associated with `hemoplasma` prevalence in either the complete dataset (FDR-adjusted p = 0.76, 0.52, and 0.52 for invertebrate-, vertebrate-, and plant-based diet, respectively) or the conservative dataset (FDR-adjusted p = 0.56 for all three variables). In all cases, adding dietary composition resulted in a higher AIC than the corresponding null model (ΔAIC < 0).

-> Interpretation: Interspecific variation in `hemoplasma` prevalence was not explained by the relative contribution of invertebrates, vertebrates, or plants in the host diet, with consistent results across both datasets.

### Visualization of association between `hemoplasma` prevalence and dietary composition variables for complete (44 species)
```
species_order <- data_hemoplasma_stat %>%
  select(species, order) %>%
  distinct()

get_predictions <- function(model, data, variable, label) {
  x <- seq(
    min(data[[variable]], na.rm = TRUE),
    max(data[[variable]], na.rm = TRUE),
    length.out = 100
  )
  newdata <- data.frame(x)
  names(newdata) <- variable
  pred <- predict(
    model,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )
  newdata$fit <- plogis(pred$fit)
  newdata$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  newdata$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  newdata$variable <- label
  newdata
}

pred_dietinv <- get_predictions(
  model_dietinv,
  data_dietinv,
  "dietinv",
  "Invertebrates"
)

pred_dietvet <- get_predictions(
  model_dietvet,
  data_dietvet,
  "dietvet",
  "Vertebrates"
)

pred_dietplant <- get_predictions(
  model_dietplant,
  data_dietplant,
  "dietplant",
  "Plants"
)

plot_data <- bind_rows(
  data_dietinv %>%
    mutate(
      diet = dietinv,
      variable = "Invertebrates"
    ),
  data_dietvet %>%
    mutate(
      diet = dietvet,
      variable = "Vertebrates"
    ),
  data_dietplant %>%
    mutate(
      diet = dietplant,
      variable = "Plants"
    )
) %>%
  left_join(
    species_order,
    by = "species"
  ) %>%
  mutate(
    species_initials = ifelse(
      n >= 15,
      paste0(
        substr(species, 1, 1),
        substr(
          sub("^[^_]+_", "", species),
          1,
          1
        )
      ),
      NA_character_
    )
  )

pred_data <- bind_rows(
  pred_dietinv %>%
    mutate(diet = dietinv),
  pred_dietvet %>%
    mutate(diet = dietvet),
  pred_dietplant %>%
    mutate(diet = dietplant)
)

label_data <- plot_data %>%
  filter(n >= 15) %>%
  group_by(variable, diet) %>%
  arrange(species, .by_group = TRUE) %>%
  mutate(
    duplicate_n = n(),
    position = row_number()
  ) %>%
  ungroup()

order_colors <- c(
  "Primates" = "#E69F00",
  "Rodentia" = "#56B4E9",
  "Pilosa" = "#009E73",
  "Didelphimorphia" = "#CC79A7",
  "Carnivora" = "#D55E00",
  "Cingulata" = "#0072B2"
)

figure_diet <- ggplot(
  plot_data,
  aes(
    x = diet,
    y = prevalence,
    colour = order
  )
) +
  geom_ribbon(
    data = pred_data,
    aes(
      x = diet,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.20
  ) +
  geom_line(
    data = pred_data,
    aes(
      x = diet,
      y = fit
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    data = subset(plot_data, n < 15),
    aes(size = n),
    alpha = 0.70
  ) +
  geom_point(
    data = subset(plot_data, n >= 15),
    aes(size = n),
    shape = 1,
    stroke = 1.1
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(
      x = diet,
      y = prevalence,
      label = species_initials
    ),
    size = 3.5,
    direction = "both",
    force = 2,
    force_pull = 0.5,
    box.padding = 0.8,
    point.padding = 0.5,
    min.segment.length = 0,
    segment.size = 0.3,
    segment.alpha = 0.5,
    max.overlaps = Inf,
    seed = 123,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ variable,
    nrow = 1,
    scales = "free_x",
    strip.position = "bottom"
  ) +
  scale_colour_manual(
    values = order_colors
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  scale_size_continuous(
    name = "Number tested"
  ) +
  labs(
    x = "Diet composition (%)",
    y = "Hemoplasma prevalence"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(
      size = 11,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 10
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.title.x = element_text(
      size = 11
    ),
    axis.title.y = element_text(
      size = 11
    ),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

print(figure_diet)

ggsave(
  "Hemoplasma_prevalence_diet_beta_binomial.png",
  figure_diet,
  width = 11,
  height = 4.8,
  dpi = 300
)
```

### Variation in `hemoplasma` infection prevalence according to the host foraging `strata` (complete dataset, 44 `species`) 
```
species_data_strata <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, strata),
    by = "species"
  ) %>%
  mutate(
    strata = factor(
      strata,
      levels = c("G", "S", "Ar")
    )
  )
cat("\nNumber of species =", nrow(species_data_strata), "\n")
cat("\nMissing values:\n")
print(colSums(is.na(species_data_strata[c("strata")])))
cat("\nStrata distribution:\n")
print(table(species_data_strata$strata))
model_strata_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = species_data_strata,
  family = betabinomial(link = "logit")
)
model_strata <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ strata,
  data = species_data_strata,
  family = betabinomial(link = "logit")
)
lrt_strata <- anova(
  model_strata_null,
  model_strata
)
cat("\n================ STRATA LRT ================\n")
print(lrt_strata)
cat("\n================ STRATA AIC ================\n")
print(AIC(model_strata_null, model_strata))
cat("\n================ STRATA MODEL ================\n")
print(summary(model_strata))
results_strata <- as.data.frame(
  summary(model_strata)$coefficients$cond
) %>%
  tibble::rownames_to_column("coefficient") %>%
  filter(coefficient != "(Intercept)") %>%
  mutate(
    variable = "strata",
    n_species = nobs(model_strata),
    estimate = Estimate,
    SE = `Std. Error`,
    z = `z value`,
    p_coefficient = `Pr(>|z|)`,
    CI_low = estimate - 1.96 * SE,
    CI_high = estimate + 1.96 * SE,
    OR = exp(estimate),
    OR_low = exp(CI_low),
    OR_high = exp(CI_high),
    AIC_null = AIC(model_strata_null),
    AIC_model = AIC(model_strata),
    delta_AIC = AIC_null - AIC_model,
    LRT_chisq = lrt_strata$Chisq[2],
    LRT_p = lrt_strata$`Pr(>Chisq)`[2]
  ) %>%
  select(
    variable,
    coefficient,
    n_species,
    estimate,
    SE,
    z,
    p_coefficient,
    LRT_chisq,
    LRT_p,
    AIC_null,
    AIC_model,
    delta_AIC,
    CI_low,
    CI_high,
    OR,
    OR_low,
    OR_high
  )
results_strata$LRT_p_FDR <- p.adjust(
  results_strata$LRT_p,
  method = "BH"
)
results_strata$p_coefficient_FDR <- p.adjust(
  results_strata$p_coefficient,
  method = "BH"
)
cat("\n================ COEFFICIENT RESULTS ================\n")
print(results_strata, row.names = FALSE)
emm_strata <- emmeans(
  model_strata,
  ~ strata,
  type = "response"
)
pairwise_strata <- pairs(
  emm_strata,
  adjust = "tukey"
)
cat("\n================ ESTIMATED PREVALENCE BY STRATA ================\n")
print(emm_strata)
cat("\n================ PAIRWISE COMPARISONS ================\n")
print(pairwise_strata)
pairwise_strata_results <- as.data.frame(
  summary(
    pairwise_strata,
    infer = TRUE
  )
)
print(pairwise_strata_results)
```
### Variation in `hemoplasma` infection prevalence according to the host’s foraging `strata` (conservative dataset, 16 `species`) 
```
species_data_strata_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  filter(
    n >= 15
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, strata),
    by = "species"
  ) %>%
  mutate(
    strata = factor(
      strata,
      levels = c("G", "S", "Ar")
    )
  )

cat(
  "\nNumber of species =",
  nrow(species_data_strata_n15),
  "\n"
)

cat("\nMissing values:\n")

print(
  colSums(
    is.na(
      species_data_strata_n15[c("strata")]
    )
  )
)

cat("\nStrata distribution:\n")

print(
  table(
    species_data_strata_n15$strata
  )
)

model_strata_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = species_data_strata_n15,
  family = betabinomial(link = "logit")
)

model_strata_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ strata,
  data = species_data_strata_n15,
  family = betabinomial(link = "logit")
)

lrt_strata_n15 <- anova(
  model_strata_null_n15,
  model_strata_n15
)

cat("\n================ STRATA LRT ================\n")

print(
  lrt_strata_n15
)

cat("\n================ STRATA AIC ================\n")

print(
  AIC(
    model_strata_null_n15,
    model_strata_n15
  )
)

cat("\n================ STRATA MODEL ================\n")

print(
  summary(
    model_strata_n15
  )
)

results_strata_n15 <- as.data.frame(
  summary(
    model_strata_n15
  )$coefficients$cond
) %>%
  tibble::rownames_to_column(
    "coefficient"
  ) %>%
  filter(
    coefficient != "(Intercept)"
  ) %>%
  mutate(
    variable = "strata",
    n_species = nobs(model_strata_n15),
    estimate = Estimate,
    SE = `Std. Error`,
    z = `z value`,
    p_coefficient = `Pr(>|z|)`,
    CI_low = estimate - 1.96 * SE,
    CI_high = estimate + 1.96 * SE,
    OR = exp(estimate),
    OR_low = exp(CI_low),
    OR_high = exp(CI_high),
    AIC_null = AIC(model_strata_null_n15),
    AIC_model = AIC(model_strata_n15),
    delta_AIC = AIC_null - AIC_model,
    LRT_chisq = lrt_strata_n15$Chisq[2],
    LRT_p = lrt_strata_n15$`Pr(>Chisq)`[2]
  ) %>%
  select(
    variable,
    coefficient,
    n_species,
    estimate,
    SE,
    z,
    p_coefficient,
    LRT_chisq,
    LRT_p,
    AIC_null,
    AIC_model,
    delta_AIC,
    CI_low,
    CI_high,
    OR,
    OR_low,
    OR_high
  )

results_strata_n15$LRT_p_FDR <- p.adjust(
  results_strata_n15$LRT_p,
  method = "BH"
)

results_strata_n15$p_coefficient_FDR <- p.adjust(
  results_strata_n15$p_coefficient,
  method = "BH"
)

cat("\n================ COEFFICIENT RESULTS ================\n")

print(
  results_strata_n15,
  row.names = FALSE
)

emm_strata_n15 <- emmeans(
  model_strata_n15,
  ~ strata,
  type = "response"
)

pairwise_strata_n15 <- pairs(
  emm_strata_n15,
  adjust = "tukey"
)

cat("\n================ ESTIMATED PREVALENCE BY STRATA ================\n")

print(
  emm_strata_n15
)

cat("\n================ PAIRWISE COMPARISONS ================\n")

print(
  pairwise_strata_n15
)

pairwise_strata_results_n15 <- as.data.frame(
  summary(
    pairwise_strata_n15,
    infer = TRUE
  )
)

print(
  pairwise_strata_results_n15
)
```
-> Results : `hemoplasma` prevalence did not differ significantly among foraging `strata` in the complete dataset (LRT: χ²₂ = 1.92, p = 0.382, ΔAIC = −2.08). In the conservative dataset, `strata` significantly improved model fit (LRT: χ²₂ = 7.72, p = 0.021; ΔAIC = 3.72). Model-estimated prevalence was 14.7% (95% CI: 5.7–32.9%) for ground-foraging `species`, 39.6% (95% CI: 14.7–71.4%) for scansorial `species`, and 60.8% (95% CI: 33.9–82.4%) for arboreal `species`. Tukey-adjusted pairwise comparisons identified a significant difference only between ground-foraging and arboreal species (p = 0.013).

-> Interpretation : Foraging `strata` was not associated with interspecific variation in hemoplasma prevalence across all 44 species, but prevalence increased from ground-foraging to arboreal species in the conservative dataset, with marginally to significantly stronger evidence for higher prevalence in arboreal species.

### Visualization of association between `hemoplasma` prevalence and foraging `strata` (conservative dataset, 44 `species`)
```
species_data_strata_plot <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, strata),
    by = "species"
  ) %>%
  left_join(
    species_order,
    by = "species"
  ) %>%
  mutate(
    strata = factor(
      strata,
      levels = c("G", "S", "Ar"),
      labels = c(
        "Ground",
        "Shrub",
        "Arboreal"
      )
    ),
    species_initials = paste0(
      substr(species, 1, 1),
      substr(
        sub("^[^_]+_", "", species),
        1,
        1
      )
    ),
    n15 = n >= 15
  ) %>%
  filter(
    !is.na(strata)
  )

label_data <- species_data_strata_plot %>%
  filter(n15) %>%
  group_by(strata) %>%
  arrange(prevalence) %>%
  mutate(
    direction = if_else(
      row_number() %% 2 == 0,
      1,
      -1
    ),
    label_x = as.numeric(strata) + direction * 0.18
  ) %>%
  ungroup()

figure_strata <- ggplot(
  species_data_strata_plot,
  aes(
    x = strata,
    y = prevalence,
    colour = order
  )
) +
  geom_boxplot(
    width = 0.50,
    outlier.shape = NA,
    colour = "black",
    fill = NA
  ) +
  geom_point(
    data = species_data_strata_plot %>%
      filter(!n15),
    aes(
      size = n
    ),
    alpha = 0.75
  ) +
  geom_point(
    data = species_data_strata_plot %>%
      filter(n15),
    aes(
      size = n
    ),
    shape = 1,
    stroke = 1.1
  ) +
  geom_segment(
    data = label_data,
    aes(
      x = as.numeric(strata),
      xend = label_x,
      y = prevalence,
      yend = prevalence,
      colour = order
    ),
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_text(
    data = label_data,
    aes(
      x = label_x,
      y = prevalence,
      label = species_initials,
      colour = order,
      hjust = ifelse(
        direction > 0,
        0,
        1
      )
    ),
    size = 3.5,
    family = "sans",
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = order_colors
  ) +
  scale_size_continuous(
    name = "Number tested",
    range = c(2.5, 8)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0.02, 0.08)
    )
  ) +
  labs(
    x = "Foraging stratum",
    y = "Hemoplasma prevalence",
    colour = "Mammalian order"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 10,
      family = "sans"
    ),
    axis.text.y = element_text(
      size = 10,
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 11,
      family = "sans"
    ),
    axis.title.y = element_text(
      size = 11,
      family = "sans"
    ),
    legend.title = element_text(
      size = 9,
      family = "sans"
    ),
    legend.text = element_text(
      size = 8,
      family = "sans"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

print(figure_strata)

ggsave(
  "Hemoplasma_prevalence_strata.png",
  figure_strata,
  width = 7,
  height = 4.8,
  dpi = 300
)
```

### Variation in `hemoplasma` infection prevalence according to the host’s activity (nocturnal, diurnal, crepuscular) (complete dataset, 44 `species`)
```
species_data_activity <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(
        species,
        activitynocturnal,
        activitycrepuscular,
        activitydiurnal
      ),
    by = "species"
  )
cat("\nNumber of species =", nrow(species_data_activity), "\n")
cat("\nMissing values:\n")
print(
  colSums(
    is.na(
      species_data_activity[
        ,
        c(
          "activitynocturnal",
          "activitycrepuscular",
          "activitydiurnal"
        )
      ]
    )
  )
)
cat("\nActivity distribution:\n")
print(table(species_data_activity$activitynocturnal, useNA = "ifany"))
print(table(species_data_activity$activitycrepuscular, useNA = "ifany"))
print(table(species_data_activity$activitydiurnal, useNA = "ifany"))
data_activitynocturnal <- species_data_activity %>%
  filter(!is.na(activitynocturnal))
data_activitycrepuscular <- species_data_activity %>%
  filter(!is.na(activitycrepuscular))
data_activitydiurnal <- species_data_activity %>%
  filter(!is.na(activitydiurnal))
model_activitynocturnal_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_activitynocturnal,
  family = betabinomial(link = "logit")
)
model_activitynocturnal <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ activitynocturnal,
  data = data_activitynocturnal,
  family = betabinomial(link = "logit")
)
model_activitycrepuscular_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_activitycrepuscular,
  family = betabinomial(link = "logit")
)
model_activitycrepuscular <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ activitycrepuscular,
  data = data_activitycrepuscular,
  family = betabinomial(link = "logit")
)
model_activitydiurnal_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_activitydiurnal,
  family = betabinomial(link = "logit")
)
model_activitydiurnal <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ activitydiurnal,
  data = data_activitydiurnal,
  family = betabinomial(link = "logit")
)
lrt_activitynocturnal <- anova(
  model_activitynocturnal_null,
  model_activitynocturnal
)
lrt_activitycrepuscular <- anova(
  model_activitycrepuscular_null,
  model_activitycrepuscular
)
lrt_activitydiurnal <- anova(
  model_activitydiurnal_null,
  model_activitydiurnal
)
cat("\n================ NOCTURNAL LRT ================\n")
print(lrt_activitynocturnal)
cat("\n================ CREPUSCULAR LRT ================\n")
print(lrt_activitycrepuscular)
cat("\n================ DIURNAL LRT ================\n")
print(lrt_activitydiurnal)
cat("\n================ AIC ================\n")
print(
  AIC(
    model_activitynocturnal_null,
    model_activitynocturnal
  )
)
print(
  AIC(
    model_activitycrepuscular_null,
    model_activitycrepuscular
  )
)
print(
  AIC(
    model_activitydiurnal_null,
    model_activitydiurnal
  )
)
extract_activity_results <- function(
  model,
  null_model,
  variable,
  lrt
) {
  coef_table <- summary(model)$coefficients$cond
  coefficient_name <- setdiff(
    rownames(coef_table),
    "(Intercept)"
  )[1]
  coef_row <- coef_table[
    coefficient_name,
    ,
    drop = FALSE
  ]
  estimate <- coef_row[1, "Estimate"]
  SE <- coef_row[1, "Std. Error"]
  z <- coef_row[1, "z value"]
  p <- coef_row[1, "Pr(>|z|)"]
  LRT_chisq <- lrt$Chisq[2]
  LRT_p <- lrt$`Pr(>Chisq)`[2]
  AIC_null <- AIC(null_model)
  AIC_model <- AIC(model)
  CI_low <- estimate - 1.96 * SE
  CI_high <- estimate + 1.96 * SE
  data.frame(
    variable = variable,
    coefficient = coefficient_name,
    n_species = nobs(model),
    estimate = estimate,
    SE = SE,
    z = z,
    p_coefficient = p,
    LRT_chisq = LRT_chisq,
    LRT_p = LRT_p,
    AIC_null = AIC_null,
    AIC_model = AIC_model,
    delta_AIC = AIC_null - AIC_model,
    CI_low = CI_low,
    CI_high = CI_high,
    OR = exp(estimate),
    OR_low = exp(CI_low),
    OR_high = exp(CI_high)
  )
}
results_activity <- bind_rows(
  extract_activity_results(
    model_activitynocturnal,
    model_activitynocturnal_null,
    "activitynocturnal",
    lrt_activitynocturnal
  ),
  extract_activity_results(
    model_activitycrepuscular,
    model_activitycrepuscular_null,
    "activitycrepuscular",
    lrt_activitycrepuscular
  ),
  extract_activity_results(
    model_activitydiurnal,
    model_activitydiurnal_null,
    "activitydiurnal",
    lrt_activitydiurnal
  )
)
results_activity$LRT_p_FDR <- p.adjust(
  results_activity$LRT_p,
  method = "BH"
)
results_activity$p_coefficient_FDR <- p.adjust(
  results_activity$p_coefficient,
  method = "BH"
)
cat("\n================ FINAL ACTIVITY RESULTS ================\n")
print(
  results_activity,
  row.names = FALSE
)
```
### Variation in `hemoplasma` infection prevalence according to the host’s activity (nocturnal, diurnal, crepuscular) (conservative dataset, 16 `species`)
```
species_data_activity_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  filter(
    n >= 15
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(
        species,
        activitynocturnal,
        activitycrepuscular,
        activitydiurnal
      ),
    by = "species"
  )

cat(
  "\nNumber of species =",
  nrow(species_data_activity_n15),
  "\n"
)

cat("\nMissing values:\n")

print(
  colSums(
    is.na(
      species_data_activity_n15[
        ,
        c(
          "activitynocturnal",
          "activitycrepuscular",
          "activitydiurnal"
        )
      ]
    )
  )
)

cat("\nActivity distribution:\n")

print(
  table(
    species_data_activity_n15$activitynocturnal,
    useNA = "ifany"
  )
)

print(
  table(
    species_data_activity_n15$activitycrepuscular,
    useNA = "ifany"
  )
)

print(
  table(
    species_data_activity_n15$activitydiurnal,
    useNA = "ifany"
  )
)

data_activitynocturnal_n15 <- species_data_activity_n15 %>%
  filter(
    !is.na(activitynocturnal)
  )

data_activitycrepuscular_n15 <- species_data_activity_n15 %>%
  filter(
    !is.na(activitycrepuscular)
  )

data_activitydiurnal_n15 <- species_data_activity_n15 %>%
  filter(
    !is.na(activitydiurnal)
  )

model_activitynocturnal_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_activitynocturnal_n15,
  family = betabinomial(link = "logit")
)

model_activitynocturnal_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ activitynocturnal,
  data = data_activitynocturnal_n15,
  family = betabinomial(link = "logit")
)

model_activitycrepuscular_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_activitycrepuscular_n15,
  family = betabinomial(link = "logit")
)

model_activitycrepuscular_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ activitycrepuscular,
  data = data_activitycrepuscular_n15,
  family = betabinomial(link = "logit")
)

model_activitydiurnal_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_activitydiurnal_n15,
  family = betabinomial(link = "logit")
)

model_activitydiurnal_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ activitydiurnal,
  data = data_activitydiurnal_n15,
  family = betabinomial(link = "logit")
)

lrt_activitynocturnal_n15 <- anova(
  model_activitynocturnal_null_n15,
  model_activitynocturnal_n15
)

lrt_activitycrepuscular_n15 <- anova(
  model_activitycrepuscular_null_n15,
  model_activitycrepuscular_n15
)

lrt_activitydiurnal_n15 <- anova(
  model_activitydiurnal_null_n15,
  model_activitydiurnal_n15
)

cat("\n================ NOCTURNAL LRT ================\n")

print(
  lrt_activitynocturnal_n15
)

cat("\n================ CREPUSCULAR LRT ================\n")

print(
  lrt_activitycrepuscular_n15
)

cat("\n================ DIURNAL LRT ================\n")

print(
  lrt_activitydiurnal_n15
)

cat("\n================ AIC ================\n")

print(
  AIC(
    model_activitynocturnal_null_n15,
    model_activitynocturnal_n15
  )
)

print(
  AIC(
    model_activitycrepuscular_null_n15,
    model_activitycrepuscular_n15
  )
)

print(
  AIC(
    model_activitydiurnal_null_n15,
    model_activitydiurnal_n15
  )
)

extract_activity_results <- function(
  model,
  null_model,
  variable,
  lrt
) {
  
  coef_table <- summary(
    model
  )$coefficients$cond
  
  coefficient_name <- setdiff(
    rownames(coef_table),
    "(Intercept)"
  )[1]
  
  coef_row <- coef_table[
    coefficient_name,
    ,
    drop = FALSE
  ]
  
  estimate <- coef_row[
    1,
    "Estimate"
  ]
  
  SE <- coef_row[
    1,
    "Std. Error"
  ]
  
  z <- coef_row[
    1,
    "z value"
  ]
  
  p <- coef_row[
    1,
    "Pr(>|z|)"
  ]
  
  LRT_chisq <- lrt$Chisq[2]
  
  LRT_p <- lrt$`Pr(>Chisq)`[2]
  
  AIC_null <- AIC(
    null_model
  )
  
  AIC_model <- AIC(
    model
  )
  
  CI_low <- estimate - 1.96 * SE
  
  CI_high <- estimate + 1.96 * SE
  
  data.frame(
    variable = variable,
    coefficient = coefficient_name,
    n_species = nobs(model),
    estimate = estimate,
    SE = SE,
    z = z,
    p_coefficient = p,
    LRT_chisq = LRT_chisq,
    LRT_p = LRT_p,
    AIC_null = AIC_null,
    AIC_model = AIC_model,
    delta_AIC = AIC_null - AIC_model,
    CI_low = CI_low,
    CI_high = CI_high,
    OR = exp(estimate),
    OR_low = exp(CI_low),
    OR_high = exp(CI_high)
  )
}

results_activity_n15 <- bind_rows(
  extract_activity_results(
    model_activitynocturnal_n15,
    model_activitynocturnal_null_n15,
    "activitynocturnal",
    lrt_activitynocturnal_n15
  ),
  extract_activity_results(
    model_activitycrepuscular_n15,
    model_activitycrepuscular_null_n15,
    "activitycrepuscular",
    lrt_activitycrepuscular_n15
  ),
  extract_activity_results(
    model_activitydiurnal_n15,
    model_activitydiurnal_null_n15,
    "activitydiurnal",
    lrt_activitydiurnal_n15
  )
)

results_activity_n15$LRT_p_FDR <- p.adjust(
  results_activity_n15$LRT_p,
  method = "BH"
)

results_activity_n15$p_coefficient_FDR <- p.adjust(
  results_activity_n15$p_coefficient,
  method = "BH"
)

cat(
  "\n================ FINAL ACTIVITY RESULTS ================\n"
)

print(
  results_activity_n15,
  row.names = FALSE
)
```
-> Results : In the complete dataset, none of the three activity categories was associated with `hemoplasma` prevalence: nocturnal activity (LRT: χ²₁ = 0.82, p = 0.364; ΔAIC = −1.18; OR = 0.45, 95% CI: 0.09–2.27), crepuscular activity (χ²₁ < 0.001, p = 0.994; ΔAIC = −2.00; OR = 1.00, 95% CI: 0.32–3.11), or diurnal activity (χ²₁ < 0.001, p = 0.994; ΔAIC = −2.00; OR = 1.00, 95% CI: 0.32–3.11). In the conservative dataset, nocturnal activity was strongly associated with lower `hemoplasma` prevalence (LRT: χ²₁ = 12.00, p < 0.001; ΔAIC = 10.00; OR = 0.022, 95% CI: 0.002–0.210), whereas neither crepuscular nor diurnal activity showed evidence of an association (χ²₁ = 1.59, p = 0.208; ΔAIC = −0.41; OR = 2.65, 95% CI: 0.60–11.70, for both).

-> Interpretation : Activity pattern was not associated with `hemoplasma` prevalence across the complete dataset, but the conservative dataset revealed a strong association with nocturnal activity, with nocturnal species showing substantially lower estimated prevalence. No comparable evidence was found for crepuscular or diurnal activity.

### Visualization of association between `hemoplasma` prevalence and nocturnal / crepuscular / diurnal activity (complete database, 44 `species`)
```
label_data <- activity_plot_data %>%
  filter(n15) %>%
  group_by(activity) %>%
  arrange(prevalence) %>%
  mutate(
    direction = if_else(row_number() %% 2 == 0, 1, -1),
    label_x = as.numeric(activity) + direction * 0.16
  ) %>%
  ungroup()

figure_activity <- ggplot(
  activity_plot_data,
  aes(
    x = activity,
    y = prevalence,
    colour = order
  )
) +
  geom_boxplot(
    width = 0.50,
    outlier.shape = NA,
    colour = "black",
    fill = NA
  ) +
  geom_point(
    data = activity_plot_data %>%
      filter(!n15),
    aes(size = n),
    alpha = 0.75
  ) +
  geom_point(
    data = activity_plot_data %>%
      filter(n15),
    aes(size = n),
    shape = 1,
    stroke = 1.1
  ) +
  geom_segment(
    data = label_data,
    aes(
      x = as.numeric(activity),
      xend = label_x,
      y = prevalence,
      yend = prevalence,
      colour = order
    ),
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_text(
    data = label_data,
    aes(
      x = label_x,
      y = prevalence,
      label = species_initials,
      colour = order,
      hjust = ifelse(direction > 0, 0, 1)
    ),
    size = 3.5,
    family = "sans",
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = order_colors
  ) +
  scale_size_continuous(
    name = "Number tested",
    range = c(2.5, 8)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0.02, 0.08)
    )
  ) +
  labs(
    x = "Activity pattern",
    y = "Hemoplasma prevalence",
    colour = "Mammalian order"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 10,
      family = "sans"
    ),
    axis.text.y = element_text(
      size = 10,
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 11,
      family = "sans"
    ),
    axis.title.y = element_text(
      size = 11,
      family = "sans"
    ),
    legend.title = element_text(
      size = 9,
      family = "sans"
    ),
    legend.text = element_text(
      size = 8,
      family = "sans"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

print(figure_activity)

ggsave(
  "Hemoplasma_prevalence_activity.png",
  figure_activity,
  width = 7,
  height = 4.8,
  dpi = 300
)
```

### Variation in `hemoplasma` infection prevalence according to the host `bodymass` (complete dataset, 44 `species`)
```
species_data_bodymass <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, bodymass) %>%
      mutate(bodymass = as.numeric(bodymass)),
    by = "species"
  ) %>%
  mutate(
    log_bodymass = log(bodymass),
    log_bodymass_scaled = as.numeric(scale(log_bodymass))
  )
cat("\nNumber of species =", nrow(species_data_bodymass), "\n")
cat("\nMissing bodymass values =", sum(is.na(species_data_bodymass$bodymass)), "\n")
data_bodymass <- species_data_bodymass %>%
  filter(
    !is.na(bodymass),
    !is.na(log_bodymass_scaled)
  )
cat("\nSpecies used for bodymass =", nrow(data_bodymass), "\n")
model_bodymass_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_bodymass,
  family = betabinomial(link = "logit")
)
model_bodymass <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ log_bodymass_scaled,
  data = data_bodymass,
  family = betabinomial(link = "logit")
)
lrt_bodymass <- anova(
  model_bodymass_null,
  model_bodymass
)
cat("\n================ BODY MASS LRT ================\n")
print(lrt_bodymass)
cat("\n================ BODY MASS AIC ================\n")
print(
  AIC(
    model_bodymass_null,
    model_bodymass
  )
)
cat("\n================ BODY MASS MODEL ================\n")
print(
  summary(model_bodymass)
)
coef_table <- summary(model_bodymass)$coefficients$cond
estimate <- coef_table["log_bodymass_scaled", "Estimate"]
SE <- coef_table["log_bodymass_scaled", "Std. Error"]
z <- coef_table["log_bodymass_scaled", "z value"]
p <- coef_table["log_bodymass_scaled", "Pr(>|z|)"]
CI_low <- estimate - 1.96 * SE
CI_high <- estimate + 1.96 * SE
results_bodymass <- data.frame(
  variable = "bodymass",
  coefficient = "log_bodymass_scaled",
  n_species = nobs(model_bodymass),
  estimate = estimate,
  SE = SE,
  z = z,
  p_coefficient = p,
  LRT_chisq = lrt_bodymass$Chisq[2],
  LRT_p = lrt_bodymass$`Pr(>Chisq)`[2],
  AIC_null = AIC(model_bodymass_null),
  AIC_model = AIC(model_bodymass),
  delta_AIC = AIC(model_bodymass_null) - AIC(model_bodymass),
  CI_low = CI_low,
  CI_high = CI_high,
  OR = exp(estimate),
  OR_low = exp(CI_low),
  OR_high = exp(CI_high)
)
results_bodymass$LRT_p_FDR <- results_bodymass$LRT_p
results_bodymass$p_coefficient_FDR <- results_bodymass$p_coefficient
cat("\n================ FINAL BODY MASS RESULTS ================\n")
print(
  results_bodymass,
  row.names = FALSE
)
```

### Variation in `hemoplasma` infection prevalence according to the host `bodymass` (conservative dataset, 16 `species`)
```
species_data_bodymass_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  filter(
    n >= 15
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, bodymass) %>%
      mutate(
        bodymass = as.numeric(bodymass)
      ),
    by = "species"
  ) %>%
  mutate(
    log_bodymass = log(bodymass),
    log_bodymass_scaled = as.numeric(
      scale(log_bodymass)
    )
  )

cat(
  "\nNumber of species =",
  nrow(species_data_bodymass_n15),
  "\n"
)

cat(
  "\nMissing bodymass values =",
  sum(is.na(species_data_bodymass_n15$bodymass)),
  "\n"
)

data_bodymass_n15 <- species_data_bodymass_n15 %>%
  filter(
    !is.na(bodymass),
    !is.na(log_bodymass_scaled)
  )

cat(
  "\nSpecies used for bodymass =",
  nrow(data_bodymass_n15),
  "\n"
)

model_bodymass_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_bodymass_n15,
  family = betabinomial(link = "logit")
)

model_bodymass_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ log_bodymass_scaled,
  data = data_bodymass_n15,
  family = betabinomial(link = "logit")
)

lrt_bodymass_n15 <- anova(
  model_bodymass_null_n15,
  model_bodymass_n15
)

cat(
  "\n================ BODY MASS LRT ================\n"
)

print(
  lrt_bodymass_n15
)

cat(
  "\n================ BODY MASS AIC ================\n"
)

print(
  AIC(
    model_bodymass_null_n15,
    model_bodymass_n15
  )
)

cat(
  "\n================ BODY MASS MODEL ================\n"
)

print(
  summary(model_bodymass_n15)
)

coef_table <- summary(
  model_bodymass_n15
)$coefficients$cond

estimate <- coef_table[
  "log_bodymass_scaled",
  "Estimate"
]

SE <- coef_table[
  "log_bodymass_scaled",
  "Std. Error"
]

z <- coef_table[
  "log_bodymass_scaled",
  "z value"
]

p <- coef_table[
  "log_bodymass_scaled",
  "Pr(>|z|)"
]

CI_low <- estimate - 1.96 * SE

CI_high <- estimate + 1.96 * SE

results_bodymass_n15 <- data.frame(
  variable = "bodymass",
  coefficient = "log_bodymass_scaled",
  n_species = nobs(
    model_bodymass_n15
  ),
  estimate = estimate,
  SE = SE,
  z = z,
  p_coefficient = p,
  LRT_chisq = lrt_bodymass_n15$Chisq[2],
  LRT_p = lrt_bodymass_n15$`Pr(>Chisq)`[2],
  AIC_null = AIC(
    model_bodymass_null_n15
  ),
  AIC_model = AIC(
    model_bodymass_n15
  ),
  delta_AIC = AIC(
    model_bodymass_null_n15
  ) - AIC(
    model_bodymass_n15
  ),
  CI_low = CI_low,
  CI_high = CI_high,
  OR = exp(estimate),
  OR_low = exp(CI_low),
  OR_high = exp(CI_high)
)

results_bodymass_n15$LRT_p_FDR <-
  results_bodymass_n15$LRT_p

results_bodymass_n15$p_coefficient_FDR <-
  results_bodymass_n15$p_coefficient

cat(
  "\n================ FINAL BODY MASS RESULTS ================\n"
)

print(
  results_bodymass_n15,
  row.names = FALSE
)
```
-> Results : In the complete dataset, `bodymass` showed a positive but marginal association with `hemoplasma` prevalence (β = 0.47 ± 0.27 SE, OR = 1.60, 95% CI: 0.95–2.70; LRT χ²₁ = 3.14, p = 0.077; ΔAIC = 1.14). In the conservative dataset , this association was stronger and statistically significant (β = 0.93 ± 0.36 SE, OR = 2.52, 95% CI: 1.26–5.07; LRT χ²₁ = 6.70, p = 0.010; ΔAIC = 4.70).

-> Interpretation : `hemoplasma` prevalence tended to increase with host `bodymass`, with stronger evidence for this association among better-sampled species.

### Visualization of association between `hemoplasma` prevalence and host `bodymass` (complete dataset, `species`)
```
species_order <- data_hemoplasma_stat %>%
  select(species, order) %>%
  distinct()

plot_data_bodymass <- data_bodymass %>%
  mutate(
    bodymass_plot = log_bodymass_scaled
  ) %>%
  left_join(
    species_order,
    by = "species"
  ) %>%
  mutate(
    species_initials = paste0(
      substr(species, 1, 1),
      substr(
        sub("^[^_]+_", "", species),
        1,
        1
      )
    )
  )

pred_data_bodymass <- pred_bodymass %>%
  mutate(
    bodymass_plot = log_bodymass_scaled
  )

label_data_bodymass <- plot_data_bodymass %>%
  filter(
    n >= 15
  )

order_colors <- c(
  "Primates" = "#E69F00",
  "Rodentia" = "#56B4E9",
  "Pilosa" = "#009E73",
  "Didelphimorphia" = "#CC79A7",
  "Carnivora" = "#D55E00",
  "Cingulata" = "#0072B2"
)

figure_bodymass <- ggplot(
  plot_data_bodymass,
  aes(
    x = bodymass_plot,
    y = prevalence,
    colour = order
  )
) +
  geom_ribbon(
    data = pred_data_bodymass,
    aes(
      x = bodymass_plot,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.20
  ) +
  geom_line(
    data = pred_data_bodymass,
    aes(
      x = bodymass_plot,
      y = fit
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    data = plot_data_bodymass %>%
      filter(n < 15),
    aes(
      size = n
    ),
    shape = 16,
    alpha = 0.70
  ) +
  geom_point(
    data = plot_data_bodymass %>%
      filter(n >= 15),
    aes(
      size = n
    ),
    shape = 21,
    fill = "white",
    alpha = 1,
    stroke = 1.1
  ) +
  ggrepel::geom_text_repel(
    data = label_data_bodymass,
    aes(
      label = species_initials,
      colour = order
    ),
    size = 3.5,
    family = "sans",
    fontface = "plain",
    box.padding = 1.0,
    point.padding = 0.8,
    force = 2.5,
    force_pull = 0.3,
    direction = "both",
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.5,
    segment.color = "black",
    segment.curvature = 0,
    segment.ncp = 1,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = order_colors
  ) +
  scale_size_continuous(
    name = "Number tested"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    x = "Log body mass (g)",
    y = "Hemoplasma prevalence",
    colour = "Mammalian order"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 10,
      family = "sans"
    ),
    axis.text.y = element_text(
      size = 10,
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 11,
      family = "sans"
    ),
    axis.title.y = element_text(
      size = 11,
      family = "sans"
    ),
    legend.title = element_text(
      size = 9,
      family = "sans"
    ),
    legend.text = element_text(
      size = 8,
      family = "sans"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

print(figure_bodymass)

ggsave(
  "Hemoplasma_prevalence_bodymass.png",
  figure_bodymass,
  width = 7,
  height = 4.8,
  dpi = 300
)
```

### Variation in `hemoplasma` infection status according to the host mean `longivity` (complete dataset, 44 `species`)
```
species_data_longevity <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, longevity),
    by = "species"
  ) %>%
  mutate(
    longevity = as.numeric(longevity),
    log_longevity = log(longevity),
    log_longevity_scaled = as.numeric(scale(log_longevity))
  )
cat("\nNumber of species =", nrow(species_data_longevity), "\n")
cat("\nMissing longevity values =", sum(is.na(species_data_longevity$longevity)), "\n")
data_longevity <- species_data_longevity %>%
  filter(
    !is.na(longevity),
    is.finite(longevity),
    !is.na(log_longevity_scaled)
  )
cat("\nSpecies used for longevity =", nrow(data_longevity), "\n")
model_longevity_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_longevity,
  family = betabinomial(link = "logit")
)
model_longevity <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ log_longevity_scaled,
  data = data_longevity,
  family = betabinomial(link = "logit")
)
lrt_longevity <- anova(
  model_longevity_null,
  model_longevity
)
cat("\n================ LONGEVITY LRT ================\n")
print(lrt_longevity)
cat("\n================ LONGEVITY AIC ================\n")
print(AIC(model_longevity_null, model_longevity))
cat("\n================ LONGEVITY MODEL ================\n")
print(summary(model_longevity))
coef_table <- summary(model_longevity)$coefficients$cond
estimate <- coef_table["log_longevity_scaled", "Estimate"]
SE <- coef_table["log_longevity_scaled", "Std. Error"]
z <- coef_table["log_longevity_scaled", "z value"]
p <- coef_table["log_longevity_scaled", "Pr(>|z|)"]
CI_low <- estimate - 1.96 * SE
CI_high <- estimate + 1.96 * SE
results_longevity <- data.frame(
  variable = "longevity",
  coefficient = "log_longevity_scaled",
  n_species = nobs(model_longevity),
  estimate = estimate,
  SE = SE,
  z = z,
  p_coefficient = p,
  LRT_chisq = lrt_longevity$Chisq[2],
  LRT_p = lrt_longevity$`Pr(>Chisq)`[2],
  AIC_null = AIC(model_longevity_null),
  AIC_model = AIC(model_longevity),
  delta_AIC = AIC(model_longevity_null) - AIC(model_longevity),
  CI_low = CI_low,
  CI_high = CI_high,
  OR = exp(estimate),
  OR_low = exp(CI_low),
  OR_high = exp(CI_high)
)
results_longevity$LRT_p_FDR <- results_longevity$LRT_p
results_longevity$p_coefficient_FDR <- results_longevity$p_coefficient
cat("\n================ FINAL LONGEVITY RESULTS ================\n")
print(results_longevity, row.names = FALSE)
```
### Variation in `hemoplasma` infection status according to the host mean `longivity` (conservative dataset, 16 `species`)
```
species_data_longevity_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  filter(
    n >= 15
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, longevity),
    by = "species"
  ) %>%
  mutate(
    longevity = as.numeric(longevity),
    log_longevity = log(longevity),
    log_longevity_scaled = as.numeric(
      scale(log_longevity)
    )
  )

cat(
  "\nNumber of species =",
  nrow(species_data_longevity_n15),
  "\n"
)

cat(
  "\nMissing longevity values =",
  sum(is.na(species_data_longevity_n15$longevity)),
  "\n"
)

data_longevity_n15 <- species_data_longevity_n15 %>%
  filter(
    !is.na(longevity),
    is.finite(longevity),
    !is.na(log_longevity_scaled)
  )

cat(
  "\nSpecies used for longevity =",
  nrow(data_longevity_n15),
  "\n"
)

model_longevity_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_longevity_n15,
  family = betabinomial(link = "logit")
)

model_longevity_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ log_longevity_scaled,
  data = data_longevity_n15,
  family = betabinomial(link = "logit")
)

lrt_longevity_n15 <- anova(
  model_longevity_null_n15,
  model_longevity_n15
)

cat(
  "\n================ LONGEVITY LRT ================\n"
)

print(
  lrt_longevity_n15
)

cat(
  "\n================ LONGEVITY AIC ================\n"
)

print(
  AIC(
    model_longevity_null_n15,
    model_longevity_n15
  )
)

cat(
  "\n================ LONGEVITY MODEL ================\n"
)

print(
  summary(model_longevity_n15)
)

coef_table <- summary(
  model_longevity_n15
)$coefficients$cond

estimate <- coef_table[
  "log_longevity_scaled",
  "Estimate"
]

SE <- coef_table[
  "log_longevity_scaled",
  "Std. Error"
]

z <- coef_table[
  "log_longevity_scaled",
  "z value"
]

p <- coef_table[
  "log_longevity_scaled",
  "Pr(>|z|)"
]

CI_low <- estimate - 1.96 * SE
CI_high <- estimate + 1.96 * SE

results_longevity_n15 <- data.frame(
  variable = "longevity",
  coefficient = "log_longevity_scaled",
  n_species = nobs(
    model_longevity_n15
  ),
  estimate = estimate,
  SE = SE,
  z = z,
  p_coefficient = p,
  LRT_chisq = lrt_longevity_n15$Chisq[2],
  LRT_p = lrt_longevity_n15$`Pr(>Chisq)`[2],
  AIC_null = AIC(
    model_longevity_null_n15
  ),
  AIC_model = AIC(
    model_longevity_n15
  ),
  delta_AIC = AIC(
    model_longevity_null_n15
  ) - AIC(
    model_longevity_n15
  ),
  CI_low = CI_low,
  CI_high = CI_high,
  OR = exp(estimate),
  OR_low = exp(CI_low),
  OR_high = exp(CI_high)
)

results_longevity_n15$LRT_p_FDR <-
  results_longevity_n15$LRT_p

results_longevity_n15$p_coefficient_FDR <-
  results_longevity_n15$p_coefficient

cat(
  "\n================ FINAL LONGEVITY RESULTS ================\n"
)

print(
  results_longevity_n15,
  row.names = FALSE
)
```
-> Results : In the complete dataset, across the 33 mammalian `species` with available `longevity` data, `longevity` was not significantly associated with `hemoplasma` prevalence (beta-binomial GLM, LRT χ²₁ = 2.13, *p* = 0.145; ΔAIC = 0.13; OR = 1.57, 95% CI: 0.86–2.88). In the conservative dataset, 13 species had available `longevity` data, and the association remained non-significant but showed a similar positive trend (LRT χ²₁ = 2.91, *p* = 0.088; ΔAIC = 0.91; OR = 1.96, 95% CI: 0.91–4.22).

-> Interpretation : `hemoplasma` prevalence tended to increase with host `longevity`, with a stronger but still marginal association among better-sampled species.

### Visualization of `hemoplasma` infection status according to the host mean `longivity`(complete dataset, 44 `species`)
```
species_order <- data_hemoplasma_stat %>%
  select(
    species,
    order
  ) %>%
  distinct()

order_colors <- c(
  "Primates" = "#E69F00",
  "Rodentia" = "#56B4E9",
  "Pilosa" = "#009E73",
  "Didelphimorphia" = "#CC79A7",
  "Carnivora" = "#D55E00",
  "Cingulata" = "#0072B2"
)

plot_longevity <- data_longevity %>%
  left_join(
    species_order,
    by = "species"
  ) %>%
  mutate(
    species_initials = paste0(
      substr(species, 1, 1),
      substr(
        sub("^[^_]+_", "", species),
        1,
        1
      )
    ),
    n15 = n >= 15
  )

label_data_longevity <- plot_longevity %>%
  filter(
    n15,
    !is.na(log_longevity),
    !is.na(prevalence)
  ) %>%
  arrange(
    log_longevity,
    prevalence
  ) %>%
  mutate(
    label_side = rep(
      c("right", "left"),
      length.out = n()
    ),
    label_x = ifelse(
      label_side == "right",
      log_longevity + 0.08,
      log_longevity - 0.08
    ),
    label_y = prevalence
  )

pred_longevity <- data.frame(
  log_longevity_scaled = seq(
    min(
      data_longevity$log_longevity_scaled,
      na.rm = TRUE
    ),
    max(
      data_longevity$log_longevity_scaled,
      na.rm = TRUE
    ),
    length.out = 100
  )
)

pred <- predict(
  model_longevity,
  newdata = pred_longevity,
  type = "link",
  se.fit = TRUE
)

pred_longevity <- pred_longevity %>%
  mutate(
    fit = plogis(pred$fit),
    lower = plogis(
      pred$fit - 1.96 * pred$se.fit
    ),
    upper = plogis(
      pred$fit + 1.96 * pred$se.fit
    )
  )

longevity_range <- range(
  data_longevity$log_longevity,
  na.rm = TRUE
)

pred_longevity <- pred_longevity %>%
  mutate(
    log_longevity = seq(
      longevity_range[1],
      longevity_range[2],
      length.out = n()
    )
  )

figure_longevity <- ggplot(
  plot_longevity,
  aes(
    x = log_longevity,
    y = prevalence
  )
) +
  geom_ribbon(
    data = pred_longevity,
    aes(
      x = log_longevity,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.20
  ) +
  geom_line(
    data = pred_longevity,
    aes(
      x = log_longevity,
      y = fit
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    data = plot_longevity %>%
      filter(!n15),
    aes(
      size = n,
      colour = order
    ),
    shape = 16,
    alpha = 0.70
  ) +
  geom_point(
    data = plot_longevity %>%
      filter(n15),
    aes(
      size = n,
      colour = order
    ),
    shape = 1,
    stroke = 1.1
  ) +
  geom_segment(
    data = label_data_longevity,
    aes(
      x = log_longevity,
      xend = label_x,
      y = prevalence,
      yend = label_y,
      colour = order
    ),
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_text(
    data = label_data_longevity,
    aes(
      x = label_x,
      y = label_y,
      label = species_initials,
      colour = order,
      hjust = ifelse(
        label_side == "right",
        0,
        1
      )
    ),
    size = 3.2,
    family = "sans",
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = order_colors,
    name = "Mammalian order"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  scale_size_continuous(
    name = "Number tested"
  ) +
  labs(
    x = "Log host longevity (years)",
    y = "Hemoplasma prevalence"
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "sans"
    ),
    axis.text.x = element_text(
      size = 10,
      family = "sans"
    ),
    axis.text.y = element_text(
      size = 10,
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 11,
      family = "sans"
    ),
    axis.title.y = element_text(
      size = 11,
      family = "sans"
    ),
    legend.title = element_text(
      size = 9,
      family = "sans"
    ),
    legend.text = element_text(
      size = 8,
      family = "sans"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

print(figure_longevity)

ggsave(
  "Hemoplasma_prevalence_longevity_beta_binomial.png",
  figure_longevity,
  width = 7,
  height = 4.8,
  dpi = 300
)
```

### Variation in `hemoplasma` infection status according to the `femalematurity` (complete dataset, 44 `species`)
```
species_data_femalematurity <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, femalematurity),
    by = "species"
  ) %>%
  mutate(
    femalematurity = as.numeric(femalematurity),
    log_femalematurity = log(femalematurity),
    log_femalematurity_scaled = as.numeric(scale(log_femalematurity))
  )
cat("\nNumber of species =", nrow(species_data_femalematurity), "\n")
cat("\nMissing femalematurity values =", sum(is.na(species_data_femalematurity$femalematurity)), "\n")
data_femalematurity <- species_data_femalematurity %>%
  filter(
    !is.na(femalematurity),
    is.finite(femalematurity),
    !is.na(log_femalematurity_scaled)
  )
cat("\nSpecies used for femalematurity =", nrow(data_femalematurity), "\n")
model_femalematurity_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_femalematurity,
  family = betabinomial(link = "logit")
)
model_femalematurity <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ log_femalematurity_scaled,
  data = data_femalematurity,
  family = betabinomial(link = "logit")
)
lrt_femalematurity <- anova(
  model_femalematurity_null,
  model_femalematurity
)
cat("\n================ FEMALE MATURITY LRT ================\n")
print(lrt_femalematurity)
cat("\n================ FEMALE MATURITY AIC ================\n")
print(AIC(model_femalematurity_null, model_femalematurity))
cat("\n================ FEMALE MATURITY MODEL ================\n")
print(summary(model_femalematurity))
coef_table <- summary(model_femalematurity)$coefficients$cond
estimate <- coef_table["log_femalematurity_scaled", "Estimate"]
SE <- coef_table["log_femalematurity_scaled", "Std. Error"]
z <- coef_table["log_femalematurity_scaled", "z value"]
p <- coef_table["log_femalematurity_scaled", "Pr(>|z|)"]
CI_low <- estimate - 1.96 * SE
CI_high <- estimate + 1.96 * SE
results_femalematurity <- data.frame(
  variable = "femalematurity",
  coefficient = "log_femalematurity_scaled",
  n_species = nobs(model_femalematurity),
  estimate = estimate,
  SE = SE,
  z = z,
  p_coefficient = p,
  LRT_chisq = lrt_femalematurity$Chisq[2],
  LRT_p = lrt_femalematurity$`Pr(>Chisq)`[2],
  AIC_null = AIC(model_femalematurity_null),
  AIC_model = AIC(model_femalematurity),
  delta_AIC = AIC(model_femalematurity_null) - AIC(model_femalematurity),
  CI_low = CI_low,
  CI_high = CI_high,
  OR = exp(estimate),
  OR_low = exp(CI_low),
  OR_high = exp(CI_high)
)
results_femalematurity$LRT_p_FDR <- results_femalematurity$LRT_p
results_femalematurity$p_coefficient_FDR <- results_femalematurity$p_coefficient
cat("\n================ FINAL FEMALE MATURITY RESULTS ================\n")
print(results_femalematurity, row.names = FALSE)
```
### Variation in `hemoplasma` infection status according to the `femalematurity` (conservative dataset, 16 `species`)
```
species_data_femalematurity_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  filter(
    n >= 15
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, femalematurity),
    by = "species"
  ) %>%
  mutate(
    femalematurity = as.numeric(femalematurity),
    log_femalematurity = log(femalematurity),
    log_femalematurity_scaled = as.numeric(
      scale(log_femalematurity)
    )
  )

cat(
  "\nNumber of species =",
  nrow(species_data_femalematurity_n15),
  "\n"
)

cat(
  "\nMissing femalematurity values =",
  sum(
    is.na(
      species_data_femalematurity_n15$femalematurity
    )
  ),
  "\n"
)

data_femalematurity_n15 <- species_data_femalematurity_n15 %>%
  filter(
    !is.na(femalematurity),
    is.finite(femalematurity),
    !is.na(log_femalematurity_scaled)
  )

cat(
  "\nSpecies used for femalematurity =",
  nrow(data_femalematurity_n15),
  "\n"
)

model_femalematurity_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_femalematurity_n15,
  family = betabinomial(link = "logit")
)

model_femalematurity_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ log_femalematurity_scaled,
  data = data_femalematurity_n15,
  family = betabinomial(link = "logit")
)

lrt_femalematurity_n15 <- anova(
  model_femalematurity_null_n15,
  model_femalematurity_n15
)

cat(
  "\n================ FEMALE MATURITY LRT ================\n"
)

print(
  lrt_femalematurity_n15
)

cat(
  "\n================ FEMALE MATURITY AIC ================\n"
)

print(
  AIC(
    model_femalematurity_null_n15,
    model_femalematurity_n15
  )
)

cat(
  "\n================ FEMALE MATURITY MODEL ================\n"
)

print(
  summary(
    model_femalematurity_n15
  )
)

coef_table <- summary(
  model_femalematurity_n15
)$coefficients$cond

estimate <- coef_table[
  "log_femalematurity_scaled",
  "Estimate"
]

SE <- coef_table[
  "log_femalematurity_scaled",
  "Std. Error"
]

z <- coef_table[
  "log_femalematurity_scaled",
  "z value"
]

p <- coef_table[
  "log_femalematurity_scaled",
  "Pr(>|z|)"
]

CI_low <- estimate - 1.96 * SE
CI_high <- estimate + 1.96 * SE

results_femalematurity_n15 <- data.frame(
  variable = "femalematurity",
  coefficient = "log_femalematurity_scaled",
  n_species = nobs(
    model_femalematurity_n15
  ),
  estimate = estimate,
  SE = SE,
  z = z,
  p_coefficient = p,
  LRT_chisq = lrt_femalematurity_n15$Chisq[2],
  LRT_p = lrt_femalematurity_n15$`Pr(>Chisq)`[2],
  AIC_null = AIC(
    model_femalematurity_null_n15
  ),
  AIC_model = AIC(
    model_femalematurity_n15
  ),
  delta_AIC = AIC(
    model_femalematurity_null_n15
  ) - AIC(
    model_femalematurity_n15
  ),
  CI_low = CI_low,
  CI_high = CI_high,
  OR = exp(estimate),
  OR_low = exp(CI_low),
  OR_high = exp(CI_high)
)

results_femalematurity_n15$LRT_p_FDR <-
  results_femalematurity_n15$LRT_p

results_femalematurity_n15$p_coefficient_FDR <-
  results_femalematurity_n15$p_coefficient

cat(
  "\n================ FINAL FEMALE MATURITY RESULTS ================\n"
)

print(
  results_femalematurity_n15,
  row.names = FALSE
)
```
-> Results: `femalematurity` was significantly associated with `hemoplasma` prevalence across the 26 mammalian `species` (beta-binomial GLM, LRT χ²₁ = 4.30, *p* = 0.038; ΔAIC = 2.30; OR = 2.03, 95% CI: 1.05–3.92). In the conservative dataset, 14 species had available `femalematurity` data, and the association was stronger (LRT χ²₁ = 10.46, *p* = 0.001; ΔAIC = 8.46; OR = 3.79, 95% CI: 1.68–8.56).

-> Interpretation: `hemoplasma` prevalence increased with later female age at maturity, with consistent evidence for this positive association in both datasets and stronger support among better-sampled species.

### Visualization of `hemoplasma` infection status according to the `femalematurity` (complete dataset, 44 `species`)
```
species_order <- data_hemoplasma_stat %>%
  select(
    species,
    order
  ) %>%
  distinct()

order_colors <- c(
  "Primates" = "#E69F00",
  "Rodentia" = "#56B4E9",
  "Pilosa" = "#009E73",
  "Didelphimorphia" = "#CC79A7",
  "Carnivora" = "#D55E00",
  "Cingulata" = "#0072B2"
)

plot_femalematurity <- data_femalematurity %>%
  left_join(
    species_order,
    by = "species"
  ) %>%
  mutate(
    species_initials = paste0(
      substr(species, 1, 1),
      substr(
        sub("^[^_]+_", "", species),
        1,
        1
      )
    ),
    n15 = n >= 15
  )

label_data_femalematurity <- plot_femalematurity %>%
  filter(
    n15,
    !is.na(log_femalematurity),
    !is.na(prevalence)
  ) %>%
  arrange(
    log_femalematurity,
    prevalence
  ) %>%
  mutate(
    label_side = rep(
      c("right", "left"),
      length.out = n()
    ),
    label_x = ifelse(
      label_side == "right",
      log_femalematurity + 0.08,
      log_femalematurity - 0.08
    ),
    label_y = prevalence
  )

pred_femalematurity <- data.frame(
  log_femalematurity_scaled = seq(
    min(
      data_femalematurity$log_femalematurity_scaled,
      na.rm = TRUE
    ),
    max(
      data_femalematurity$log_femalematurity_scaled,
      na.rm = TRUE
    ),
    length.out = 100
  )
)

pred <- predict(
  model_femalematurity,
  newdata = pred_femalematurity,
  type = "link",
  se.fit = TRUE
)

pred_femalematurity <- pred_femalematurity %>%
  mutate(
    fit = plogis(pred$fit),
    lower = plogis(
      pred$fit - 1.96 * pred$se.fit
    ),
    upper = plogis(
      pred$fit + 1.96 * pred$se.fit
    )
  )

maturity_range <- range(
  data_femalematurity$log_femalematurity,
  na.rm = TRUE
)

pred_femalematurity <- pred_femalematurity %>%
  mutate(
    log_femalematurity = seq(
      maturity_range[1],
      maturity_range[2],
      length.out = n()
    )
  )

figure_femalematurity <- ggplot(
  plot_femalematurity,
  aes(
    x = log_femalematurity,
    y = prevalence
  )
) +
  geom_ribbon(
    data = pred_femalematurity,
    aes(
      x = log_femalematurity,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.20
  ) +
  geom_line(
    data = pred_femalematurity,
    aes(
      x = log_femalematurity,
      y = fit
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    data = plot_femalematurity %>%
      filter(!n15),
    aes(
      size = n,
      colour = order
    ),
    shape = 16,
    alpha = 0.70
  ) +
  geom_point(
    data = plot_femalematurity %>%
      filter(n15),
    aes(
      size = n,
      colour = order
    ),
    shape = 1,
    stroke = 1.1
  ) +
  geom_segment(
    data = label_data_femalematurity,
    aes(
      x = log_femalematurity,
      xend = label_x,
      y = prevalence,
      yend = label_y,
      colour = order
    ),
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_text(
    data = label_data_femalematurity,
    aes(
      x = label_x,
      y = label_y,
      label = species_initials,
      colour = order,
      hjust = ifelse(
        label_side == "right",
        0,
        1
      )
    ),
    size = 3.2,
    family = "sans",
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = order_colors,
    name = "Mammalian order"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  scale_size_continuous(
    name = "Number tested"
  ) +
  labs(
    x = "Log female age at maturity (days)",
    y = "Hemoplasma prevalence"
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "sans"
    ),
    axis.text.x = element_text(
      size = 10,
      family = "sans"
    ),
    axis.text.y = element_text(
      size = 10,
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 11,
      family = "sans"
    ),
    axis.title.y = element_text(
      size = 11,
      family = "sans"
    ),
    legend.title = element_text(
      size = 9,
      family = "sans"
    ),
    legend.text = element_text(
      size = 8,
      family = "sans"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

print(figure_femalematurity)

ggsave(
  "Hemoplasma_prevalence_female_maturity_beta_binomial.png",
  figure_femalematurity,
  width = 7,
  height = 4.8,
  dpi = 300
)
```

### Variation of `hemoplasma` infection status according to the host `littersize` (complete dataset, 44 `species`)
```
species_data_littersize <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, littersize),
    by = "species"
  ) %>%
  mutate(
    littersize = as.numeric(littersize),
    log_littersize = log(littersize),
    log_littersize_scaled = as.numeric(scale(log_littersize))
  )
cat("\nNumber of species =", nrow(species_data_littersize), "\n")
cat("\nMissing littersize values =", sum(is.na(species_data_littersize$littersize)), "\n")
data_littersize <- species_data_littersize %>%
  filter(
    !is.na(littersize),
    is.finite(littersize),
    !is.na(log_littersize_scaled)
  )
cat("\nSpecies used for littersize =", nrow(data_littersize), "\n")
model_littersize_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_littersize,
  family = betabinomial(link = "logit")
)
model_littersize <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ log_littersize_scaled,
  data = data_littersize,
  family = betabinomial(link = "logit")
)
lrt_littersize <- anova(
  model_littersize_null,
  model_littersize
)
cat("\n================ LITTER SIZE LRT ================\n")
print(lrt_littersize)
cat("\n================ LITTER SIZE AIC ================\n")
print(AIC(model_littersize_null, model_littersize))
cat("\n================ LITTER SIZE MODEL ================\n")
print(summary(model_littersize))
coef_table <- summary(model_littersize)$coefficients$cond
estimate <- coef_table["log_littersize_scaled", "Estimate"]
SE <- coef_table["log_littersize_scaled", "Std. Error"]
z <- coef_table["log_littersize_scaled", "z value"]
p <- coef_table["log_littersize_scaled", "Pr(>|z|)"]
CI_low <- estimate - 1.96 * SE
CI_high <- estimate + 1.96 * SE
results_littersize <- data.frame(
  variable = "littersize",
  coefficient = "log_littersize_scaled",
  n_species = nobs(model_littersize),
  estimate = estimate,
  SE = SE,
  z = z,
  p_coefficient = p,
  LRT_chisq = lrt_littersize$Chisq[2],
  LRT_p = lrt_littersize$`Pr(>Chisq)`[2],
  AIC_null = AIC(model_littersize_null),
  AIC_model = AIC(model_littersize),
  delta_AIC = AIC(model_littersize_null) - AIC(model_littersize),
  CI_low = CI_low,
  CI_high = CI_high,
  OR = exp(estimate),
  OR_low = exp(CI_low),
  OR_high = exp(CI_high)
)
results_littersize$LRT_p_FDR <- results_littersize$LRT_p
results_littersize$p_coefficient_FDR <- results_littersize$p_coefficient
cat("\n================ FINAL LITTER SIZE RESULTS ================\n")
print(results_littersize, row.names = FALSE)
```
### Variation of `hemoplasma` infection status according to the host `littersize` (conservative dataset, 16 `species`)
```
species_data_littersize_n15 <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  filter(
    n >= 15
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, littersize),
    by = "species"
  ) %>%
  mutate(
    littersize = as.numeric(littersize),
    log_littersize = log(littersize),
    log_littersize_scaled = as.numeric(
      scale(log_littersize)
    )
  )

cat(
  "\nNumber of species =",
  nrow(species_data_littersize_n15),
  "\n"
)

cat(
  "\nMissing littersize values =",
  sum(
    is.na(
      species_data_littersize_n15$littersize
    )
  ),
  "\n"
)

data_littersize_n15 <- species_data_littersize_n15 %>%
  filter(
    !is.na(littersize),
    is.finite(littersize),
    !is.na(log_littersize_scaled)
  )

cat(
  "\nSpecies used for littersize =",
  nrow(data_littersize_n15),
  "\n"
)

model_littersize_null_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_littersize_n15,
  family = betabinomial(link = "logit")
)

model_littersize_n15 <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ log_littersize_scaled,
  data = data_littersize_n15,
  family = betabinomial(link = "logit")
)

lrt_littersize_n15 <- anova(
  model_littersize_null_n15,
  model_littersize_n15
)

cat(
  "\n================ LITTER SIZE LRT ================\n"
)

print(
  lrt_littersize_n15
)

cat(
  "\n================ LITTER SIZE AIC ================\n"
)

print(
  AIC(
    model_littersize_null_n15,
    model_littersize_n15
  )
)

cat(
  "\n================ LITTER SIZE MODEL ================\n"
)

print(
  summary(
    model_littersize_n15
  )
)

coef_table <- summary(
  model_littersize_n15
)$coefficients$cond

estimate <- coef_table[
  "log_littersize_scaled",
  "Estimate"
]

SE <- coef_table[
  "log_littersize_scaled",
  "Std. Error"
]

z <- coef_table[
  "log_littersize_scaled",
  "z value"
]

p <- coef_table[
  "log_littersize_scaled",
  "Pr(>|z|)"
]

CI_low <- estimate - 1.96 * SE
CI_high <- estimate + 1.96 * SE

results_littersize_n15 <- data.frame(
  variable = "littersize",
  coefficient = "log_littersize_scaled",
  n_species = nobs(
    model_littersize_n15
  ),
  estimate = estimate,
  SE = SE,
  z = z,
  p_coefficient = p,
  LRT_chisq = lrt_littersize_n15$Chisq[2],
  LRT_p = lrt_littersize_n15$`Pr(>Chisq)`[2],
  AIC_null = AIC(
    model_littersize_null_n15
  ),
  AIC_model = AIC(
    model_littersize_n15
  ),
  delta_AIC = AIC(
    model_littersize_null_n15
  ) - AIC(
    model_littersize_n15
  ),
  CI_low = CI_low,
  CI_high = CI_high,
  OR = exp(estimate),
  OR_low = exp(CI_low),
  OR_high = exp(CI_high)
)

results_littersize_n15$LRT_p_FDR <-
  results_littersize_n15$LRT_p

results_littersize_n15$p_coefficient_FDR <-
  results_littersize_n15$p_coefficient

cat(
  "\n================ FINAL LITTER SIZE RESULTS ================\n"
)

print(
  results_littersize_n15,
  row.names = FALSE
)
```
-> Results: `littersize` was not associated with `hemoplasma` prevalence among the 39 mammalian `species` with available data (beta-binomial GLM, LRT χ²₁ = 0.10, *p* = 0.749; ΔAIC = −1.90; OR = 1.09, 95% CI: 0.64–1.85). This result was consistent in the conservative dataset, where `littersize` showed no significant association with `hemoplasma` prevalence among the 14 species with available data (LRT χ²₁ = 1.78, *p* = 0.182; ΔAIC = −0.22; OR = 0.60, 95% CI: 0.29–1.27).

-> Interpretation: There was no evidence that litter size explained interspecific variation in `hemoplasma` prevalence, with no consistent association in either dataset.

### Visualization of `hemoplasma` infection status according to host `littersize` (complete dataset, 44 `species`)
```
species_order <- data_hemoplasma_stat %>%
  select(
    species,
    order
  ) %>%
  distinct()

order_colors <- c(
  "Primates" = "#E69F00",
  "Rodentia" = "#56B4E9",
  "Pilosa" = "#009E73",
  "Didelphimorphia" = "#CC79A7",
  "Carnivora" = "#D55E00",
  "Cingulata" = "#0072B2"
)

plot_littersize <- data_littersize %>%
  left_join(
    species_order,
    by = "species"
  ) %>%
  mutate(
    species_initials = paste0(
      substr(species, 1, 1),
      substr(
        sub("^[^_]+_", "", species),
        1,
        1
      )
    ),
    conservative = n >= 15
  )

x_seq <- seq(
  min(
    data_littersize$littersize,
    na.rm = TRUE
  ),
  max(
    data_littersize$littersize,
    na.rm = TRUE
  ),
  length.out = 100
)

newdata_littersize <- data.frame(
  littersize = x_seq
)

newdata_littersize$log_littersize_scaled <- (
  log(newdata_littersize$littersize) -
    mean(
      log(data_littersize$littersize),
      na.rm = TRUE
    )
) / sd(
  log(data_littersize$littersize),
  na.rm = TRUE
)

pred_littersize <- predict(
  model_littersize,
  newdata = newdata_littersize,
  type = "link",
  se.fit = TRUE
)

newdata_littersize$fit <- plogis(
  pred_littersize$fit
)

newdata_littersize$lower <- plogis(
  pred_littersize$fit -
    1.96 * pred_littersize$se.fit
)

newdata_littersize$upper <- plogis(
  pred_littersize$fit +
    1.96 * pred_littersize$se.fit
)

label_data_littersize <- plot_littersize %>%
  filter(
    conservative,
    !is.na(littersize),
    !is.na(prevalence),
    !is.na(order)
  )

figure_littersize <- ggplot(
  plot_littersize,
  aes(
    x = littersize,
    y = prevalence
  )
) +
  geom_ribbon(
    data = newdata_littersize,
    aes(
      x = littersize,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.20
  ) +
  geom_line(
    data = newdata_littersize,
    aes(
      x = littersize,
      y = fit
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    data = plot_littersize %>%
      filter(!conservative),
    aes(
      size = n,
      colour = order
    ),
    alpha = 0.70
  ) +
  geom_point(
    data = plot_littersize %>%
      filter(conservative),
    aes(
      size = n,
      colour = order
    ),
    shape = 1,
    stroke = 1.1
  ) +
  ggrepel::geom_text_repel(
    data = label_data_littersize,
    aes(
      label = species_initials,
      colour = order
    ),
    size = 3.2,
    family = "sans",
    fontface = "plain",
    box.padding = 0.45,
    point.padding = 0.35,
    force = 1.5,
    force_pull = 0.4,
    direction = "both",
    min.segment.length = 0,
    segment.size = 0.5,
    segment.curvature = 0,
    segment.ncp = 1,
    segment.angle = 90,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = order_colors,
    name = "Mammalian order"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0.02, 0.08)
    )
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0.04, 0.08)
    )
  ) +
  scale_size_continuous(
    name = "Number tested"
  ) +
  labs(
    x = "Litter size",
    y = "Hemoplasma prevalence"
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "sans"
    ),
    axis.text.x = element_text(
      size = 10
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.title.x = element_text(
      size = 11
    ),
    axis.title.y = element_text(
      size = 11
    ),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

print(figure_littersize)

ggsave(
  "Hemoplasma_prevalence_littersize_beta_binomial.png",
  figure_littersize,
  width = 7,
  height = 4.8,
  dpi = 300
)
```
