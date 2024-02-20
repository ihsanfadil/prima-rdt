
# Project     : Double-band positivity of HRP2/pLDH rapid diagnostic tests with
#               high-density Plasmodium falciparum parasitaemia in a low
#               transmission setting (Christian et al.)
#
# Code author : Ihsan Fadilah
# Email       : ifadilah@oucru.org

# Setup -------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(ggbeeswarm)
library(broom)
library(scales)
library(here)
library(ggExtra)
library(xtable); library(knitr)
library(extrafont); loadfonts()

theme_set(theme_bw())
theme_update(
  text = element_text(size = 9.5, family = "Fira Code"), # Font
  plot.title = element_text(hjust = 0),      # Centre-align title
  plot.subtitle = element_text(hjust = 0),   # Centre-align subtitle
  legend.title = element_blank(),            # Remove legend title
  # legend.position = c(0.80, 0.15),         # Move legend to bottom right
  legend.background = element_blank(),       # Remove legend background
  legend.box.background = element_blank(),   # Remove lengend-box background
  legend.spacing.y = unit(0.01, 'mm'),       # Make legend closer
  legend.key.height = unit(0.4, "cm"),       # Make legend closer
  # panel.grid.minor = element_blank(),      # Remove minor lines
  panel.grid.minor.x = element_blank(),      # Remove minor lines on the x axis
  axis.title.x = element_text(hjust = 1),    # Move title for x-axis
  axis.title.y = element_text(hjust = 0.5)   # Move title for y-axis
)

# Import data -------------------------------------------------------------

prima_raw <- read_excel('ACROSS For RDT_20240127_MC.xlsx')
mal_dens_raw <- read_excel('PRIMA-ACROSS_Microcopy and RDT_20230223_MC.xlsx')

# Clean data --------------------------------------------------------------

# Malaria density
mal_dens <- mal_dens_raw |> 
  select(CODE, MalDens) |> 
  janitor::clean_names() |> 
  arrange(code)

# All data (n = 317)
prima <- prima_raw |> 
  select(CODE, AGE, SEX, WEIGHT, TEMP, MRDT, SPC, APC) |>
  janitor::clean_names() |> 
  arrange(code) |> 
  mutate(mrdt = case_when(mrdt == 0 ~ 'Negative',
                          mrdt == 1 ~ 'Pf (HRP-2)',
                          mrdt == 2 ~ 'Pan (pLDH)',
                          mrdt == 3 ~ 'Pf/Pan (HRP-2/pLDH)',
                          mrdt == 5 ~ 'Not done', # None found
                          is.na(mrdt) ~ NA_character_,
                          TRUE ~ 'Check me!') |>
           factor(levels = c('Pf/Pan (HRP-2/pLDH)',
                             'Pf (HRP-2)',
                             'Pan (pLDH)',
                             'Negative',
                             'ND')),
         spc = case_when(spc == 0 ~ 'Negative',
                         spc == 1 ~ 'Pf',
                         spc == 2 ~ 'Pv',
                         spc == 3 ~ 'Pm',
                         spc == 5 ~ 'Pf and Pv',
                         spc == 6 ~ 'Pf and Pm',
                         is.na(spc) ~ NA_character_,
                         TRUE ~ 'Check me!') |> factor(),
         sex = if_else(sex == 1, 'Male', 'Female') |> factor(),
         apc_per1000 = apc / 1000,
         log2_apc = if_else(apc == 0, 0, log2(apc)))  # 0 not transformed
         
# Pf-only data (n = 149, n = 148 excluding one patient with missing `apc`)
prima_pf <- prima |> 
  drop_na(apc) |> 
  
  # Exclude negative and Pan-only MRDT results
  # i.e., include patients with Pf
  filter(mrdt != 'Negative', mrdt != 'Pan (pLDH)') |> 
  
  # Outcome for regression
  mutate(double = if_else(mrdt == 'Pf/Pan (HRP-2/pLDH)', 1L, 0L))

# Descriptives ------------------------------------------------------------

(means <- prima_pf |> 
   group_by(mrdt) |> 
   summarise(mean_log2 = mean(log2_apc),
             geo_mean = 2^(mean_log2),
             median_ori = median(apc),
             mean_ori = mean(apc)))

(apc_histogram_ori <- prima_pf |> 
  ggplot(aes(x = apc)) +
  geom_histogram(aes(fill = mrdt), alpha = 0.8, bins = 13) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(expand = c(0, 0, 0, 5),
                     breaks = seq(0, 120, by = 20)) +
  scale_fill_manual(values = c('gray30', 'gray60')) +
  theme(legend.position = c(0.8, 0.15),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = '\nParasite density (per \U00B5L)',
       y = 'Number of patients\n'))

ggsave(plot = apc_histogram_ori,
       filename = "apc_histogram_ori.png",
       height = 3, width = 3,
       dpi = 1200)

(apc_histogram_log <- prima_pf |> 
    ggplot(aes(x = log2_apc)) +
    geom_histogram(aes(fill = mrdt), alpha = 0.8, bins = 13) +
    # scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(expand = c(0, 0, 0, 5)) +
    scale_fill_manual(values = c('gray30', 'gray60')) +
    theme(legend.position = 'right',
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank()) +
    labs(x = '\nParasite density (per \U00B5L)',
         y = 'Number of patients\n'))

ggsave(plot = apc_histogram_ori,
       filename = "apc_histogram_ori.png",
       height = 3, width = 5,
       dpi = 1200)

# Modelling ---------------------------------------------------------------

  # negative spc, 0 mal_dens. exclude?
  # adjust for 4 sites, try also age.
  # show confidence intervals along the expected probability (ggpredict)














