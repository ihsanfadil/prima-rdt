
# Project     : Double-band positivity of HRP2/pLDH rapid diagnostic tests with
#               high-density Plasmodium falciparum parasitaemia in a low
#               transmission setting (Christian et al.)
#
# Code author : Ihsan Fadilah
# Email       : ifadilah@oucru.org

# Setup -------------------------------------------------------------------

library(tidyverse)
library(ggbeeswarm)
library(broom)
library(scales)
library(ggExtra)
library(ggeffects)
library(cowplot)
library(rms)
library(xtable); library(knitr)
library(extrafont); loadfonts()

options(scipen = 999)
theme_set(theme_bw())
theme_update(
  text = element_text(size = 8.5, family = "Fira Code"), # Font
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

prima_raw <- rio::import('ACROSS For RDT_20240222_revised.xlsx')

# Clean data --------------------------------------------------------------

# All data (n = 330)
prima <- prima_raw |> 
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
         film = case_when(film == 1 ~ 'Thick',
                          film == 2 ~ 'Thin',
                          is.na(film) ~ NA_character_,
                          TRUE ~ 'Check me!') |> factor(),
         bcell_counted = case_when(per == 1 ~ 'RBC',
                                   per %in% c(2, 3, 4) ~ 'WBC',
                                   is.na(per) ~ NA_character_,
                                   TRUE ~ 'Check me!') |>
           factor(levels = c('WBC', 'RBC')),
         bcell_mult = case_when(
           bcell_counted == 'RBC' ~ 5e6,
           bcell_counted == 'WBC' ~ 8000,
           is.na(bcell_counted) ~ NA_real_,
           TRUE ~ 9999999),
         bcell_denom = case_when(
           bcell_counted == 'RBC' & per == 1 ~ thin_denom,
           bcell_counted == 'WBC' & per == 2 ~ 200,
           bcell_counted == 'WBC' & per == 3 ~ 500,
           bcell_counted == 'WBC' & per == 4 ~ exper,
           is.na(bcell_counted) | is.na(per) ~ NA_real_,
           TRUE ~ 9999999),
         par_dens = round((apc * bcell_mult) / bcell_denom, digits = 0),
         log10_par_dens = if_else(par_dens == 0, 0, log10(par_dens))
         ) 

rio::export(prima, "prima_clean.xlsx")

# # Check
# select(.data = prima,
#        apc, per, exper,
#        bcell_counted, bcell_mult, bcell_denom,
#        par_dens, md) |> view()
         
# Pf-only data
prima_pf <- prima |> 
  
  # Exclude negative and Pan-only MRDT results
  # i.e., include patients with Pf
  filter(mrdt != 'Negative', mrdt != 'Pan (pLDH)') |> 
  
  # Outcome for regression
  mutate(double_band = if_else(mrdt == 'Pf/Pan (HRP-2/pLDH)', 1L, 0L) |>
           factor()) |>
  filter(spc != 'Negative') # Exclude false-positive MRDT

rio::export(prima_pf, "prima_pf.xlsx")

# Descriptives ------------------------------------------------------------

(means <- prima_pf |> 
   group_by(mrdt) |> 
   summarise(mean_log10 = mean(log10_par_dens),
             geo_mean = 10^(mean_log10),
             median_ori = median(par_dens),
             mean_ori = mean(par_dens)))

dens_histogram_ori <- prima_pf |> 
  ggplot(aes(x = par_dens / 1000)) +
  geom_histogram(aes(fill = mrdt),
                 alpha = 0.8, bins = 14, position = 'identity') +
  scale_x_continuous(labels = scales::comma,
                     breaks = seq(0, 2000, by = 500)) +
  scale_y_continuous(expand = c(0, 0, 0, 5),
                     breaks = seq(0, 120, by = 20),
                     limits = c(0, 121)) +
  scale_fill_manual(values = c('gray30', 'gray75')) +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = 5.5, family = "Fira Code")) +
  labs(x = '',
       y = '')

dens_histogram_ori2 <- prima_pf |> 
    filter(par_dens < 0.5e6) |> 
    ggplot(aes(x = par_dens / 1000)) +
    geom_histogram(aes(fill = mrdt),
                   alpha = 0.8, bins = 15, position = 'identity') +
    scale_x_continuous(labels = scales::comma,
                       breaks = seq(0, 500, by = 50)) +
    scale_y_continuous(expand = c(0, 0, 0, 5),
                       breaks = seq(0, 50, by = 10)) +
    scale_fill_manual(values = c('gray30', 'gray75')) +
    theme(legend.position = c(0.8, 0.2),
          panel.grid = element_blank(),
          axis.ticks = element_blank()) +
    labs(x = '\nParasite density (in thousand per \U00B5L)',
         y = 'Number of patients\n')

combined_histogram <- ggdraw() + draw_plot(dens_histogram_ori2)
(combined_histogram <- combined_histogram +
  draw_plot(dens_histogram_ori, x = 0.45, y = 0.45, width = 0.5, height = 0.5))

ggsave(plot = combined_histogram,
       filename = "combined_histogram.png",
       height = 5, width = 5.5,
       dpi = 1200)

(dens_histogram_log <- prima_pf |> 
    ggplot(aes(x = log10_par_dens)) +
    geom_histogram(aes(fill = mrdt),
                   alpha = 0.8, bins = 17, position = 'identity') +
    scale_x_continuous(labels = function(x) parse(text = sprintf("10^%d", x)),
                       breaks = seq(1, 6)) +
    scale_y_continuous(expand = c(0, 0, 0, 5),
                       breaks = seq(0, 50, by = 5)) +
    scale_fill_manual(values = c('gray30', 'gray75')) +
    theme(legend.position = c(0.2, 0.5),
          panel.grid = element_blank(),
          axis.ticks = element_blank()) +
    labs(x = '\nParasite density (per \U00B5L)',
         y = 'Number of patients\n')
    # + facet_wrap(~bcell_counted)
    )

ggsave(plot = dens_histogram_log,
       filename = "dens_histogram_log.png",
       height = 5, width = 5.5,
       dpi = 1200)

(dens_beeswarm <- prima_pf |> 
  ggplot(aes(x = mrdt, y = log10_par_dens)) +
  geom_beeswarm(aes(group = mrdt), alpha = 0.2, size = 1.5, cex = 1.8) +
  geom_boxplot(alpha = 0, width = 0.5/2) +
  # geom_point(data = means, shape = 18, size = 3, colour = '#D60B00',
  #            aes(x = mrdt, y = mean_log10)) +
  scale_y_continuous(limits = c(0.9, 6.5),
                     labels = function(x) parse(text = sprintf("10^%d", x)),
                     breaks = seq(1, 6)) +
  # scale_colour_manual(values = c('gray60', 'gray30')) +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(10, 50, 10, 50)) +
  labs(x = '\nMalaria RDT positive for Plasmodium falciparum',
       y = 'Parasite density (per \U00B5L)\n'))

ggsave(plot = dens_beeswarm,
       filename = "dens_beeswarm.png",
       height = 4.5, width = 5.5,
       dpi = 1200)

# Modelling ---------------------------------------------------------------

linear_model <- lm(formula = log10_par_dens ~ double_band + film,
                   data = prima_pf)
(out_linear <- tidy(linear_model, conf.int = TRUE))
10^out_linear[2, 2] # Interpretation: higher/lower by a factor of ...
ggpredict(linear_model, terms = c("double_band", "film")) |> plot()

logit_model <- glm(formula = double_band ~ log10_par_dens,
                   family = "binomial", data = prima_pf)
(out_logit <- tidy(logit_model, conf.int = TRUE, exponentiate = TRUE))

(marginal_plot <- ggpredict(logit_model,
                           terms = c("log10_par_dens [0:6, by = 0.01]")) |>
  plot() +
  theme_bw() +
  scale_x_continuous(labels = function(x) parse(text = sprintf("10^%d", x)),
                     breaks = seq(0, 6),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.025, 0),
                     breaks = seq(0, 1, by = 0.2),
                     labels = label_number(scale = 100, suffix = '%')) +
  theme(panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(10, 50, 10, 50),
        text = element_text(size = 9.5, family = "Fira Code"),
        axis.title.x = element_text(hjust = 1),   
        axis.title.y = element_text(hjust = 0.5)) +
  labs(y = 'Estimated probability\n',
       x = '\nParasite density (per \U00B5L)',
       title = ''))

ggsave(plot = marginal_plot,
       filename = "marginal_plot.png",
       height = 5.5, width = 7.5,
       dpi = 1200)

# Linearity check
prima_check <- select(prima_pf, c(log10_par_dens, film, double_band))
ddist <- datadist(prima_check)
options(datadist = 'ddist')

lrm_linear <- lrm(double_band ~ log10_par_dens, data = prima_pf)
lrm_rcs <- lrm(double_band ~ rcs(log10_par_dens, 3), data = prima_pf)

## Approximately linear
(linearity_plot <- Predict(lrm_rcs, name = 'log10_par_dens') |>
  ggplot() +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          plot.margin = margin(10, 50, 10, 50),
          text = element_text(size = 9.5, family = "Fira Code"),
          axis.title.x = element_text(hjust = 1),   
          axis.title.y = element_text(hjust = 0.5),
          panel.grid = element_blank()) +
    scale_x_continuous(labels = function(x) parse(text = sprintf("10^%d", x)),
                       breaks = seq(0, 6),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(-4, 12, by = 4)) +
    labs(y = 'Log-odds\n',
         x = '\nParasite density (per \U00B5L)',
         title = ''))

ggsave(plot = linearity_plot,
       filename = "linearity_plot.png",
       height = 5.5, width = 7.5,
       dpi = 1200)









