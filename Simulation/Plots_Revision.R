
# Power Plots
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggpubr)
library(data.table)

createPowerPlot <- function(combined_data, Scenario, Eff, MAF) {
  
  Pwr_plot_single <- combined_data %>% 
    dplyr::filter(Scen==Scenario, EffectSize==Eff, maxMAF==MAF) %>% 
    pivot_longer(cols = c(crSKAT_pwr, crBurden_pwr, ICSKAT_pwr, ICBurden_pwr, WVIC_pwr), names_to = "Test", values_to = "Power") %>% 
    mutate(Test = factor(Test, levels = c("crSKAT_pwr", "crBurden_pwr","ICSKAT_pwr", "ICBurden_pwr", "WVIC_pwr"), 
                         labels = c("crSKAT", "crBurden","ICSKAT", "ICBurden", "WVIC"))) %>% 
    rename("CausalSNPs" = "Q")
  
  ggplot(Pwr_plot_single, aes(x = CausalSNPs, y = Power, group = Test)) +
    scale_x_continuous(breaks = seq(2, 14, 1)) +
   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
    scale_color_brewer(palette = "Set1") +  # Custom color palette
    theme_minimal(base_size = 15) +  # Using theme_minimal
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = "bold"),
          strip.text = element_text(size = 15),
          axis.line = element_line(linetype = "solid", color = "gray0"),
          panel.grid.major = element_blank(),  # Removing major grid lines
          panel.grid.minor = element_blank(),  # Removing minor grid lines
          panel.background = element_rect(fill = NA, color = NA),
          axis.ticks = element_line(color = "gray0", linewidth = 0.5),
          axis.ticks.length = unit(0.2, "cm")) +  # Adjusting legend position
    labs(x = "Number of Causal SNPs", y = "Power") +
    stat_smooth(aes(color = Test, linetype = Test), method = "lm",
                formula = y ~ poly(x, 3), se = FALSE) + 
    expand_limits(y = c(0, 1))
}

createPowerPlot2 <- function(combined_data, Scenario, Eff, MAF) {
  
  Pwr_plot_single <- combined_data %>% 
    dplyr::filter(Scen==Scenario, EffectSize==Eff, maxMAF==MAF) %>% 
    pivot_longer(cols = c(crSKAT_pwr, crBurden_pwr, ICSKAT_pwr, ICBurden_pwr, WVIC_pwr), names_to = "Test", values_to = "Power") %>% 
    mutate(Test = factor(Test, levels = c("crSKAT_pwr", "crBurden_pwr","ICSKAT_pwr", "ICBurden_pwr", "WVIC_pwr"), 
                         labels = c("crSKAT", "crBurden","ICSKAT", "ICBurden", "WVIC"))) %>% 
    rename("CausalSNPs" = "Q")
  
  ggplot(Pwr_plot_single, aes(x = CausalSNPs, y = Power, group = Test)) +
    scale_x_continuous(breaks = seq(2, 11, 1), limits = c(2, 11)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
    scale_color_brewer(palette = "Set1") +  # Custom color palette
    theme_minimal(base_size = 15) +  # Using theme_minimal
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = "bold"),
          strip.text = element_text(size = 15),
          axis.line = element_line(linetype = "solid", color = "gray0"),
          panel.grid.major = element_blank(),  # Removing major grid lines
          panel.grid.minor = element_blank(),  # Removing minor grid lines
          panel.background = element_rect(fill = NA, color = NA),
          axis.ticks = element_line(color = "gray0", linewidth = 0.5),
          axis.ticks.length = unit(0.2, "cm")) +  # Adjusting legend position
    labs(x = "Number of Causal SNPs", y = "Power") +
    stat_smooth(aes(color = Test, linetype = Test), method = "lm",
                formula = y ~ poly(x, 3), se = FALSE) + 
    expand_limits(y = c(0, 1))
}

 
setwd("/YourPathToResults/")
filelist = list.files(path = getwd(), 
                      pattern = "P2_40_.*\\SimG.txt", 
                      full.names = TRUE)


# Read each file into a list of data frames
data_list <- lapply(filelist, read.table, header = TRUE) # Adjust read.table parameters as necessary
# Combine all data frames into one
combined_data <- do.call(rbind, data_list)

Pwr_plot_ALL <- combined_data %>% 
  dplyr::filter(Scen==1) %>% 
  pivot_longer(cols = c(crSKAT_pwr, crBurden_pwr, ICSKAT_pwr, ICBurden_pwr, WVIC_pwr), names_to = "Test", values_to = "Power") %>% 
  mutate(Test = factor(Test, levels = c("crSKAT_pwr", "crBurden_pwr","ICSKAT_pwr", "ICBurden_pwr", "WVIC_pwr"), 
                       labels = c("crSKAT", "crBurden","ICSKAT", "ICBurden", "WVIC"))) %>% 
  rename("CausalSNPs" = "Q")

Pwr_plot_ALL %>% 
  ggplot(aes(x = CausalSNPs, y = Power, group = Test)) + 
  scale_x_continuous(breaks = seq(2, 15, 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))  + 
  theme(strip.text = element_text(size = 12, face = "bold")) + 
  facet_wrap(~ EffectSize + maxMAF, ncol = 6)+
  # facet_wrap(~ EffectSize, ncol = 6)+
  # geom_line(aes(colour = Test, linetype = Test), linewidth = 1) +
  stat_smooth(aes(x = CausalSNPs, y = Power, colour = Test, linetype = Test), method = "lm",
              formula = y ~ poly(x,3), se = FALSE)+
  expand_limits(y = c(0, 1))


# ------- Combine all figures
# Example usage with your data frame (Pwr_plot)
# -------------- MAF = 0.03
Figure2 <- createPowerPlot(combined_data_50000, Scenario=1, Eff=0.02, MAF=0.03)
Figure4 <- createPowerPlot(combined_data_50000, Scenario=2, Eff=0.1, MAF=0.03)
Figure6 <- createPowerPlot(combined_data_50000, Scenario=3, Eff=0.05, MAF=0.03)

Figure1 <- createPowerPlot(combined_data_2500, Scenario=1, Eff=0.08, MAF=0.03)
Figure3 <- createPowerPlot(combined_data_2500, Scenario=2, Eff=2.5, MAF=0.03)
Figure5 <- createPowerPlot(combined_data_2500, Scenario=3, Eff=0.5, MAF=0.03)

library(ggpubr)
ggpubr::ggarrange(Figure1, Figure2, Figure3, Figure4, Figure5, Figure6, 
                  ncol = 2, nrow = 3, labels = c("a", "b", "c", "d", "e", "f"), 
                  common.legend = T, legend = "right", font.label = list(size=20))
ggpubr::ggarrange(Figure3, Figure4, Figure5, Figure6, 
                  ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), 
                  common.legend = T, legend = "right", font.label = list(size=20))

# -------------- MAF = 0.05
Figure4 <- createPowerPlot(combined_data_50000, Scenario=1, Eff=0.02, MAF=0.05)
Figure5 <- createPowerPlot(combined_data_50000, Scenario=2, Eff=0.075, MAF=0.05)
Figure6 <- createPowerPlot(combined_data_50000, Scenario=3, Eff=0.04, MAF=0.05)

Figure1 <- createPowerPlot(combined_data_2500, Scenario=1, Eff=0.08, MAF=0.05)
Figure2 <- createPowerPlot(combined_data_2500, Scenario=2, Eff=1, MAF=0.05)
Figure3 <- createPowerPlot(combined_data_2500, Scenario=3, Eff=0.8, MAF=0.05)
ggpubr::ggarrange(Figure1, Figure2, Figure3, Figure4, Figure5, Figure6, 
                  ncol = 3, nrow = 2, labels = c("a", "b", "c", "d", "e", "f"), 
                  common.legend = T, legend = "bottom", font.label = list(size=20))

# -------------- MAF = 0.15
Figure4 <- createPowerPlot(combined_data_50000, Scenario=1, Eff=0.01, MAF=0.15)
Figure5 <- createPowerPlot(combined_data_50000, Scenario=2, Eff=0.075, MAF=0.15)
Figure6 <- createPowerPlot(combined_data_50000, Scenario=3, Eff=0.03, MAF=0.15)

Figure1 <- createPowerPlot(combined_data_2500, Scenario=1, Eff=0.05, MAF=0.15)
Figure2 <- createPowerPlot(combined_data_2500, Scenario=2, Eff=1.5, MAF=0.15)
Figure3 <- createPowerPlot(combined_data_2500, Scenario=3, Eff=0.4, MAF=0.15)
ggpubr::ggarrange(Figure1, Figure2, Figure3, Figure4, Figure5, Figure6, 
                  ncol = 3, nrow = 2, labels = c("a", "b", "c", "d", "e", "f"), 
                  common.legend = T, legend = "bottom", font.label = list(size=20))

# ----------------- Real Genotype --------------
# -------------- MAF = 0.15
Figure1 <- createPowerPlot2(combined_data_1000_Sim, Scenario=1, Eff=0.07, MAF=0.05)
Figure2 <- createPowerPlot2(combined_data_1000_Sim, Scenario=2, Eff=1, MAF=0.05)
Figure3 <- createPowerPlot2(combined_data_1000_Sim, Scenario=3, Eff=0.6, MAF=0.05)
Figure4 <- createPowerPlot2(combined_data_1000_Real, Scenario=1, Eff=0.1, MAF=0.15)
Figure5 <- createPowerPlot2(combined_data_1000_Real, Scenario=2, Eff=2.5, MAF=0.15)
Figure6 <- createPowerPlot2(combined_data_1000_Real, Scenario=3, Eff=0.4, MAF=0.15)
ggpubr::ggarrange(Figure1, Figure2, Figure3, Figure4, Figure5, Figure6, 
                  ncol = 3, nrow = 2, labels = c("a", "b", "c", "d", "e", "f"), 
                  common.legend = T, legend = "bottom", font.label = list(size=20))

# ----------------- WVIC best --------------
Figure1 <- createPowerPlot(combined_data_1000_Sim, Scenario=1, Eff=0.07, MAF=0.15)
Figure2 <- createPowerPlot(combined_data_1000_Sim, Scenario=2, Eff=2.5, MAF=0.15)
Figure3 <- createPowerPlot(combined_data_1000_Sim, Scenario=3, Eff=0.8, MAF=0.15)

ggpubr::ggarrange(Figure1, Figure2, Figure3, 
                  ncol = 3, nrow = 1, labels = c("a", "b", "c"), 
                  common.legend = T, legend = "bottom", font.label = list(size=20))






























