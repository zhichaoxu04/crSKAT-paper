# Real Data Analysis - Diagnosis Plots and Table

library(tidyverse)
library(ggplot2)
library(data.table)


# ---- Read in all results
filelist = list.files(path = "your/path/to/results", 
                      recursive = FALSE,
                      pattern = ".*.txt",
                      full.names = TRUE)
# assuming tab separated values with a header    
datalist <-  lapply(filelist, function(x) data.table::fread(x)) 
#assuming the same header/columns for all files
datafrcrSKAT <-  do.call("rbind", datalist)

# Merge all tests
FinalResult <- datafrcrSKAT %>% 
  drop_na(crskatp, crburdenp) %>% 
  filter(crskatp != 0 & crburdenp != 0) %>% 
  dplyr::select(-skatp, -burdenp, -SKATOp) %>% 
  left_join(datafrICSKAT %>% 
              drop_na(skatp, burdenp, SKATOp) %>% 
              filter(skatp != 0 & burdenp != 0 & SKATOp != 0) %>% 
              dplyr::select(gene, chr, start, skatp, burdenp, SKATOp),
            by = c("gene", "chr", "start")) %>% 
  drop_na(crskatp, crburdenp, skatp, burdenp, SKATOp)

# ---- Q-Q Plots For ICSKAT & crSKAT
library(lattice)
QQplotDF <- FinalResult %>% 
  mutate(ex_crs = rank(crskatp)/(n()+1), ex_crb = rank(crburdenp)/(n()+1),
         ex_ics = rank(skatp)/(n()+1), ex_icb = rank(burdenp)/(n()+1), ex_icso = rank(SKATOp)/(n()+1)) %>% 
  rename(crSKAT=crskatp, crBurden=crburdenp, ICSKAT=skatp, ICBurden=burdenp, ICSKATO=SKATOp) %>% 
  pivot_longer(cols = c("crSKAT", "crBurden", "ICSKAT", "ICBurden", "ICSKATO"),
               names_to = "Tests") %>% 
  mutate(Expected = case_when(Tests == "crSKAT" ~ ex_crs,
                              Tests == "crBurden" ~ ex_crb,
                              Tests == "ICSKAT" ~ ex_ics,
                              Tests == "ICBurden" ~ ex_icb,
                              Tests == "ICSKATO" ~ ex_icso,
                              TRUE ~ as.numeric(NA)),
         Tpvalue = -log10(value),
         TExpected = -log10(Expected))

QQplotDF %>% ggplot(aes(x = TExpected, y = Tpvalue, color = Tests)) +
  geom_point(alpha = 0.9) + 
  geom_abline(intercept = 0, slope = 1)+
  labs(x = "Expected", y = "Observed", colour = "Test")+ 
  theme(axis.line = element_line(linetype = "solid", linewidth = 1), 
        panel.grid.major = element_line(colour = "azure4", linetype = "blank"), 
        panel.grid.minor = element_line(colour = NA, linetype = "blank"), 
        axis.text = element_text(size = 13, colour = "gray0"), 
        panel.background = element_rect(fill = NA)) +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))  

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# --- Manhattan Plot
don <- FinalResult %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len = max(end)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot = lag(cumsum(as.numeric(chr_len)), default = 0)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(FinalResult, ., by = c("chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, end) %>%
  mutate(BPcum = end + tot)

axisdf <- don %>% group_by(chr) %>% summarize(center = (max(BPcum) + min(BPcum) ) / 2 )

ManhCRICSKAT <- ggplot(don, aes(x=BPcum, y=-log10(crskatp))) +
  # Show all points
  geom_point(aes(color = as.factor(chr)), alpha = 0.7, size = 1.3) +
  scale_color_manual(values = rep(c("#999999", "#006400"), 22 )) + 
  # custom X axis:
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), breaks = c(1:10) ) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())+
  labs(x = "Chromosome", y = "-log10(p)")  + 
  theme(axis.line = element_line(linetype = "solid", linewidth = 1),
        panel.background = element_rect(fill = NA)) + 
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank")) + 
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 


# ---- Table for real data results
Table4 <- FinalResult %>% 
  rename(CRICSKAT=crskatp, CRBurden=crburdenp, ICSKAT=skatp, Burden=burdenp, ICSKATO=SKATOp) %>% 
  pivot_longer(cols = c("CRICSKAT", "CRBurden", "ICSKAT", "Burden", "ICSKATO"),
               names_to = "Tests") %>% 
  slice_min(value, n = 50) %>% 
  pivot_wider(names_from = Tests, values_from = "value") %>% 
  left_join(FinalResult, by = c("gene", "chr")) %>% 
  dplyr::select(gene, chr, crskatp, crburdenp, skatp, burdenp, SKATOp) %>% 
  rename(Gene = gene, Chr=chr, CRICSKAT=crskatp, CRBurden=crburdenp, ICSKAT=skatp, ICBurden=burdenp, ISCKATO=SKATOp)



