######################################################################
### The script is used to generate delayed profiles as the boundary 
## condition for diffusion ----
# y.guo@lacdr.leidenuniv.nl - July 2025
######################################################################

rm(list = ls(all = TRUE))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(rxode2)
library(patchwork)
library(RColorBrewer) # color-blind people friendly

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Source files   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function
source("functions/get_D_coeff.R") ## get the diffusion coefficient
source("functions/pde.R") ## PDE function 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Plasma conc profiles (No.1 plasma half-lives)   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clearance <- c(0.05775,  # L/min, corresponding to half-life 3h
               0.0289,   # L/min, corresponding to half-life 6h
               0.01444,  # L/min, corresponding to half-life 12h
               0.00722)  # L/min, corresponding to half-life 24h

mod_lung_delay <- RxODE({
  
  d/dt(blood) ~ - CL* blood/Vbl
  d/dt(C1) ~  kt * blood/Vbl - kt * C1
  d/dt(C2) ~ kt * C1 - kt * C2
  d/dt(C3) ~ kt * C2 - kt * C3
  
  ####
  C_blood = blood/Vbl
  Con1 = C3
  
})

kt_values <- c(1, 0.02, 0.002) # 1/60, 0.5/60
distribution <- c("rapid distribution", "intermediate distribution", "slow distribution")

dat.pk.plasma_clct_delay <- list()
dat.pk.plasma_clct <- list()

for (m in 1: length(kt_values)) {
  
  kt_value = kt_values[m]
  
  for (j in 1:4) {
    
    cl_ok = clearance[j]
    # parameters for mod_lung_delay
    
    par <- c(CL = cl_ok,
             Vbl = 15,
             kt = kt_value) # ktr = (n + 1)/MTT; MTT = 8 h corrsponds to 0.5 h-1, /60 to min-1
    
    ## Create event table; in the time unit of min 
    amt <- 100 # mg
    dur <- 0.5 # hour
    ii <- 8 # hour
    event.pk <- et() %>%
      # dosing events
      et(id = 1,
         amt = amt,
         cmt = 1,
         dur = dur*60,
         ii = ii*60,  
         addl = 24*3/ii - 1
      ) %>% 
      # sampling events
      et(seq(from = 0, to = 24*3*60, by = 1)) %>% 
      as_tibble()
    
    ## Simulate
    set.seed(2021)
    dat.pk.plasma <-
      rxSolve(mod_lung_delay, 
              par, event = event.pk, addDosing = F) %>%
      distinct() %>%
      as.data.frame() %>%
      mutate(distribution = distribution[m],
             halflife = halflife[j])
    
    dat.pk.plasma_clct_delay[[j]] <- dat.pk.plasma
    
    j = j + 1
  }
  
  dat.pk.plasma_clct[[m]] <- dat.pk.plasma_clct_delay
  
  m = m + 1
}

df_long <- bind_rows(dat.pk.plasma_clct) %>%
  pivot_longer(cols = c(C_blood, Con1), 
               names_to = "compartment",
               values_to = "concentration") %>%
  mutate(compartment = recode(compartment,
                              C_blood = "Plasma",
                              Con1 = "ELF"),
         compartment = factor(compartment, levels = c("Plasma", "ELF")))

cb_palette2 <- brewer.pal(n = 3, name = "Dark2")
unique(df_long$distribution)

df_long$distribution <-  factor(df_long$distribution, 
                                levels = unique(df_long$distribution))
figure_PBPK <- 
  ggplot(df_long, 
         aes(x = time / 60, y = concentration, color = factor(distribution), linetype = compartment)) +
  geom_line(linewidth = 0.7, alpha = 0.5) +
  scale_color_manual(values = cb_palette2) +
  scale_linetype_manual(values = c("Plasma" = "solid","ELF" = "dashed"),name = "Compartment") +
  labs(x = "Time (hour)", y = "Concentration (mg/L)", color = "Distribution") +
  facet_grid(cols = vars(halflife)) +
  theme_bw() +
  theme(text = element_text(family = "sans", color = "black", size = 10),
        axis.text = element_text(family = "sans", color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.width = unit(1, "cm"),
        strip.text = element_text(color = "black", size =8, face = "plain"))
figure_PBPK
