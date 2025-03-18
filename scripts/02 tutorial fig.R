############################################################
### The script is used as tutorial figure of                ----
### the spatial PK model with small molecule as an example  ----
# y.guo@lacdr.leidenuniv.nl - March 2024
############################################################
rm(list = ls(all = TRUE))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(rxode2)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Source files   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function
source("functions/pde.R") ## PDE function (elimination)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Plasma conc profile   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plasma Pk model for the hypothetical drug
mod_H <- RxODE({
  
  CL ~ 0.05775; # L/min, corresponding to half-life 3h
  V1 ~ 15;# L , corresponding to small molecule
  
  d/dt(C1) ~ -(CL/V1)*C1;# Infusion
  
  Con1 = C1/V1;
})

## Dosing
amt <- 100 # mg
dur <- 0.5 # hour
ii <- 8 # hour

## Create event table
event.pk <- et() %>%
  # dosing events
  et(id = 1,
     amt = amt,
     cmt = 1,
     dur = dur*60,
     ii = ii*60,  # ii= ii*60 for small molecule, ii=0 for antibody
     addl = 24/ii - 1) %>% 
  # sampling events
  et(seq(from = 0, to = 24*60, by = 1)) %>% 
  as_tibble()

## Simulate!
set.seed(2021)

dat.pk.plasma <- rxSolve(mod_H, event = event.pk, addDosing = F) %>% 
  distinct()

# draw the plasma PK
dat.pk.plasma <- dat.pk.plasma %>% 
  mutate(loc = "plasma")

PlasmaPK <- dat.pk.plasma %>% 
  ggplot() +
  geom_line(aes(time/60, Con1)) +
  scale_x_continuous(name = "Time (hour)", breaks = c(0:100)*4) +
  scale_y_continuous(name = "Concentration \n(mg/L)") +
  #labs(fill = "Plasma concentration") +
  theme(text = element_text(family = "sans", color = "black", size = 12),
        axis.text = element_text(family = "sans", color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.width = unit(0.01, "npc"),
        strip.text = element_text(color = "black", size = 10, face = "plain"))+
  facet_grid(cols = vars(loc)) +
  theme_bw()
PlasmaPK

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - Mucus conc profile   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat.pk.mucus <- pde(db.pk = dat.pk.plasma$Con1, # plasma PK
                    coef.diff = 114 , # diffusion coefficient um^2/sec
                    kon = 100 , # binding rate /(M*sec)
                    koff = 0.1, # unbinding rate /sec
                    Tmucin = 2e-5,# concentration of total mucin, mol/L
                    dt = 0.5, # step size time in seconds
                    span.t = 24, # observation time span in hours 
                    dx = 25, # step size depth in um
                    span.x = 500, # observation depth span in um
                    clr = 0.00001) # elimination rate in 1/s

dat.pk.mucus.tp <- dat.pk.mucus %>% 
  rename(Conc_F = CONC_FREE) %>% 
  rename(Conc_B = CONC_BIND) %>%
  left_join(dat.pk.plasma, by = c("TIME" = "time")) %>% 
  mutate(Ratio = Conc_F/Con1) %>%
  mutate(ConcTT = Conc_F+Conc_B)  %>% 
  mutate(BRatio = Conc_B/ConcTT)  %>% 
  mutate(tRatio = Conc_B*0.01)%>%  # unit change: Conc_B(mol/L)=Conc_B(mg/L)*0.001/500; Ratio: Conc_B(mol/L)/(2e-4)=0.01
  mutate(loc = "mucus")

pCon_f <- ggplot(dat.pk.mucus.tp, aes(TIME/60, DEPTH, fill = Conc_F)) +
  geom_tile() +
  guides(fill = guide_colorbar(
    title.position = "top",
    direction = "vertical")) +
  scale_fill_distiller("Conc_F \n(mg/L)", palette = 4,
                       direction = 1) +
  scale_x_continuous(name = "Time (hour)", breaks = c(0:100)*10) + 
  scale_y_continuous(name = "Distance (μm)") +
  theme_bw() +
  theme(text = element_text(family = "sans", color = "black", size = 12),
        axis.text = element_text(family = "sans", color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.width = unit(0.01, "npc"),
        strip.text = element_text(color = "black", size = 10, face = "plain"))+
  facet_grid(cols = vars(loc)) 
pCon_f


pRatio <- ggplot(dat.pk.mucus.tp, aes(TIME/60, DEPTH, fill = Ratio)) +
  geom_tile() +
  guides(fill = guide_colorbar(
    title.position = "top",
    direction = "vertical")) +
  scale_fill_distiller(palette = 4,
                       direction = 1) +
  scale_x_continuous(name = "Time (hour)", breaks = c(0:100)*10) + 
  scale_y_continuous(name = "Distance (μm)") +
  theme_bw() +
  theme(text = element_text(family = "sans", color = "black", size = 12),
        axis.text = element_text(family = "sans", color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.width = unit(0.01, "npc"),
        strip.text = element_text(color = "black", size = 10, face = "plain"))+
  facet_grid(cols = vars(loc)) 
pRatio

p_blank <- ggplot() +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

figure2_1 <- ggarrange(PlasmaPK,p_blank, 
                     labels = c("A",""),
                     nrow = 1, ncol = 2, 
                     widths = c(7.5, 1.1))

figure2_2 <- ggarrange(pCon_f,pRatio,
                     labels = c("B","C"),
                     nrow = 2, ncol = 1)

figure2 <- ggarrange(figure2_1, figure2_2,
                     # labels = c("A","", "B","C"),
                     nrow = 2, ncol = 1, 
                     heights=c(7.5, 16.8))
figure2
ggsave("figure/figure2.pdf", plot = figure2, width = 14, height = 14, units = "cm", device = cairo_pdf) 
ggsave("figure/figure2.tiff", plot = figure2, width = 14, height = 14, units = "cm",dpi=600)
save(file = "result/f2_spatialPK.Rdata", dat.pk.mucus.tp)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - calculation of the coefficient of variance   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# coefficient of variance, CV
# CV = (σ / μ) * 100%

# load the data
str(dat.pk.mucus.tp)
input_data <- dat.pk.mucus.tp[-1,]

# calculate the CV for each time
# cv_x <- (sd_x / mean_x) * 100

input_data <- input_data %>%
  group_by(TIME) %>%
  summarise(CV = 100*sd(Conc_F,na.rm = TRUE)/mean(Conc_F,na.rm = TRUE))

# plot the CV over time
CV_plot <- input_data %>% 
  mutate(loc = "mucus") %>% 
  ggplot() +
  geom_line(aes(TIME/60, CV)) +
  scale_x_continuous(name = "Time (hour)", breaks = c(0:100)*4) +
  scale_y_continuous(name = "Coefficient \nof variation (%)") +
  #labs(fill = "Plasma concentration") +
  theme_bw() +
  theme(text = element_text(family = "sans", color = "black", size = 12),
        axis.text = element_text(family = "sans", color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.width = unit(0.01, "npc"),
        strip.text = element_text(color = "black", size = 10, face = "plain"))+
  annotate('text', x = 3.3, y = 133, label = '145%, 0.02 hour',size = 8/.pt) +
  annotate('text', x = 11.5, y =33, label = '26.7%, 8.18 hour',size = 8/.pt) +
  annotate('text', x = 19.8, y =33, label = '24.7%, 16.20 hour',size = 8/.pt) +
  facet_grid(cols = vars(loc)) 


CV_plot
ggsave("figure/figure2_cv.pdf", plot = CV_plot, width = 9, height = 7, units = "cm") 
ggsave("figure/figure2_cv.tiff", plot = CV_plot, width = 9, height = 7, units = "cm",dpi=600)

# numerical analysis
# 145, 1st peak, 0.02 hour
# 26.7, 2nd peak, 8.18 hour
# 24.7, 3rd peak, 16.20 hour

p_blank <- ggplot() +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

figure2_1 <- ggarrange(PlasmaPK,p_blank, 
                       labels = c("A",""),
                       nrow = 1, ncol = 2, 
                       widths = c(7.5, 1.1))

figure2_4 <- ggarrange(CV_plot,p_blank, 
                       labels = c("D",""),
                       nrow = 1, ncol = 2, 
                       widths = c(7.5, 0.8))

figure2_4figs <- ggarrange(figure2_1, figure2_2,figure2_4,
                     # labels = c("A","", "B","C"),
                     nrow = 3, ncol = 1, 
                     heights=c(7.5, 16.8, 7.5))
figure2_4figs
ggsave("figure/figure2_4figs.pdf", plot = figure2_4figs, width = 14, height = 18, units = "cm", device = cairo_pdf) 
ggsave("figure/figure2_4figs.tiff", plot = figure2_4figs, width = 14, height = 18, units = "cm",dpi=600)
       