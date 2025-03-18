######################################################################
### The script is used to investigate the impact of dosing schedules  ----
# y.guo@lacdr.leidenuniv.nl - March 2024
######################################################################
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
source("functions/get_D_coeff.R") ## get the diffusion coefficient
source("functions/pde.R") ## PDE function (elimination)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Group the impact factors   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No.1 dosing schedules
ii <- c(8, 24, 8) # unit: hour
dur <- c(0.5, 0.5, 8) # unit: hour
amt <- c(100, 300, 100)
dosing_name <- c(
  "intermittent infusion\n dose interval = 8 h",
  "intermittent infusion\n dose interval = 24 h",
  "continuous infusion"
)
dosing <- cbind(ii, dur, amt, dosing_name)

# No.2 molecular size
radius <- c(1, 50) # unit: nm
radius_lb <- c("radius = 1 nm","radius = 50 nm")

# merge the parameters
dosing <- cbind(rbind(dosing, dosing),
                rep(radius, each = 3),
                rep(radius_lb, each = 3)) %>% as.data.frame()
colnames(dosing) <- c("interval", "duration", "amount", "description","radius", "feature")

for(j in c(1,2,3,5)){
  dosing[,j] <- as.numeric(dosing[,j])
}

dosing <- dosing %>% 
  mutate(sum_concf = 0,
         sum_concb = 0)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - Plasma conc profiles (dosing schedules)   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Plasma Pk model for the hypothetical model
mod_H <- RxODE({
  CL ~ 0.05775
  # L/min, corresponding to half-life 3h
  V1 ~ 15
  # L , corresponding to small molecule
  
  d / dt(C1) ~ -(CL / V1) * C1
  # Infusion
  
  Con1 = C1 / V1
  
})

dat.pk.plasma_clct <- list()

for (i in 1:nrow(dosing)) {
  
    # assign parameter
    ii <- dosing[i, 1]
    dur <- dosing[i, 2]
    amt <- dosing[i, 3]
  
    # create event table
    event.pk <- et() %>%
      # dosing events
      et(
        id = 1,
        amt = amt,
        cmt = 1,
        dur = dur * 60,
        ii = ii * 60,
        # ii= ii*60 for small molecule, ii=0 for antibody
        addl = 72 / ii - 1
      ) %>% # addl = 72/ii - 1 for small molecule, addl =0 for antibody ; #previously observe to = 72*60, now change it to 24hour
      # sampling events
      et(seq(from = 0, to = 72 * 60, by = 1)) %>% #previously observe to = 72*60, now change it to 24hour
      as_tibble()
  
  ## Simulate!
  set.seed(2021)
  
  dat.pk.plasma <-
    rxSolve(mod_H, event = event.pk, addDosing = F) %>%
    mutate(description = dosing[i, 4],
           feature = dosing[i, 6]) %>%
    distinct()
  
  dat.pk.plasma_clct[[i]] <- dat.pk.plasma
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - Mucus conc profiles   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate a blank list
listdata <- list()
listdata_CV <- list()
# test
i = 1

# run in loop
for (i in 1:nrow(dosing)) {
  
  start_time <- Sys.time()
  
  # assign parameter
  dat.pk.plasma_ok <- dat.pk.plasma_clct[[i]]
  r_ok = dosing[i, "radius"]
  
  # run pde model    
  dat.pk.mucus <- pde(
    db.pk = dat.pk.plasma_ok$Con1,
    # plasma PK
    coef.diff = Diffusion_coeff(r_ok),
    # diffusion coefficient um^2/sec,113 for 1nm, 15.13 for 5nm
    kon = 100 ,
    # binding rate /(M*sec)
    koff = 0.1,
    # unbinding rate /sec
    Tmucin = 2e-5,
    # concentration of total mucin mol/L
    dt = 0.5,
    # step size time in seconds
    span.t = 72,
    # observation time span in hours
    dx = 25,
    # step size depth in um
    span.x = 500,
    # observation depth span in um
    clr = 0.00001 
    # elimination rate in 1/sec
  ) 
  
  dat.pk.mucus.ok <- dat.pk.mucus %>% 
    rename(Conc_F = CONC_FREE) %>% 
    rename(Conc_B = CONC_BIND) %>%
    left_join(dat.pk.plasma, by = c("TIME" = "time")) %>% 
    mutate(feature = dosing[i, "feature"]) %>%  # outcome
    mutate(description = dosing[i, "description"]) %>%  # outcome
    mutate(Ratio = Conc_F/Con1) %>%  # outcome
    mutate(ConcTT = Conc_F+Conc_B)  %>%  # outcome
    mutate(BRatio = Conc_B/ConcTT)  %>%  # outcome
    mutate(tRatio = Conc_B*0.01)
  
  listdata[[i]] <- dat.pk.mucus.ok
  
  # to calculate the cv
  input_data <- dat.pk.mucus.ok[-1,] %>%
    group_by(TIME) %>%
    summarise(CV_F = 100*sd(Conc_F,na.rm = TRUE)/mean(Conc_F,na.rm = TRUE),
              CV_B = 100*sd(Conc_B,na.rm = TRUE)/mean(Conc_B,na.rm = TRUE)) %>% 
    mutate(description = dosing[i, 4],
           feature = dosing[i, 6])  
  
  listdata_CV[[i]] <- input_data
  
  
  dosing[i, c("sum_concf","sum_concb")] <-
    c(sum(dat.pk.mucus.ok$Conc_F),
      sum(dat.pk.mucus.ok$Conc_B))
  
  end_time <- Sys.time()
  
  cat("time of simulation round: ", end_time - start_time, i, "\n")
  i = i + 1
} 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 5 - plot for the spatial distribution ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_data <- do.call(rbind, listdata)
# organize the order of different values for each factor
all_data$description <- factor(all_data$description, levels = unique(all_data$description))

# plot
pCon_f <- ggplot(all_data, aes(TIME/60, DEPTH, fill = Conc_F)) +
  geom_tile() +
  guides(fill = guide_colorbar(
    title.position = "top",
    direction = "vertical")) +
  scale_fill_distiller("Conc_F \n(mg/L)", palette = 4,
                       direction = 1) +
  scale_x_continuous(name = "Time (hour)", breaks = c(0:100)*10) + 
  scale_y_continuous(name = "Distance (μm)") +
  labs(fill = "Conc_F\n (mg/L)") +
  theme_bw() +
  theme(text = element_text(family = "sans", color = "black", size = 9),
        axis.text = element_text(family = "sans", color = "black"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7.5),
        legend.key.width = unit(0.005, "npc"),
        legend.key.height = unit(0.3, 'cm'),
        strip.text = element_text(color = "black", size = 7, face = "plain"))+
  facet_grid(feature~description)
pCon_f

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 6 - plot for the summary exposure ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dosing <- dosing %>%
  mutate(AUCf = round(sum_concf/(60 * 21))) %>% # during calculation: 60 steps in one hour; 20 steps in space
  mutate(AUCb = round(sum_concb/(60 * 21)))     # the unit is mg*hour/L

dosing$description <- factor(dosing$description, 
                             levels = unique(dosing$description))

exposure <- ggplot(dosing, aes(x = description,y = AUCf,fill = feature)) +
  geom_bar(stat = "identity",position = "dodge") +
  scale_x_discrete() +
  scale_y_continuous(name = "AUCf (mg*h/L)",limits = c(0,300)) +
  scale_fill_manual(values=c("#7e97ad","#a6cba9"), labels = c("1 nm", "50 nm"))+ 
  labs(fill = "Radius") +
  geom_text(aes(label = AUCf), hjust=0.5, vjust=-0.5, position = position_dodge(.9),size = 2.75)+ # make the lable right above
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        axis.text.x = element_text(family = "sans", color = "black",size = 6.7),
        axis.text.y = element_text(family = "sans", color = "black",size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.width = unit(0.02, "npc"),
        legend.key.height = unit(0.3, "cm"))
exposure 

figure4 <- ggarrange(pCon_f,exposure,
                   labels = c("A","B"),
                   nrow = 1, ncol = 2, widths = c(1.15, 1), heights=1,
                   font.label = list(size = 12))
figure4 

# save the plot
ggsave("figure/figure4.pdf", plot = figure4, width = 22, height = 7, units = "cm", device = cairo_pdf) # stick with this version


# save the data
save(file = "result/f4_impactDosing.Rdata", listdata)
save(file = "result/f4_dosing_sum.Rdata", dosing)
save(file = "result/f4_impactDosing_cv.Rdata", listdata_CV)
