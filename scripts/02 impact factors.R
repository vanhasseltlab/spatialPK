######################################################################
### The script is used to investigate the impact of relevant factors  ----
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Source files   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function
source("functions/get_D_coeff.R") ## get the diffusion coefficient
source("functions/pde.R") ## PDE function (elimination)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Group the impact factors   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No.1 plasma halflives 
halflife <- c(3, 6, 12, 24) # unit: hour

# No.2 molecular size 
radius <- c(1, 5, 20, 50) # unit: nm
# https://www.fidabio.com/molecular-weight-to-size-protein-radius-calculator
# weight <- c(500, 153000, 6020000, 68058000)
weight <- c(500, 150000, 1000000, 10000000)

# No.3 binding affinity 
ka <- c(100, 1000, 10000, 100000) # unit: 1/M
kon <- c(10, 100, 1000, 10000)
koff = 0.1

# No.4 mucin concentrations 
Mc <- c(5e-6, 1e-5, 5e-5, 1e-4) # unit: M

# No.5 muco-cilliary clearance 
rate <- c(0, 1e-6, 1e-5, 1e-4) # unit: 1/s, 5e-06, 1e-05, 5e-05 update: 1e-06 1e-05 1e-041e-05
rateh <- rate * 3600 # unit: 1/h

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - Plasma conc profiles (No.1 plasma half-lives)   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clearance <- c(0.05775,  # L/min, corresponding to half-life 3h
               0.0289,   # L/min, corresponding to half-life 6h
               0.01444,  # L/min, corresponding to half-life 12h
               0.00722)  # L/min, corresponding to half-life 24h
## Create function
mod_plasma <- RxODE({
  d / dt(C1) ~ -(CL / V1) * C1
  Con1 = C1 / V1
  
})


dat.pk.plasma_clct <- list()

for (i in 1:length(clearance)) {
  cl = clearance[i]
  par <- c(CL = cl,
           V1 = 15)
  
  ## Create event table
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
       addl = 72/ii - 1) %>% 
    # sampling events
    et(seq(from = 0, to = 72*60, by = 1)) %>% 
    as_tibble()
  
  ## Simulate
  set.seed(2021)
  dat.pk.plasma <-
    rxSolve(mod_plasma, par, event = event.pk, addDosing = F) %>%
    distinct() %>%
    as.data.frame() %>%
    mutate(halflife = clearance[i])
  
  dat.pk.plasma_clct[[i]] <- dat.pk.plasma
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - Parameter list  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set the reference drug
col_names <-c("halflife",
              "radius",
              "kon",
              "Mc",
              "rate",
              "feature",
              "sum_concf",
              "sum_concb"
)
values <- c(1,   # "halflife" # 1 indicates the 3-hour halflife
            1,   # "radius"
            100, # "kon"
            2e-5, # "Mc"
            0.00001, # "rate"
            NA, # "feature"
            0, #"sum_concf"
            0) # "sum_concb"

num_rows <- 20
para <- as.data.frame(matrix(rep(values, num_rows), nrow = num_rows, byrow = TRUE))
colnames(para) <- col_names


# Efficiently populate parameter data using indexing
para[1:4, c("halflife", "feature")] <- cbind(c(1:4), "Plasma\n Half-life") # use index because the change is in the plasma PK profile
para[5:8, c("radius", "feature")] <- cbind(radius, "Molecule/particle\n size")
para[9:12, c("kon", "feature")] <- cbind(kon, "Binding\n affinity")
para[13:16, c("Mc", "feature")] <- cbind(Mc, "Mucin\n concentration")
para[17:20, c("rate", "feature")] <- cbind(rate, "Muco-ciliary\n clearance")


para$des <- text1 <- c(
  "t1/2 = 3 h","t1/2 = 6 h","t1/2 = 12 h","t1/2 = 24 h",                           # half-life
  "r = 1 nm","r = 5 nm","r = 20 nm","r = 50 nm",                                   # radius
  "Ka = 10 (1/M)","Ka = 100 (1/M)","Ka = 1000 (1/M)","Ka = 10000 (1/M)",           # binding
  "M = 1e-6 mol/L","M = 5e-6 mol/L", "M = 1e-5 mol/L","M = 5e-5 mol/L",            # mucin conc
  "Ke,m = 0 (1/h)","Ke,m = 0.0036 (1/h)","Ke,m = 0.036 (1/h)","Ke,m = 0.36 (1/h)") # muco-ciliary clearance


# Convert relevant columns to numeric
num_cols <- c("halflife", "radius", "kon", "Mc", "rate")
para[num_cols] <- lapply(para[num_cols], as.numeric)

save(file = "result/para_list.Rdata", para)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 5 - Mucus conc profiles  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate a blank list(20)
listdata <- list()
listdata_CV <- list()

# test
k = 1

# run in loop
for (k in 1:nrow(para)) {
  
    start_time <- Sys.time()
    
    # assign parameter
    dat.pk.plasma_ok <- dat.pk.plasma_clct[[para[k, "halflife"]]]
    r_ok = para[k, "radius"]
    kon_ok  = para[k, "kon"]
    Mc_ok  = para[k, "Mc"]
    clr_ok  = para[k, "rate"]
    
    # run pde model    
    dat.pk.mucus <- pde(
      db.pk = dat.pk.plasma_ok$Con1,
      # plasma PK
      coef.diff = Diffusion_coeff(r_ok,vis(Mc_ok)),
      # diffusion coefficient um^2/sec,113 for 1nm, 15.13 for 5nm
      kon = kon_ok ,
      # binding rate /(M*sec)
      koff = 0.1,
      # unbinding rate /sec
      Tmucin = Mc_ok,
      # concentration of total mucin mol/L
      dt = 0.5,
      # step size time in seconds
      span.t = 72,
      # observation time span in hours
      dx = 25,
      # step size depth in um
      span.x = 500,
      # observation depth span in um
      clr = clr_ok 
      # elimination rate in 1/sec
    ) 
    
    dat.pk.mucus.ok <- dat.pk.mucus %>% 
      rename(Conc_F = CONC_FREE) %>% 
      rename(Conc_B = CONC_BIND) %>%
      left_join(dat.pk.plasma, by = c("TIME" = "time")) %>% 
      mutate(feature = para[k,6]) %>%  
      mutate(description = para[k,9]) %>% 
      mutate(Ratio = Conc_F/Con1) %>% 
      mutate(ConcTT = Conc_F+Conc_B)  %>% 
      mutate(BRatio = Conc_B/ConcTT)  %>% 
      mutate(tRatio = Conc_B*0.01)
    
    listdata[[k]] <- dat.pk.mucus.ok
    
    para[k, 7:8] <-
      c(sum(dat.pk.mucus.ok$Conc_F),
        sum(dat.pk.mucus.ok$Conc_B))
    
    input_data <- dat.pk.mucus.ok[-1,] %>%
      group_by(TIME) %>%
      summarise(CV_F = 100*sd(Conc_F,na.rm = TRUE)/mean(Conc_F,na.rm = TRUE),
                CV_B = 100*sd(Conc_B,na.rm = TRUE)/mean(Conc_B,na.rm = TRUE)) %>% 
      mutate(description = para[k,9]) 
    
    listdata_CV[[k]] <- input_data

    
    end_time <- Sys.time()
    
    cat("time of simulation round: ", end_time - start_time, "\n")
    
    k = k + 1
} 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 6 - Data arrangement for plots  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_list <- list()

# put every 4 units of list data in one dataframe  
split_data <- lapply(seq(1, length(listdata), by = 4), function(j) {
  do.call(rbind, listdata[j:min(j+3, length(listdata))])
})

for (j in 1:5) {
  set_data <- split_data[[j]]
 
  # organize the order of different values for each factor
  set_data$description <- factor(set_data$description, levels=unique(set_data$description))
  
  # plot
  pCon_f <- ggplot(set_data, aes(TIME/60, DEPTH, fill = Conc_F)) +
    geom_tile() +
    guides(fill = guide_colorbar(
      title.position = "top",
      direction = "vertical")) +
    scale_fill_distiller("Conc_F \n(mg/L)", palette = 4,
                         direction = 1) +
    scale_x_continuous(name = "Time (hour)", breaks = c(0:100)*10) + 
    scale_y_continuous(name = "Distance (μm)") +
    theme_bw() +
    theme(text = element_text(family = "sans", color = "black", size = 9),
          axis.text = element_text(family = "sans", color = "black"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7.5),
          legend.key.width = unit(0.005, "npc"),
          legend.key.height = unit(0.3, 'cm'),
          strip.text = element_text(color = "black", size = 7.5, face = "plain"))+
    facet_grid(feature~description)
  pCon_f
  
  plot_list[[j]] <- pCon_f
}


# save the data
save(file = "result/f3_impactFactor.Rdata", listdata)
save(file = "result/para_list.Rdata", para)
save(file = "result/f3_impactFactor_CV.Rdata", listdata_CV)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 7 - Numeric analysis of AUC  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
factor_AUC <- para %>%
  mutate(AUCf = round(sum_concf/(60 * 21))) %>% # during calculation: 60 steps in one hour; 20 steps in space
  mutate(AUCb = round(sum_concb/(60 * 21),4))     # the unit is mg*hour/L

factor_AUC$des <- factor(factor_AUC$des, 
                             levels = unique(factor_AUC$des))
factor_AUC$type <- c(rep("Plasma half-life",4),
                     rep("Molecule/particle size",4),
                     rep("Binding affinity",4),
                     rep("Mucin concentration",4),
                     rep("Muco-ciliary clearance",4))
factor_AUC$type <- factor(factor_AUC$type, 
                      levels = c("Molecule/particle size",
                                 "Binding affinity",
                                 "Plasma half-life",
                                 "Mucin concentration",
                                 "Muco-ciliary clearance"))
factor_AUC$Num <- rep(1:4,5)
factor_AUC <- left_join(factor_AUC,plasma_AUC, by = "halflife") %>% 
  mutate(Ratio_AUCf_plasma = AUCf/AUC_plasma,
         Ratio_AUCb_plasma = AUCb/AUC_plasma)


AUC_plot_full <- factor_AUC %>% 
  filter(type != "Molecule/particle size") %>% 
  ggplot( aes(y = AUCf, x = des)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = round(AUCf, 3)), vjust = -0.5,size = 2.75) +
  facet_wrap(~type, nrow = 1, scales = "free_x") +
  scale_y_continuous(name = "AUCf (mg*h/L)", limits = c(0, 1500)) + # to replace # AUCf (mg*h/L)
  theme_bw() +
  theme(text = element_text(family = "sans", color = "black", size = 10),
        axis.title =  element_text(family = "sans", color = "black", size = 10),
        axis.text = element_text(family = "sans", color = "black", size = 9),
        axis.title.x = element_blank(),
        # text = element_text(family = "sans", color = "black", size = 10),
        # axis.title.x = element_blank(),
        # axis.title.y = element_text(size = 9),
        axis.text.x = element_text(
          size = 8,
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ),
        strip.text = element_text(size = 9, face = "plain"),
        panel.spacing.x = unit(0.8, "lines") # Adjust space between panels
  )
AUC_plot_full

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 8 - Numeric analysis of CV  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_list_cvf <- list()
plot_list_cvb <- list()

# put every 4 units of list data in one dataframe  
split_data <- lapply(seq(1, length(listdata_CV), by = 4), function(j) {
  do.call(rbind, listdata_CV[j:min(j+3, length(listdata_CV))])
})

alldata_cv <- do.call(rbind,split_data)
alldata_cv$Type <- rep(c(
  "Plasma half-life",
  "Molecule/particle size",
  "Binding affinity",
  "Mucin concentration",
  "Muco-ciliary clearance"), each = 17280) 

alldata_cv$description <- factor(alldata_cv$description, levels=
                                   c("r = 1 nm","r = 5 nm","r = 20 nm","r = 50 nm",                                   # radius
                                     "Ka = 10 (1/M)","Ka = 100 (1/M)","Ka = 1000 (1/M)","Ka = 10000 (1/M)",      
                                     "t1/2 = 3 h","t1/2 = 6 h","t1/2 = 12 h","t1/2 = 24 h",   # binding
                                     "M = 1e-6 mol/L","M = 5e-6 mol/L", "M = 1e-5 mol/L","M = 5e-5 mol/L",            # mucin conc
                                     "Ke,m = 0 (1/h)","Ke,m = 0.0036 (1/h)","Ke,m = 0.036 (1/h)","Ke,m = 0.36 (1/h)") # muco-ciliary clearance
)
alldata_cv$Type <- factor(alldata_cv$Type, levels=c("Molecule/particle size",
                                                    "Binding affinity",
                                                    "Plasma half-life",
                                                    "Mucin concentration",
                                                    "Muco-ciliary clearance"))

Cate <- c("Molecule/particle size",
          "Binding affinity",
          "Plasma half-life",
          "Mucin concentration",
          "Muco-ciliary clearance")
# plot
plot_list_cvf_legend <- list()

for(i in 1:5) {
  CVf_plot <- alldata_cv %>% 
    filter(Type == Cate[i]) %>% 
    ggplot() +
    geom_line(aes(TIME/60, CV_F, color = description)) +
    scale_x_continuous(name = "Time (hour)", breaks = c(0:100)*10) +
    scale_y_continuous(name = "Coefficient of variation (%)") +
    scale_color_manual(values = rep(c("#1F77B4", "#FF7F0E", "#2CA02C","#D62728"))) +
    guides(color = guide_legend(nrow = 4)) +
    labs(color = "") +
    facet_grid(~Type) +
    coord_cartesian(ylim = c(0,400)) +
    theme_bw() +
    theme(text = element_text(family = "sans", color = "black", size = 10),
          axis.text = element_text(family = "sans", color = "black"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          legend.position = c(0.5,0.8),
          legend.background = element_rect(fill = alpha("white", 0)),
          legend.key = element_rect(fill = NA),
          strip.text = element_text(color = "black", size = 9, face = "plain"))
  plot_list_cvf_legend[[i]] <- CVf_plot
  i = i + 1
}

CVf_plot <- cowplot::plot_grid(plot_list_cvf_legend[[1]], 
                   plot_list_cvf_legend[[3]]+ 
                     theme(axis.text.y = element_blank(),
                           axis.line.y = element_blank(),
                           axis.title.y= element_blank(),
                           axis.ticks.y= element_blank()),
                   plot_list_cvf_legend[[2]]+ 
                     theme(axis.text.y = element_blank(),
                           axis.line.y = element_blank(),
                           axis.title.y= element_blank(),
                           axis.ticks.y= element_blank()),
                   plot_list_cvf_legend[[4]]+ 
                     theme(axis.text.y = element_blank(),
                           axis.line.y = element_blank(),
                           axis.title.y= element_blank(),
                           axis.ticks.y= element_blank()),
                   plot_list_cvf_legend[[5]]+ 
                     theme(axis.text.y = element_blank(),
                           axis.line.y = element_blank(),
                           axis.title.y= element_blank(),
                           axis.ticks.y= element_blank()),
                   nrow = 1,
                   rel_widths = c(1.15,1,1,1,1),
                   align = 'h', axis = 'tb')
CVf_plot
