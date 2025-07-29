######################################################################
### The script is used to compare the diffusion under 
## thin source and infinite source  ----
# y.guo@lacdr.leidenuniv.nl - July 2025
######################################################################


rm(list = ls(all = TRUE))
source("functions/pde_bound.R") ## PDE function (elimination)
source_values <- c("thin", "infinite") #  10, 100, 1000,
list_data <- list()

for (i in 1:length(source_values)) {
  
  start_time <- Sys.time()
  
  source_value <- source_values[i]
  
  if(!is.na(as.numeric(source_value))){
    source_value = as.numeric(source_value)
  }
  
  source_value
  
  dat.pk.mucus <- pde(
    db.pk = 10,
    # plasma PK
    coef.diff = 78.6,
    # diffusion coefficient um^2/sec,113 for 1nm, 15.13 for 5nm
    kon = 100 ,
    # binding rate /(M*sec)
    koff = 0.5,
    # unbinding rate /sec
    Tmucin = 2e-5,
    # concentration of total mucin mol/L
    dt = 0.2, # adjust to meet the calculation criterion
    # step size time in seconds
    span.t = 24*3,
    # observation time span in hours
    dx = 10,
    # step size depth in um
    span.x = 500,
    # observation depth span in um
    clr = 0,
    # elimination rate in 1/sec
    MW = 500,
    source = source_value
  ) 
  
  
  dat.pk.mucus.ok <- dat.pk.mucus %>% 
    rename(Conc_F = CONC_FREE) %>% 
    rename(Conc_B = CONC_BIND) %>% 
    mutate(source = source_value)
  
  list_data[[i]] = dat.pk.mucus.ok
  
  end_time <- Sys.time()
  
  cat("time of simulation round: ", end_time - start_time, "\n")
}

plot_list <- list()

all_data <- do.call(rbind,list_data)

pCon_f_thin <- ggplot(list_data[[1]] %>% mutate(source = "thin source"),
                      aes(TIME/60, DEPTH, fill = Conc_F)) +
  geom_tile() +
  guides(fill = guide_colorbar(
    title.position = "top",
    direction = "vertical")) +
  coord_cartesian(xlim = c(0,5)) + 
  scale_fill_distiller("Conc_F \n(mg/L)", palette = 4,
                       direction = 1) +
  scale_x_continuous(name = "Time (hour)", 
                     breaks = seq(0, 150, by = 2)) + 
  scale_y_continuous(name = "Distance (μm)") +
  theme_bw() +
  theme(text = element_text(family = "sans", color = "black", size = 9),
        axis.text = element_text(family = "sans", color = "black"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7.5),
        legend.key.width = unit(0.005, "npc"),
        legend.key.height = unit(0.3, 'cm'),
        strip.text = element_text(color = "black", size =9, face = "plain"),
        strip.background = element_rect(colour = "black", fill = "white"))+
  facet_grid(.~source)
pCon_f_thin
pCon_f_infinite <- ggplot(list_data[[2]] %>% mutate(source = "infinite source"),#%>% filter(source %in% c(10, 100, 1000)), 
                          aes(TIME/60, DEPTH, fill = Conc_F)) +
  geom_tile() +
  guides(fill = guide_colorbar(
    title.position = "top",
    direction = "vertical")) +
  scale_fill_distiller("Conc_F \n(mg/L)", palette = 4,
                       direction = 1) +
  scale_x_continuous(name = "Time (hour)", 
                     breaks = seq(0, 150, by = 2)) + 
  scale_y_continuous(name = "Distance (μm)") +
  coord_cartesian(xlim = c(0,5)) + 
  theme_bw() +
  theme(text = element_text(family = "sans", color = "black", size = 9),
        axis.text = element_text(family = "sans", color = "black"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7.5),
        legend.key.width = unit(0.005, "npc"),
        legend.key.height = unit(0.3, 'cm'),
        strip.text = element_text(color = "black", size =9, face = "plain"),
        strip.background = element_rect(colour = "black", fill = "white"))+
  facet_grid(.~source)
pCon_f_infinite

figure_source <-ggarrange(pCon_f_thin, pCon_f_infinite, labels = c("", ""), ncol = 2, nrow = 1)
figure_source