# PDE solved by finite difference method
# modified by Yuchen Guo
# DATE 26/06/2025
# Main addition: comparison with the thin source and infinite source

pde <- function(db.pk, coef.diff, kon, koff, Tmucin, dt, span.t, dx, span.x, clr, MW=500, source ="infinite") {
  require(dplyr)
  require(tibble)
  
  # # # sanity check ----
  # db.pk = 100 # plamsa conc, mg/L, per minite
  # coef.diff =104 # #um^2/s
  # kon = 100 # binding rate /(M*sec)
  # koff = 0.1 # unbinding rate /sec
  # Tmucin = 2e-5 # concentration of total mucin mol/L
  # dt = 0.2 #  s, # adjust to meet the calculation criterion
  # span.t = 24*3 # observation time span in hours
  # dx = 20 # step size depth in um
  # span.x = 500 # observation depth span in um
  # clr = 0 # elimination rate in 1/sec
  # MW = 500 # molecular weight (g/mol), to calculate the concentration with mol/L for binding kinetics
  # source = "thin" # if character("infinite"/"thin"), specify the source to be either "thin" or "infinite"; if numeric, specify the volume ratio of donor-receiver compartments
  
  Lnode <- span.x/dx
  Tnode <- span.t*60*60/dt
  
  
  if (coef.diff * dt / dx ^ 2 > 0.5) {
    stop("stability criterion unmet\n", coef.diff * dt / dx ^ 2)
  }
  
  # Initialization ----
  A <- rep(0, Lnode + 1)                   # free drug concentration
  B <- rep(0, Lnode + 1)                   # bound drug concentration
  A[1] <- db.pk
  An <- rep(NA_real_, Lnode)
  Bn <- rep(NA_real_, Lnode)
  df_C <- matrix(0, ncol = 4, nrow = 0)    # results storage
  
  # Coefficients ----
  a <- ((coef.diff * dt) / dx ^ 2)
  b <- 1 - 2 * a
  c <- a
  d <- kon
  e <- koff
  f <- clr
  
  # boundary condition ----
  
  # Time loop ----
  for (i in 1:Tnode) {
    
    
    if (i %% (60/dt) == 0) { 
      df_C_i <- matrix(data = c(rep(i/(60/dt), Lnode+1),
                                seq(from = 0, by = dx, length.out = Lnode+1),
                                A,
                                B), ncol = 4)
      
      df_C <- rbind(df_C, df_C_i)
      
    } 
    
    if(source == "infinite"){
      # infinite source
      A[1] <- db.pk
    } else if (source == "thin") { # thin source
      A[1] <- (1-a)*A[1] + c*A[2] - d*(Tmucin-(B[1]*0.001/MW))*A[1]*dt + e*B[1]*dt - f*A[1]*dt
    } 

    
    if(is.na(as.numeric(source))){
      Bn[1] <- d*(Tmucin-(B[1]*0.001/MW))*A[1]*dt - e*B[1]*dt + B[1] - f*B[1]*dt 
      B[1] <- Bn[1] 
    } 

    
    for (j in 2:Lnode) {
      An[j] <- a*A[j-1] + b*A[j] + c*A[j+1] - d*(Tmucin-(B[j]*0.001/MW))*A[j]*dt + e*B[j]*dt - f*A[j]*dt
      Bn[j] <- d*(Tmucin-(B[j]*0.001/MW))*A[j]*dt - e*B[j]*dt + B[j] - f*B[j]*dt 
      A[j] <- An[j]
      B[j] <- Bn[j]
    }
    A[Lnode+1] <- A[Lnode]
    B[Lnode+1] <- B[Lnode]
  }
  
  df_C %>%
    as_tibble() %>%
    rename(TIME = V1, DEPTH = V2, CONC_FREE = V3, CONC_BIND = V4) %>% 
    return()
}