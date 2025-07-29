# Final: Imipenem-Planktonic

pde <- function(db.pk, coef.diff, kon, koff, Tmucin, dt, span.t, dx, span.x, clr,
                pd.bl ,kgs ,kgr,ksr, slopeS, slopeR, Bmax) {
  require(dplyr)
  require(tibble)
  
  # back up parameters: Emaxs, EC50s, hills
  
  # sanity check ----
  # db.pk = rep(100, 4321) # plamsa conc, mg/L, per minite
  # coef.diff =104 # #um^2/s
  # kon = 100 # binding rate /(M*sec)
  # koff = 0.1 # unbinding rate /sec
  # Tmucin = 2e-5 # concentration of total mucin mol/L
  # dt = 1 #  s, # adjust to meet the calculation criterion
  # span.t = 24*3 # observation time span in hours
  # dx = 10 # step size depth in um
  # span.x = 500 # observation depth span in um
  # clr = 0 # elimination rate in 1/sec
  # MW = 500 # molecular weight (g/mol), to calculate the concentration with mol/L for binding kinetics
  # bl: baseline of bacteria (log10 CFU/mL), pd.bl is the baseline dataset for different depths
  # kgs :natural growth rate of sensitive species (/s)
  # kgr :natural growth rate of resistant species (/s)
  # ksr : transfer rate
  # slopeS: killing rate of antibiotics on sensitive species (L/mg/s)
  # slopeR: killing rate of antibiotics on resistant species (L/mg/s)
  # Bmax: a parameter, maximum bacteria count (CFU/mL) 
  
  Lnode <- span.x/dx
  Tnode <- span.t*60*60/dt
  
  if (length(db.pk) < span.t * 60) {
    stop("Error: length of db.pk is shorter than required number of time steps (span.t * 60).")
  }
  
  if (coef.diff * dt / dx ^ 2 > 0.5) {
    stop("stability criterion unmet\n", coef.diff * dt / dx ^ 2)
  }
  
  # Defining nodal points at different intervals (1000 um100 intervals101 nodal points)
  # x <- 0:Lnode
  # Initializing drug concentration at the mid nodes
  A <- rep(0, Lnode+1)
  B <- rep(0, Lnode+1)
  S <- rep(0, Lnode+1) # sensitive species
  R <- rep(0, Lnode+1) # resistant species 
  
  # Constants used in calculation
  a <- ((coef.diff * dt) / dx ^ 2)
  b <- 1 - (2 * coef.diff * dt / dx ^ 2)
  c <- ((coef.diff * dt) / dx ^ 2)
  d <- kon # 1/second linear binding rate
  e <- koff # 1/second linear unbinding rate
  f <- clr # 1/second linear elimination rate
  
  sumA <- rep(0, Tnode)
  # Increasing the value of the time step or the end time (i=1:timestep:end
  # time) will result in convergence
  
  ## Initialization
  k <- 2
  An <- rep(NA_real_, Lnode) # Free concentration
  Bn <- rep(NA_real_, Lnode) # Bound concentration
  Sn <- rep(NA_real_, Lnode) # Bacterial concentration (sensitive)
  Rn <- rep(NA_real_, Lnode) # Bacterial concentration (resistant)
  
  # df_C <- tibble(TIME = numeric(0), DEPTH = numeric(0), CONC_FREE = numeric(0), CONC_BIND = numeric(0))
  df_C <- matrix(rep(0, 6), ncol = 6) # times of rep and ncol used to be 4
  
  for (i in 1:Tnode) { # Time progression 
    
    # Concentration at the Boundary nodes(Boundary Condition)
    if (i %% (60/dt) == 0) { # Drug concentration data is for every minute (120 time step=1 min). 
      
      # cat("Time unit = ", k, "\n", sep = "")
      
      df_C_i <- matrix(data = c(rep(i/(60/dt), Lnode+1),
                                seq(from = 0, by = dx, length.out = Lnode+1),
                                A,
                                B,
                                S,
                                R), ncol = 6) # add a S to represent bacteria, change ncol to be 5
      
      # df_C_i <- tibble(TIME = rep(i/(60/dt), Lnode+1),
      #                  DEPTH = seq(from = 0, by = dx, length.out = Lnode+1),
      #                  CONC_FREE = A,
      #                  CONC_BIND = B)
      df_C <- rbind(df_C, df_C_i)
      
      k <- k + 1 # So for each minute we update the boundary condition
      
      An[1] <- db.pk[k] # Exported drug concentration on tissue-mucus interface
      A[1] <- An[1]
      
      # } else if (k < 480) {
    } else if (k < span.t*60) {
      An[1] <- db.pk[k] # Drug concentration data will stay constant in time steps within every minute (120 time step=1 min).
      A[1] <- An[1]
      
    } 
    # else {
    #   An[1] <- 0 # For the time steps when plasma concentraion=0.
    #   A[1] <- An[1]  
    # }
    
    # Bound concentration
    Bn[1] <- d*(Tmucin-(B[1]*0.001/500))*A[1]*dt - e*B[1]*dt + B[1] - f*B[1]*dt 
    B[1] <- Bn[1]
    
    # Sensitive population of bacteria; Resistant population of bacteria
    if (i == 1) {
      Sn[1] <- pd.bl[1]
      Rn[1] <- 0
    }else { 
      growths <- kgs*(1-((S[1]+R[1])/Bmax))
      # effects <- Emaxs*(A[1])^hills /(EC50s^hills + (A[1])^hills)
      
      Sn[1] <- S[1]*(1 + growths*dt - slopeS*A[1]*dt-ksr*dt) - f*S[1]*dt 
      growthr <- kgr*(1-((S[1]+R[1])/Bmax))
      Rn[1] <- R[1]*(1 + growthr*dt - slopeR*A[1]*dt) + ksr*dt*S[1] - f*R[1]*dt 
    }
    
    S[1] <- Sn[1]
    R[1] <- Rn[1]
    # if (Sn[1] >= 1) {Sn[1] <- Sn[1]} else {Sn[1] <- 0}
    # if (Rn[1] >= 1) {Rn[1] <- Rn[1]} else {Rn[1] <- 0}
    
    for (j in 2:Lnode) { # Solution of PDE in space for each time step 
      # a, b, c diffusion term
      # d binding term
      # e unbinding term
      # by finite difference method
      An[j] <- a*A[j-1] + b*A[j] + c*A[j+1] - d*(Tmucin-(B[j]*0.001/500))*A[j]*dt + e*B[j]*dt - f*A[j]*dt ##500 is the assumed molecular weight of The drug, help to adjust the unit of Kon*M
      Bn[j] <- d*(Tmucin-(B[j]*0.001/500))*A[j]*dt - e*B[j]*dt + B[j] - f*B[j]*dt 
      
      A[j] <- An[j]
      B[j] <- Bn[j]
      
      if (i == 1) {
        Sn[j] <- pd.bl[j]
        Rn[j] <- 0
      }else { 
        growths <- kgs*(1-((S[j]+R[j])/Bmax))
        # effects <- Emaxs*(A[j])^hills /(EC50s^hills + (A[j])^hills)
        
        Sn[j] <- S[j]*(1 + growths*dt - slopeS*A[j]*dt-ksr*dt) - f*S[j]*dt 
        growthr <- kgr*(1-((S[j]+R[j])/Bmax))
        Rn[j] <- R[j]*(1 + growthr*dt - slopeR*A[j]*dt) + ksr*dt*S[j] - f*R[j]*dt 
      }
      
      S[j] <- Sn[j]
      R[j] <- Rn[j]
      # S[j] <- if (Sn[j] >= 1) Sn[j] else 0
      # R[j] <- if (Rn[j] >= 1) Rn[j] else 0
    }
    
    
    A[Lnode+1] <- A[Lnode] # Reflective boundary at upper surface of mucus
    B[Lnode+1] <- B[Lnode] # Reflective boundary at upper surface of mucus
    S[Lnode+1] <- S[Lnode] # Reflective boundary at upper surface of mucus
    R[Lnode+1] <- R[Lnode] # Reflective boundary at upper surface of mucus
  }
  
  df_C %>%
    as_tibble() %>%
    rename(TIME = V1, DEPTH = V2, CONC_FREE = V3, CONC_BIND = V4, Sensitive = V5, Resistant = V6) %>% 
    return()
}
