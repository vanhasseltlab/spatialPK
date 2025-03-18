

pde <- function(db.pk, coef.diff, kon, koff, Tmucin, dt, span.t, dx, span.x, clr) {
  require(dplyr)
  require(tibble)
  
  # Diffusivity coefficient um^2/s
  # alpha <- 10 ^ 2 #um^2/s
  # Time step
  # dt <- 0.5 #s
  # Distance between nodal points in um
  # dx <- 10 #um
  # Number of time nodes (total simulation time = dt*Tnode)
  # Tnode <- 57598
  Tnode <- span.t*60*60/dt
  # Number of space nodes (total mucus lenght = dx*Lnode)
  # Lnode <- 100
  Lnode <- span.x/dx
  # Note:Satisfying the stability criterion (alpha*dt/dx^2 <= 0.5)
  if (coef.diff * dt / dx ^ 2 > 0.5) {
    stop("stability criterion unmet\n", coef.diff * dt / dx ^ 2)
  }
  
  
  # Defining nodal points at different intervals (1000 um100 intervals101 nodal points)
  # x <- 0:Lnode
  # Initializing drug concentration at the mid nodes
  A <- rep(0, Lnode+1)
  B <- rep(0, Lnode+1)
  
  
  # Constants used in calculation
  a <- ((coef.diff * dt) / dx ^ 2)
  b <- 1 - (2 * coef.diff * dt / dx ^ 2)
  c <- ((coef.diff * dt) / dx ^ 2)
  # d <- 1.12 * 10 ^ -4 #1/second linear binding rate
  # e <- 0 * 10 ^ -3 #1/second linear unbinding rate
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
  
  # df_C <- tibble(TIME = numeric(0), DEPTH = numeric(0), CONC_FREE = numeric(0), CONC_BIND = numeric(0))
  df_C <- matrix(rep(0, 4), ncol = 4)
  
  for (i in 1:Tnode) { # Time progression 
    
    # Concentration at the Boundary nodes(Boundary Condition)
    if (i %% (60/dt) == 0) { # Drug concentration data is for every minute (120 time step=1 min). 
      
      # cat("Time unit = ", k, "\n", sep = "")
      
      df_C_i <- matrix(data = c(rep(i/(60/dt), Lnode+1),
                                seq(from = 0, by = dx, length.out = Lnode+1),
                                A,
                                B), ncol = 4)

      df_C <- rbind(df_C, df_C_i)
      
      k <- k + 1 # So for each minute we update the boundary condition
      
      An[1] <- db.pk[k] # Exported drug concentration on tissue-mucus interface
      A[1] <- An[1]
      
      # } else if (k < 480) {
    } else if (k < span.t*60) {
      An[1] <- db.pk[k] # Drug concentration data will stay constant in time steps within every minute (120 time step=1 min).
      A[1] <- An[1]
      
    } 
    
    
    
    # Bound concentration
    Bn[1] <- d*(Tmucin-(B[1]*0.001/500))*A[1]*dt - e*B[1]*dt + B[1] - f*B[1]*dt 
    B[1] <- Bn[1]
    
    for (j in 2:Lnode) { # Solution of PDE in space for each time step
      # a, b, c diffusion term
      # d binding term
      # e unbinding term
      # by finite difference method
      An[j] <- a*A[j-1] + b*A[j] + c*A[j+1] - d*(Tmucin-(B[j]*0.001/500))*A[j]*dt + e*B[j]*dt - f*A[j]*dt ##500 is the assumed molecular weight of The drug, help to adjust the unit of Kon*M
      Bn[j] <- d*(Tmucin-(B[j]*0.001/500))*A[j]*dt - e*B[j]*dt + B[j] - f*B[j]*dt 
      
      A[j] <- An[j]
      B[j] <- Bn[j]
      
    }
    A[Lnode+1] <- A[Lnode] # Reflective boundary at upper surface of mucus
    B[Lnode+1] <- B[Lnode] # Reflective boundary at upper surface of mucus
    
  }
  
  df_C %>%
    as_tibble() %>%
    rename(TIME = V1, DEPTH = V2, CONC_FREE = V3, CONC_BIND = V4) %>% 
    return()
}
