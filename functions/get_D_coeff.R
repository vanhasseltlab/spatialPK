##################################################################################
### The script is used for calculation of diffusion coefficients based on radius  ----
# y.guo@lacdr.leidenuniv.nl - March 2024
##################################################################################
rm(list = ls(all = TRUE))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Mechanistic knowledge   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reference link: https://faseb.onlinelibrary.wiley.com/doi/abs/10.1096/fasebj.2020.34.s1.04366
# small molecule (<1000 nm): microviscosity + modified ES equation
# large molecule (>1000 nm): macroviscosity + ES equation

### Diffusion efficient
### Equation is 
### D = k * t / (6*pi*nu*r) 
# nu is viscosity of liquid
# r is radius of molecule

### Boltzmann's constant
### 1.380649e-23
### unit: J/K or N*m/K
k_value = 1.380649e-23

### Absolute temperature
### unit: K
### 273 correspond to 0 celsius
t_value = 273 + 37

### microviscosity (https://doi.org/10.1096/fasebj.2020.34.s1.04366)
### mucus in CF patient: 3 cP
### unit: N*s/m^2
### 1000 cP = 1 N*s/m^2
nu1 = 3 * 0.001 

### the dependency between mucin conc and micro-viscosity (η)
### source: https://pubs-acs-org.ezproxy.leidenuniv.nl/doi/10.1021/acsnano.4c14927
### info: 
# mucin conc 10-20 mg/mL, η being 3-7 cP
# mucin conc 200 mg/mL, η being 45 cP
# molecular weight: MUC5B ~ 596 kDa (596*10^3 g/mol)
# x <- c(15/(596*10^3), 200/(596*10^3)) # mol/L
# y <- c(5, 45) # cP
# fit <- lm(y ~ x)
# coef(fit)
# coef(fit)[1]
# linear regression
# unit of vis is N*s/m^2
vis <- function(mucin) 0.001*(1.756757 + 1.288649e+05 * mucin)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Calculation of diffusion coefficient [classical equations]  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate the viscosity
vis <- function(mucin) 0.001*(1.756757 + 1.288649e+05 * mucin)
# unit
# mucin: mol/L
# vis: N*s/m^2


# Function to calculate the diffusion coefficients
Diffusion_coeff <- function(r, nu1 = vis_value){
  if (r <= 1){
    D <- 1e+21*k_value*t_value/(4*pi*nu1*r)
  }
  else if (r > 1){
    D <- 1e+21*k_value*t_value/(6*pi*nu1*r)
  }
  return(D)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - Find out the diffusion coefficient for specific size molecule  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# find the diffusion coefficient value for radius=1nm
# Diffusion_coeff(1, nu1 = vis(2e-5))
# 78.6