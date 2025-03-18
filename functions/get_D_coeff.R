##################################################################################
### The script is used for calculation of diffusion coefficients based on radius  ----
# y.guo@lacdr.leidenuniv.nl - March 2024
##################################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Mechanistic knowledge   ----
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
nu1 = 3 * 0.001 

### macroviscosity (https://doi.org/10.1096/fasebj.2020.34.s1.04366)
### mucus in CF patient: 14-110*10^3 cP
### take the average value: 62*10^3 cP
### unit: N*s/m^2
nu2 = 62*10^3 * 0.001 

### Radius
### 1 nm
### unit: m
### range : 0.1 - 1000 nm

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Calculation of diffusion coefficient [modified equations]  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to calculate the diffusion coefficients with modified version
Diffusion_coeff <- function(r){
    if (r <= 1){
      D <- 1e+21*k_value*t_value/(4*pi*nu1*r)
    }
    else if (r > 1){
      D <- 1e+21*k_value*t_value/(6*pi*nu1*r)
    }
    return(D)
}
