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
k = 1.380649e-23

### Absolute temperature
### unit: K
### 273 correspond to 0 celsius
t = 273 + 37

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


# since current radius of drugs is within 1000 nm
# just use the microviscosity to calculate the diffusion coefficient
### Diffusion coefficient
D = 1e+21*k*t/(6*pi*nu1*r) # um^2/s


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Calculation of diffusion coefficient [classical equations]  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to calculate the diffusion coefficients
Diffusion_coeff <- function(r){
  D <- 1e+21*k*t/(6*pi*nu1*r)
  data <- data.frame(x=r,y=D)
  return(data)
}

# calculation and plot
r<-seq(0.1,50,0.01)
d<-Diffusion_coeff(r)
P1 <- ggplot(d,aes(x,y))+
  geom_line(color="blue",linetype=1,alpha=0.3,size = 1)+
  scale_x_log10(name = "Radius (nm)") + 
  scale_y_continuous(name = expression(paste("Diffusion coefficient (",mu,"m^2/s)"))) +
  theme_bw() 
P1

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - Calculation of diffusion coefficient [modified equations]  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to calculate the diffusion coefficients with modified version
Diffusion_coeff <- function(r){
  sapply(r, function(r){
    if (r <= 1){
      D <- 1e+21*k*t/(4*pi*nu1*r)
    }
    else if (r > 1){
      D <- 1e+21*k*t/(6*pi*nu1*r)
    }
    data <- data.frame(x=r,y=D)
    return(data)
  })
}

r<-seq(0.1,50,0.01)
d<-Diffusion_coeff(r)
d <-t(d)
d <- as.data.frame(d)
d$x <- as.numeric(d$x)
d$y <- as.numeric(d$y)
P1M <- ggplot(d,aes(x,y))+
  geom_line(color="blue",linetype=1,alpha=0.3,size = 1)+
  scale_x_log10(name = "Radius (nm)") + 
  scale_y_continuous(name = expression(paste("Diffusion coefficient (",mu,"m^2/s)"))) +
  theme_bw() 
P1M

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - Find out the diffusion coefficient for specific size molecule  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# find the diffusion coefficient value for radius=1nm
Diffusion_coeff(1)
# 113.5308