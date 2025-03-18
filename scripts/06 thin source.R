# literature:
# https://ocw.mit.edu/courses/3-205-thermodynamics-and-kinetics-of-materials-fall-2006/90b4b7376505ab14ba82d95309083961_lecture03_slides.pdf
# https://ebrary.net/184562/engineering/fick_s_laws_diffusion#359705

# Load libraries
library(ggplot2)
library(viridis)

# Define parameters
D <- 200             # Diffusion coefficient
M0 <- 1              # Initial mass at the source
x_vals <- seq(0, 50, by = 1)   # Distance values
time_vals <- seq(1, 50, by = 1)  # Time values (start from 1 to avoid division by zero)

# Calculate thin-film solution
thin_film_solution <- function(x, t, M0, D) {
  (M0 / sqrt(4 * pi * D * t)) * exp(-x^2 / (4 * D * t))
}

# Create a data frame for all time and distance combinations
data <- expand.grid(Distance = x_vals, Time = time_vals)
data$Concentration <- mapply(thin_film_solution, x = data$Distance, t = data$Time, MoreArgs = list(M0 = M0, D = D))

# Plot heatmap
p_thin <- ggplot(data, aes(x = Time, y = Distance, fill = Concentration)) +
  geom_tile() +
  scale_fill_distiller("Conc_F \n(mg/L)", palette = 4,
                       direction = 1) +
  labs(
    title = "thin-source diffusion",
    x = "Time",
    y = "Distance from source"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )



# infinite-constant film diffusion
# Load libraries
library(ggplot2)
library(viridis)

# Define parameters
D <- 200            # Diffusion coefficient
C0 <- 1              # Constant concentration at the source
x_vals <- seq(0, 50, by = 1)   # Distance values
time_vals <- seq(1, 50, by = 1)  # Time values (start from 1 to avoid division by zero)

# Define the infinite-film solution using the complementary error function (erfc)
infinite_film_solution <- function(x, t, C0, D) {
  C0 * pracma::erfc(x / (2 * sqrt(D * t)))  # erfc from the pracma package
}

# Create a data frame for all time and distance combinations
data <- expand.grid(Distance = x_vals, Time = time_vals)
data$Concentration <- mapply(infinite_film_solution, x = data$Distance, t = data$Time, MoreArgs = list(C0 = C0, D = D))

# Plot heatmap
p_infinite <- ggplot(data, aes(x = Time, y = Distance, fill = Concentration)) +
  geom_tile() +
  scale_fill_distiller("Conc_F \n(mg/L)", palette = 4,
                       direction = 1,limits = c(0, 1)) +
  labs(
    title = "infinite-source diffusion",
    x = "Time",
    y = "Distance from source"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

figure_source <-ggarrange(p_thin, p_infinite, labels = c("A", "B"), ncol = 2, nrow = 1)
figure_source
ggsave(figure_source, file = "figure/source.pdf", width = 18, height = 10, units = "cm", device = cairo_pdf)
ggsave(figure_source, file = "figure/source.tiff", width = 22, height = 10, units = "cm", dpi = 600)
