setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

rm(list = ls())  
library(tidyverse)
theme_set(theme_bw())
library(latex2exp)
library(patchwork)
library(grid)

ggplot() + theme_void() + theme(plot.background=element_rect(fill="lightblue1"))
select <- dplyr::select
filter <- dplyr::filter
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                               Rosenbrock 2d                                 #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# ----------------------------------------- #
###               Data prep               ###
# ----------------------------------------- #
data_list <- list("mu" = 1, "a" = 0.05, "b" = 0.05)

# ----------------------------------------- #
###               Fit model               ###
# ----------------------------------------- #
nuts_controls <- list(max_treedepth = 15, adapt_delta = 0.99)

fit <- stan(file = "stan/5.rosenbrock_2d.stan",
            data = data_list,
            seed = 42,
            chains = 4,
            iter = 2e3, # warmup = iter/2 = 1e3
            refresh = 500,
            control = nuts_controls)

fit

util$check_all_diagnostics(fit)

draws <- as_draws_array(fit)
hmc_diagnostics <- nuts_params(fit)
color_scheme_set("darkgray")
div_style <- scatter_style_np(div_color = "green", div_size = 2.5, div_alpha = 0.75)
mcmc_scatter(draws,
             pars = c("x1", "x2"),
             np = hmc_diagnostics,
             np_style = div_style)


# --------------------------------------------------- #
#                  Rosenbrock ncp                     #
# --------------------------------------------------- #
fit <- stan(file = "stan/5.rosenbrock_2d_ncp.stan",
            data = data_list,
            seed = 42,
            chains = 4,
            iter = 2e3, # warmup = iter/2 = 1e3
            refresh = 500,
            control = nuts_controls)

fit

util$check_all_diagnostics(fit)

draws <- as_draws_array(fit)
hmc_diagnostics <- nuts_params(fit)
color_scheme_set("darkgray")
div_style <- scatter_style_np(div_color = "green", div_size = 2.5, div_alpha = 0.75)
mcmc_scatter(draws,
             pars = c("x1", "x2"),
             np = hmc_diagnostics,
             np_style = div_style)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                           Plot Rosenbrock 2d                                #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

# Rosenbrock log-density (unnormalized)
log_pi <- function(x1, x2) {
  - (100*(x2 - x1^2)^2 + (1 - x1)^2) / 20
}

# Grid (roughly like the paper’s axes)
x1_rng <- c(-10, 15)
x2_rng <- c(-5, 150)

# Choose a resolution that renders fast but looks smooth
dx1 <- 0.05
dx2 <- 0.5
x1 <- seq(x1_rng[1], x1_rng[2], by = dx1)
x2 <- seq(x2_rng[1], x2_rng[2], by = dx2)

grid <- expand_grid(x1 = x1, x2 = x2) |>
  mutate(lpdf = log_pi(x1, x2),
         pdf  = exp(lpdf))

# Numerically normalize to get a proper pdf on this grid
cell_area <- dx1 * dx2
Z <- sum(grid$pdf) * cell_area
grid <- grid |> mutate(pdf_n = pdf / Z)

# Target probabilities to be enclosed by contours (labels like the figure)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999)

# Compute HPD cutoffs: the density threshold c_p such that
# integral{ pdf_n(x) * I(pdf_n(x) >= c_p) dx } = p
hpd_cutoff <- function(p, dens, area) {
  ord <- order(dens, decreasing = TRUE)
  dens_sorted <- dens[ord]
  cumprob <- cumsum(dens_sorted) * area
  # smallest cutoff with cumprob >= p
  idx <- which(cumprob >= p)[1]
  dens_sorted[idx]
}

cutoffs <- map_dbl(probs, ~ hpd_cutoff(.x, grid$pdf_n, cell_area))

# Plot: contour levels at those HPD cutoffs, colored by enclosed probability
ggplot(grid, aes(x = x1, y = x2, z = pdf_n)) +
  geom_contour(aes(color = after_stat(level)),
               breaks = cutoffs, linewidth = 0.6) +
  # relabel the legend in terms of probability contained
  scale_color_continuous(
    name = "Probability Contained",
    breaks = cutoffs,
    labels = as.character(probs)
  ) +
  labs(x = expression(x[1]), y = expression(x[2]),
       title = "Contour plot of the 2-D Rosenbrock density") +
  coord_cartesian(xlim = x1_rng, ylim = x2_rng) +
  theme_minimal(base_size = 12)
