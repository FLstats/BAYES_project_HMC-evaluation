setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

rm(list = ls())  
library(tidyverse)
theme_set(theme_bw())
library(latex2exp)
library(patchwork)
library(grid)
library(glue)

ggplot() + theme_void() + theme(plot.background=element_rect(fill="lightblue1"))
select <- dplyr::select
filter <- dplyr::filter

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
"#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(rstan)
rstan_options(auto_write = TRUE)              # Cache compiled Stan programs
options(mc.cores = parallel::detectCores())   # Parallelize chains

util <- new.env()                             # Creates clean new environment
source('stan_utility.R', local=util)          # Source the file located in WD
ls(util)                                      # Check the utils functionality

library(loo)
library(posterior)
options(posterior.num_args=list(digits=2))
library(bayesplot)
ggplot2::theme_set(bayesplot::theme_default(base_family="sans", base_size=14))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                               Neal's funnel test1                           #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# ----------------------------------------- #
###               Data prep               ###
# ----------------------------------------- #
data_list <- list("P" = 9)

# ----------------------------------------- #
###               Fit model               ###
# ----------------------------------------- #
nuts_controls <- list(max_treedepth = 15, adapt_delta = 0.99)

fit <- stan(file = "stan/4.neals_funnel.stan",
            data = data_list,
            seed = 42,
            chains = 4,
            iter = 2e3, # warmup = iter/2 = 1e3
            refresh = 500,
            control = nuts_controls)

fit

# Divergence scatter plot
draws <- as_draws_array(fit)
hmc_diagnostics <- nuts_params(fit)
color_scheme_set("darkgray")
div_style <- scatter_style_np(div_color = "green", div_size = 2.5, div_alpha = 0.75)
mcmc_scatter(draws,
             pars = c("x[1]", "v"),
             np = hmc_diagnostics,
             transformations = list(v = "log"),
             np_style = div_style) +
  labs(x = TeX("$x_1$"), y = TeX("$\\log(\\nu)$"))


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                             Hybrid Rosenbrock R                             #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

# --------------------------------------------------- #
#               DATA-SIMULATING FUNCTION              #
# --------------------------------------------------- #
hybrid_rosenbrock <- function(N_samp, ni, nj, mu=1, a=1/20, b=100/20, seed=42) {
  
  set.seed(seed)
  
  ### CALCULATE INPUTS
  sd1 <- 1/sqrt(2*a)      # sd x1
  sd_b <- 1/sqrt(2*b)   # sd next variables
  n_dim <- (ni-1)*nj+1    # total variable dim
  
  ### DEFINE RESULT MATRIX
  X <- matrix(NA, nrow = N_samp, ncol = n_dim)
  
  ### DEFINE MATRIX COLNAMES
  make_colnames <- function(ni, nj) {
    cols <- "x1"
    for (j in 1:nj) {
      for (i in 2:ni) {
        cols <- c(cols, paste0("x", j, i))
      }
    }
    cols
  }
  
  colnames(X) <- make_colnames(ni, nj)

  ### SIMULATE VARIABLES
  for (k in seq_len(N_samp)) {
    x1 <- rnorm(1, mu, sd1)
    
    sample_vec <- numeric(n_dim)
    sample_vec[1] <- x1
    idx <- 2
    
    for (j in 1:nj) {
      x_prev <- x1
      
      for (i in 2:ni) {
        x_new <- rnorm(1, x_prev^2, sd_b)
        sample_vec[idx] <- x_new
        x_prev <- x_new
        idx <- idx + 1
      }
    }
    
    X[k, ] <- sample_vec
  }
  
  return(as.data.frame(X))
}

# SIMULATE
hrb <- hybrid_rosenbrock(N_samp = 1e4, ni=3, nj=2)

# --------------------------------------------------- #
#                      Plots                          #
# --------------------------------------------------- #
library(GGally)

# ----------------------------------------- #
###             BIVARIATE PLOTS           ###
# ----------------------------------------- #
### CONTOUR PLOT
ggplot(data = hrb, aes(x=x1, y=x12)) +
  geom_point(alpha = 0.1, size = 0.6) +
  stat_density2d(aes(fill = after_stat(level)),
                 geom = "polygon",
                 contour_var = "ndensity",
                 bins = 9, # number of contours
                 n = 300, # higher grid resolution
                 h = c(bw.nrd(hrb$x1), bw.nrd(hrb$x12)) # gentler smoothing
                 ) +
  scale_fill_viridis_c(guide = "none")

### HISTOGRAM
ggplot(data = hrb, aes(x=x1)) +
  geom_histogram(color="black", fill=cbbPalette[1],
                 bins = 80)

# ----------------------------------------- #
###               PAIRS PLOT              ###
# ----------------------------------------- #
# Diagonal: histogram
panel_hist <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_histogram(color="black", fill=cbbPalette[1],
                   bins = 80) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank())
}

# Lower triangle: filled density contours
panel_contour <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point(alpha = 0.1, size = 0.6) +
    stat_density2d(aes(fill = after_stat(level)),
                   geom = "polygon",
                   contour_var = "ndensity",
                   bins = 9, # number of contours
                   n = 300, # higher grid resolution
                   h = c(bw.nrd(hrb$x1), bw.nrd(hrb$x12)) # gentler smoothing
    ) +
    scale_fill_viridis_c(guide = "none") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank())
}

# Build the plot
pairs_plot <- ggpairs(hrb,
             columns = 1:ncol(hrb),
             upper = "blank",
             diag  = list(continuous = panel_hist),
             lower = list(continuous = panel_contour),
             progress = FALSE) + 
  theme(legend.position = "none",
        strip.text = element_text(size = 9),
        panel.grid = element_blank())

pairs_plot







