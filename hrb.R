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
#                                 GENERATE DATA                               #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# n_i and n_j specify the structure of the target Rosenbrock distribution. We
# want to test the same model against different data sets, so we will keep 
# these variable dimensions fixed. We will only vary mu, a, and b.

# --------------------------------------------------- #
#           Specify variables & parameters            #
# --------------------------------------------------- #
# Mod 1
N <- 1e4; ni <- 3; nj <- 2; mu <- 1; a <- 0.05; b <- 5
# Mod 2
N <- 1e4; ni <- 3; nj <- 2; mu <- 1; a <- 0.5; b <- 5

# --------------------------------------------------- #
#             RUN STAN DATA GENERATION                #
# --------------------------------------------------- #
data_list <- list(
  "N" = N,
  "ni" = ni,
  "nj" = nj,
  "mu"= mu,
  "a" = a,
  "b" = b
)

mod_sim_hbr <- stan_model(file = "stan/rosenbrock/hrb_simulate-data.stan")

sim <- sampling(mod_sim_hbr,
                data = data_list,
                algorithm = "Fixed_param", # Stan only runs generated quantities
                seed = 42,
                iter = 1, # generates one data set (iter*chains)
                chains = 1
)

# --------------------------------------------------- #
#               EXTRACT SIMULATED DATA                #
# --------------------------------------------------- #
# 10 000 samples, 5 variables per sample = 50 000 variables + 1 "lp__" variable.
# iter = 1 gives only 1 draw.
dim(sim) # 1 draw, 1 chain, 50001 variables

# Extract generated variables X as array. dim = (iter, N, n_dim)
X_array <- extract(sim)$X 
# Drop one dimension to make it into a matrix
hrb_sim_data <- drop(X_array)
dim(hrb_sim_data)

# Define colnames
make_colnames <- function(ni, nj) {
  cols <- "x1"
  for (j in 1:nj) {
    for (i in 2:ni) {
      cols <- c(cols, paste0("x", j, i))
    }
  }
  cols
}

colnames(hrb_sim_data) <- make_colnames(ni, nj)
head(hrb_sim_data)

# --------------------------------------------------- #
#             SAVE AND PLOT SIMULATED DATA            #
# --------------------------------------------------- #
hrb_sim_df <- as.data.frame(hrb_sim_data)

# saveRDS(hrb_sim_df, file = "data/hrb_gen_data1.RData")
# saveRDS(hrb_sim_df, file = "data/hrb_gen_data2.RData")

# ----------------------------------------- #
###               PAIRS PLOT              ###
# ----------------------------------------- #
library(GGally)
library(rlang)

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
                   h = c(bw.nrd(data[[as_name(mapping$x)]]), 
                         bw.nrd(data[[as_name(mapping$y)]])) # gentler smoothing
    ) +
    scale_fill_viridis_c(guide = "none") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank())
}

# Build the plot
pairs_plot <- ggpairs(hrb_sim_df,
                      columns = 1:ncol(hrb_sim_df),
                      upper = "blank",
                      diag  = list(continuous = panel_hist),
                      lower = list(continuous = panel_contour),
                      progress = FALSE) + 
  theme(legend.position = "none",
        strip.text = element_text(size = 9),
        panel.grid = element_blank())

pairs_plot


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                     FIT HMC MODEL FOR PARAMETER SAMPLING                    #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# --------------------------------------------------- #
#           LOAD SAVED GENERATED DATA SETS            #
# --------------------------------------------------- #
# n_i and n_j define the model structure so they are kept fixed.
make_data_list <- function(path, ni=3, nj=2) {
  dat <- readRDS(path)
  N <- nrow(dat)
  n_dim <- ncol(dat)
  data_list <- list(N = N, ni = ni, nj = nj, n_dim = n_dim, X = dat)
}

# DATA SET 1
data_list1 <- make_data_list("data/hrb_gen_data1.RData")
# DATA SET 2
data_list2 <- make_data_list("data/hrb_gen_data2.RData")

# --------------------------------------------------- #
#                     FIT STAN MODEL                  #
# --------------------------------------------------- #
# Compile Stan program
hrb_mod <- stan_model("stan/rosenbrock/hrb.stan")

# Fit the models with Stan
nuts_controls <- list(max_treedepth = 10, adapt_delta = 0.80)

stan_sampling <- function(data_list) {
  sampling(hrb_mod, data = data_list, seed = 42, iter = 2e3,
           chains = 4, refresh = 500, control = nuts_controls)
}

mod1 <- stan_sampling(data_list1)
mod2 <- stan_sampling(data_list2)

# --------------------------------------------------- #
#               SAVE FITTED STAN MODELS               #
# --------------------------------------------------- #
save(mod1, file = "results/stanfit_mod1.RData")
save(mod2, file = "results/stanfit_mod2.RData")

# --------------------------------------------------- #
#                     DIAGNOSTICS                     #
# --------------------------------------------------- #
util$check_all_diagnostics(mod2)
mod2

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                               EVALUATION                                    #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# --------------------------------------------------- #
#                         RMSE                        #
# --------------------------------------------------- #
draws <- as_draws_df(mod2)
truth <- c(mu = 1, a = 0.5, b = 5)
sq_err <- sweep(draws[, c("mu", "a", "b")], 2, truth, "-")^2
rmse_point <- colMeans(sq_err) %>% sqrt()

# Find (2.5%, 97.5%) quantiles for each parameter.
# Take sqrt() to get back squared errors to parameter scale.
rmse_ci <- apply(sq_err, 2, function(se) {
  quantile(sqrt(se), probs = c(0.025, 0.975))
})

# df for plotting
df <- data.frame(
  rmse = rmse_point,
  param = names(truth),
  lower = rmse_ci[1, ],
  upper = rmse_ci[2, ]
)

# RMSE scatter with error bars
ggplot(data = df, aes(x= param, y = rmse)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.2)) +
  coord_cartesian(ylim = c(0, 1e-1)) +
  labs(y="RMSE", x="Parameter")






