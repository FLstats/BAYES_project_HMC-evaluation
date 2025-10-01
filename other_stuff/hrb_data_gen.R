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

# --------------------------------------------------- #
#                   DEFINE COLNAMES                   #
# --------------------------------------------------- #
make_colnames <- function(ni, nj) {
  cols <- "x1"
  for (j in 1:nj) {
    for (i in 2:ni) {
      cols <- c(cols, paste0("x", j, i))
    }
  }
  cols
}
# --------------------------------------------------- #
#                 COMPILE STAN MODEL                  #
# --------------------------------------------------- #
mod_sim_hbr <- stan_model(file = "stan/rosenbrock/hrb_simulate-data.stan")


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                           SIMULATION FUNCTION                               #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
data_gen_fun <- function(configs) {
  
  nj <- 2
  ni <- 3
  K <- nrow(configs)
  
  for (i in seq_len(K)) {
    
    pars <- configs[i, ]
    
    # --------------------------------------------------- #
    #             RUN STAN DATA GENERATION                #
    # --------------------------------------------------- #
    data_list <- list(
      "N" = 1e4,
      "ni" = ni,
      "nj" = nj,
      "mu"= pars$mu,
      "a" = pars$a,
      "b" = pars$b
    )
    
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
    X_array <- rstan::extract(sim)$X # dim: [1, N, n_dim]
    X_mat <- drop(X_array) # dim: [N, n_dim]
    colnames(X_mat) <- make_colnames(ni, nj)
    hrb_sim_df <- as.data.frame(X_mat)
    
    # --------------------------------------------------- #
    #                 SAVE SIMULATED DATA                 #
    # --------------------------------------------------- #
    save_name <- paste0("data/hrb_gen_data", i, ".rds")
    saveRDS(hrb_sim_df, file = save_name)
  }
  
  # End of function
}

# Define parameters
configs <- data.frame(
  mu = c(-4, 0, 10, 1, 1, 1, 1, 1, 1),
  a = c(0.05, 0.05, 0.05, 0.05, 0.05, 5e-3, 0.05, 0.05, 0.05),
  b = c(5, 5, 5, 5, 5, 5, 50, 0.05, 5e-4)
)

# Run function
data_gen_fun(configs)

