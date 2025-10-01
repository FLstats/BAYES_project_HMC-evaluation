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
#                           SIMULATION FUNCTION                               #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
data_gen_fun <- function(configs, model, stan_file) {
  
  if (model == "hrb") {
    nj <- 2
    ni <- 3
    make_colnames <- function(ni, nj) {
      cols <- "x1"
      for (j in 1:nj) {
        for (i in 2:ni) {
          cols <- c(cols, paste0("x", j, i))
        }
      }
      return(cols)
    }
  }
  
  K <- nrow(configs)
  mod_sim <- stan_model(file = stan_file)
  
  # SIMULATE ONE DATA SET AT A TIME
  for (i in seq_len(K)) {
    
    pars <- configs[i, ]
    
    # --------------------------------------------------- #
    #             RUN STAN DATA GENERATION                #
    # --------------------------------------------------- #
    if (model == "hrb") {
      data_list <- list(
        "N" = 1e4,
        "ni" = ni,
        "nj" = nj,
        "mu"= pars$mu,
        "a" = pars$a,
        "b" = pars$b
      )
    }
    
    if (model == "funnel") {
      data_list <- list(
        "N" = 1e4,
        "P" = 9,
        "sigma_v" = pars
      )
    }
    
    sim <- sampling(mod_sim,
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
    if (model == "hrb") {
      colnames(X_mat) <- make_colnames(ni, nj)
    }
    sim_df <- as.data.frame(X_mat)
    
    # --------------------------------------------------- #
    #                 SAVE SIMULATED DATA                 #
    # --------------------------------------------------- #
    save_name <- paste0("data", "/", model, "_gen_data", i, ".rds")
    saveRDS(sim_df, file = save_name)
  }
  
  # End of function
}

# ----------------------------------------- #
###                  hrb                  ###
# ----------------------------------------- #
cnfgs_hrb <- data.frame(
  mu = c(-4, 0, 10, 1, 1, 1, 1, 1, 1),
  a = c(0.05, 0.05, 0.05, 0.05, 0.05, 5e-3, 0.05, 0.05, 0.05),
  b = c(5, 5, 5, 5, 5, 5, 50, 0.05, 5e-4)
)

data_gen_fun(configs = cnfgs_hrb, 
             model = "hrb", 
             stan_file = "stan/rosenbrock/hrb_simulate-data.stan")

# ----------------------------------------- #
###                 funnel                ###
# ----------------------------------------- #
cnfgs_funnel <- data.frame(
  sigma_v = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
)

data_gen_fun(configs = cnfgs_funnel,
             model = "funnel",
             stan_file = "stan/funnel/neals-funnel_simulate-data.stan")


