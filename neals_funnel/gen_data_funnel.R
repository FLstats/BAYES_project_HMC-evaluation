setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

rm(list = ls())  
library(tidyverse)
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

library(tidybayes)
library(loo)
library(posterior)
options(posterior.num_args=list(digits=2))
library(bayesplot)
theme_set(theme_bw(base_family = "serif", base_size = 14))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                           SIMULATION FUNCTION                               #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
data_gen_fun <- function(stan_file, configs) {
  
  K <- nrow(configs)
  mod_sim <- stan_model(file = stan_file)
  
  # SIMULATE ONE DATA SET AT A TIME
  for (i in seq_len(K)) {
    
    P <- configs[i, ]
    
    # --------------------------------------------------- #
    #             RUN STAN DATA GENERATION                #
    # --------------------------------------------------- #
    data_list <- list(
      "N" = 1e4,
      "P" = P
      )
    
    sim <- sampling(mod_sim,
                    data = data_list,
                    algorithm = "Fixed_param", # Stan only runs generated quantities
                    seed = 42+i,
                    iter = 1, # generates one data set (iter*chains)
                    chains = 1
    )
    
    # --------------------------------------------------- #
    #               EXTRACT SIMULATED DATA                #
    # --------------------------------------------------- #
    true_v <- rstan::extract(sim)$v
    X_array <- rstan::extract(sim)$X # dim: [1, N, P]
    X_mat <- drop(X_array) # dim: [N, P]
    sim_df <- as.data.frame(X_mat)
    colnames(sim_df) <- paste0("x", rep(1:P))
    
    out <- list(
      "data" = sim_df,
      "true_v" = true_v,
      "P" = P
    )
    
    save_name <- paste0("data/funnel_gen_data_", "P=", P, ".rds")
    saveRDS(out, file = save_name)
  }
  # End of function
}


# --------------------------------------------------- #
#                         RUN                         #
# --------------------------------------------------- #
cnfgs <- data.frame(
  # P = c(3, 6, 9, 12, 15, 18, 21, 24, 27)
  P = c(1, 3, 9, 20, 50, 75, 100, 150, 200)
)

data_gen_fun(stan_file = "stan/simulate-data_funnel.stan",
             configs = cnfgs)



###
d <- read_rds(file.path(getwd(), "data", "funnel_gen_data_P=9.rds"))
d$data











