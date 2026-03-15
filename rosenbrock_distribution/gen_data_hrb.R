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
data_gen_fun <- function(stan_file, nj, ni, configs) {
  
  # COLNAMES FUNCTION
  make_colnames <- function(nj, ni) {
    cols <- "x1"
    for (j in 1:nj) {
      for (i in 2:ni) {
        cols <- c(cols, paste0("x", j, i))
      }
    }
    return(cols)
  }
  
  K <- nrow(configs)
  mod_sim <- stan_model(file = stan_file)
  
  # SIMULATE ONE DATA SET AT A TIME
  for (i in seq_len(K)) {
    
    pars <- configs[i, ]
    
    # --------------------------------------------------- #
    #             RUN STAN DATA GENERATION                #
    # --------------------------------------------------- #
    data_list <- list(
      "N" = 1e4,
      "nj" = nj,
      "ni" = ni,
      "mu"= pars$mu,
      "a" = pars$a,
      "b" = pars$b
    )
    
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
    sim_df <- as.data.frame(X_mat)
    colnames(sim_df) <- make_colnames(nj, ni)
    
    out <- list(
      data = sim_df,
      config = pars
    )
    
    save_name <- paste0("data/hrb_gen_data_", "nj=", nj, "_ni=", ni, "_set", i, ".rds")
    saveRDS(out, file = save_name)
  }
  # End of function
}


# --------------------------------------------------- #
#                         RUN                         #
# --------------------------------------------------- #
cnfgs <- data.frame(
  mu = c(-4, 0, 10, 1, 1, 1, 1, 1, 1),
  a = c(0.05, 0.05, 0.05, 0.5, 0.05, 5e-3, 0.05, 0.05, 0.05),
  b = c(5, 5, 5, 5, 5, 5, 50, 0.05, 5e-4)
)

data_gen_fun(stan_file = "stan/simulate-data_hrb.stan",
             nj = 2,
             ni = 3,
             configs = cnfgs)


###
d <- read_rds(file.path(getwd(), "data", "hrb_gen_data_nj=2_ni=3_set2.rds"))
d$data

