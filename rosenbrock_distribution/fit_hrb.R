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

fit_fun <- function(stan_file, nj, ni) {
  
  # READ ALL GENERATED DATA SETS INTO A LIST
  mypath <- getwd()
  data_path <- file.path(mypath, "data")
  data_files <- list.files(data_path,
                           # .* = "any number (*) of any characters (.)
                           # \\. = "a literal dot (.) escaped with two "\\"
                           # $ = "end of string"
                           pattern = paste0("nj=", nj, "_ni=", ni, ".*\\.rds$"),
                           full.names = TRUE)
  datasets <- lapply(data_files, readRDS)
  K <- length(datasets)
  
  # Compile Stan program
  mod <- stan_model(stan_file)
  
  # FIT EACH DATA SET ONE AT A TIME
  for (i in seq_len(K)) {
    
    X <- datasets[[i]]$data
    cnfg <- datasets[[i]]$config
    
    N <- nrow(X)
    n_dim <- ncol(X)
    data_list <- list(N=N, ni=3, nj=2, n_dim=n_dim, X=X)
    nuts_controls <- list(max_treedepth = 10, adapt_delta = 0.80)
    
    fit <- sampling(mod, data = data_list, seed = 42, iter = 2e3,
                    chains = 4, refresh = 500, control = nuts_controls)
    
    out <- list(
      fit = fit,
      config = cnfg
    )
    
    # Save fitted data file
    save_name <- paste0("stanfits/hrb_fit_nj=", nj, "_ni=", ni,
                        "_set", i, ".rds")
    saveRDS(out, file = save_name)
    
  }
  
  # End of function
}

fit_fun(stan_file = "stan/hrb.stan", nj = 2, ni = 3)





