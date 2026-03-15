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

fit_fun <- function(stan_file) {
  
  # READ ALL GENERATED DATA SETS INTO A LIST
  mypath <- getwd()
  data_path <- file.path(mypath, "data")
  data_files <- list.files(data_path,
                           # .* = "any number (*) of any characters (.)
                           # \\. = "a literal dot (.) escaped with two "\\"
                           # $ = "end of string"
                           pattern = paste0("funnel_gen_data", ".*\\.rds$"),
                           full.names = TRUE)
  datasets <- lapply(data_files, readRDS)
  K <- length(datasets)
  
  # Compile Stan program
  mod <- stan_model(stan_file)
  
  # FIT EACH DATA SET ONE AT A TIME
  for (i in seq_len(K)) {
    
    X <- datasets[[i]]$data
    true_v <- datasets[[i]]$true_v
    P <- datasets[[i]]$P
    
    N <- nrow(X)
    data_list <- list(N=N, P=P, X=X)
    nuts_controls <- list(max_treedepth = 10, adapt_delta = 0.99)
    
    fit <- sampling(mod, data = data_list, seed = 42, iter = 2e3,
                    chains = 4, refresh = 500, control = nuts_controls)
    
    out <- list(
      "fit" = fit,
      "true_v" = true_v,
      "P" = P
    )
    
    # Save fitted data file
    save_name <- paste0("stanfits/funnel_fit_P=", P, ".rds")
    saveRDS(out, file = save_name)
    
  }
  
  # End of function
}

fit_fun(stan_file = "stan/funnel.stan")



