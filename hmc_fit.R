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
#               GENERAL FITTING FUNCTION              #
# --------------------------------------------------- #
fit_fun <- function(model, stan_file) {
  
  # READ ALL GENERATED DATA SETS INTO A LIST
  mypath <- getwd()
  model <- model
  data_path <- file.path(mypath, "data")
  data_files <- list.files(data_path,
                           pattern = model,
                           full.names = TRUE)
  datasets <- lapply(data_files, readRDS)
  K <- length(datasets)
  # K <- length(datasets)
  
  # Compile Stan program
  mod <- stan_model(stan_file)
  
  # FIT EACH DATA SET ONE AT A TIME
  for (i in seq_len(K)) {
    
    X <- datasets[[i]]
    
    if (model == "hrb") {
      N <- nrow(X)
      n_dim <- ncol(X)
      data_list <- list(N=N, ni=3, nj=2, n_dim=n_dim, X=X)
    }
    
    if (model == "funnel") {
      N <- nrow(X)
      P <- ncol(X)
      data_list <- list(N=N, P=P, X=X)
    }
    
    # Fit the models with Stan
    nuts_controls <- list(max_treedepth = 10, adapt_delta = 0.80)
    fit <- sampling(mod, data = data_list, seed = 42, iter = 2e3,
                    chains = 4, refresh = 500, control = nuts_controls)
    
    # Save the fits
    save_name <- paste0("fits/", model,"_fit_mod", i, ".rds")
    saveRDS(fit, file = save_name)
    
  }
  
  # End of function
}

# --------------------------------------------------- #
#                     RUN FUNCTION                    #
# --------------------------------------------------- #
fit_fun(model = "funnel", stan_file = "stan/funnel/neals_funnel.stan")





