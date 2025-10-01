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
#                             Neal's funnel                                   #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# --------------------------------------------------- #
#                   GENERATE DATA                     #
# --------------------------------------------------- #
data_list <- list("N" = 1e4, "P" = 9)

mod_nf <- stan_model(file = "stan/funnel/neals-funnel_simulate-data.stan")

fit <- sampling(mod_nf,
                data = data_list,
                algorithm = "Fixed_param",
                seed = 42,
                iter = 1,
                chains = 1)

# --------------------------------------------------- #
#               EXTRACT SIMULATED DATA                #
# --------------------------------------------------- #
X_array <- extract(fit)$X
dim(X_array)
X_mat <- drop(X_array)
dim(X_mat)
head(X_mat)

nf_sim_df <- as.data.frame(X_mat)
colnames(nf_sim_df) <- c("v", paste0("x", seq(1,9,1)))
head(nf_sim_df)

# --------------------------------------------------- #
#                     SAVE DATA                       #
# --------------------------------------------------- #
saveRDS(nf_sim_df, file = "data/nf_gen_data1.RData")

ggplot(data = nf_sim_df, aes(x=x1, y=v)) +
  geom_point()

# --------------------------------------------------- #
#                         *                           #
# --------------------------------------------------- #
X <- readRDS("data/funnel_gen_data1.rds")
data_list <- list(N=nrow(X), P=ncol(X), X=X)
hrb_mod <- stan_model("stan/funnel/neals_funnel.stan")
nuts_controls <- list(max_treedepth = 10, adapt_delta = 0.80)
fit <- sampling(hrb_mod, data = data_list, seed = 42, iter = 2e3,
                chains = 4, refresh = 500, control = nuts_controls)





