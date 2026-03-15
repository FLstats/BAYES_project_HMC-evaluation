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
#                               DATA PREP                                     #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
dat <- read.csv(file = "data/german_raw.csv", sep = ",", header = T)

# Define response
y <- dat[, length(dat)] %>% {ifelse(.=="good", 1, 0)}

# Define predictor matrix
X <- dat[, -c(1, length(dat))] 

# Index vector of complete case observations
cc <- complete.cases(X)

# Build model matrix. Factors become dummies with k-1 columns 
# (the alphabetically first level gets dropped and set as baseline).
# The function drops NAs.
# ----------------------------------------- #
###           only main effects           ###
# ----------------------------------------- #
# X <- model.matrix(~ . - 1, data = X)

# ----------------------------------------- #
###     and all two-way interactions      ###
# ----------------------------------------- #
X <- model.matrix(~ (.)^2 - 1, data = X)

# ----------------------------------------- #
###           Grouping variables          ###
# ----------------------------------------- #
g_levels <- sort(unique(X[, "Job"]))
f_g <- factor(X[, "Job"], levels = g_levels)
g_id <- as.integer(f_g)
N_g <- length(unique(g_id))

# Remove all variables containing Job
X <- as_tibble(X) %>% select(-contains("Job")) %>% as.matrix()

# Remove NAs
y <- y[cc]


# ----------------------------------------- #
###          Center and scale             ###
# ----------------------------------------- #
is_binary <- apply(X, 2, function(col) all(col %in% c(0,1)))
X_num <- X[, !is_binary]
X_bin <- X[,  is_binary]

X_num <- scale(X_num)
X <- cbind(X_num, X_bin)


# ----------------------------------------- #
###           Stan data list              ###
# ----------------------------------------- #
N <- nrow(X)
P <- ncol(X)

data_list <- list(N = N, P = P,
                  N_g = N_g, g_id = g_id,
                  y = as.integer(y), X = X)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                 Fit multilevel hlr, alpha(cp), beta(cp)                     #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
nuts_controls <- list(max_treedepth = 10, adapt_delta = 0.99)
# cp
# (a_delta = 0.99)
# P=146 --> 73 div, 4 chains low EFMI, low bESS
fit <- stan(file = "stan/hlr_multilevel_cp.stan",
            data = data_list,
            seed = 42,
            chains = 4,
            iter = 2e3,
            refresh = 500,
            control = nuts_controls)

util$check_all_diagnostics(fit)

# saveRDS(fit, file.path(getwd(), "stanfits", "hlr_fit_mlvl_alpha=cp_beta=cp_P=19.rds"))
saveRDS(fit, file.path(getwd(), "stanfits", "hlr_fit_mlvl_alpha=cp_beta=cp_P=146.rds"))


# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# #                                                                             #
# #                 Fit multilevel hlr, alpha(ncp), beta(cp)                    #
# #                                                                             #
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# nuts_controls <- list(max_treedepth = 10, adapt_delta = 0.99)
# # (a_delta = 0.80)
# # P=19 --> 77 div
# # P=146 --> 7 div
# 
# # (a_delta = 0.99)
# # P=19 --> 0 div
# # P=146 --> 0 div, low EFMI
# fit <- stan(file = "stan/hlr_multilevel_alpha=ncp_beta=cp.stan",
#             data = data_list,
#             seed = 42,
#             chains = 4,
#             iter = 2e3,
#             refresh = 500,
#             control = nuts_controls)
# 
# util$check_all_diagnostics(fit)
# 
# # saveRDS(fit, file.path(getwd(), "stanfits", "hlr_fit_mlvl_alpha=ncp_beta=cp_P=19.rds"))
# # saveRDS(fit, file.path(getwd(), "stanfits", "hlr_fit_mlvl_alpha=ncp_beta=cp_P=146.rds"))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                 Fit multilevel hlr, alpha(ncp), beta(ncp)                   #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
nuts_controls <- list(max_treedepth = 10, adapt_delta = 0.99)
# (a_delta = 0.80)
# P=146 --> 9 divergences

# (a_delta = 0.99)
# P=146 --> 1 divergences
fit <- stan(file = "stan/hlr_multilevel_alpha=ncp_beta=ncp.stan",
            data = data_list,
            seed = 42,
            chains = 4,
            iter = 2e3,
            refresh = 500,
            control = nuts_controls)

util$check_all_diagnostics(fit)

# saveRDS(fit, file.path(getwd(), "stanfits", "hlr_fit_mlvl_alpha=ncp_beta=ncp_P=19.rds"))
saveRDS(fit, file.path(getwd(), "stanfits", "hlr_fit_mlvl_alpha=ncp_beta=ncp_P=146.rds"))





