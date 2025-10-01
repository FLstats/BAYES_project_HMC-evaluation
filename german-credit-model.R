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
#                   Hierarchical LR (German credit model)                     #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# ----------------------------------------- #
###               Data prep               ###
# ----------------------------------------- #
dat <- read.csv(file = "data/german.csv", sep = ",", header = T)

# Define response
y <- dat[, length(dat)] %>% {ifelse(.=="good", 1, 0)}

# Define predictor matrix
X <- dat[, -c(1, length(dat))] 

# Index vector of complete case observations
cc <- complete.cases(X)

# Build model matrix. Factors become dummies with k-1 columns 
# (the alphabetically first level gets dropped and set as baseline).
# The function drops NAs.
X <- model.matrix(~ . - 1, data = X)

# Check NAs in columns
apply(X, 2, function(col) sum(is.na(col)))

# ----------------------------------------- #
###           Grouping variables          ###
# ----------------------------------------- #
g1_levels <- sort(unique(X[, "Job"]))
f_g1 <- factor(X[, "Job"], levels = g1_levels)
g1_id <- as.integer(f_g1)
N_g1 <- length(unique(g1_id))

X <- as_tibble(X) %>% select(-Job) %>% as.matrix()

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
                  N_g1 = N_g1, g1_id = g1_id,
                  y = as.integer(y), X = X)

nuts_controls <- list(max_treedepth = 10, adapt_delta = 0.99)

# --------------------------------------------------- #
#                   multilevel cp                     #
# --------------------------------------------------- #
fit <- stan(file = "stan/HLR/HLR_multilevel_cp.stan",
            data = data_list,
            seed = 42,
            chains = 4,
            iter = 2e3,
            refresh = 500,
            control = nuts_controls)

util$check_all_diagnostics(fit)

# --------------------------------------------------- #
#           multilevel ncp-alpha, cp-beta             #
# --------------------------------------------------- #
fit <- stan(file = "stan/HLR/HLR_multilevel_ncp-alpha.stan",
            data = data_list,
            seed = 42,
            chains = 4,
            iter = 2e3,
            refresh = 500,
            control = nuts_controls)

util$check_all_diagnostics(fit)

# --------------------------------------------------- #
#           multilevel ncp-alpha, ncp-beta            #
# --------------------------------------------------- #
fit <- stan(file = "stan/HLR/HLR_multilevel_ncp-alpha-beta.stan",
            data = data_list,
            seed = 42,
            chains = 4,
            iter = 2e3,
            refresh = 500,
            control = nuts_controls)

util$check_all_diagnostics(fit)

# ----------------------------------------- #
###               Diagnostics             ###
# ----------------------------------------- #
# Quick summary
fit

# See fitted parameters
names(fit)
fit@sim$fnames_oi

# ----------------------------------------- #
###                 Plots                 ###
# ----------------------------------------- #
mcmc_rank_overlay(fit, pars = c("tau", "sigma_alpha_g1"))

posterior_df <- as_draws_df(fit)
posterior_array <- as_draws_array(fit)

# Get nuts parameters
np <- nuts_params(fit)
lp <- log_posterior(fit)

mcmc_nuts_acceptance(np, lp)
mcmc_nuts_divergence(np, lp)

# Scatter plot
color_scheme_set("darkgray")
div_style <- scatter_style_np(div_color = "green", div_size = 3)

# Set sd param on vertical axis to potentially see funnel shape.
# Actual param and sd should be correlated.
# Raw param and sd should ideally be uncorrelated.
# If we see funnel shape in Raw + sd, something needs to be changed.
mcmc_scatter(posterior_array, 
             pars = c("alpha_g1[1]", "sigma_alpha_g1"), 
             np = np,
             transform = list(sigma_alpha_g1 = "log"), # more interpretable axis
             size=1,
             np_style = div_style)



#







#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                                 SLASK                                       #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# recode to numeric
X <- X %>%
  mutate(
    sex = ifelse(Sex == "male", 1, 0),
    housing = recode(Housing, "own"=1, "rent"=2, "free"=3),
    saving_accounts = recode(Saving.accounts, "little"=1, "moderate"=2,
                             "quite rich"=3, "rich"=4),
    checking_account = recode(Checking.account, "little"=1, "moderate"=2,
                              "rich"=3),
    purpose = recode(Purpose, "business"=1, "car"=2, "domestic appliances"=3,
                     "education"=4, "furniture/equipment"=5, "radio/TV"=6,
                     "repairs"=7, "vacation/others"=8)
  ) %>%
  select(-c(Sex, Housing, Saving.accounts, Checking.account, Purpose)) %>%
  rename("age" = "Age", "job" = "Job")

# ----------------------------------------- #
###       All two way interactions        ###
# ----------------------------------------- #
# Build design matrix, including all two-way interactions, removing intercept
X <- model.matrix(~ (.)^2 - 1, data = X)
# # Center and scale
m <- colMeans(X)
s <- apply(X, 2, sd)
s[s == 0] <- 1  # prevent division by 0
X2 <- X %>%
  sweep(2, m, "-") %>%
  sweep(2, s, "/")
