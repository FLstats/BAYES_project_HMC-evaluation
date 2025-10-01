setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

rm(list = ls())  
library(tidyverse)
theme_set(theme_bw())
library(latex2exp)
library(patchwork)
library(grid)

ggplot() + theme_void() + theme(plot.background=element_rect(fill="lightblue1"))
select <- dplyr::select
filter <- dplyr::filter
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                     3. Banana posterior Stan (ncp model)                    #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# ----------------------------------------- #
###               Data prep               ###
# ----------------------------------------- #
data_list <- list("b" = -0.06)

# ----------------------------------------- #
###               Fit model               ###
# ----------------------------------------- #
nuts_controls <- list(max_treedepth = 15, adapt_delta = 0.99)

fit <- stan(file = "stan/3.banana_ncp.stan",
            data = data_list,
            seed = 42,
            chains = 4,
            iter = 2e3, # warmup = iter/2 = 1e3
            refresh = 500,
            control = nuts_controls)

# ----------------------------------------- #
###              Diagnostics              ###
# ----------------------------------------- #
util$check_all_diagnostics(fit)

# ----------------------------------------- #
###               Plots                   ###
# ----------------------------------------- #
draws <- as_draws_array(fit)
summary(draws)

mcmc_scatter(draws, pars = c("x1", "x2"))
mcmc_scatter(draws, pars = c("z1", "z2"))


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                           Banana posterior in R                             #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
library(MASS)
library(mvtnorm)

# --------------------------------------------------- #
#                Uncorrelated normal                  #
# --------------------------------------------------- #
mu <- c(0, 0)
C1 <- matrix(c(100, 0,
               0, 1), nrow = 2, byrow = T)

x1 <- seq(-20, 20, length.out=100)
x2 <- seq(-20, 20, length.out=100)
grid <- expand.grid(x = x1, y = x2)

# contour plot
z <- dmvnorm(grid, mean = mu, sigma = C1)
z <- matrix(z, nrow=length(x1), ncol=length(x2))

contour(x=x1, y=x2, z=z, nlevels = 3, drawlabels = F)


# --------------------------------------------------- #
#               Banana shaped normal                  #
# --------------------------------------------------- #
mu <- c(0, 0)
C1 <- matrix(c(100, 0,
               0, 1), nrow = 2, byrow = T)

x1 <- seq(-20, 20, length.out=100)
x2 <- seq(-20, 20, length.out=100)
grid <- expand.grid(x = x1, y = x2)

# Twisting function
b <- 0.06
grid_banana <- cbind(grid$x, grid$y + b*grid$x^2 - 100*b)

z <- dmvnorm(grid_banana, mean = mu, sigma = C1) %>%
  matrix(nrow = length(x1), ncol = length(x2))

contour(x=x1, y=x2, z=z, nlevels = 3, drawlabels = F)

# --------------------------------------------------- #
#      Banana shaped normal, rotated -90 degrees      #
# --------------------------------------------------- #
mu <- c(0, 0)
C1 <- matrix(c(5, 0,
               0, 100), nrow = 2, byrow = T)

x1 <- seq(-20, 20, length.out=100)
x2 <- seq(-20, 20, length.out=100)
grid <- expand.grid(x = x1, y = x2)

# Twisting function
b <- 0.1
grid_banana <- cbind(grid$x + b*grid$y^2 - 100*b, grid$y)

z <- dmvnorm(grid_banana, mean = mu, sigma = C1) %>%
  matrix(nrow = length(x1), ncol = length(x2))

contour(x=x1, y=x2, z=z, nlevels = 3, drawlabels = F)

# --------------------------------------------------- #
#                 Dripping normal                     #
# --------------------------------------------------- #
mu <- c(0, 0)
C1 <- matrix(c(100, 0,
               0, 1), nrow = 2, byrow = T)

x1 <- seq(-20, 20, length.out=100)
x2 <- seq(-20, 20, length.out=100)
grid <- expand.grid(x = x1, y = x2)

# Twisting function
b <- 0.06
grid_banana <- cbind(grid$x, grid$y + b*exp(grid$x) - 100*b)

z <- dmvnorm(grid_banana, mean = mu, sigma = C1) %>%
  matrix(nrow = length(x1), ncol = length(x2))

contour(x=x1, y=x2, z=z, nlevels = 3, drawlabels = F)