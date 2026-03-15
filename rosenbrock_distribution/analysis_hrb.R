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

calc_metrics <- function(nj, ni) {
  
  mypath <- getwd()
  fit_path <- file.path(mypath, "stanfits")
  fit_files <- list.files(fit_path,
                          pattern = paste0("nj=", nj, "_ni=", ni, ".*\\.rds$"),
                          full.names = TRUE)
  stanfits <- lapply(fit_files, readRDS)
  K <- length(stanfits)
  
  # empty result-list of length K
  results <- vector("list", K)
  
  for (i in seq_len(K)) {
    
    draws <- stanfits[[i]]$fit %>% as_draws_df()
    truth <- stanfits[[i]]$config
    params <- names(truth)
    
    # RMSE
    sq_err <- sweep(draws[, params], 2, as.numeric(truth), "-")^2
    rmse_point <- sqrt(colMeans(sq_err))
    # Take sqrt() to get back squared errors to "error scale".
    rmse_abs_q <- apply(sq_err, 2, function(se) {
      quantile(sqrt(se), probs = c(0.025, 0.975))
    })
    
    # RHAT
    rhats <- summary(stanfits[[i]]$fit)$summary[, "Rhat"]
    
    # POSTERIOR QUANTITIES
    mcse_mu <- draws %>% pull(mu) %>% mcse_mean()
    mcse_a <- draws %>% pull(a) %>% mcse_mean()
    mcse_b <- draws %>% pull(b) %>% mcse_mean()
    
    post_mean_mu <- draws %>% pull(mu) %>% mean()
    post_mean_a <- draws %>% pull(a) %>% mean()
    post_mean_b <- draws %>% pull(b) %>% mean()
    
    post_sd_mu <- draws %>% pull(mu) %>% sd()
    post_sd_a <- draws %>% pull(a) %>% sd()
    post_sd_b <- draws %>% pull(b) %>% sd()
    
    # Rounding
    count_leading_zeros <- function(x) {
      s <- format(x, scientific = FALSE, trim = TRUE)
      if (!grepl("\\.", s)) return(0)  # no decimal point
      digits <- sub("^[^.]*\\.", "", s)           # remove everything before "."
      zeros <- regexpr("[1-9]", digits) - 1       # count zeros before first nonzero
      ifelse(zeros < 0, 0, zeros)
    }
    digits_mu <- count_leading_zeros(mcse_mu)
    post_mean_mu <- round(post_mean_mu, digits_mu)
    post_sd_mu <- round(post_sd_mu, digits_mu)
    
    digits_a <- count_leading_zeros(mcse_a)
    post_mean_a <- round(post_mean_a, digits_a)
    post_sd_a <- round(post_sd_a, digits_a)
    
    digits_b <- count_leading_zeros(mcse_b)
    post_mean_b <- round(post_mean_b, digits_b)
    post_sd_b <- round(post_sd_b, digits_b)
    
    # Place in df
    mcse <- data.frame(
      mu = mcse_mu, a = mcse_a, b = mcse_b
    )
    
    posterior_means <- data.frame(
      mu = post_mean_mu, a = post_mean_a, b = post_mean_b
    )
    
    posterior_sd <- data.frame(
      mu = post_sd_mu, a = post_sd_a, b = post_sd_b
    )
    
    # OUTPUT LIST OF RESULTS
    results[[i]] <- list(
      rmse_point = rmse_point,
      rmse_abs_q = rmse_abs_q,
      rhats = rhats,
      true_values = truth,
      mcse = mcse,
      post_mean = posterior_means,
      post_sd = posterior_sd
    )
    # End of for loop
  }
  
  return(results)
  # End of function
}


metrics <- calc_metrics(nj = 2, ni = 3)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                             RMSE ERRORBAR PLOT                              #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
plot_df_mu <- matrix(NA, nrow = length(metrics), ncol = 4) %>% as.data.frame()
colnames(plot_df_mu) <- c("rmse_point", "rmse_lower", "rmse_upper", "set")

for (i in 1:length(metrics)) {
  plot_df_mu[i, 1] <- metrics[[i]]$rmse_point[1]
  plot_df_mu[i, 2] <- metrics[[i]]$rmse_ci[1, 1]
  plot_df_mu[i, 3] <- metrics[[i]]$rmse_ci[2, 1]
  plot_df_mu[i, 4] <- i
}

ggplot(data = plot_df_mu, aes(x = rmse_point, y = set)) +
  geom_point() +
  geom_errorbar(aes(xmin = rmse_lower, xmax = rmse_upper), width = 0.2) +
  scale_y_continuous(breaks = 1:9) +
  scale_x_log10() +
  labs(x = "RMSE", y = "Set") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                           ***********                                       #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
rhats <- lapply(metrics, function(x) x$rhats)
imap(rhats, ~ enframe(.x, name = "parameter", value = "rhat") %>%
       mutate(dataset = .y)) %>%
  bind_rows() %>%
  arrange(desc(rhat)) %>%
  mutate(rhat = format(rhat, digits = 3)) %>%
  print(n = Inf)











