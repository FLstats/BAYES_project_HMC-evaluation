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

calc_metrics <- function(P) {
  
  mypath <- getwd()
  fit_path <- file.path(mypath, "stanfits")
  fit_files <- list.files(fit_path, pattern = paste0("P=", P, ".*\\.rds$"),
                          full.names = TRUE)
  stanfits <- lapply(fit_files, readRDS)
  K <- length(stanfits)
  
  # empty result-list of length K
  results <- vector("list", K)
  
  for (i in seq_len(K)) {
    
    draws <- stanfits[[i]]$fit %>% as_draws_df()
    truth <- stanfits[[i]]$config
    params <- "sigma_v"
    
    # RMSE
    sq_err <- sweep(draws[, params], 2, as.numeric(truth), "-")^2
    rmse_point <- sqrt(colMeans(sq_err))
    # Take sqrt() to get back squared errors to "error scale".
    rmse_ci <- apply(sq_err, 2, function(se) {
      quantile(sqrt(se), probs = c(0.025, 0.975))
    })
    
    # RHAT
    rhats <- summary(stanfits[[i]]$fit)$summary[, "Rhat"]
    
    # OUTPUT LIST OF RESULTS
    results[[i]] <- list(
      rmse_point = rmse_point,
      rmse_ci = rmse_ci,
      rhats = rhats,
      config = truth
    )
    # End of for loop
  }
  
  return(results)
  # End of function
}


metrics <- calc_metrics(P = 9)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                             RMSE ERRORBAR PLOT                              #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
plot_df <- matrix(NA, nrow = length(metrics), ncol = 4) %>% as.data.frame()
colnames(plot_df) <- c("rmse_point", "rmse_lower", "rmse_upper", "sigma_v")

for (i in 1:length(metrics)) {
  plot_df[i, 1] <- metrics[[i]]$rmse_point
  plot_df[i, 2] <- metrics[[i]]$rmse_ci[1, ]
  plot_df[i, 3] <- metrics[[i]]$rmse_ci[2, ]
  plot_df[i, 4] <- metrics[[i]]$config
}

ggplot(data = plot_df, aes(x = rmse_point, y = sigma_v)) +
  geom_point(size = 2, shape = 15, color = cbbPalette[6]) +
  geom_errorbar(aes(xmin = rmse_lower, xmax = rmse_upper), width = 0.3) +
  scale_y_continuous(breaks = 1:9) +
  scale_x_log10(breaks = c(0.1, 1, 10, 30)) +
  labs(x = "RMSE", y = bquote(sigma[nu])) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("plots/rmse_errorbar.pdf", width = 6, height = 5)





#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                                   SLASK                                     #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
### Make global sd (does not work)
dfs <- lapply(stanfits, function(x) as_draws_df(x$fit))
all_draws <- do.call(rbind, dfs)
sd_global <- sd(all_draws$sigma_v)

### Largest Rhats
# imap() allplies a function all elements of the rhat-list
# where each element is a vector
#   .x = contents of each list element
#   .y = name of each list element
# enframe() converts the rhat-vector to a data frame
rhats <- lapply(metrics, function(x) x$rhats)
imap(rhats, ~ enframe(.x, name = "parameter", value = "rhat") %>%
       mutate(dataset = .y)) %>%
  bind_rows() %>%
  arrange(desc(rhat)) %>%
  mutate(rhat = format(rhat, digits = 3)) %>%
  print(n = Inf)













