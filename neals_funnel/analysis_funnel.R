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

calc_metrics <- function() {
  
  mypath <- getwd()
  fit_path <- file.path(mypath, "stanfits")
  fit_files <- list.files(fit_path, pattern = paste0("funnel_fit", ".*\\.rds$"),
                          full.names = TRUE)
  stanfits <- lapply(fit_files, readRDS)
  K <- length(stanfits)
  
  # empty result-list of length K
  results <- vector("list", K)
  
  for (i in seq_len(K)) {
    
    draws <- stanfits[[i]]$fit %>% as_draws_df()
    true_v <- stanfits[[i]]$true_v
    params <- "v"
    
    # RMSE
    sq_err <- sweep(draws[, params], 2, true_v, "-")^2
    rmse_point <- sqrt(colMeans(sq_err))
    # Take sqrt() to get back squared errors to "error scale".
    rmse_abs_q <- apply(sq_err, 2, function(se) {
      quantile(sqrt(se), probs = c(0.025, 0.975))
    })
    
    # RHAT
    rhats <- summary(stanfits[[i]]$fit)$summary[, "Rhat"]
    
    # POSTERIOR QUANTITIES
    post_v <- draws$v
    post_mcse <- mcse_mean(post_v)
    post_mean <- mean(post_v)
    post_sd <- sd(post_v)
    
    # Rounding
    count_leading_zeros <- function(x) {
      s <- format(x, scientific = FALSE, trim = TRUE)
      if (!grepl("\\.", s)) return(0)  # no decimal point
      digits <- sub("^[^.]*\\.", "", s)           # remove everything before "."
      zeros <- regexpr("[1-9]", digits) - 1       # count zeros before first nonzero
      ifelse(zeros < 0, 0, zeros)
    }
    digits <- count_leading_zeros(post_mcse)
    post_mean <- round(post_mean, digits)
    post_sd <- round(post_sd, digits)
    
    # OUTPUT LIST OF RESULTS
    results[[i]] <- list(
      rmse_point = rmse_point,
      rmse_abs_q = rmse_abs_q,
      rhats = rhats,
      true_v = stanfits[[i]]$true_v,
      P = stanfits[[i]]$P,
      post_mean = post_mean,
      post_sd = post_sd,
      post_mcse = post_mcse
    )
    # End of for loop
  }
  
  return(results)
  # End of function
}


metrics <- calc_metrics()


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                                                        #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
dfit1 <- readRDS(file.path(getwd(), "stanfits", "funnel_fit_P=1.rds"))
draws <- dfit1$fit %>% as_draws_df()
dfit1$fit

mcse_mean(draws$v)
mcse_quantile(draws$v, probs = c(0.05, 0.95))

summ <- data.frame(
  mcse_mean = mcse_mean(draws$v),
  q05 = quantile(draws$v, 0.05),
  mcse_q05 = mcse_quantile(draws$v, probs = 0.05),
  q95 = quantile(draws$v, 0.95),
  mcse_q95 = mcse_quantile(draws$v, probs = 0.95)
)



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                           ***********                                       #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
d3 <- readRDS(file.path(getwd(), "data", "funnel_gen_data_P=3.rds"))
d9 <- readRDS(file.path(getwd(), "data", "funnel_gen_data_P=9.rds"))
d200 <- readRDS(file.path(getwd(), "data", "funnel_gen_data_P=200.rds"))
d100 <- readRDS(file.path(getwd(), "data", "funnel_gen_data_P=100.rds"))



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


###
# X is your N x P matrix from the dataset you fit
d <- read_rds(file.path(getwd(), "data", "funnel_gen_data_P=12.rds"))
X <- d$data

s2 <- mean(X^2)
v_mle <- log(s2)

post <- draws$v
c(mean_post = mean(post),
  sd_post   = sd(post),
  v_mle     = v_mle,
  approx_sd = sqrt(2 / (nrow(X)*ncol(X))))

# ----------------------------------------- #
###            RMSE errorbar              ###
# ----------------------------------------- #
plot_df <- matrix(NA, nrow = length(metrics), ncol = 4) %>% as.data.frame()
colnames(plot_df) <- c("rmse_point", "rmse_lower", "rmse_upper", "P")

for (i in 1:length(metrics)) {
  plot_df[i, 1] <- metrics[[i]]$rmse_point
  plot_df[i, 2] <- metrics[[i]]$rmse_abs_q[1, ]
  plot_df[i, 3] <- metrics[[i]]$rmse_abs_q[2, ]
  plot_df[i, 4] <- metrics[[i]]$P
}

# Correct order of P
plot_df <- plot_df %>%
  arrange(P)

ggplot(data = plot_df, aes(x = rmse_point, y = P)) +
  geom_point(size = 2, shape = 15, color = cbbPalette[6]) +
  geom_errorbar(aes(xmin = rmse_lower, xmax = rmse_upper), width = 0.3) +
  scale_y_continuous(breaks = seq(3,27,3)) +
  # scale_x_log10() +
  labs(x = "RMSE (point)", y = "Dimensionality of X",
       subtitle = "Bars show 5–95% posterior interval of |v − v_true| within each datasetr") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("plots/rmse_errorbar.pdf", width = 6, height = 5)








