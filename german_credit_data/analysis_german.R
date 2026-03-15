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

fit <- readRDS(file.path(getwd(), "stanfits", "hlr_fit_mlvl_alpha=cp_beta=cp_P=146.rds"))
fit <- readRDS(file.path(getwd(), "stanfits", "hlr_fit_mlvl_alpha=ncp_beta=ncp_P=146.rds"))

draws <- posterior::as_draws_df(fit)



# --------------------------------------------------- #
#       Caterpillar plot of OR median for betas       #
# --------------------------------------------------- #
# (top 20 distances from 1)

### Summarize betas (assumes Stan name beta[1], beta[2], ...)
# Long-format posterior draws.
# 4 chains x 1000 iterations x 146 beta[i] = 584 000 rows.
# Compute with odds ratios per coef and draw ("given this sampled beta[i], 
# the odds multiply by this factor per +1 increase in x").
# Summarise the draws per coef --> end up with P rows.
# e.g. ORmed is the median of the OR for the draws
betas <- draws %>%
  gather_draws(beta[i]) %>%
  mutate(OR = exp(.value)) %>%
  group_by(i) %>%
  summarise(
    med   = median(.value),
    lo95  = quantile(.value, .025),
    hi95  = quantile(.value, .975),
    ORmed = median(OR),
    ORlo  = quantile(OR, .025),
    ORhi  = quantile(OR, .975),
    pr_pos = mean(.value > 0)
  ) %>%
  mutate(delta = abs(ORmed - 1)) %>% # absolute distance from 1
  arrange(desc(delta)) # sort descending

# Keep top 20 for the main figure
betas_top <- betas %>% 
  slice_head(n = 20) %>%
  mutate(i = reorder(as.factor(i), delta))

# Caterpillar plot
ggplot(betas_top, aes(x = i, y = ORmed, ymin = ORlo, ymax = ORhi)) +
  geom_hline(yintercept = 1, linetype = 2, lwd = 0.5, color = "darkgrey") +
  # geom_pointrange(size = 0.5, color = cbbPalette[6]) +
  geom_errorbar(width = 0.3) +
  geom_point(size = 2, color = cbbPalette[6]) +
  coord_flip() +
  labs(x = "Beta coefficient index",
       y = "Median odds ratio") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 12)
        )

ggsave("plots/cpplot_betas_med-or.pdf", width = 6, height = 5)

### INSPECT ONE PARAMETER
draws %>% pull("beta[3]") %>% exp() %>% median()

# --------------------------------------------------- #
#         Caterpillar plot of alpha_g effects         #
# --------------------------------------------------- #
### Group effects (random intercepts), e.g. alpha[g]
# 1000 iter * 4 chains * 4 groups = 16 000
# Interpretation: baseline odds of Y=1 for each group when all predictors
# are at their mean levels
alpha_groups <- draws %>%
  gather_draws(alpha_g[g]) %>%
  mutate(OR = exp(.value)) %>%
  group_by(g) %>%
  median_qi(OR, .width = 0.95)

alpha_groups %>%
  mutate(g = reorder(g, OR)) %>%
  ggplot(aes(x=OR, y=g)) +
  geom_point(size = 2, color = cbbPalette[6]) +
  geom_errorbar(aes(xmin=.lower, xmax=.upper), width = 0.1) +
  labs(x = "Group-intercept effect (OR)", y = "Group")



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                               SCATTER PLOTS                                 #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# store draws
draws <- posterior::as_draws_df(fit)
posterior_array <- as_draws_array(fit)

# check variables names in the array
variables(posterior_array)

# Get nuts parameters
np <- nuts_params(fit)
lp <- log_posterior(fit)

color_scheme_set("darkgray")
div_style <- scatter_style_np(div_color = cbbPalette[7],
                              div_size = 1.5)
# ----------------------------------------- #
###           SCATTER PLOT GRID           ###
# ----------------------------------------- #
library(bayesplot)
library(cowplot)

# parameters to plot
g_idx <- 1:4

# COMMON LIMITS AND BREAKS BETWEEN PLOTS
draws %>% select(starts_with("alpha_g[")) %>% unlist(use.names = F) %>% range()
draws %>% select("sigma2_alpha_g") %>% log() %>% range()

xlim <- c(-2.2, 2)
xbreaks <- seq(-2, 2, by=1)
ylim <- c(-10.1, 4.5)
ybreaks <- seq(-10, 4, by=2)

# create one scatter plot per group
plots <- lapply(g_idx, function(i) {
  mcmc_scatter(
    posterior_array,
    pars = c(glue("alpha_g[{i}]"), "sigma2_alpha_g"),
    np = np,
    transform = list(sigma2_alpha_g = "log"),
    size = 0.7,
    np_style = div_style
  ) +
    labs(
      # x = bquote(alpha[g[.(i)]]),
      # x = bquote(alpha[g*.(i)]),
      # x = bquote(alpha[g[.(as.character(i))]]),
      # x = bquote(alpha * plain("_g[")(.(i)) * plain("]")),
      x = bquote(alpha[.(i)]),
      y = expression(log(sigma[alpha]^2))
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    scale_x_continuous(breaks = xbreaks) +
    scale_y_continuous(breaks = ybreaks) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 12)
    )
})

# combine into 2x2 grid (ncp)
plot_grid(plotlist = plots, ncol = 2, align = "hv")

ggsave("plots/scatter_alpha_g_ncp.pdf", width = 6, height = 5)
# ggsave("plots/scatter_alpha_g_cp.pdf", width = 6, height = 5)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                                   RHAT                                      #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
rhats <- rhat(fit)
summary(rhats)
rhats[order(rhats, decreasing = T)[1:10]] # 10 largest rhats
rhats[rhats > 1.01]





#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                                 SLASK                                       #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# ----------------------------------------- #
###     Histogram of log_sigma_alpha_g    ###
# ----------------------------------------- #
posterior_df %>%
  mutate(log_sigma_alpha_g = log(sigma_alpha_g)) %>%
  mcmc_hist(., pars = "log_sigma_alpha_g")

# ----------------------------------------- #
###             SCATTER PLOT              ###
# ----------------------------------------- #
# Set sd param on vertical axis to potentially see funnel shape.
# Actual param and sd should be correlated.
# Raw param and sd should ideally be uncorrelated.
# If we see funnel shape in Raw + sd, something needs to be changed.
mcmc_scatter(posterior_array, 
             pars = c("alpha_g[4]", "sigma_alpha_g"), 
             np = np,
             transform = list(sigma_alpha_g = "log"), # more interpretable axis
             size=1,
             np_style = div_style)
# ----------------------------------------- #
###               RANDOM PLOTS            ###
# ----------------------------------------- #
mcmc_nuts_acceptance(np, lp)
mcmc_nuts_divergence(np, lp)


# --------------------------------------------------- #
#       Marginal effect for one covariate xj          #
# --------------------------------------------------- #
# Construct a small design grid Xnew varying xj, others at typical values
# (adapt to your design naming)
x_seq <- seq(from = -2, to = 2, length.out = 100)
# example: use posterior means for other covariates = 0; add a group if needed
# eta = alpha + beta_j * x_seq  (extend to full X as appropriate)
eta_draws <- draws %>%
  spread_draws(alpha_g[g], beta[i]) %>%
  filter(g == 1, i == 1) %>%
  slice_sample(n = 1000) %>%
  tidyr::crossing(x = x_seq) %>%
  mutate(eta = alpha_g + beta * x,  # pick your j
         p = plogis(eta))

marg <- eta_draws %>%
  group_by(x) %>%
  summarise(med = median(p), lo = quantile(p, .025), hi = quantile(p, .975))

ggplot(marg, aes(x, med, ymin = lo, ymax = hi)) +
  geom_ribbon(alpha = .15) + geom_line() +
  labs(y = "Predicted probability", x = "x_j (standardized)")


