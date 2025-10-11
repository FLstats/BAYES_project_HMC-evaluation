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

# ONLY P = 9

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#               3x3 GRID OF ALL X-VARIABLES FOR ONE V-VALUE                   #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
read_model <- function(P, sigma_v) {
  readRDS(paste0(getwd(), "/data/", "funnel_gen_data_","P=", P, "_", "v=",
                 sigma_v, ".rds"))
}

df <- read_model(9, 1)$data
colnames(df) <- c("v", paste0("x", 1:9))
head(df)

# --------------------------------------------------- #
#               GRID OF SCATTER PLOTS                 #
# --------------------------------------------------- #
plots <- lapply(paste0("x", 1:9), function(xi) {
  ggplot(df, aes(x = .data[[xi]], y = v)) +
    geom_point(alpha = 0.3, size = 0.6) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom = "polygon",
      contour_var = "ndensity",
      bins = 9,
      n = 300,
      h = c(bw.nrd(df[[xi]]), bw.nrd(df$v))
    ) +
    scale_fill_viridis_c(guide = "none") +
    labs(x = xi, y = "v") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 9, hjust = 0.5))
})

wrap_plots(plots, ncol = 3)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#               3x3 GRID OF ONE X-VARIABLE FOR ALL V-VALUES                   #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
mypath <- getwd()
data_path <- file.path(mypath, "data")
data_files <- list.files(data_path, pattern = "P=9_.*\\.rds$", full.names = TRUE)
datasets <- lapply(data_files, readRDS)
K <- length(datasets)


plot_grid <- function(xi, sigma_vals) {
  
  plots <- list()
  
  for (i in seq_len(K)) {
    df <- datasets[[i]]$data
    colnames(df) <- c("v", paste0("x", 1:9))
    
    plot <- ggplot(df, aes(x = .data[[xi]], y = v)) +
      geom_point(alpha = 0.3, size = 0.6) +
      stat_density_2d(
        aes(fill = after_stat(level)),
        geom = "polygon",
        contour_var = "ndensity",
        bins = 9,
        n = 300,
        h = c(bw.nrd(df[[xi]]), bw.nrd(df$v))
      ) +
      scale_fill_viridis_c(guide = "none") +
      labs(x = xi, y = "v",
           # title = paste("sigma_v =", sigma_v[i]),
           title = bquote(sigma[v] == .(sigma_vals[i]))) +
      theme_bw(base_size = 9) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = 9, hjust = 0.5))
    
    plots[[i]] <- plot
    
  }
  
  wrap_plots(plots, ncol = 3)
  
}

plot_grid(xi = "x1", sigma_vals = 1:9)






