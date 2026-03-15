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

fit <- readRDS(file.path(getwd(), "stanfits/funnel_fit_P=3.rds"))
draws <- as_draws_df(fit[["fit"]])

# --------------------------------------------------- #
#                   2d funnel plot                    #
# --------------------------------------------------- #
set.seed(0)
S <- 1e4
v <- rnorm(S, 0, 3)
x1 <- rnorm(S, 0, sqrt(exp(v)))
x2 <- rnorm(S, 0, sqrt(exp(v)))

df <- data.frame(
  v = v,
  x1 = x1,
  x2 = x2
)

# Funnel plot
plot_funnel <- function(xi, idx) {
  ggplot(data = df) +
    geom_point(aes(x=xi, y=v), color = "#320000cc") +
    labs(x=bquote(x[.(idx)]),
         y=expression(nu)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 10),
    )  
}

plot_funnel(x1, 1)

ggsave("plots/2d-funnel.png", width = 6, height = 5)

# # Remove outliers if needed
# range(v)
# range(x1)
# which(x1 == min(x1))
# df <- df[-which(x1 == min(x1)), ]

# --------------------------------------------------- #
#                   3d funnel plot                    #
# --------------------------------------------------- #
library(plotly)

plot_ly(
  data = df,
  x = ~x1, y = ~x2, z = ~v,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size=4,
    opacity=0.55,
    color = "rgba(50, 0, 0, 0.8)"
  )
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = list(
        text = "x\u2081",
        font = list(size = 20, family = "Arial", color = "black")
        )
      ),
      yaxis = list(
        title = list(
          text = "x\u2082",
          font = list(size = 20, family = "Arial", color = "black")
          )
      ),
      zaxis = list(
        title = list(
          text = "\u03BD",
          font = list(size = 20, family = "Arial", color = "black")
          )
      )
    )
  )

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#                                                                             #
#                                 SLASK                                       #
#                                                                             #
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

# --------------------------------------------------- #
#    3x3 GRID OF ALL X-VARIABLES FOR ONE V-VALUE      #
# --------------------------------------------------- #
read_model <- function(P, sigma_v) {
  readRDS(paste0(getwd(), "/data/", "funnel_gen_data_","P=", P, "_", "v=",
                 sigma_v, ".rds"))
}

df <- read_model(9, 1)$data
colnames(df) <- c("v", paste0("x", 1:9))
head(df)


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

# --------------------------------------------------- #
#     3x3 GRID OF ONE X-VARIABLE FOR ALL V-VALUES     #
# --------------------------------------------------- #
mypath <- getwd()
data_path <- file.path(mypath, "data")
data_files <- list.files(data_path, pattern = "funnel_gen_data.*\\.rds$", full.names = TRUE)
datasets <- lapply(data_files, readRDS)
K <- length(datasets)


plot_grid <- function(xi) {
  
  plots <- list()
  
  for (i in seq_len(K)) {
    df <- datasets[[i]]$data
    
    plot <- ggplot(df, aes(x = .data[[xi]], y = V)) +
      geom_point(alpha = 0.3, size = 0.6) +
      stat_density_2d(
        aes(fill = after_stat(level)),
        geom = "polygon",
        contour_var = "ndensity",
        bins = 9,
        n = 300,
        h = c(bw.nrd(df[[xi]]), bw.nrd(df$V))
      ) +
      scale_fill_viridis_c(guide = "none") +
      labs(x = xi,
           # title = paste("sigma_v =", sigma_v[i]),
           # title = bquote(sigma[v] == .(sigma_vals[i])),
           y = "v"
           ) +
      theme_bw(base_size = 9) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = 9, hjust = 0.5))
    
    plots[[i]] <- plot
    
  }
  
  wrap_plots(plots, ncol = 3)
  
}

plot_grid(xi = "X1")






