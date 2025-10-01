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

# --------------------------------------------------- #
#                 READ GENERATED DATA                 #
# --------------------------------------------------- #
read_model <- function(model) {
  readRDS(paste0(getwd(), "/data", "/hrb_gen_data", model, ".rds"))
}

df <- read_model(6)
head(df)

# --------------------------------------------------- #
#                     PAIRS PLOT                      #
# --------------------------------------------------- #
pairs_plot <- function(data) {
  
  library(GGally)
  library(rlang)
  
  # Diagonal: histogram
  panel_hist <- function(data, mapping, ...) {
    ggplot(data, mapping) +
      geom_histogram(color="black", fill=cbbPalette[1],
                     bins = 80) +
      theme_bw(base_size = 9) +
      theme(panel.grid = element_blank())
  }
  
  # Lower triangle: filled density contours
  panel_contour <- function(data, mapping, ...) {
    ggplot(data, mapping) +
      geom_point(alpha = 0.1, size = 0.6) +
      stat_density2d(aes(fill = after_stat(level)),
                     geom = "polygon",
                     contour_var = "ndensity",
                     bins = 9, # number of contours
                     n = 300, # higher grid resolution
                     h = c(bw.nrd(data[[as_name(mapping$x)]]), 
                           bw.nrd(data[[as_name(mapping$y)]])) # gentler smoothing
      ) +
      scale_fill_viridis_c(guide = "none") +
      theme_bw(base_size = 9) +
      theme(panel.grid = element_blank())
  }
  
  # Build the plot
  pairs_plot <- ggpairs(data,
                        columns = 1:ncol(data),
                        upper = "blank",
                        diag  = list(continuous = panel_hist),
                        lower = list(continuous = panel_contour),
                        progress = FALSE) + 
    theme(legend.position = "none",
          strip.text = element_text(size = 9),
          panel.grid = element_blank())
  
  return(pairs_plot)
}

# CREATE PAIRS PLOT
pairs_plot(data = df)

