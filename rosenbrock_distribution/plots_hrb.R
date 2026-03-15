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
read_data <- function(nj, ni, set) {
  readRDS(paste0(getwd(), "/data", "/hrb_gen_data_nj=", nj,
                 "_ni=", ni, "_set", set, ".rds"))
}

nj <- 2; ni <- 3; set <- 1
df <- read_data(nj, ni, set)$data
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
      geom_point(alpha = 0.5, size = 0.3) +
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

# Save plot
ggsave(paste0("plots/pairsplot_", "nj=", nj, "_ni=", ni, "_set", set, ".pdf"),
       width = 6, height = 5)


# --------------------------------------------------- #
#                       2D PLOT                       #
# --------------------------------------------------- #
ggplot(df, aes(x=x1, y=x12)) +
  geom_point(alpha = 0.5, size = 0.3) +
  # stat_density2d(aes(fill = after_stat(level)),
  #                geom = "polygon",
  #                contour_var = "ndensity",
  #                bins = 9, # number of contours
  #                n = 300, # higher grid resolution
  #                h = c(bw.nrd(df$x1), 
  #                      bw.nrd(df$x12)) # gentler smoothing
  # ) +
  scale_fill_viridis_c(guide = "none") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank())

ggsave("plots/2d_(x1,x12)_set1.pdf", width = 6, height = 5)

# --------------------------------------------------- #
#                       3D PLOT                       #
# --------------------------------------------------- #
library(plotly)

plot_ly(
  data = df,
  x = ~x1, y = ~x12, z = ~x13,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size=4,
    opacity=0.55,
    color = "rgba(50, 0, 0, 0.8)"
  )
)

ggsave("plots/3d_(x1,x12,x13)_set1.pdf", width = 6, height = 5)











