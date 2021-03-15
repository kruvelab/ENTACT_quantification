library(tidyverse)
library(extrafont)
loadfonts(device = "win")

font <- choose_font("Raleway")
fontsize <- 16
basecolor <- "#14213d" 

my_theme <-   theme(
  plot.background = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(size = 0.5,
                           color = basecolor),
  text = element_text(family = font,
                      size = fontsize,
                      color = basecolor),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text = element_blank(),
  legend.position = "none",
  legend.title = element_blank(), 
  legend.text = element_text(family = font,
                             size = fontsize,
                             color = basecolor),
  axis.text = element_text(family = font,
                           size = fontsize,
                           color = basecolor),
  axis.ticks = element_blank(),
  aspect.ratio = 1,
  axis.title.x = element_text(hjust = c(1), vjust = c(0)),
  axis.title.y = element_text(hjust = c(1), vjust = c(1))
)
