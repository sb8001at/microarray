library(ggplot2)
library(dplyr)
library(tidyr)

mytheme <- function(){
	theme_bw() +
	theme(
		plot.title = element_text(size = 25, vjust=1.5),
		axis.title.x = element_text(size = 22.5),
		axis.title.y = element_text(size = 22.5, vjust=0.25),
		axis.text = element_text(size = 15, color="black"),
		strip.text = element_text(size = 15),
		legend.title = element_text(size = 15),
		legend.text = element_text(size = 15),
		legend.position = "bottom", 
		legend.background = element_rect(color="black"),
		panel.grid.major = element_line(size=1),
		panel.grid.minor = element_line(size=0.25)
	)
}
