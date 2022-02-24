rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/Pancreas/dynamical.o.qs"))
subset.o <- qread(here::here("data/Pancreas/sce_subset.o.qs"))

g <- "Gnas"
p1 <- plot_dynamicalmums2(sce.o, rownm = g, color_by = "clusters", colors.v = metadata(sce.o)$clusters_colors, point.size = 4.01,
													title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),")")) + 
	theme(legend.position = "none")

p2 <- plot_dynamicalmums2(subset.o, rownm = g, color_by = "clusters", colors.v = metadata(subset.o)$clusters_colors[-4], point.size = 4.01,
													title = str_c(g, " (LL:", sprintf("%.3f", rowData(subset.o[g, ])$fit_likelihood),")")) + 
	theme(legend.position = "none")

mp <- plot_grid(p1, p2,
								nrow = 1, ncol = 2, label_size = 10, labels = c("a", "b"))

save_plot(here::here("figs", "pancreas", "pancreas_gnas.pdf"), mp,
					base_height = 2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 1*1.2, ncol = 2, device = cairo_pdf)



