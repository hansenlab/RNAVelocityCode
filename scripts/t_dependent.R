rm(list=ls())
source(here::here("scripts/utils.R"))

scv <- reticulate::import("scvelo")
np <- reticulate::import("numpy")

adata <- scv$datasets$simulation(n_obs=500L, n_vars=10L, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=0.8)
scv$pp$pca(adata)
scv$tl$velocity(adata)


# simulate data with time-dep rate
n_obs <- 500L
alpha <- as.list(rep(5, 10))
gamma <- as.list(rep(0.5, 10))
# three genes with time-dependent rates
tnew <- np$linspace(.1, 5, num=n_obs)
gamma_t <- 5 * (1 - np$exp(-.3*tnew))
for (i in c(1, 2, 6)) {
	gamma[[i]] <- gamma_t
}

adata2 <- scv$datasets$simulation(n_obs=n_obs, n_vars=10L, t_max=25, alpha=alpha, beta=.3, gamma=gamma, noise_level=0.8)

scv$pp$pca(adata2)
scv$tl$velocity(adata2)

adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)

adata2$var$velocity_genes <- TRUE
scv$tl$velocity_graph(adata2)

sce.o <- zellkonverter::AnnData2SCE(adata)
sce2.o <- zellkonverter::AnnData2SCE(adata2)



### plot
embedded <- embedVelocity(reducedDim(sce.o, "X_pca"), sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded)
# p1 <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t, title = "Constant \u03B3 and \u03B1", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 3.01) +
# 	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
p1 <- plotVelocityStream(sce.o, embedded[, 1:2], color_by = "true_t", color.alpha = 0.7, grid.resolution = 15,
												 stream.L = 5, stream.min.L = 0.1, stream.res = 8, stream.width = 0.3,
												 arrow.angle = 18, arrow.length = 0.2) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487",
																		"#35608D", "#2A788E", "#21908C",
																		"#22A884", "#43BF71", "#7AD151",
																		"#BBDF27", "#FDE725"),name = "t", limits = range(sce.o$true_t, na.rm = TRUE)) +
	labs( y = "PCA-2", x = "PCA-1") +
	ggtitle("\u03B1 = 5, constant \u03B3") +
	.theme_noframe + theme(legend.position = c(1, 1),
												 legend.justification = c(1, 1))
	
# p1_inset  <- plotxy(x = sce.o$true_t, y = rep(0.5, 10), x_lab = "t", y_lab = "\u03B3", point.size = 2.01) +
# 	theme_half_open(6) +
# 	theme(
# 		axis.line = element_line(colour="black", linetype = 1, size= 0.25),
# 		legend.background = element_blank(),
# 		panel.grid.minor = element_blank(),
# 		panel.grid.major.x = element_blank(),
# 		panel.grid.major.y = element_blank(),
# 		axis.ticks = element_blank(),
# 		axis.text.x = element_blank(),
# 		axis.text.y = element_blank())
# 
# p1 <- ggdraw(p1 ) +
# 	draw_plot(p1_inset, .5, .4, .3, .3) 

p2 <- plotxy(assay(sce.o, "spliced")[6, ], assay(sce.o, "unspliced")[6, ], title = "\u03B1 = 5, constant \u03B3", x_lab = "Spliced",
						 y_lab = "Unspliced", color_by = "true_t", point.size = 5.01,
						 colors.v = sce.o$true_t,
						 hue.colors = c("#440154", "#482576", "#414487","#35608D",
						 							 "#2A788E", "#21908C","#22A884", "#43BF71",
						 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_abline(slope = rowData(sce.o)$velocity_gamma[6], intercept = 0, linetype = "dashed", size = 0.5) +
	theme(legend.position = "none")

de.v <- assay(sce.o, "Ms")[, 80] - assay(sce.o, "Ms")[, 50]
p3 <- plotxy(x = de.v, assay(sce.o, "velocity")[, 50], title = "\u03B1 = 5, constant \u03B3", x_lab = "Difference of Ms between two cells",
						 y_lab = "Estimated velocity", point.size = 8.01) +
	annotate(geom = "text", x = .percent_range(de.v, 0.5),
					 y = .percent_range(assay(sce.o, "velocity")[, 50],0.6), size = 3, hjust = 0.5, vjust = 1,
					 label = str_c("Cosine similarity= ", sprintf("%.2f", lsa::cosine(de.v, assay(sce.o, "velocity")[, 50])[1, 1]))) +
	theme(legend.position = "none")


embedded2 <- embedVelocity(reducedDim(sce2.o, "X_pca"), sce2.o)
grid2.df <- gridVectors(reducedDim(sce2.o, "X_pca"), embedded2)
# p3 <- plotEmbExp(sce2.o, "X_pca", color.v = sce2.o$true_t, title = "Time-dependent \u03B3", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 3.01) +
# 	geom_segment(data=grid2.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
p4 <- plotVelocityStream(sce2.o, embedded2[, 1:2], color_by = "true_t", color.alpha = 0.7, grid.resolution = 15,
												 stream.L = 5, stream.min.L = 0.1, stream.res = 8, stream.width = 0.3,
												 arrow.angle = 18, arrow.length = 0.2) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487",
																		"#35608D", "#2A788E", "#21908C",
																		"#22A884", "#43BF71", "#7AD151",
																		"#BBDF27", "#FDE725"),name = "t", limits = range(sce2.o$true_t, na.rm = TRUE)) +
	labs( y = "PCA-2", x = "PCA-1") +
	ggtitle("\u03B1 = 5, time-dep. \u03B3") +
	.theme_noframe + theme(legend.position = "none",
												 legend.justification = c(1, 1))
p4_inset  <- plotxy(x = sce2.o$true_t, y = gamma_t, x_lab = "t", y_lab = "\u03B3", point.size = 2.01) +
	theme_half_open(6) +
	theme(
		axis.line = element_line(colour="black", linetype = 1, size= 0.25),
		legend.background = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_blank(),
		axis.ticks = element_blank(),
		axis.text.x = element_blank(),
		axis.text.y = element_blank())

p4 <- ggdraw(p4 ) +
	draw_plot(p4_inset, .65, .55, .3, .3) 
	
	
p5 <- plotxy(assay(sce2.o, "spliced")[6, ], assay(sce2.o, "unspliced")[6, ], title = "\u03B1 = 5, time-dep. \u03B3", x_lab = "Spliced",
						 y_lab = "Unspliced", color_by = "true_t", point.size = 5.01,
						 colors.v = sce.o$true_t,
						 hue.colors = c("#440154", "#482576", "#414487","#35608D",
						 							 "#2A788E", "#21908C","#22A884", "#43BF71",
						 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_abline(slope = rowData(sce2.o)$velocity_gamma[6], intercept = 0, linetype = "dashed", size = 0.5) +
	theme(legend.position = "none")

de2.v <- assay(sce2.o, "Ms")[, 80] - assay(sce2.o, "Ms")[, 50]
p6 <- plotxy2(x = de2.v, assay(sce2.o, "velocity")[, 50], title = "\u03B1 = 5, time-dep. \u03B3", x_lab = "Difference of Ms between two cells",
			 y_lab = "Estimated velocity", point.size = 8.01, color_var = factor(c(1, 1, 0, 0, 0, 1, rep(0, 4))),
				colors.v = c("black", "#377EB8"), labels = c(0, 1)) +
	annotate(geom = "text", x = .percent_range(de2.v, 0.5),
					 y = .percent_range(assay(sce2.o, "velocity")[, 50], 0.6), size = 3, hjust = 0.5, vjust = 1,
					 label = str_c("Cosine similarity= ", sprintf("%.2f", lsa::cosine(de2.v, assay(sce2.o, "velocity")[, 50])[1, 1]))) +
	theme(legend.position = "none")
### simulate Fig2D


# simulate data with time-dep rate
alpha <- as.list(rep(5, 10))
gamma <- as.list(rep(0.5, 10))
for (i in c(1, 2, 6)) {
	alpha[[i]] <- 30
	gamma[[i]] <- gamma_t
}

adata3 <- scv$datasets$simulation(n_obs=n_obs, n_vars=10L, t_max=25, alpha=alpha, beta=.3, gamma=gamma, noise_level=0.8)
scv$pp$pca(adata3)
scv$tl$velocity(adata3)
adata3$var$velocity_genes <- TRUE
scv$tl$velocity_graph(adata3)
sce3.o <- zellkonverter::AnnData2SCE(adata3)


### plot
embedded3 <- embedVelocity(reducedDim(sce3.o, "X_pca"), sce3.o)
grid3.df <- gridVectors(reducedDim(sce3.o, "X_pca"), embedded3, resolution = 10)
# p5 <- plotEmbExp(sce3.o, "X_pca", color.v = sce3.o$true_t, title = "Different \u03B1", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 3.01) +
# 	geom_segment(data=grid3.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.08) 
p7 <- plotVelocityStream(sce3.o, embedded3[, 1:2], color_by = "true_t", color.alpha = 0.7, grid.resolution = 15,
									 stream.L = 5, stream.min.L = 0.1, stream.res = 8, stream.width = 0.3,
									 arrow.angle = 18, arrow.length = 0.2) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487",
																		"#35608D", "#2A788E", "#21908C",
																		"#22A884", "#43BF71", "#7AD151",
																		"#BBDF27", "#FDE725"),name = "t", limits = range(sce3.o$true_t, na.rm = TRUE)) +
	labs( y = "PCA-2", x = "PCA-1") +
	ggtitle("\u03B1 = 30, time-dep. \u03B3") +
	.theme_noframe + theme(legend.position = "none",
												 legend.justification = c(1, 1))

p8 <- plotxy(assay(sce3.o, "spliced")[6, ], assay(sce3.o, "unspliced")[6, ], title = "\u03B1 = 30, time-dep. \u03B3", x_lab = "Spliced",
						 y_lab = "Unspliced", color_by = "true_t", point.size = 5.01,
						 colors.v = sce.o$true_t,
						 hue.colors = c("#440154", "#482576", "#414487","#35608D",
						 							 "#2A788E", "#21908C","#22A884", "#43BF71",
						 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_abline(slope = rowData(sce3.o)$velocity_gamma[6], intercept = 0, linetype = "dashed", size = 0.5) +
	theme(legend.position = "none")

de3.v <- assay(sce3.o, "Ms")[, 80] - assay(sce3.o, "Ms")[, 50]
p9 <- plotxy2(x = de3.v, assay(sce3.o, "velocity")[, 50], title = "\u03B1 = 30, time-dep. \u03B3", x_lab = "Difference of Ms between two cells",
						 y_lab = "Estimated velocity", point.size = 8.01, color_var = factor(c(1, 1, 0, 0, 0, 1, rep(0, 4))),
							colors.v = c("black", "#E41A1C"), labels = c(0, 1)) +
	annotate(geom = "text", x = .percent_range(de3.v, 0.5),
					 y = .percent_range(assay(sce3.o, "velocity")[, 50], 0.6), size = 3, hjust = 0.5, vjust = 1,
					 label = str_c("Cosine similarity= ", sprintf("%.2f", lsa::cosine(de3.v, assay(sce3.o, "velocity")[, 50])[1, 1]))) +
	theme(legend.position = "none")


mp <- plot_grid(p1+ theme(legend.position = c(1, 1),
													legend.justification = c(1, 1)),
								p4 , p7,
								p2 + theme(legend.position = "none"), 
								p5, p8,
								p3, p6, p9,
								nrow = 3, ncol = 3, label_size = 10, labels = c("a", "b", "c"))

save_plot(here::here("figs", "t_dependent", "t_dependent.pdf"), mp,
					base_height = 2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 3*1.2, ncol = 3, device = cairo_pdf)




