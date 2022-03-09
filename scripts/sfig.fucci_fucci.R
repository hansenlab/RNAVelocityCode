rm(list=ls())
source(here::here("scripts/utils.R"))

steady.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))
dynamical.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))

### flip UMAP-1 
reducedDim(steady.o, "UMAP")[, 1] <-  - reducedDim(steady.o, "UMAP")[, 1]
reducedDim(dynamical.o, "UMAP")[, 1] <-  - reducedDim(dynamical.o, "UMAP")[, 1]


library(org.Hs.eg.db)
go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
															columns = "SYMBOL")[, "SYMBOL"]

cc.df <- rowData(dynamical.o)[rownames(dynamical.o) %in% go.v, ] %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

g <- rownames(cc.df)[1]

### embeddings vectors
embedded <- embedVelocity(reducedDim(steady.o, "UMAP"), steady.o)
# 
steady_vec.p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "fucci_time", color.alpha = 0.7, grid.resolution = 13,
																	 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																	 arrow.angle = 18, arrow.length = 0.45) +
	scale_color_gradientn(colours = c("#2E22EA","#9E3DFB","#F86BE2",
																		"#FCCE7B","#C4E416","#4BBA0F",
																		"#447D87","#2C24E9"),name = "FUCCI", limits = range(steady.o$fucci_time, na.rm = TRUE)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle("Steady model") +
	.theme_noframe + theme(legend.position = "none")

embedded <- embedVelocity(reducedDim(dynamical.o, "UMAP"), dynamical.o)
grid.df <- gridVectors(reducedDim(dynamical.o, "UMAP"), embedded)

dynamical_vec.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "fucci_time", color.alpha = 0.7, grid.resolution = 13,
																			stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																			arrow.angle = 18, arrow.length = 0.45) +
	scale_color_gradientn(colours = c("#2E22EA","#9E3DFB","#F86BE2",
																		"#FCCE7B","#C4E416","#4BBA0F",
																		"#447D87","#2C24E9"),name = "FUCCI", limits = range(dynamical.o$fucci_time, na.rm = TRUE)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle("Dynamical model") +
	.theme_noframe + theme(legend.position = c(1, 0),
												 legend.justification = c(1, 0),
												 legend.key.size = unit(6, "pt"))

mums.p <- plot_dynamicalmums(dynamical.o, gene = g, color_by = "fucci_time",
														 title = as.expression(substitute(italic(g)~" phase portrait", list(g=g))), color.name = "FUCCI") +
	geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
	annotate(geom = "text", x = .percent_range(assay(dynamical.o, "Ms")[g, ], 0.63),
					 y = .percent_range(assay(dynamical.o, "Mu")[g, ], 0.25), size = 2, hjust = 0, vjust = 1,
					 label = str_c("italic(R) ^ 2: ", sprintf("%.3f", rowData(dynamical.o[g, ])$fit_r2)), parse = TRUE) +
	annotate(geom = "text", x = .percent_range(assay(dynamical.o, "Ms")[g, ], 0.63),
					 y = .percent_range(assay(dynamical.o, "Mu")[g, ], 0.15), size = 2, hjust = 0, vjust = 1,
					 label = str_c("Steady \u03B3: ", sprintf("%.3f", rowData(steady.o[g, ])$velocity_gamma)), parse = FALSE) +
	annotate(geom = "text", x = .percent_range(assay(dynamical.o, "Ms")[g, ], 0.63),
					 y = .percent_range(assay(dynamical.o, "Mu")[g, ], 0.05), size = 2, hjust = 0, vjust = 1,
					 label = str_c("Dyn. LL: ", sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)), parse = FALSE) +
	theme(legend.position = "none")


### fit periodic loess and choose cutoffs
loess.l <- tricycle::fit_periodic_loess(theta.v = dynamical.o$fucci_time * 2 * pi, y = assay(dynamical.o, "Ms")[g, ])
cutoff1.v <- loess.l$pred.df$x[which.min(loess.l$pred.df$y)] / (2 * pi)
cutoff2.v <- loess.l$pred.df$x[which.max(loess.l$pred.df$y)] / (2 * pi)

if (cutoff1.v < cutoff2.v) {
	g_color.v <- ifelse((dynamical.o$fucci_time > cutoff1.v) & (dynamical.o$fucci_time < cutoff2.v), "up", "down")
} else {
	g_color.v <- ifelse((dynamical.o$fucci_time < cutoff2.v) | (dynamical.o$fucci_time > cutoff1.v), "up", "down")
}
ms_t.p <- plotFucciLoess2(dynamical.o, assay(dynamical.o, "Ms")[g, ], y_lab = "Ms", color_var = g_color.v,
													colors.v = c("#80B1D3", "#FB8072"), labels.v = c("Down", "Up"), color.name = "Direction",
													title = as.expression(substitute(italic(g)~" Ms over FUCCI", list(g=g)))) +
	annotate("segment", x = 0.125, xend = 0.125, y = 60, yend = 15, alpha = 0.4,
					 colour = "#80B1D3", size = 2, arrow = arrow(length = unit(0.15, "inches"))) +
	annotate("segment", x = 0.55, xend = 0.55, y = 15, yend = 60, alpha = 0.4,
					 colour = "#FB8072", size = 2, arrow = arrow(length = unit(0.15, "inches"))) +
	annotate("segment", x = 0.875, xend = 0.875, y = 70, yend = 25, alpha = 0.4,
					 colour = "#80B1D3", size = 2, arrow = arrow(length = unit(0.15, "inches"))) +
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key.size = unit(6, "pt"))



gc_vcolor2.v <- as.numeric(factor(g_color.v)) * 4 - 6
gc_vcolor2.v[(gc_vcolor2.v == -2) & (assay(dynamical.o, "velocity")[g, ] > 0)] <- -1
gc_vcolor2.v[(gc_vcolor2.v == 2) & (assay(dynamical.o, "velocity")[g, ] < 0)] <- 1
gc_vcolor2.v <- factor(gc_vcolor2.v)
dynamical.o$color <- gc_vcolor2.v


dynamical_velo.p <-  plotFucciLoess2(dynamical.o, assay(dynamical.o, "velocity")[g, ], y_lab = "Velocity",color_var = gc_vcolor2.v,
																		 colors.v = c("#80B1D3", "black", "yellow","#FB8072"), point.size = 4.01,
																		 labels.v =  c("Down", "Down incons.", "Up incons.", "Up"), color.name = "Direction",
																		 title = as.expression(substitute(italic(g)~" dynamical velocity", list(g=g))), x.pos = 0.42) +
	geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key.size = unit(6, "pt"))


fucci_t.p <- plotFucciLoess(dynamical.o, assay(dynamical.o, "fit_t")[g, ], y_lab = "Velocity latent time t",
																point.size = 4.01,
																title = as.expression(substitute(italic(g)~" dynamical velocity", list(g=g))), x.pos = 0.45) +
	theme(legend.position = "none")





### Companion SFig
gc_vcolor1.v <- as.numeric(factor(g_color.v)) * 4 - 6
gc_vcolor1.v[(gc_vcolor1.v == -2) & (assay(steady.o, "velocity")[g, ] > 0)] <- -1
gc_vcolor1.v[(gc_vcolor1.v == 2) & (assay(steady.o, "velocity")[g, ] < 0)] <- 1

gc_vcolor1.v <- factor(gc_vcolor1.v)
# dynamical.o$color <- gc_vcolor1.v

steady_velo.p <-  plotFucciLoess2(steady.o, assay(steady.o, "velocity")[g, ], y_lab = "Velocity", color_var = gc_vcolor1.v,
																	colors.v = c("#80B1D3", "black", "yellow","#FB8072"), point.size = 4.01,
																	labels.v =  c("Down", "Down incons.", "Up incons.", "Up"), color.name = "Direction",
																	title = as.expression(substitute(italic(g)~" steady velocity", list(g=g))), x.pos = 0.45) +
	geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
	theme(legend.position = "none")


mp <- plot_grid(dynamical_vec.p,  mums.p, ggplot() + theme_nothing(),
								ms_t.p, dynamical_velo.p, fucci_t.p,
								steady_vec.p, steady_velo.p,
								nrow = 3, ncol = 3, label_size = 10, labels = c(letters[1:2], "", letters[3:7]))

save_plot(here::here("figs", "sfigs", "sfig.fucci_fucci.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 3, ncol = 3, device = cairo_pdf)




### another Sfig for top 10 cc genes
top10.lp <- lapply(rownames(cc.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums(dynamical.o, gene = g, color_by = "fucci_time",
															 title = as.expression(substitute(italic(g)~" phase portrait", list(g=g))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		annotate(geom = "text", x = .percent_range(assay(dynamical.o, "Ms")[g, ], 0.63),
						 y = .percent_range(assay(dynamical.o, "Mu")[g, ], 0.25), size = 2, hjust = 0, vjust = 1,
						 label = str_c("italic(R) ^ 2: ", sprintf("%.3f", rowData(dynamical.o[g, ])$fit_r2)), parse = TRUE) +
		annotate(geom = "text", x = .percent_range(assay(dynamical.o, "Ms")[g, ], 0.63),
						 y = .percent_range(assay(dynamical.o, "Mu")[g, ], 0.15), size = 2, hjust = 0, vjust = 1,
						 label = str_c("Steady \u03B3: ", sprintf("%.3f", rowData(steady.o[g, ])$velocity_gamma)), parse = FALSE) +
		annotate(geom = "text", x = .percent_range(assay(dynamical.o, "Ms")[g, ], 0.63),
						 y = .percent_range(assay(dynamical.o, "Mu")[g, ], 0.05), size = 2, hjust = 0, vjust = 1,
						 label = str_c("Dyn. LL: ", sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)), parse = FALSE) +
		theme(legend.position = "none")
	
	
	
	## ms t
	loess.l <- tricycle::fit_periodic_loess(theta.v = dynamical.o$fucci_time * 2 * pi, y = assay(dynamical.o, "Ms")[g, ])
	cutoff1.v <- loess.l$pred.df$x[which.min(loess.l$pred.df$y)] / (2 * pi)
	cutoff2.v <- loess.l$pred.df$x[which.max(loess.l$pred.df$y)] / (2 * pi)
	
	if (cutoff1.v < cutoff2.v) {
		g_color.v <- ifelse((dynamical.o$fucci_time > cutoff1.v) & (dynamical.o$fucci_time < cutoff2.v), "up", "down")
	} else {
		g_color.v <- ifelse((dynamical.o$fucci_time < cutoff2.v) | (dynamical.o$fucci_time > cutoff1.v), "up", "down")
	}
	
	
	
	
	
	ms_t.p <- plotFucciLoess2(dynamical.o, assay(dynamical.o, "Ms")[g, ], y_lab = "Ms", color_var = g_color.v,
														colors.v = c("#80B1D3", "#FB8072"), labels.v = c("Down", "Up"), color.name = "Direction",
														title = as.expression(substitute(italic(g)~" Ms over tricycle", list(g=g)))) +
		theme(legend.position = "none") 
	
	##
	gc_vcolor2.v <- as.numeric(factor(g_color.v)) * 4 - 6
	gc_vcolor2.v[(gc_vcolor2.v == -2) & (assay(dynamical.o, "velocity")[g, ] > 0)] <- -1
	gc_vcolor2.v[(gc_vcolor2.v == 2) & (assay(dynamical.o, "velocity")[g, ] < 0)] <- 1
	gc_vcolor2.v <- factor(gc_vcolor2.v)
	dynamical.o$color <- gc_vcolor2.v
	
	dynamical.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "color", colors.v = c("#80B1D3", "black", "yellow","#FB8072"),
																		 point.size = 4.01, labels.v = c("Down", "Down incons.", "Up incons.", "Up"),
																		 title = as.expression(substitute(italic(g)~" dynamical velocity", list(g=g))), color.name = "Direction") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = c(0, 1),
					legend.justification = c(0, 1))
	
	dynamical_velo.p <-  plotFucciLoess2(dynamical.o, assay(dynamical.o, "velocity")[g, ], y_lab = "Velocity",color_var = gc_vcolor2.v,
																			 colors.v = c("#80B1D3", "black", "yellow","#FB8072"), labels.v =c("Down", "Down incons.", "Up incons.", "Up"), color.name = "Direction",
																			 title = as.expression(substitute(italic(g)~" dynamical velocity", list(g=g))), x.pos = 0.45) +
		geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
		theme(legend.position = "none") 
	
	gc_vcolor1.v <- as.numeric(factor(g_color.v)) * 4 - 6
	gc_vcolor1.v[(gc_vcolor1.v == -2) & (assay(steady.o, "velocity")[g, ] > 0)] <- -1
	gc_vcolor1.v[(gc_vcolor1.v == 2) & (assay(steady.o, "velocity")[g, ] < 0)] <- 1
	
	gc_vcolor1.v <- factor(gc_vcolor1.v)
	dynamical.o$color <- gc_vcolor1.v
	
	steady.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "color", colors.v = c("#80B1D3", "black", "yellow","#FB8072"),
																	point.size = 4.01, labels.v = c("Down", "Down incons.", "Up incons.", "Up"),
																	title = as.expression(substitute(italic(g)~" steady velocity", list(g=g))), color.name = "Direction") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	
	steady_velo.p <-  plotFucciLoess2(steady.o, assay(steady.o, "velocity")[g, ], y_lab = "Velocity", color_var = gc_vcolor1.v,
																		colors.v = c("#80B1D3", "black", "yellow","#FB8072"), labels.v =c("Down", "Down incons.", "Up incons.", "Up"), color.name = "Direction",
																		title = as.expression(substitute(italic(g)~" steady velocity", list(g=g))), x.pos = 0.45) +
		geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
		theme(legend.position = "none") 
	
	mp <- plot_grid(mums.p, ms_t.p, dynamical.p, steady.p, dynamical_velo.p, steady_velo.p, 
									nrow = 1, ncol = 6, label_size = 10, labels = NULL)
	return(mp)
})
mp <- plot_grid(plotlist = top10.lp,
								nrow = 10, ncol = 1, label_size = 10, labels = NULL)
save_plot(here::here("figs", "sfigs", "sfig.fucci_top10_fucci.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 10, ncol = 6, device = cairo_pdf)


