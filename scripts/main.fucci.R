rm(list=ls())
source(here::here("scripts/utils.R"))

steady.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))
dynamical.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))

### flip UMAP-1 
reducedDim(steady.o, "UMAP")[, 1] <-  - reducedDim(steady.o, "UMAP")[, 1]
reducedDim(dynamical.o, "UMAP")[, 1] <-  - reducedDim(dynamical.o, "UMAP")[, 1]


### use tricycle position instead of fucci time
steady.o <- tricycle::estimate_cycle_position(steady.o, gname.type = "SYMBOL", species = "human", exprs_values = "X")
dynamical.o <- tricycle::estimate_cycle_position(dynamical.o, gname.type = "SYMBOL", species = "human", exprs_values = "X")





library(org.Hs.eg.db)
go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
															columns = "SYMBOL")[, "SYMBOL"]

cc.df <- rowData(dynamical.o)[rownames(dynamical.o) %in% go.v, ] %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

g <- rownames(cc.df)[1]

### embeddings vectors
embedded <- embedVelocity(reducedDim(steady.o, "UMAP"), steady.o)
# 
steady_vec.p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "tricyclePosition", color.alpha = 0.7, grid.resolution = 13,
									 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
									 arrow.angle = 18, arrow.length = 0.45) +
	scale_color_gradientn(colours = c("#2E22EA","#9E3DFB","#F86BE2",
																		"#FCCE7B","#C4E416","#4BBA0F",
																		"#447D87","#2C24E9"),name = "\u03b8", limits = c(0, 2 * pi)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle("Steady model") +
	.theme_noframe + theme(legend.position = "none")

embedded <- embedVelocity(reducedDim(dynamical.o, "UMAP"), dynamical.o)

dynamical_vec.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "tricyclePosition", color.alpha = 0.7, grid.resolution = 13,
																			stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2, point.size = 4.01,
																			arrow.angle = 18, arrow.length = 0.45) +
	scale_color_gradientn(colours = c("#2E22EA","#9E3DFB","#F86BE2",
																		"#FCCE7B","#C4E416","#4BBA0F",
																		"#447D87","#2C24E9"),name = "\u03b8", limits = c(0, 2 * pi)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle("Dynamical model") +
	.theme_noframe + theme(legend.position = "none")

mums.p <- plot_dynamicalmums(dynamical.o, gene = g, color_by = "tricyclePosition", point.size = 4.01,
														 title = as.expression(substitute(italic(g)~" phase portrait", list(g=g))), color.name = "\u03b8") +
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
loess.l <- tricycle::fit_periodic_loess(theta.v = dynamical.o$tricyclePosition, y = assay(dynamical.o, "Ms")[g, ])
cutoff1.v <- loess.l$pred.df$x[which.min(loess.l$pred.df$y)] 
cutoff2.v <- loess.l$pred.df$x[which.max(loess.l$pred.df$y)] 

if (cutoff1.v < cutoff2.v) {
	g_color.v <- ifelse((dynamical.o$tricyclePosition > cutoff1.v) & (dynamical.o$tricyclePosition < cutoff2.v), "up", "down")
} else {
	g_color.v <- ifelse((dynamical.o$tricyclePosition < cutoff2.v) | (dynamical.o$tricyclePosition > cutoff1.v), "up", "down")
}

ms_t.p <- plotFucciLoess3(dynamical.o, assay(dynamical.o, "Ms")[g, ], y_lab = "Ms", color_var = g_color.v, point.size = 4.01,
													colors.v = c("#80B1D3", "#FB8072"), labels.v = c("Down", "Up"), color.name = "Direction",
													title = as.expression(substitute(italic(g)~" Ms over \u03b8", list(g=g)))) +
	annotate("segment", x = 1.5 * pi, xend = 1.5 * pi, y = 80, yend = 20, alpha = 0.4,
					 colour = "#80B1D3", size = 2, arrow = arrow(length = unit(0.18, "inches"))) +
	annotate("segment", x = 0.3 * pi, xend = 0.3 * pi, y = 20, yend = 80, alpha = 0.4,
					 colour = "#FB8072", size = 2, arrow = arrow(length = unit(0.18, "inches"))) +
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key.size = unit(6, "pt"))



gc_vcolor2.v <- as.numeric(factor(g_color.v)) * 4 - 6
gc_vcolor2.v[(gc_vcolor2.v == -2) & (assay(dynamical.o, "velocity")[g, ] > 0)] <- -1
gc_vcolor2.v[(gc_vcolor2.v == 2) & (assay(dynamical.o, "velocity")[g, ] < 0)] <- 1
gc_vcolor2.v <- factor(gc_vcolor2.v)
dynamical.o$color <- gc_vcolor2.v

dynamical_velo.p <-  plotFucciLoess3(dynamical.o, assay(dynamical.o, "velocity")[g, ], y_lab = "Velocity",color_var = gc_vcolor2.v,
																		 colors.v = c("#80B1D3", "black", "yellow","#FB8072"), point.size = 4.01,
																		 labels.v =  c("Down", "Down incons.", "Up incons.", "Up"), color.name = "Direction",
																		title = as.expression(substitute(italic(g)~" dynamical velocity", list(g=g))), x.pos = 0.45) +
	geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
	theme(legend.position = c(0.1, 0),
				legend.justification = c(0.1, 0),
				legend.key.size = unit(6, "pt"))


tricycle_t.p <- plotFucciLoess3(dynamical.o, assay(dynamical.o, "fit_t")[g, ], y_lab = "Velocity latent time t",color_var = dynamical.o$tricyclePosition,
																point.size = 4.01,
																title = as.expression(substitute(italic(g)~" dynamical velocity", list(g=g))), x.pos = 0.25) +
	scale_color_gradientn(colours = c("#2E22EA","#9E3DFB","#F86BE2",
																		"#FCCE7B","#C4E416","#4BBA0F",
																		"#447D87","#2C24E9"),name = "\u03b8", limits = c(0, 2 * pi)) +
	theme(legend.position = "none")

	# plotxy(dynamical.o$tricyclePosition, assay(dynamical.o, "fit_t")[g, ], title = as.expression(substitute(italic(g)~" dynamical velocity", list(g=g))), x_lab = "Cell cycle position \u03b8",
	# 										y_lab = "Velocity latent time t", color_by = "true_t", point.size = 4.01,
	# 										colors.v = dynamical.o$tricyclePosition,
	# 										hue.colors = c("#2E22EA","#9E3DFB","#F86BE2",
	# 																	 "#FCCE7B","#C4E416","#4BBA0F",
	# 																	 "#447D87","#2C24E9"), color.name = "\u03b8") 

cir_legend.p <- tricycle::circle_scale_legend(text.size = 3, y.inner =  0.7, ymax = 4.5, y.outer = 2.1, y.text = 2.8, addStageLabel = TRUE)


mp <- plot_grid(plot_grid(dynamical_vec.p, ggplot() + theme_nothing(), mums.p, nrow = 1, ncol = 3, label_size = 10, labels = c("a", "", "b")),
								# ggplot() + xlim(c(0, 1)) + ylim(c(0, 1)) +
								# 	annotate("text", x = 0.5, y = 0.5, label = "Placeholder") +
								# 	theme_nothing(),
								plot_grid(# ggplot() + theme_nothing(), 
													ms_t.p, dynamical_velo.p,
													tricycle_t.p,
													nrow = 1, ncol = 3,
													label_size = 10, labels = c("c", "d", "e")),
								# ggplot() + xlim(c(0, 1)) + ylim(c(0, 1)) +
								# 	annotate("text", x = 0.5, y = 0.5, label = "Placeholder") +
								# 	theme_nothing(),
								nrow = 2, ncol = 1, label_size = 10, labels = NULL)

mp2 <- ggdraw(mp) +
	draw_plot(cir_legend.p,
						0.095, 0.37, 0.8, 0.8, hjust = 0, vjust = 0, halign = 0, valign = 0)

save_plot(here::here("figs", "main", "main.fucci.pdf"), mp2,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 3, device = cairo_pdf)









### Companion SFig
gc_vcolor1.v <- as.numeric(factor(g_color.v)) * 4 - 6
gc_vcolor1.v[(gc_vcolor1.v == -2) & (assay(steady.o, "velocity")[g, ] > 0)] <- -1
gc_vcolor1.v[(gc_vcolor1.v == 2) & (assay(steady.o, "velocity")[g, ] < 0)] <- 1

gc_vcolor1.v <- factor(gc_vcolor1.v)
# dynamical.o$color <- gc_vcolor1.v

steady_velo.p <-  plotFucciLoess3(steady.o, assay(steady.o, "velocity")[g, ], y_lab = "Velocity", color_var = gc_vcolor1.v,
																	colors.v = c("#80B1D3", "black", "yellow","#FB8072"), point.size = 4.01,
																	labels.v =  c("Down", "Down incons.", "Up incons.", "Up"), color.name = "Direction",
																	title = as.expression(substitute(italic(g)~" steady velocity", list(g=g)))) +
	geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
	theme(legend.position = c(0.1, 0),
				legend.justification = c(0.1, 0),
				legend.key.size = unit(6, "pt"))


compare.p <- plot_comV(steady.o, dynamical.o, gene = g, point.size = 3.01, title = as.expression(substitute(italic(g)~" velocity comparison", list(g=g))),
											 x.pos = 0.65) +
	theme(legend.position = "none")


### idx for velocity genes for both
gene.v <- intersect(rownames(steady.o)[rowData(steady.o)$velocity_genes],
										rownames(dynamical.o)[rowData(dynamical.o)$velocity_genes])

v_steady.m <- assay(steady.o[gene.v, ], "velocity")
v_dynamical.m <- assay(dynamical.o[gene.v, ], "velocity")

inconsistency.v <- sapply(seq_len(nrow(v_steady.m)), function(i) mean(sign(v_steady.m[i, ] * v_dynamical.m[i, ]) == -1))  ### median 0.3060686
pcc.v <- sapply(seq_len(nrow(v_steady.m)), function(i) cor(v_steady.m[i, ], v_dynamical.m[i, ]))   ### median 0.4027132

### scatter plot
tmp.df <- data.frame(x = inconsistency.v, y = pcc.v)

incon_pcc.scat.p <- ggplot(data = tmp.df, aes(x = x, y = y)) +
	geom_scattermore(pointsize = 3.01,  alpha = 0.8) +
	geom_vline(xintercept = median(inconsistency.v), color = "blue", alpha = 0.5) +
	geom_hline(yintercept = median(pcc.v), color = "blue", alpha = 0.5) +
	scale_x_continuous(breaks = seq(0, 0.8, 0.2), labels = c(0, 20, 40, 60, 80), name = "Perct. of incon. directions(%)") +
	labs(y = "PCC", title = "Comparison of velocities")



mp <- plot_grid(steady_vec.p, steady_velo.p, 
								 compare.p,  incon_pcc.scat.p,
								nrow = 2, ncol = 2, label_size = 10, labels = c("auto"))

save_plot(here::here("figs", "sfigs", "sfig.fucci_steady.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 2, device = cairo_pdf)




### another Sfig for top 10 cc genes
top10.lp <- lapply(rownames(cc.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums(dynamical.o, gene = g, color_by = "tricyclePosition",
															 title = as.expression(substitute(italic(g)~" phase portrait", list(g=g))), color.name = "\u03b8") +
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
	loess.l <- tricycle::fit_periodic_loess(theta.v = dynamical.o$tricyclePosition, y = assay(dynamical.o, "Ms")[g, ])
	cutoff1.v <- loess.l$pred.df$x[which.min(loess.l$pred.df$y)] 
	cutoff2.v <- loess.l$pred.df$x[which.max(loess.l$pred.df$y)] 
	
	if (cutoff1.v < cutoff2.v) {
		g_color.v <- ifelse((dynamical.o$tricyclePosition > cutoff1.v) & (dynamical.o$tricyclePosition < cutoff2.v), "up", "down")
	} else {
		g_color.v <- ifelse((dynamical.o$tricyclePosition < cutoff2.v) | (dynamical.o$tricyclePosition > cutoff1.v), "up", "down")
	}
	
	
	
	ms_t.p <- plotFucciLoess3(dynamical.o, assay(dynamical.o, "Ms")[g, ], y_lab = "Ms", color_var = g_color.v,
														colors.v = c("#80B1D3", "#FB8072"), labels.v = c("Down", "Up"), color.name = "Direction",
														title = as.expression(substitute(italic(g)~" Ms over \u03b8", list(g=g)))) +
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
	
	dynamical_velo.p <-  plotFucciLoess3(dynamical.o, assay(dynamical.o, "velocity")[g, ], y_lab = "Velocity",color_var = gc_vcolor2.v,
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
	
	steady_velo.p <-  plotFucciLoess3(steady.o, assay(steady.o, "velocity")[g, ], y_lab = "Velocity", color_var = gc_vcolor1.v,
																		colors.v = c("#80B1D3", "black", "yellow","#FB8072"), labels.v =c("Down", "Down incons.", "Up incons.", "Up"), color.name = "Direction",
																		title = as.expression(substitute(italic(g)~" steady velocity", list(g=g))), x.pos = 0.45) +
		geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
		theme(legend.position = "none") 
	
	
	compare.p <- plot_comV(steady.o, dynamical.o, gene = g, point.size = 3.01, title = as.expression(substitute(italic(g)~" velocity comparison", list(g=g))),
												 x.pos = 0.65) +
		theme(legend.position = "none")
	
	
	
	mp <- plot_grid(mums.p, ms_t.p, dynamical.p, steady.p, dynamical_velo.p, steady_velo.p, compare.p, 
									nrow = 1, ncol = 7, label_size = 10, labels = NULL)
	return(mp)
})
mp <- plot_grid(plotlist = top10.lp,
								nrow = 10, ncol = 1, label_size = 10, labels = NULL)
save_plot(here::here("figs", "sfigs", "sfig.fucci_top10.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 10, ncol = 7, device = cairo_pdf)




