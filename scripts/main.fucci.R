rm(list=ls())
source(here::here("scripts/utils.R"))

steady.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))
dynamical.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))

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
	.theme_noframe + theme(legend.position = c(0, 0),
												 legend.justification = c(0, 0),
												 legend.key.size = unit(6, "pt"))

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
	.theme_noframe + theme(legend.position = c(0, 0),
												 legend.justification = c(0, 0),
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

g_color.v <- ifelse((dynamical.o$fucci_time < cutoff1.v) | (dynamical.o$fucci_time > cutoff2.v), "down", "up")

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



gc_vcolor2.v <- g_color.v
gc_vcolor2.v[(c(-1, 1)[as.numeric(factor(g_color.v))] * assay(dynamical.o, "velocity")[g, ]) < 0] <- "non"

dynamical_velo.p <-  plotFucciLoess2(dynamical.o, assay(dynamical.o, "velocity")[g, ], y_lab = "Velocity",color_var = gc_vcolor2.v,
																		 colors.v = c("#80B1D3", "black","#FB8072"), labels.v = c("Down", "Incons.", "Up"), color.name = "Direction",
																		title = as.expression(substitute(italic(g)~" dynamical velocity", list(g=g))), x.pos = 0.45) +
	geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key.size = unit(6, "pt"))






mp <- plot_grid(dynamical_vec.p, mums.p,
								# ggplot() + xlim(c(0, 1)) + ylim(c(0, 1)) +
								# 	annotate("text", x = 0.5, y = 0.5, label = "Placeholder") +
								# 	theme_nothing(),
								ms_t.p, dynamical_velo.p,
								# ggplot() + xlim(c(0, 1)) + ylim(c(0, 1)) +
								# 	annotate("text", x = 0.5, y = 0.5, label = "Placeholder") +
								# 	theme_nothing(),
								nrow = 2, ncol = 2, label_size = 10, labels = c("a", "b", "c", "d"))

save_plot(here::here("figs", "main", "main.fucci.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 2, device = cairo_pdf)









### Companion SFig
gc_vcolor1.v <- g_color.v
gc_vcolor1.v[(c(-1, 1)[as.numeric(factor(g_color.v))] * assay(steady.o, "velocity")[g, ]) < 0] <- "non"

steady_velo.p <-  plotFucciLoess2(steady.o, assay(steady.o, "velocity")[g, ], y_lab = "Velocity", color_var = gc_vcolor1.v,
																	colors.v = c("#80B1D3", "black","#FB8072"), labels.v = c("Down", "Incons.", "Up"), color.name = "Direction",
																	title = as.expression(substitute(italic(g)~" steady velocity", list(g=g)))) +
	geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1))


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

save_plot(here::here("figs", "sfigs", "sfig.fucci.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 2, device = cairo_pdf)




### another Sfig for top 10 cc genes
top10.lp <- lapply(rownames(cc.df)[seq_len(10)[-1]], function(g) {
	mums.p <- plot_dynamicalmums(dynamical.o, gene = g, color_by = "fucci_time",
										 title = as.expression(substitute(italic(g)~" phase portrait", list(g=g))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	
	
	## ms t
	loess.l <- tricycle::fit_periodic_loess(theta.v = dynamical.o$fucci_time * 2 * pi, y = assay(dynamical.o, "Ms")[g, ])
	cutoff1.v <- loess.l$pred.df$x[which.min(loess.l$pred.df$y)] / (2 * pi)
	cutoff2.v <- loess.l$pred.df$x[which.max(loess.l$pred.df$y)] / (2 * pi)
	
	g_color.v <- ifelse((dynamical.o$fucci_time < cutoff1.v) | (dynamical.o$fucci_time > cutoff2.v), "down", "up")
	
	ms_t.p <- plotFucciLoess2(dynamical.o, assay(dynamical.o, "Ms")[g, ], y_lab = "Ms", color_var = g_color.v,
														colors.v = c("#80B1D3", "#FB8072"), labels.v = c("Down", "Up"), color.name = "Direction",
														title = as.expression(substitute(italic(g)~" Ms over FUCCI", list(g=g)))) +
		theme(legend.position = "none")
	
	##
	gc_vcolor2.v <- g_color.v
	gc_vcolor2.v[(c(-1, 1)[as.numeric(factor(g_color.v))] * assay(dynamical.o, "velocity")[g, ]) < 0] <- "non"
	
	dynamical_velo.p <-  plotFucciLoess2(dynamical.o, assay(dynamical.o, "velocity")[g, ], y_lab = "Velocity",color_var = gc_vcolor2.v,
																			 colors.v = c("#80B1D3", "black","#FB8072"), labels.v = c("Down", "Incons.", "Up"), color.name = "Direction",
																			 title = as.expression(substitute(italic(g)~" dynamical velocity", list(g=g))), x.pos = 0.7, y.pos = 0.1) +
		geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
		theme(legend.position = "none")
	
	gc_vcolor1.v <- g_color.v
	gc_vcolor1.v[(c(-1, 1)[as.numeric(factor(g_color.v))] * assay(steady.o, "velocity")[g, ]) < 0] <- "non"
	
	steady_velo.p <-  plotFucciLoess2(steady.o, assay(steady.o, "velocity")[g, ], y_lab = "Velocity", color_var = gc_vcolor1.v,
																		colors.v = c("#80B1D3", "black","#FB8072"), labels.v = c("Down", "Incons.", "Up"), color.name = "Direction",
																		title = as.expression(substitute(italic(g)~" steady velocity", list(g=g))), x.pos = 0.7, y.pos = 0.1) +
		geom_hline(yintercept = 0, linetype = "solid", alpha = 0.7) +
		theme(legend.position = "none")
	
	
	compare.p <- plot_comV(steady.o, dynamical.o, gene = g, point.size = 3.01, title = as.expression(substitute(italic(g)~" velocity comparison", list(g=g))),
												 x.pos = 0.65) +
		theme(legend.position = "none")
	
	mp <- plot_grid(mums.p, ms_t.p, dynamical_velo.p, steady_velo.p, compare.p, 
									nrow = 1, ncol = 5, label_size = 10, labels = NULL)
	return(mp)
})
mp <- plot_grid(plotlist = top10.lp,
								nrow = 9, ncol = 1, label_size = 10, labels = NULL)
save_plot(here::here("figs", "sfigs", "sfig.fucci_top10.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 9, ncol = 5, device = cairo_pdf)




