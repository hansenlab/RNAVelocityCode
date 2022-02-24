rm(list=ls())
source(here::here("scripts/utils.R"))


sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

fucci.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "fucci_time", color.alpha = 0.7, grid.resolution = 13,
															stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.4,
															arrow.angle = 18, arrow.length = 0.45) +
	scale_color_gradientn(colours = c("#2E22EA","#9E3DFB","#F86BE2",
																		"#FCCE7B","#C4E416","#4BBA0F",
																		"#447D87","#2C24E9"),name = "FUCCI", limits = range(sce.o$fucci_time, na.rm = TRUE)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle("FUCCI (dynamical)") +
	.theme_noframe + theme(legend.position = c(0, 0),
												 legend.justification = c(0, 0),
												 legend.key.size = unit(6, "pt"))

sce.o <- qread(here::here("data/Pancreas/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

pancreas.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = 15,
																			stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.4,
																			arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key.size = unit(6, "pt"))



sce.o <-qread(here::here("data/chromaffin/dynamical.o.qs"))
colors.v <- levels(factor(sce.o$cell_type))
sce.o$cell_type <- fct_recode(factor(sce.o$cell_type), "SCPs" = "#0066FFFF", "chromaffin" = "#00FF66FF",
															"sympathoblasts" = "#CC00FFFF", "non-labelled" = "#CCFF00FF", "bridge" = "#FF0000FF")
# reducedDim(sce.o, "PCA") <- reducedDim(sce.o, "PCA")[, 1:2]
embedded <- embedVelocity(reducedDim(sce.o, "PCA"), sce.o)


chromaffin.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "PCA", color_by = "cell_type", color.alpha = 0.7, grid.resolution = 20,
																	 stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4, point.size = 5.01,
																	 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = colors.v, name = "Cell type", labels = levels(sce.o$cell_type), limits = levels(sce.o$cell_type)) +
	labs( y = "PCA-2", x = "PCA-1") +
	ggtitle(str_c("Chromaffin (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/bonemarrow/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "X_tsne"), sce.o)

Bonemarrow.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "X_tsne", color_by = "clusters", color.alpha = 0.7, grid.resolution = 20,
																 stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4,
																 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "TSNE-2", x = "TSNE-1") +
	ggtitle(str_c("Bonemarrow (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/dentategyrus_lamanno/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "TSNE"), sce.o)
sce.o$ClusterName <- factor(sce.o$ClusterName, levels = c("RadialGlia", "RadialGlia2", "ImmAstro", 
																													"GlialProg", "OPC", "nIPC", "Nbl1", "Nbl2",
																													"ImmGranule1", "ImmGranule2", "Granule", "CA", 
																													"CA1-Sub", "CA2-3-4"))
Dentategyrus_lamanno.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "TSNE", color_by = "ClusterName", color.alpha = 0.7, grid.resolution = 20,
																	 stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4, point.size = 2.01,
																	 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$cluster_colors, name = "Cell type", labels = levels(sce.o$ClusterName), limits = levels(sce.o$ClusterName)) +
	labs( y = "TSNE-2", x = "TSNE-1") +
	ggtitle(str_c("Dentategyrus Lamanno (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/dentategyrus/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

Dentategyrus_hochgerner.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.6, grid.resolution = 25,
																						 stream.L = 3, stream.min.L = 1, stream.res = 10, stream.width = 0.4, point.size = 3.01,
																						 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Dentategyrus Hochgerner (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/erythroid/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

erythroid.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "celltype", color.alpha = 0.7, grid.resolution = 20,
																 stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4,
																 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$celltype_colors, name = "Cell type", labels = levels(sce.o$celltype), limits = levels(sce.o$celltype)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Gastrulation erythroid (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/forebrain/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "PCA"), sce.o)
sce.o$Clusters <- factor(sce.o$Clusters)
forebrain.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "PCA", color_by = "Clusters", color.alpha = 0.7, grid.resolution = 20,
																	stream.L = 4, stream.min.L = 1, stream.res = 10, stream.width = 0.4,
																	arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = brewer.pal(nlevels(factor(sce.o$Clusters)), "Set3"), name = "Cell type",
										 labels = c("Radial glia", "Radial glia",
										 					 "Neuroblast", "Neuroblast",
										 					 "Immature neuron", #"Immature neuron",
										 					 "Neuron", "Neuron"), limits = levels(sce.o$Clusters)) +
	labs( y = "PCA-2", x = "PCA-1") +
	ggtitle(str_c("Forebrain (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/gastrulation_e75/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

gastrulation_e75.p <-  plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "celltype", color.alpha = 0.7, grid.resolution = 20,
																	stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4, point.size = 2.01,
																	arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = names(table(sce.o$colour)), name = "Cell type", labels = levels(sce.o$celltype), limits = levels(sce.o$celltype)) +
	guides(color=guide_legend(ncol=1)) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Gastrulation E7.5 (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key.size = unit(5, "pt"))


sce.o <- qread(here::here("data/pbmc68k/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "TSNE"), sce.o)

pbmc68k.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "TSNE", color_by = "celltype", color.alpha = 0.8, grid.resolution = 20,
																	stream.L = 4, stream.min.L = 1, stream.res = 10, stream.width = 0.4, point.size = 2.01,
																	arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = brewer.pal(nlevels(factor(sce.o$celltype)), "Set3"), name = "Cell type",
										 labels = levels(factor(sce.o$celltype)), limits = levels(factor(sce.o$celltype))) +
	guides(color=guide_legend(nrow = 6)) + 
	labs( y = "TSNE-2", x = "TSNE-1") +
	ggtitle(str_c("PBMC68k (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.key.size = unit(6, "pt"))

gc()

### steady
sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

fucci_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "fucci_time", color.alpha = 0.7, grid.resolution = 13,
																		 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.4,
																		 arrow.angle = 18, arrow.length = 0.45) +
	scale_color_gradientn(colours = c("#2E22EA","#9E3DFB","#F86BE2",
																		"#FCCE7B","#C4E416","#4BBA0F",
																		"#447D87","#2C24E9"),name = "FUCCI", limits = range(sce.o$fucci_time, na.rm = TRUE)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle("FUCCI (steady)") +
	.theme_noframe + theme(legend.position = c(0, 0),
												 legend.justification = c(0, 0),
												 legend.key.size = unit(6, "pt"))

sce.o <- qread(here::here("data/Pancreas/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

pancreas_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = 15,
																				stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.4,
																				arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key.size = unit(6, "pt"))



sce.o <-qread(here::here("data/chromaffin/steady.o.qs"))
colors.v <- levels(factor(sce.o$cell_type))
sce.o$cell_type <- fct_recode(factor(sce.o$cell_type), "SCPs" = "#0066FFFF", "chromaffin" = "#00FF66FF",
															"sympathoblasts" = "#CC00FFFF", "non-labelled" = "#CCFF00FF", "bridge" = "#FF0000FF")
# reducedDim(sce.o, "PCA") <- reducedDim(sce.o, "PCA")[, 1:2]
embedded <- embedVelocity(reducedDim(sce.o, "PCA"), sce.o)


chromaffin_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "PCA", color_by = "cell_type", color.alpha = 0.7, grid.resolution = 20,
																					stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4, point.size = 5.01,
																					arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = colors.v, name = "Cell type", labels = levels(sce.o$cell_type), limits = levels(sce.o$cell_type)) +
	labs( y = "PCA-2", x = "PCA-1") +
	ggtitle(str_c("Chromaffin (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/bonemarrow/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "X_tsne"), sce.o)

Bonemarrow_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "X_tsne", color_by = "clusters", color.alpha = 0.7, grid.resolution = 20,
																					stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4,
																					arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "TSNE-2", x = "TSNE-1") +
	ggtitle(str_c("Bonemarrow (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/dentategyrus_lamanno/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "TSNE"), sce.o)
sce.o$ClusterName <- factor(sce.o$ClusterName, levels = c("RadialGlia", "RadialGlia2", "ImmAstro", 
																													"GlialProg", "OPC", "nIPC", "Nbl1", "Nbl2",
																													"ImmGranule1", "ImmGranule2", "Granule", "CA", 
																													"CA1-Sub", "CA2-3-4"))
Dentategyrus_lamanno_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "TSNE", color_by = "ClusterName", color.alpha = 0.7, grid.resolution = 20,
																										stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4, point.size = 2.01,
																										arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$cluster_colors, name = "Cell type", labels = levels(sce.o$ClusterName), limits = levels(sce.o$ClusterName)) +
	labs( y = "TSNE-2", x = "TSNE-1") +
	ggtitle(str_c("Dentategyrus Lamanno (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/dentategyrus/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

Dentategyrus_hochgerner_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.6, grid.resolution = 25,
																											 stream.L = 3, stream.min.L = 1, stream.res = 10, stream.width = 0.4, point.size = 3.01,
																											 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Dentategyrus Hochgerner (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/erythroid/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

erythroid_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "celltype", color.alpha = 0.7, grid.resolution = 20,
																				 stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4,
																				 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$celltype_colors, name = "Cell type", labels = levels(sce.o$celltype), limits = levels(sce.o$celltype)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Gastrulation erythroid (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/forebrain/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "PCA"), sce.o)
sce.o$Clusters <- factor(sce.o$Clusters)
forebrain_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "PCA", color_by = "Clusters", color.alpha = 0.7, grid.resolution = 20,
																				 stream.L = 4, stream.min.L = 1, stream.res = 10, stream.width = 0.4,
																				 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = brewer.pal(nlevels(factor(sce.o$Clusters)), "Set3"), name = "Cell type",
										 labels = c("Radial glia", "Radial glia",
										 					 "Neuroblast", "Neuroblast",
										 					 "Immature neuron", #"Immature neuron",
										 					 "Neuron", "Neuron"), limits = levels(sce.o$Clusters)) +
	labs( y = "PCA-2", x = "PCA-1") +
	ggtitle(str_c("Forebrain (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/gastrulation_e75/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

gastrulation_e75_steady.p <-  plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "celltype", color.alpha = 0.7, grid.resolution = 20,
																								 stream.L = 4, stream.min.L = 0.1, stream.res = 10, stream.width = 0.4, point.size = 2.01,
																								 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = names(table(sce.o$colour)), name = "Cell type", labels = levels(sce.o$celltype), limits = levels(sce.o$celltype)) +
	guides(color=guide_legend(ncol=1)) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Gastrulation E7.5 (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key.size = unit(5, "pt"))


sce.o <- qread(here::here("data/pbmc68k/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "TSNE"), sce.o)

pbmc68k_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "TSNE", color_by = "celltype", color.alpha = 0.8, grid.resolution = 20,
																			 stream.L = 4, stream.min.L = 1, stream.res = 10, stream.width = 0.4, point.size = 2.01,
																			 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = brewer.pal(nlevels(factor(sce.o$celltype)), "Set3"), name = "Cell type",
										 labels = levels(factor(sce.o$celltype)), limits = levels(factor(sce.o$celltype))) +
	guides(color=guide_legend(nrow = 6)) + 
	labs( y = "TSNE-2", x = "TSNE-1") +
	ggtitle(str_c("PBMC68k (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.key.size = unit(6, "pt"))








mp <- plot_grid( forebrain.p, forebrain_steady.p, chromaffin.p, chromaffin_steady.p,
								 fucci.p, fucci_steady.p, Bonemarrow.p, Bonemarrow_steady.p, 
								 pancreas.p, pancreas_steady.p, Dentategyrus_hochgerner.p, Dentategyrus_hochgerner_steady.p,
								 Dentategyrus_lamanno.p, Dentategyrus_lamanno_steady.p, erythroid.p, erythroid_steady.p,
								 gastrulation_e75.p, gastrulation_e75_steady.p, pbmc68k.p, pbmc68k_steady.p,
								 byrow = TRUE, 
								 nrow = 5, ncol = 4, label_size = 10, labels = as.vector(rbind(letters[1:10], rep("", 10))))
save_plot(here::here("figs", "sfigs", "sfig.realdata_steamline.pdf"), mp,
					base_height = 2 *1.2 , base_width = 2*1.4 , nrow = 5 * 1.5, ncol = 4 * 1.5, device = cairo_pdf)





