rm(list=ls())
source(here::here("scripts/utils.R"))


sce1.o <- qread(here::here("data/neurosphere_wt1/dynamical.o.qs"))
sce2.o <- qread(here::here("data/neurosphere_wt2/dynamical.o.qs"))

genes.v <- intersect(rownames(sce1.o), rownames(sce2.o))
data.m <- cbind(assay(sce1.o, "spliced")[genes.v, ], assay(sce2.o, "spliced")[genes.v, ])
merged.m <- seuratIntegrate(count.m = data.m, batch = c(rep("WT1", ncol(sce1.o)), rep("WT2", ncol(sce2.o))))
pca.m <- calculatePCA(merged.m, ncomponents = 30)
library(umap)
umap.o <- umap::umap(pca.m, min_dist = 0.5)
plot(umap.o$layout, col = factor(c(rep("WT1", ncol(sce1.o)), rep("WT2", ncol(sce2.o)))))
umap.m <- umap.o$layout

reducedDim(sce1.o, "matched.UMAP") <- umap.m[seq_len(ncol(sce1.o)), ]
reducedDim(sce2.o, "matched.UMAP") <- umap.m[seq_len(ncol(sce2.o)) + ncol(sce1.o), ]

res <- 15
### plot
embedded1 <- embedVelocity(reducedDim(sce1.o, "matched.UMAP"), sce1.o)

wt1.p <- plotVelocityStream(sce1.o, embedded1[, 1:2], use.dimred = "matched.UMAP", color_by = "sample", color.alpha = 0.7, grid.resolution = res,
									 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
									 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = "#66C2A5", name = "Cell_type", guide = "none") +
	# annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.3, color= "#E41A1C", fill= NA, linetype = "solid", size = 0.7) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Neurosphere WT1")) +
	.theme_noframe + theme(legend.position = "none")

embedded2 <- embedVelocity(reducedDim(sce2.o, "matched.UMAP"), sce2.o)

wt2.p <- plotVelocityStream(sce2.o, embedded2[, 1:2], use.dimred = "matched.UMAP", color_by = "sample", color.alpha = 0.7, grid.resolution = res,
														stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
														arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = "#FC8D62", name = "Cell_type", guide = "none") +
	# annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.3, color= "#E41A1C", fill= NA, linetype = "solid", size = 0.7) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Neurosphere WT2")) +
	.theme_noframe + theme(legend.position = "none")


mp <- plot_grid(wt1.p, wt2.p,
								nrow = 1, ncol = 2, label_size = 10, labels = c("a", "b"))

save_plot(here::here("figs", "main", "main.neurosphere.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 1.5, ncol = 3, device = cairo_pdf)



