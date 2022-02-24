rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))

sce.o <- sce.o[, which(!is.na(sce.o$fucci_time))]



rowData.df <- rowData(sce.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

t.m <- assay(sce.o, "fit_t")[rownames(rowData.df)[1:500], ]
min.v <- rowMins(t.m)
t.m <- sweep(t.m, MARGIN = 1, min.v, FUN = "-")
max.v <- rowMaxs(t.m)
t.m <- sweep(t.m, MARGIN = 1, max.v, FUN = "/")

order.idx <- order(sce.o$fucci_time)


library(ComplexHeatmap)
library(seriation)
o1 = seriate(dist(t.m[, order.idx]), method = "TSP")

fucci.p <- Heatmap(matrix(sort(sce.o$fucci_time), nrow = 1),
									 cluster_columns = FALSE,
									 cluster_rows = FALSE,
									 circlize::colorRamp2(breaks = seq(from = 0, to = 1, length.out = 8), c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9")),
									 name = "FUCCI",
									 column_title = "Scaled latent time of cells sorted by FUCCI pseudotime",
									 height = 0.3
									 )
heat.p <- Heatmap(t.m[, order.idx], name = "Latent\ntime", col = circlize::colorRamp2(breaks = seq(from = 0, to = 1, length.out = 11), c("#440154", "#482576", "#414487",
																																																																		 "#35608D", "#2A788E", "#21908C",
																																																																		 "#22A884", "#43BF71", "#7AD151",
																																																																		 "#BBDF27", "#FDE725"), space = "RGB"),
				column_title = NULL,
				row_title = "Top 500 genes (high LL)",
				clustering_distance_rows = "pearson",
				row_order = get_order(o1),
				cluster_columns = FALSE,
				cluster_rows = FALSE,
				row_dend_reorder = TRUE,
				show_row_dend = FALSE,
				show_row_names = FALSE,
				show_column_names = FALSE,
				height = 10
				)

ht_list <- fucci.p %v% heat.p

Cairo::CairoPDF(here::here("figs", "sfigs", "sfig.fucci_t.pdf"), width = 8, height = 4)
draw(ht_list)
dev.off()

Cairo::CairoJPEG(here::here("figs", "sfigs", "sfig.fucci_t.jpeg"), width = 960, height = 480, quality = 200)
draw(ht_list)
dev.off()

# p <- grid::grid.grabExpr(draw(ht_list), width = 8, height = 7)
# 
# 
# save_plot(here::here("figs", "sfigs", "sfig.fucci_t.pdf"), p,
# 					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 3, ncol = 3.2, device = cairo_pdf)
# 



# 
# pheatmap::pheatmap(t.m[, order.idx], cluster_rows = TRUE, cluster_cols = FALSE, filename = "test.jpg",
# 									 clustering_method = "ward.D2",
# 									 show_rownames = FALSE, show_colnames = FALSE,
# 									 width = 9, height = 6,
# 									 color = colorRampPalette(c("#440154", "#482576", "#414487",
# 									 													 "#35608D", "#2A788E", "#21908C",
# 									 													 "#22A884", "#43BF71", "#7AD151",
# 									 													 "#BBDF27", "#FDE725"))(100))
# 
# pheatmap(data.m[order(rank.v), order.idx], cluster_rows = FALSE, cluster_cols = FALSE, #breaks = seq(from = -.85, to = .85, length.out = 101),
# 				 main = str_c("mNeurosphere (n.genes = 500; n.cells=", ncol(neurosphere.o), ")"), show_rownames = FALSE, show_colnames = FALSE, fontsize = 8,
# 				 clustering_method = "ward.D2", #clustering_distance_rows = dist.o,  clustering_distance_cols = dist.o,
# 				 annotation_col = pheno.df, annotation_colors = list(SchwabeCC = setNames(ccColors.v, ccLabels.v),
# 				 																										"PCA\u03B8.\u03C0" = setNames(brewer.pal(9, "Blues")[c(2, 4, 6, 8)], as.character(c(0.5, 1, 1.5, 2)))),
# 				 breaks = seq(-2, 2, length.out = 100),
# 				 silent = TRUE,
# 				 width = 9, height = 6)


