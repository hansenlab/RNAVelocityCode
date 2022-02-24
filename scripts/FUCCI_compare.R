rm(list=ls())
source(here::here("scripts/utils.R"))

steady.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))
dynamical.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))

### embeddings vectors
embedded <- embedVelocity(reducedDim(steady.o, "UMAP"), steady.o)
grid.df <- gridVectors(reducedDim(steady.o, "UMAP"), embedded)

steady_vec.p <- plotEmbFucci(steady.o, "UMAP", title = "Steady model", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)

embedded <- embedVelocity(reducedDim(dynamical.o, "UMAP"), dynamical.o)
grid.df <- gridVectors(reducedDim(dynamical.o, "UMAP"), embedded)

dynamical_vec.p <- plotEmbFucci(dynamical.o, "UMAP", title = "Dynamical model", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)



### idx for velocity genes for both
gene.v <- intersect(rownames(steady.o)[rowData(steady.o)$velocity_genes],
									rownames(dynamical.o)[rowData(dynamical.o)$velocity_genes])

v_steady.m <- assay(steady.o[gene.v, ], "velocity")
v_dynamical.m <- assay(dynamical.o[gene.v, ], "velocity")

inconsistency.v <- sapply(seq_len(nrow(v_steady.m)), function(i) mean(sign(v_steady.m[i, ] * v_dynamical.m[i, ]) == -1))
pcc.v <- sapply(seq_len(nrow(v_steady.m)), function(i) cor(v_steady.m[i, ], v_dynamical.m[i, ]))

names(inconsistency.v) <- gene.v
names(pcc.v) <- gene.v


### scatter plot
tmp.df <- data.frame(x = inconsistency.v, y = pcc.v)

incon_pcc.scat.p <- ggplot(data = tmp.df, aes(x = x, y = y)) + 
	geom_scattermore(pointsize = 2.01,  alpha = 0.8) +
	geom_vline(xintercept = median(inconsistency.v), color = "blue", alpha = 0.5) +
	geom_hline(yintercept = median(pcc.v), color = "blue", alpha = 0.5) +
	scale_x_continuous(breaks = seq(0, 0.8, 0.2), labels = c(0, 20, 40, 60, 80), name = "Percatage of \ninconsistent directions(%)") +
	labs(y = "PCC", title = "Comparison of velocities from two models")



#### plot 3 genes
ll_max.gene <- gene.v[which.max(rowData(dynamical.o[gene.v, ])$fit_likelihood)]
incon_small.gene <- gene.v[which.min(inconsistency.v)]
pcc_max.gene <- gene.v[which.max(pcc.v)]


g <- ll_max.gene
g1_com.p <- plot_comV(steady.o, dynamical.o, gene = g)
g1_steady.p <- plot_steadymums(steady.o, gene = g, color_by = "fucci_time", title = str_c(g, " - steady model"), color.name = "FUCCI")
g1_dynamical.p <- plot_dynamicalmums(dynamical.o, gene = g, color_by = "fucci_time",
									 title = str_c(g, " dynamical (LL:", sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood),
									 							")"), color.name = "FUCCI")

g <- pcc_max.gene
g2_com.p <- plot_comV(steady.o, dynamical.o, gene = g, x.pos = 0.75)
g2_steady.p <- plot_steadymums(steady.o, gene = g, color_by = "fucci_time", title = str_c(g, " - steady model"), color.name = "FUCCI")
g2_dynamical.p <- plot_dynamicalmums(dynamical.o, gene = g, color_by = "fucci_time",
																		 title = str_c(g, " dynamical (LL:", sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood),
																		 							")"), color.name = "FUCCI")

g <- incon_small.gene
g3_com.p <- plot_comV(steady.o, dynamical.o, gene = g, x.pos = 0.7)
g3_steady.p <- plot_steadymums(steady.o, gene = g, color_by = "fucci_time", title = str_c(g, " - steady model"), color.name = "FUCCI")
g3_dynamical.p <- plot_dynamicalmums(dynamical.o, gene = g, color_by = "fucci_time",
																		 title = str_c(g, " dynamical (LL:", sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood),
																		 							")"), color.name = "FUCCI")

mp <- wrap_plots(list(steady_vec.p, dynamical_vec.p, incon_pcc.scat.p,
											g1_com.p, g1_steady.p, g1_dynamical.p,
											g2_com.p, g2_steady.p, g2_dynamical.p,
											g3_com.p, g3_steady.p, g3_dynamical.p), ncol = 3)

save_plot(here::here("figs", "FUCCI", str_c( "compare_model.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 4, ncol = 3, device = cairo_pdf)



