rm(list=ls())
source(here::here("scripts/utils.R"))

noise_level = 3
n_obs = 500L
n_vars = 10L

scv <- reticulate::import("scvelo")
adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
ncomp <- as.integer( min(30L, n_vars - 1L))


scv$pp$pca(adata, n_comps = ncomp)
scv$pp$neighbors(adata)
scv$pp$moments(adata)
rotation.m <- adata$varm["PCs"]



### get dynamical
scv$tl$recover_dynamics(adata)
scv$tl$velocity(adata, mode = "dynamical")
dynamical.o <- zellkonverter::AnnData2SCE(adata)

### get steady
scv$tl$velocity(adata, mode = "deterministic")
steady.o <- zellkonverter::AnnData2SCE(adata)


###
i <- 1
tmp.m <- dynamical_thet(alpha = rowData(dynamical.o)$true_alpha[i], beta = rowData(dynamical.o)$true_beta[i],
												gamma = rowData(dynamical.o)$true_gamma[i], scaling = rowData(dynamical.o)$true_scaling[i],
												t_ = rowData(dynamical.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = dynamical.o$true_t)
tmp.df <- data.frame(x = dynamical.o$true_t, s = tmp.m[, 2], v = tmp.m[, 4], u =tmp.m[, 1],
										 est = assay(dynamical.o, "fit_t")[i, ],
										 ms = assay(dynamical.o, "Ms")[i, ],
										 mu = assay(dynamical.o, "Mu")[i, ],
										 spliced = assay(dynamical.o, "spliced")[i, ],
										 unspliced = assay(dynamical.o, "unspliced")[i, ],
										 steady_v = assay(steady.o, "velocity")[i,],
										 dynamical_v = assay(dynamical.o, "velocity")[i,])

mums.p <- plot_dynamicalmums(dynamical.o, rownm = as.character(i - 1), color_by = "true_t", hue.colors = c("#440154", "#482576", "#414487",
																																																					 "#35608D", "#2A788E", "#21908C",
																																																					 "#22A884", "#43BF71", "#7AD151",
																																																					 "#BBDF27", "#FDE725"),
														 title = str_c("Using learned KNN"), color.name = "True t",
														 point.size = 5.01) +
	geom_abline(slope = rowData(steady.o[i, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
	theme(legend.position = "none")


ms_t.p <- plotxy(tmp.df$x, tmp.df$ms, title = "Using learned KNN", x_lab = "True t.",
								 y_lab = "Ms", color_by = "true_t", point.size = 4.01,
								 colors.v = tmp.df$x,
								 hue.colors = c("#440154", "#482576", "#414487","#35608D",
								 							 "#2A788E", "#21908C","#22A884", "#43BF71",
								 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = s), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")

mu_t.p <- plotxy(tmp.df$x, tmp.df$mu, title = "Using learned KNN", x_lab = "True t.",
								 y_lab = "Mu", color_by = "true_t", point.size = 4.01,
								 colors.v = tmp.df$x,
								 hue.colors = c("#440154", "#482576", "#414487","#35608D",
								 							 "#2A788E", "#21908C","#22A884", "#43BF71",
								 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = u), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")

spliced_t.p <- plotxy(tmp.df$x, tmp.df$spliced, title = "Using learned KNN", x_lab = "True t.",
											y_lab = "Spliced", color_by = "true_t", point.size = 4.01,
											colors.v = tmp.df$x,
											hue.colors = c("#440154", "#482576", "#414487","#35608D",
																		 "#2A788E", "#21908C","#22A884", "#43BF71",
																		 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = s), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key.size = unit(6, "pt"))

unspliced_t.p <- plotxy(tmp.df$x, tmp.df$unspliced, title = "Using learned KNN", x_lab = "True t.",
												y_lab = "Unspliced", color_by = "true_t", point.size = 4.01,
												colors.v = tmp.df$x,
												hue.colors = c("#440154", "#482576", "#414487","#35608D",
																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																			 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = u), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")



dynamical_t.p <- plotxy(tmp.df$x, tmp.df$dynamical_v, title = "Using learned KNN", x_lab = "True t.",
												y_lab = "Dynamical v.", color_by = "true_t", point.size = 4.01,
												colors.v = tmp.df$x,
												hue.colors = c("#440154", "#482576", "#414487","#35608D",
																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																			 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = v), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")


# steady_v.p <- plotxy(tmp.df$v, tmp.df$steady_v, title = "Steady v.", x_lab = "True v.",
# 										 y_lab = "Steady v.", color_by = "true_t", point.size = 4.01,
# 										 colors.v = tmp.df$x,
# 										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
# 										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
# 										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
# 	.theme_noframe +
# 	theme(plot.title = element_blank(), legend.position = "none",
# 				plot.margin = unit(c(2, 2, 3, 3), "pt"))

dynamical_v.p <- plotxy(tmp.df$v, tmp.df$dynamical_v, title = "Using learned KNN", x_lab = "True v.",
												y_lab = "Dynamical v.", color_by = "true_t", point.size = 4.01,
												colors.v = tmp.df$x,
												hue.colors = c("#440154", "#482576", "#414487","#35608D",
																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																			 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	theme(legend.position = "none") +
	annotate(geom = "text", x = .percent_range(tmp.df$v, 0.5),
					 y = .percent_range(tmp.df$dynamical_v, 0.15), size = 2, hjust = 0, vjust = 1,
					 label = str_c("PCC: ", sprintf("%.3f", cor(tmp.df$v, tmp.df$dynamical_v))), parse = FALSE) +
	annotate(geom = "text", x = .percent_range(tmp.df$v, 0.5),
					 y = .percent_range(tmp.df$dynamical_v, 0.05), size = 2, hjust = 0, vjust = 1,
					 label = str_c("NRMSE: ", sprintf("%.3f", nrmse(tmp.df$v, tmp.df$dynamical_v))), parse = FALSE) 

est_t.p <- plotxy(tmp.df$x, tmp.df$est, title = "Using learned KNN", x_lab = "True t.",
									y_lab = "Estimated t.", color_by = "true_t", point.size = 5.01,
									colors.v = tmp.df$x,
									hue.colors = c("#440154", "#482576", "#414487","#35608D",
																 "#2A788E", "#21908C","#22A884", "#43BF71",
																 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	theme(legend.position = "none")




panel_a_h.p <- plot_grid(spliced_t.p, ms_t.p, unspliced_t.p,  mu_t.p, 
								mums.p, dynamical_t.p, dynamical_v.p,  est_t.p,
								nrow = 2, ncol = 4, label_size = 10, labels = "auto")

panel_efg.p <- plot_grid(
												 mums.p, dynamical_t.p, dynamical_v.p,  
												 nrow = 1, ncol = 3, label_size = 10, labels = "auto")

### use true KNN
adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
ncomp <- as.integer( min(30L, n_vars - 1L))

### get true s
sce.o <- zellkonverter::AnnData2SCE(adata)
s.m <- t(sapply(seq_len(nrow(sce.o)), function(i) {
	dynamical_thet(alpha = rowData(sce.o)$true_alpha[i], beta = rowData(sce.o)$true_beta[i],
								 gamma = rowData(sce.o)$true_gamma[i], scaling = rowData(sce.o)$true_scaling[i],
								 t_ = rowData(sce.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = sce.o$true_t)[, 2]
}))
dimnames(s.m) <- dimnames(sce.o)
assay(sce.o, "X") <- s.m
adata <- zellkonverter::SCE2AnnData(sce.o)

### learn KNN graph using true S
scv$pp$pca(adata, n_comps = ncomp)
scv$pp$neighbors(adata)
scv$pp$moments(adata)

### switch back and rerun PCA
sce.o <- zellkonverter::AnnData2SCE(adata)
reducedDim(sce.o, "true_pca") <- reducedDim(sce.o, "X_pca")
assay(sce.o, "X") <- assay(sce.o, "spliced")
adata <- zellkonverter::SCE2AnnData(sce.o)

scv$pp$pca(adata, n_comps = ncomp)
rotation.m <- adata$varm["PCs"]



### get dynamical
scv$tl$recover_dynamics(adata)
scv$tl$velocity(adata, mode = "dynamical")
dynamical.o <- zellkonverter::AnnData2SCE(adata)

### get steady
scv$tl$velocity(adata, mode = "deterministic")
steady.o <- zellkonverter::AnnData2SCE(adata)


###
i <- 1
tmp.m <- dynamical_thet(alpha = rowData(dynamical.o)$true_alpha[i], beta = rowData(dynamical.o)$true_beta[i],
												gamma = rowData(dynamical.o)$true_gamma[i], scaling = rowData(dynamical.o)$true_scaling[i],
												t_ = rowData(dynamical.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = dynamical.o$true_t)
tmp.df <- data.frame(x = dynamical.o$true_t, s = tmp.m[, 2], v = tmp.m[, 4], u =tmp.m[, 1],
										 est = assay(dynamical.o, "fit_t")[i, ],
										 ms = assay(dynamical.o, "Ms")[i, ],
										 mu = assay(dynamical.o, "Mu")[i, ],
										 spliced = assay(dynamical.o, "spliced")[i, ],
										 unspliced = assay(dynamical.o, "unspliced")[i, ],
										 steady_v = assay(steady.o, "velocity")[i,],
										 dynamical_v = assay(dynamical.o, "velocity")[i,])

mums.p <- plot_dynamicalmums(dynamical.o, rownm = as.character(i - 1), color_by = "true_t", hue.colors = c("#440154", "#482576", "#414487",
																																																					 "#35608D", "#2A788E", "#21908C",
																																																					 "#22A884", "#43BF71", "#7AD151",
																																																					 "#BBDF27", "#FDE725"),
														 title = str_c("Using true KNN"), color.name = "True t",
														 point.size = 5.01) +
	geom_abline(slope = rowData(steady.o[i, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
	theme(legend.position = "none")


ms_t.p <- plotxy(tmp.df$x, tmp.df$ms, title = "Using true KNN", x_lab = "True t.",
								 y_lab = "Ms", color_by = "true_t", point.size = 4.01,
								 colors.v = tmp.df$x,
								 hue.colors = c("#440154", "#482576", "#414487","#35608D",
								 							 "#2A788E", "#21908C","#22A884", "#43BF71",
								 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = s), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")

mu_t.p <- plotxy(tmp.df$x, tmp.df$mu, title = "Using true KNN", x_lab = "True t.",
								 y_lab = "Mu", color_by = "true_t", point.size = 4.01,
								 colors.v = tmp.df$x,
								 hue.colors = c("#440154", "#482576", "#414487","#35608D",
								 							 "#2A788E", "#21908C","#22A884", "#43BF71",
								 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = u), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")

spliced_t.p <- plotxy(tmp.df$x, tmp.df$spliced, title = "Using true KNN", x_lab = "True t.",
											y_lab = "Spliced", color_by = "true_t", point.size = 4.01,
											colors.v = tmp.df$x,
											hue.colors = c("#440154", "#482576", "#414487","#35608D",
																		 "#2A788E", "#21908C","#22A884", "#43BF71",
																		 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = s), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
				legend.key.size = unit(6, "pt"))

unspliced_t.p <- plotxy(tmp.df$x, tmp.df$unspliced, title = "Using true KNN", x_lab = "True t.",
												y_lab = "Unspliced", color_by = "true_t", point.size = 4.01,
												colors.v = tmp.df$x,
												hue.colors = c("#440154", "#482576", "#414487","#35608D",
																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																			 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = u), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")



dynamical_t.p <- plotxy(tmp.df$x, tmp.df$dynamical_v, title = "Using true KNN", x_lab = "True t.",
												y_lab = "Dynamical v.", color_by = "true_t", point.size = 4.01,
												colors.v = tmp.df$x,
												hue.colors = c("#440154", "#482576", "#414487","#35608D",
																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																			 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = v), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")


# steady_v.p <- plotxy(tmp.df$v, tmp.df$steady_v, title = "Steady v.", x_lab = "True v.",
# 										 y_lab = "Steady v.", color_by = "true_t", point.size = 4.01,
# 										 colors.v = tmp.df$x,
# 										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
# 										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
# 										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
# 	.theme_noframe +
# 	theme(plot.title = element_blank(), legend.position = "none",
# 				plot.margin = unit(c(2, 2, 3, 3), "pt"))

dynamical_v.p <- plotxy(tmp.df$v, tmp.df$dynamical_v, title = "Using true KNN", x_lab = "True v.",
												y_lab = "Dynamical v.", color_by = "true_t", point.size = 4.01,
												colors.v = tmp.df$x,
												hue.colors = c("#440154", "#482576", "#414487","#35608D",
																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																			 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	theme(legend.position = "none") +
	annotate(geom = "text", x = .percent_range(tmp.df$v, 0.5),
					 y = .percent_range(tmp.df$dynamical_v, 0.15), size = 2, hjust = 0, vjust = 1,
					 label = str_c("PCC: ", sprintf("%.3f", cor(tmp.df$v, tmp.df$dynamical_v))), parse = FALSE) +
	annotate(geom = "text", x = .percent_range(tmp.df$v, 0.5),
					 y = .percent_range(tmp.df$dynamical_v, 0.05), size = 2, hjust = 0, vjust = 1,
					 label = str_c("NRMSE: ", sprintf("%.3f", nrmse(tmp.df$v, tmp.df$dynamical_v))), parse = FALSE) 


est_t.p <- plotxy(tmp.df$x, tmp.df$est, title = "Using true KNN", x_lab = "True t.",
									y_lab = "Estimated t.", color_by = "true_t", point.size = 5.01,
									colors.v = tmp.df$x,
									hue.colors = c("#440154", "#482576", "#414487","#35608D",
																 "#2A788E", "#21908C","#22A884", "#43BF71",
																 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	theme(legend.position = "none")




panel_i_p.p <- plot_grid(spliced_t.p, ms_t.p, unspliced_t.p,  mu_t.p, 
												 mums.p, dynamical_t.p, dynamical_v.p,  est_t.p,
												 nrow = 2, ncol = 4, label_size = 10, labels = letters[9:16])

panel_mno.p <- plot_grid(
												 mums.p, dynamical_t.p, dynamical_v.p,  
												 nrow = 1, ncol = 3, label_size = 10, labels = letters[4:6])


mp <- plot_grid(panel_a_h.p, panel_i_p.p,
								nrow = 2, ncol = 1, label_size = 10, labels = NULL)

save_plot(here::here("figs", "sfigs", "sfig.noise_genes.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 4 * 0.8, ncol = 3, device = cairo_pdf)

mp <- plot_grid(panel_efg.p, panel_mno.p,
								nrow = 2, ncol = 1, label_size = 10, labels = NULL)

save_plot(here::here("figs", "main", "main.noise_genes.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2 * 0.8, ncol = 3, device = cairo_pdf)









