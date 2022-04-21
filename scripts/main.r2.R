rm(list=ls())
source(here::here("scripts/utils.R"))

make_sim <- function(noise_level = 1, n_obs = 500L, n_vars = 10L) {
	scv <- reticulate::import("scvelo")
	adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
	ncomp <- as.integer( min(30L, n_vars - 1L))
	scv$pp$pca(adata, n_comps = ncomp)
	scv$pp$neighbors(adata)
	scv$pp$moments(adata)
	
	### get dynamical
	scv$tl$recover_dynamics(adata)
	scv$tl$velocity(adata, mode = "dynamical")
	sce.o <- zellkonverter::AnnData2SCE(adata)
	return(sce.o)
}

new_nvar.v <- 10L
noise_level <- 5

sce.o <- make_sim(noise_level = noise_level, n_obs = 500L, n_vars = new_nvar.v)
i <- which(!is.na(rowData(sce.o)$fit_t_))[1]
tmp.m <- dynamical_thet(alpha = rowData(sce.o)$true_alpha[i], beta = rowData(sce.o)$true_beta[i],
												gamma = rowData(sce.o)$true_gamma[i], scaling = rowData(sce.o)$true_scaling[i],
												t_ = rowData(sce.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = sce.o$true_t)
tmp.df <- data.frame(x = sce.o$true_t, s = tmp.m[, 2], v = tmp.m[, 4],
										 est = assay(sce.o, "fit_t")[i, ],
										 ms = assay(sce.o, "Ms")[i, ],
										 mu = assay(sce.o, "Mu")[i, ],
										 dynamical_v = assay(sce.o, "velocity")[i,])

true_t.p <- plotTLoess(x = tmp.df$x, y = tmp.df$ms, x_lab = "True t", x.pos = 0.55, y.pos = 0.9, point.size = 5.01,
											 color_var = tmp.df$x, y_lab = "Ms", title = str_c("Noise level ", noise_level),
											 scale_color = scale_color_gradientn(name = "t", limits = range(tmp.df$x), 
											 																		colors = c("#440154", "#482576", "#414487","#35608D",
											 																							 "#2A788E", "#21908C","#22A884", "#43BF71",
											 																							 "#7AD151", "#BBDF27", "#FDE725"))) +
	theme(legend.position = "none",
				plot.margin = unit(c(3, 3, 4, 4), "pt"))
true_ut.p <- plotTLoess(x = tmp.df$x, y = tmp.df$mu, x_lab = "True t", x.pos = 0.55, y.pos = 0.9, point.size = 5.01,
												color_var = tmp.df$x, y_lab = "Mu", title = str_c("Noise level ", noise_level),
												scale_color = scale_color_gradientn(name = "t", limits = range(tmp.df$x), 
																														colors = c("#440154", "#482576", "#414487","#35608D",
																																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																																			 "#7AD151", "#BBDF27", "#FDE725"))) +
	theme(legend.position = "none",
				plot.margin = unit(c(3, 3, 4, 4), "pt"))

est_t.p <- plotTLoess(x = tmp.df$est, y = tmp.df$ms, x_lab = "Estimated latent time", x.pos = 0.15, y.pos = 0.9,point.size = 5.01,
											color_var = tmp.df$x, y_lab = "Ms", title = str_c("Noise level ", noise_level),
											scale_color = scale_color_gradientn(name = "t", limits = range(tmp.df$x), 
																													colors = c("#440154", "#482576", "#414487","#35608D",
																																		 "#2A788E", "#21908C","#22A884", "#43BF71",
																																		 "#7AD151", "#BBDF27", "#FDE725"))) +
	theme(legend.position = "none",
				plot.margin = unit(c(3, 3, 4, 4), "pt"))
est_ut.p <- plotTLoess(x = tmp.df$est, y = tmp.df$mu, x_lab = "Estimated latent time", x.pos = 0.15, y.pos = 0.9,point.size = 5.01,
											 color_var = tmp.df$x, y_lab = "Mu", title = str_c("Noise level ", noise_level),
											 scale_color = scale_color_gradientn(name = "t", limits = range(tmp.df$x), 
											 																		colors = c("#440154", "#482576", "#414487","#35608D",
											 																							 "#2A788E", "#21908C","#22A884", "#43BF71",
											 																							 "#7AD151", "#BBDF27", "#FDE725"))) +
	theme(legend.position = "none",
				plot.margin = unit(c(3, 3, 4, 4), "pt"))



panel_abcd.p <- plot_grid(true_t.p, est_t.p, true_ut.p,  est_ut.p,
												nrow = 1, ncol = 4, label_size = 10, labels = c("auto"))





### add real data r2
tmp.df <- qread(here::here("data/realdata_r2.df.qs"))

tmp.df %<>% pivot_longer(
	cols = starts_with("r2"),
	names_to = "su",
	names_prefix = "r2",
	values_to = "r2"
)


## get orders
tmp.df$name <- factor(tmp.df$name, levels = (tmp.df %>% filter(su == "u") %>% group_by(name) %>%
																						 	summarize(medians = median(r2, na.rm = TRUE)) %>% arrange(desc(medians)))$name)

r2_real.p <- ggplot(tmp.df, aes(x = name, y = r2, dodge = su, fill = su)) +
	geom_boxplot(color = "#377EB8",outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.5) +
	scale_fill_manual(values = c("white", "grey80"), name = "S/U") + 
	scale_x_discrete(labels = c("Forebrain\n(10x)", "Chromaffin\n(SMART-seq2)", "FUCCI\n(SMART-seq2)", "Bonemarrow\n(10x)",
															"Pancreas\n(10x)", "Dentategyrus\nHochgerner (10x)", "Dentategyrus\nLaManno(10x)", 
															"Gastrulation\nerythroid (10x)", 
															"Gastrulation\nE7.5 (10x)", "PBMC68k\n(10x)")) +
	annotate("rect", xmin= 9.6, xmax= 10.4, ymin= 0 , ymax= 1, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	labs( x = "", 
				y =  expression(italic(R)^2), 
				title = "Explained variations in real data") +
	theme(legend.position = "none", 
				axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 0.5, angle = 30, color = c(rep('black', 9), "red")),
				axis.title.x = element_blank())






panel_abcde.p <- plot_grid(panel_abcd.p, r2_real.p,
													rel_heights = c(1 / 1.4, 1.1),
												nrow = 2, ncol = 1, label_size = 10, labels = c("", "e"))


save_plot(here::here("figs", "main", "main.r2.pdf"), panel_abcde.p,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 3, device = cairo_pdf)






