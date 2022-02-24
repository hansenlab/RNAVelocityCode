rm(list=ls())
source(here::here("scripts/utils.R"))

make_sim <- function(noise_level = 1, n_obs = 500L, n_vars = 10L) {
	scv <- reticulate::import("scvelo")
	adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
	ncomp <- as.integer( min(30L, n_vars - 1L))
	scv$pp$pca(adata, n_comps = ncomp)
	
	### get dynamical
	scv$tl$recover_dynamics(adata)
	scv$tl$velocity(adata, mode = "dynamical")
	sce.o <- zellkonverter::AnnData2SCE(adata)
	return(sce.o)
}

panel_abc.lp <- lapply(c(1, 5, 10), function(noise_level) {
	sce.o <- make_sim(noise_level = noise_level, n_obs = 500L, n_vars = 10L)
	i <- 1
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
	
	est_t.p <- plotTLoess(x = tmp.df$est, y = tmp.df$ms, x_lab = "Estimated latent time", x.pos = 0.55, y.pos = 0.9,point.size = 5.01,
												color_var = tmp.df$x, y_lab = "Ms", title = str_c("Noise level ", noise_level),
												scale_color = scale_color_gradientn(name = "t", limits = range(tmp.df$x), 
																														colors = c("#440154", "#482576", "#414487","#35608D",
																																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																																			 "#7AD151", "#BBDF27", "#FDE725"))) +
		theme(legend.position = "none",
					plot.margin = unit(c(3, 3, 4, 4), "pt"))
	est_ut.p <- plotTLoess(x = tmp.df$est, y = tmp.df$mu, x_lab = "Estimated latent time", x.pos = 0.55, y.pos = 0.9,point.size = 5.01,
												color_var = tmp.df$x, y_lab = "Mu", title = str_c("Noise level ", noise_level),
												scale_color = scale_color_gradientn(name = "t", limits = range(tmp.df$x), 
																														colors = c("#440154", "#482576", "#414487","#35608D",
																																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																																			 "#7AD151", "#BBDF27", "#FDE725"))) +
		theme(legend.position = "none",
					plot.margin = unit(c(3, 3, 4, 4), "pt"))
	
	return(list(true_t.p, est_t.p, true_ut.p,  est_ut.p))
})

panel_abc.p <- plot_grid(plotlist = do.call(c, panel_abc.lp),
												nrow = 3, ncol = 4, label_size = 10, labels = c("a", "", "", "",
																																				"b", "", "","",
																																				"c", "", "",""))


### get all R2
r2.lm <- lapply(seq_len(10), function(noise_level) {
	sce.o <- make_sim(noise_level = noise_level, n_obs = 500L, n_vars = 10L)
	r2_true.v <- apply(assay(sce.o, "Ms"), MARGIN = 1, function(ms) fit_t_loess(x = sce.o$true_t, y = ms)$rsquared)
	r2_est.v <- sapply(seq_len(nrow(sce.o)), function(i) fit_t_loess(x = assay(sce.o, "fit_t")[i, ], y = assay(sce.o, "Ms")[i, ])$rsquared)
	r2u_true.v <- apply(assay(sce.o, "Mu"), MARGIN = 1, function(ms) fit_t_loess(x = sce.o$true_t, y = ms)$rsquared)
	r2u_est.v <- sapply(seq_len(nrow(sce.o)), function(i) fit_t_loess(x = assay(sce.o, "fit_t")[i, ], y = assay(sce.o, "Mu")[i, ])$rsquared)
	return(list(r2_true = r2_true.v, r2_est = r2_est.v, r2u_true = r2u_true.v, r2u_est = r2u_est.v))
})

r2_true.m <- do.call(cbind, lapply(r2.lm, "[[", "r2_true"))
r2_est.m <- do.call(cbind, lapply(r2.lm, "[[", "r2_est"))
r2u_true.m <- do.call(cbind, lapply(r2.lm, "[[", "r2u_true"))
r2u_est.m <- do.call(cbind, lapply(r2.lm, "[[", "r2u_est"))


# r2_300.lm <- lapply(seq_len(10), function(noise_level) {
# 	sce.o <- make_sim(noise_level = noise_level, n_obs = 800L, n_vars = 300L)
# 	r2_true.v <- apply(assay(sce.o, "Ms"), MARGIN = 1, function(ms) fit_t_loess(x = sce.o$true_t, y = ms)$rsquared)
# 	r2_est.v <- sapply(seq_len(nrow(sce.o)), function(i) fit_t_loess(x = assay(sce.o, "fit_t")[i, ], y = assay(sce.o, "Ms")[i, ])$rsquared)
# 	r2u_true.v <- apply(assay(sce.o, "Mu"), MARGIN = 1, function(ms) fit_t_loess(x = sce.o$true_t, y = ms)$rsquared)
# 	r2u_est.v <- sapply(seq_len(nrow(sce.o)), function(i) fit_t_loess(x = assay(sce.o, "fit_t")[i, ], y = assay(sce.o, "Mu")[i, ])$rsquared)
# 	return(list(r2_true = r2_true.v, r2_est = r2_est.v, r2u_true = r2u_true.v, r2u_est = r2u_est.v))
# })
# qsave(r2_300.lm, file = here::here("data/r2_300.lm.qs"))

r2_300.lm <- qread(here::here("data/r2_300.lm.qs"))


### panel d
colnames(r2_true.m) <- str_c("noise_", seq_len(10))
r2_true.df <- r2_true.m %>% as.data.frame() %>% add_column(type = "true", su = "s")
colnames(r2_est.m) <- str_c("noise_", seq_len(10))
r2_est.df <- r2_est.m %>% as.data.frame() %>% add_column(type = "est", su = "s")
colnames(r2u_true.m) <- str_c("noise_", seq_len(10))
r2u_true.df <- r2u_true.m %>% as.data.frame() %>% add_column(type = "true", su = "u")
colnames(r2u_est.m) <- str_c("noise_", seq_len(10))
r2u_est.df <- r2u_est.m %>% as.data.frame() %>% add_column(type = "est", su = "u")

tmp.df <- rbind(r2_true.df, r2_est.df, r2u_true.df, r2u_est.df) %>%
	pivot_longer(
		cols = starts_with("noise"),
		names_to = "noise",
		names_prefix = "noise_",
		values_to = "r2"
	)
tmp.df$noise <- factor(as.numeric(tmp.df$noise))

panel_d.p <- ggplot(tmp.df, aes(x = noise, y = r2, color = type,  dodge = type, fill = su)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.8) +
	scale_color_manual(values = c("#377EB8", "#FF7F00"), name = "Type") +
	scale_fill_manual(values = c("white", "grey80"), name = "S/U") + 
	labs( x = "Noise level", 
				y =  expression(italic(R)^2), 
				title = "10 genes and 500 cells") +
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0))



r2_300_true.m <- do.call(cbind, lapply(r2_300.lm, "[[", "r2_true"))
r2_300_est.m <- do.call(cbind, lapply(r2_300.lm, "[[", "r2_est"))
r2u_300_true.m <- do.call(cbind, lapply(r2_300.lm, "[[", "r2u_true"))
r2u_300_est.m <- do.call(cbind, lapply(r2_300.lm, "[[", "r2u_est"))

### panel d
colnames(r2_300_true.m) <- str_c("noise_", seq_len(10))
r2_300_true.df <- r2_300_true.m %>% as.data.frame() %>% add_column(type = "true", su = "s")
colnames(r2_300_est.m) <- str_c("noise_", seq_len(10))
r2_300_est.df <- r2_300_est.m %>% as.data.frame() %>% add_column(type = "est", su = "s")
colnames(r2u_300_true.m) <- str_c("noise_", seq_len(10))
r2u_300_true.df <- r2u_300_true.m %>% as.data.frame() %>% add_column(type = "true", su = "u")
colnames(r2u_300_est.m) <- str_c("noise_", seq_len(10))
r2u_300_est.df <- r2u_300_est.m %>% as.data.frame() %>% add_column(type = "est", su = "u")

tmp.df <- rbind(r2_300_true.df, r2_300_est.df, r2u_300_true.df, r2u_300_est.df) %>%
	pivot_longer(
		cols = starts_with("noise"),
		names_to = "noise",
		names_prefix = "noise_",
		values_to = "r2"
	)
tmp.df$noise <- factor(as.numeric(tmp.df$noise))

panel_e.p <- ggplot(tmp.df, aes(x = noise, y = r2, color = type,  dodge = type, fill = su)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.8) +
	scale_color_manual(values = c("#377EB8", "#FF7F00"), name = "Type") +
	scale_fill_manual(values = c("white", "grey80"), name = "S/U") + 
	labs( x = "Noise level", 
				y =  expression(italic(R)^2), 
				title = "300 genes and 800 cells") +
	theme(legend.position = "none")


panel_de.p <- plot_grid(panel_d.p, panel_e.p,
												nrow = 1, ncol = 2, label_size = 10, labels = c("d", "e"))


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






panel_abcdef.p <- plot_grid(panel_abc.p, panel_de.p, r2_real.p,
													rel_heights = c(3 / 1.4, 1, 1.1),
												nrow = 3, ncol = 1, label_size = 10, labels = c("", "", "f"))


save_plot(here::here("figs", "main", "main.r2.pdf"), panel_abcdef.p,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 4.5, ncol = 3, device = cairo_pdf)






