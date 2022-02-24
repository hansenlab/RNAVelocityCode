rm(list=ls())
source(here::here("scripts/utils.R"))

make_sim <- function(noise_level = 1, n_obs = 500L, n_vars = 10L) {
	### steady model
	scv <- reticulate::import("scvelo")
	adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
	ncomp <- as.integer( min(30L, n_vars - 1L))
	scv$pp$pca(adata, n_comps = ncomp)
	scv$tl$recover_dynamics(adata)
	scv$tl$velocity(adata, mode = "dynamical")
	# adata$var["velocity_genes"] <- TRUE
	# scv$tl$velocity_graph(adata)
	sce.o <- zellkonverter::AnnData2SCE(adata)

	
	
	### true v 2
	
	plot.lp <- lapply(seq_len(min(10, nrow(sce.o))), function(i) {
		
		tmp.m <- dynamical_thet(alpha = rowData(sce.o)$true_alpha[i], beta = rowData(sce.o)$true_beta[i],
														gamma = rowData(sce.o)$true_gamma[i], scaling = rowData(sce.o)$true_scaling[i],
														t_ = rowData(sce.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = sce.o$true_t)
		p1 <- plot_dynamicalmums(sce.o, rownm = as.character(i - 1), color_by = "true_t", hue.colors = c("#440154", "#482576", "#414487",
																																																		 "#35608D", "#2A788E", "#21908C",
																																																		 "#22A884", "#43BF71", "#7AD151",
																																																		 "#BBDF27", "#FDE725"),
														 title = str_c(i, " (LL:", sprintf("%.3f", rowData(sce.o[i, ])$fit_likelihood),")"), color.name = "t",
														 point.size = 3.01) +
			theme(legend.position = c(1, 0),
						legend.justification = c(1, 0))
		p2 <-  plotxy(sce.o$true_t, assay(sce.o, "Ms")[i, ], title = "Ms ~ t", x_lab = "t",
									y_lab = "Ms", color_by = "true_t", point.size = 3.01,
									colors.v = sce.o$true_t,
									hue.colors = c("#440154", "#482576", "#414487","#35608D",
																 "#2A788E", "#21908C","#22A884", "#43BF71",
																 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
			theme(legend.position = "none")
		p3 <-  plotxy(sce.o$true_t, tmp.m[, 2], title = "True s ~ t", x_lab = "t",
									y_lab = "True s", color_by = "true_t", point.size = 3.01,
									colors.v = sce.o$true_t,
									hue.colors = c("#440154", "#482576", "#414487","#35608D",
																 "#2A788E", "#21908C","#22A884", "#43BF71",
																 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
			theme(legend.position = "none")
		p4 <-  plotxy(tmp.m[, 2], assay(sce.o, "Ms")[i, ], title = "Ms ~ True s", x_lab = "True s",
									y_lab = "Ms", color_by = "true_t", point.size = 3.01,
									colors.v = sce.o$true_t,
									hue.colors = c("#440154", "#482576", "#414487","#35608D",
																 "#2A788E", "#21908C","#22A884", "#43BF71",
																 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
			theme(legend.position = "none")
		
		p5 <-  plotxy(sce.o$true_t, assay(sce.o, "velocity")[i,], title = "Dynamical v ~ t", x_lab = "t",
									y_lab = "Dynamical v", color_by = "true_t", point.size = 3.01,
									colors.v = sce.o$true_t,
									hue.colors = c("#440154", "#482576", "#414487","#35608D",
																 "#2A788E", "#21908C","#22A884", "#43BF71",
																 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
			theme(legend.position = "none")
		p6 <-  plotxy(sce.o$true_t, tmp.m[, 4], title = "True v ~ t", x_lab = "t",
									y_lab = "True v", color_by = "true_t", point.size = 3.01,
									colors.v = sce.o$true_t,
									hue.colors = c("#440154", "#482576", "#414487","#35608D",
																 "#2A788E", "#21908C","#22A884", "#43BF71",
																 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
			theme(legend.position = "none")
		p7 <-  plotxy(tmp.m[, 4], assay(sce.o, "velocity")[i,], title = "Dynamical v ~ True v", x_lab = "True v",
									y_lab = "Dynamical v", color_by = "true_t", point.size = 3.01,
									colors.v = sce.o$true_t,
									hue.colors = c("#440154", "#482576", "#414487","#35608D",
																 "#2A788E", "#21908C","#22A884", "#43BF71",
																 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
			theme(legend.position = "none")
		mp <- wrap_plots(list(p1, p2, p3, p4, p5, p6, p7), ncol = 7)
		
		
		return(mp)
	})
	
	
	mp <- wrap_plots(plot.lp, ncol =1)
	
	save_plot(here::here("figs", "sim_genes", str_c("nvars_", n_vars, "_nobs_", n_obs, "_noiselevel_",noise_level, ".pdf")), mp,
						base_height = 2, base_width = 2*1.4, nrow = 10, ncol = 7 , device = cairo_pdf)
	
	
}

make_sim(1)
make_sim(2)
make_sim(3)
make_sim(4)
make_sim(5)
make_sim(6)
make_sim(10)
make_sim(20)


make_sim(1, n_obs = 800L,  n_vars = 300L)
make_sim(2, n_obs = 800L,  n_vars = 300L)
make_sim(3, n_obs = 800L,  n_vars = 300L)
make_sim(4, n_obs = 800L,  n_vars = 300L)
make_sim(5, n_obs = 800L,  n_vars = 300L)
make_sim(6, n_obs = 800L,  n_vars = 300L)
make_sim(10, n_obs = 800L,  n_vars = 300L)
make_sim(20, n_obs = 800L,  n_vars = 300L)