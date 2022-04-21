rm(list=ls())
source(here::here("scripts/utils.R"))



scv <- reticulate::import("scvelo")

get_umap <- function(noise_level = 1, n_obs = 500L, n_vars = 10L) {
	
	adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
	sce.o <- zellkonverter::AnnData2SCE(adata)
	ncomp <- as.integer( min(30L, n_vars - 1L))
	scv$pp$pca(adata, n_comps = ncomp)
	set.seed(100)
	umap.o <- umap::umap(adata$obsm["X_pca"][, seq_len(ncomp)], min_dist = 0.5, random_state = 123)
	reducedDim(sce.o, "UMAP") <- umap.o$layout
	
	p <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = str_c("Noise level ", noise_level), x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "True t", point.size = 4.01, point.alpha = 0.6) + 
		.theme_noframe +
		theme(legend.position = "none",
					plot.margin = unit(c(1, 1, 2, 2), "pt"))
	

	return(p)
}

plot.lp <- lapply(seq_len(10), get_umap)
plot.lp[[1]] <- plot.lp[[1]] + theme(legend.position = c(0, 0), legend.justification = c(0, 0))
	
mp <- plot_grid(plotlist = plot.lp,
								nrow = 4, ncol = 3, label_size = 10, labels = c("auto"))

save_plot(here::here("figs", "sfigs", "sfig.sim_umap.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 4, ncol = 3, device = cairo_pdf)




