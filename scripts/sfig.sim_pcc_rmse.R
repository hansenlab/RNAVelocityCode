rm(list=ls())
source(here::here("scripts/utils.R"))

scv <- reticulate::import("scvelo")

make_sim_eva <- function(noise_level = 1, n_obs = 500L, n_vars = 10L) {
	
	adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
	sce.o <- zellkonverter::AnnData2SCE(adata)
	ncomp <- as.integer( min(30L, n_vars - 1L))
	scv$pp$pca(adata, n_comps = ncomp)
	scv$pp$neighbors(adata)
	scv$pp$moments(adata)
	
	
	### get dynamical
	scv$tl$recover_dynamics(adata)
	scv$tl$velocity(adata, mode = "dynamical")
	dynamical.o <- zellkonverter::AnnData2SCE(adata)
	
	### get steady
	scv$tl$velocity(adata, mode = "deterministic")
	steady.o <- zellkonverter::AnnData2SCE(adata)
	
	
	true_v.m <- t(sapply(seq_len(nrow(dynamical.o)), function(i) {
		tmp.m <- dynamical_thet(alpha = rowData(dynamical.o)$true_alpha[i], beta = rowData(dynamical.o)$true_beta[i],
														gamma = rowData(dynamical.o)$true_gamma[i], scaling = rowData(dynamical.o)$true_scaling[i],
														t_ = rowData(dynamical.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = dynamical.o$true_t)
		return(tmp.m[, "vt"])
	}))
	dimnames(true_v.m) <- dimnames(dynamical.o)
	
	true_s.m <- t(sapply(seq_len(nrow(dynamical.o)), function(i) {
		tmp.m <- dynamical_thet(alpha = rowData(dynamical.o)$true_alpha[i], beta = rowData(dynamical.o)$true_beta[i],
														gamma = rowData(dynamical.o)$true_gamma[i], scaling = rowData(dynamical.o)$true_scaling[i],
														t_ = rowData(dynamical.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = dynamical.o$true_t)
		return(tmp.m[, "st"])
	}))
	dimnames(true_s.m) <- dimnames(dynamical.o)
	
	### get PCC
	# pcc_ms.v <- sapply(seq_len(nrow(dynamical.o)), function(i) cor(assay(steady.o, "Ms")[i,], true_s.m[i,]))
	# rmse_ms.v <- sapply(seq_len(nrow(dynamical.o)), function(i) nrmse(true_s.m[i,], assay(steady.o, "Ms")[i,]) )
	pcc_steady.v <- sapply(seq_len(nrow(dynamical.o)), function(i) cor(assay(steady.o, "velocity")[i,], true_v.m[i,]))
	rmse_steady.v <- sapply(seq_len(nrow(dynamical.o)), function(i) nrmse(true_v.m[i,], assay(steady.o, "velocity")[i,]) )
	pcc_dynamical.v <- sapply(seq_len(nrow(dynamical.o)), function(i) cor(assay(dynamical.o, "velocity")[i,], true_v.m[i,]))
	rmse_dynamical.v <- sapply(seq_len(nrow(dynamical.o)), function(i) nrmse(true_v.m[i,], assay(dynamical.o, "velocity")[i,]) )
	dynamical_pct.v <- sum(!is.na(rmse_dynamical.v)) / nrow(dynamical.o)
	
	get_len_diff <- function(est.m, true.m) {
		row.idx <- which(!rowAnyNAs(est.m))
		est_len.v <- sqrt(colSums(est.m[row.idx, ] ^ 2))
		true_len.v <- sqrt(colSums(true.m[row.idx, ] ^ 2))
		return(abs(est_len.v - true_len.v))
	}
	get_len_pcc <- function(est.m, true.m) {
		row.idx <- which(!rowAnyNAs(est.m))
		est_len.v <- sqrt(colSums(est.m[row.idx, ] ^ 2))
		true_len.v <- sqrt(colSums(true.m[row.idx, ] ^ 2))
		return(cor(est_len.v, true_len.v))
	}
	get_cos <- function(est.m, true.m) {
		row.idx <- which(!rowAnyNAs(est.m))
		return(sapply(seq_len(ncol(est.m)), function(j) lsa::cosine(est.m[row.idx, j], true.m[row.idx, j])))
	}
	cos_steady.v <- get_cos(assay(steady.o, "velocity"), true_v.m)
	len_diff_steady.v <- get_len_diff(assay(steady.o, "velocity"), true_v.m)
	len_pcc.v <- get_len_pcc(assay(steady.o, "velocity"), true_v.m)
	cos_dynamical.v <- get_cos(assay(dynamical.o, "velocity"), true_v.m)
	len_diff_dynamical.v <- get_len_diff(assay(dynamical.o, "velocity"), true_v.m)
	len_pcc_dynamical.v <- get_len_pcc(assay(dynamical.o, "velocity"), true_v.m)
	
	### use true KNN
	assay(sce.o, "X") <- true_s.m
	adata <- zellkonverter::SCE2AnnData(sce.o)
	
	### learn KNN graph using true S
	scv$pp$pca(adata, n_comps = ncomp)
	scv$pp$neighbors(adata)
	scv$pp$moments(adata)
	
	### switch back 
	sce.o <- zellkonverter::AnnData2SCE(adata)
	reducedDim(sce.o, "true_pca") <- reducedDim(sce.o, "X_pca")
	assay(sce.o, "X") <- assay(sce.o, "spliced")
	adata <- zellkonverter::SCE2AnnData(sce.o)
	
	
	### get steady
	scv$tl$velocity(adata, mode = "deterministic")
	steady.o <- zellkonverter::AnnData2SCE(adata)
	
	### get dynamical
	scv$tl$recover_dynamics(adata)
	scv$tl$velocity(adata, mode = "dynamical")
	dynamical.o <- zellkonverter::AnnData2SCE(adata)
	
	### 
	pcc_steady_trueknn.v <- sapply(seq_len(nrow(dynamical.o)), function(i) cor(assay(steady.o, "velocity")[i,], true_v.m[i,]))
	rmse_steady_trueknn.v <- sapply(seq_len(nrow(dynamical.o)), function(i) nrmse(true_v.m[i,], assay(steady.o, "velocity")[i,]) )
	pcc_dynamical_trueknn.v <- sapply(seq_len(nrow(dynamical.o)), function(i) cor(assay(dynamical.o, "velocity")[i,], true_v.m[i,]))
	rmse_dynamical_trueknn.v <- sapply(seq_len(nrow(dynamical.o)), function(i) nrmse(true_v.m[i,], assay(dynamical.o, "velocity")[i,]) )
	dynamical_pct_trueknn.v <- sum(!is.na(rmse_dynamical_trueknn.v)) / nrow(dynamical.o)
	
	cos_steady_trueknn.v <- get_cos(assay(steady.o, "velocity"), true_v.m)
	len_diff_steady_trueknn.v <- get_len_diff(assay(steady.o, "velocity"), true_v.m)
	len_pcc_trueknn.v <- get_len_pcc(assay(steady.o, "velocity"), true_v.m)
	cos_dynamical_trueknn.v <- get_cos(assay(dynamical.o, "velocity"), true_v.m)
	len_diff_dynamical_trueknn.v <- get_len_diff(assay(dynamical.o, "velocity"), true_v.m)
	len_pcc_dynamical_trueknn.v <- get_len_pcc(assay(dynamical.o, "velocity"), true_v.m)
	
	pcc.df <- rbind(data.frame(steady = pcc_steady.v, dynamical = pcc_dynamical.v, noise = noise_level, type = "learned"),
									data.frame(steady = pcc_steady_trueknn.v, dynamical = pcc_dynamical_trueknn.v, noise = noise_level, type = "true"))
	rmse.df <- rbind(data.frame(steady = rmse_steady.v, dynamical = rmse_dynamical.v, noise = noise_level, type = "learned"),
									data.frame(steady = rmse_steady_trueknn.v, dynamical = rmse_dynamical_trueknn.v, noise = noise_level, type = "true"))
	dynamical_pct.df <- data.frame(pct = c(dynamical_pct.v, dynamical_pct_trueknn.v), noise = noise_level, type = c("learned", "true"))
	cos.df <- rbind(data.frame(steady = cos_steady.v, dynamical = cos_dynamical.v, noise = noise_level, type = "learned"),
									data.frame(steady = cos_steady_trueknn.v, dynamical = cos_dynamical_trueknn.v, noise = noise_level, type = "true"))
	len_diff.df <- rbind(data.frame(steady = len_diff_steady.v, dynamical = len_diff_dynamical.v, noise = noise_level, type = "learned"),
									data.frame(steady = len_diff_steady_trueknn.v, dynamical = len_diff_dynamical_trueknn.v, noise = noise_level, type = "true"))
	len_pcc.df <- data.frame(pcc = c(len_pcc.v, len_pcc_trueknn.v,
																	 len_pcc_dynamical.v, len_pcc_dynamical_trueknn.v),
													 noise = noise_level, type = c("learned", "true", "learned", "true"),
													 name = c("steady", "steady", "dynamical", "dynamical"))
	return(list(pcc = pcc.df, rmse = rmse.df, dynamical_pct = dynamical_pct.df, cos = cos.df, len_diff = len_diff.df, len_pcc = len_pcc.df))
}

### get all R2
tmp.l <- lapply(seq_len(10), function(noise_level) {
	out.l <- make_sim_eva(noise_level = noise_level, n_obs = 500L, n_vars = 10L)
	return(out.l)
})

dynamical_pct.df <- do.call(rbind, lapply(tmp.l, "[[", "dynamical_pct"))
dynamical_pct.df$noise <- as.numeric(factor(as.numeric(dynamical_pct.df$noise)))


pcc.df <- do.call(rbind, lapply(tmp.l, "[[", "pcc"))
tmp.df <- pcc.df %>% pivot_longer(
	cols = 1:2,
	values_to = "pcc"
)
tmp.df$noise <- factor(as.numeric(tmp.df$noise))
tmp.df$name <- factor(tmp.df$name, levels = c("steady", "dynamical"))
ymin.v <- min(tmp.df$pcc, na.rm = TRUE)
ymax.v <- max(tmp.df$pcc, na.rm = TRUE)
yrange.v <- ymax.v - ymin.v

pcc.p <- ggplot(tmp.df, aes(x = noise, y = pcc, color = name,  dodge = name, fill = type)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.7) +
	scale_color_manual(values = c("#377EB8", "#4DAF4A"), name = "Model", limits = c("steady", "dynamical"), labels = c("Steady", "Dynamical")) +
	scale_fill_manual(values = c("white", "grey80"), name = "k-NN", labels = c("Learned", "True")) + 
	geom_segment(data = dynamical_pct.df %>% dplyr::filter(type == "learned"), linetype = "solid", size = 0.5, alpha = 0.7,
							 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange.v + ymin.v, yend = pct * yrange.v + ymin.v), inherit.aes = FALSE, color = "#4DAF4A") +
	geom_segment(data = dynamical_pct.df %>% dplyr::filter(type == "true"), linetype = "solid", size = 0.5, alpha = 0.9,
							 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange.v + ymin.v, yend = pct * yrange.v + ymin.v), inherit.aes = FALSE, color = "gray") +
	scale_y_continuous(name = "PCC", 
										 sec.axis = sec_axis(~ (. - ymin.v) / yrange.v , name = "Percentage of recovered genes\nby dynamical model", 
										 										labels = function(b) { paste0(round(b * 100, 0), "%")})) +
	labs( x = "Noise level", 
				title = Tex(r'(PCC between estimated $v_g$ and true $v_g$ of each gene)') ) +
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0))
																 
				

	

rmse.df <- do.call(rbind, lapply(tmp.l, "[[", "rmse"))
tmp.df <- rmse.df %>% pivot_longer(
	cols = 1:2,
	values_to = "rmse"
)
tmp.df$noise <- factor(as.numeric(tmp.df$noise))
tmp.df$name <- factor(tmp.df$name, levels = c("steady", "dynamical"))
ymin2.v <- min(tmp.df$rmse, na.rm = TRUE)
ymax2.v <- max(tmp.df$rmse, na.rm = TRUE)
yrange2.v <- ymax2.v - ymin2.v

rmse.p <- ggplot(tmp.df, aes(x = noise, y = rmse, color = name,  dodge = name, fill = type)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.7) +
	scale_color_manual(values = c("#377EB8", "#4DAF4A"), name = "Type", limits = c( "steady", "dynamical"), labels = c("Steady v", "Dynamical V")) +
	scale_fill_manual(values = c("white", "grey80"), name = "kNN", labels = c("Learned", "True")) + 
	geom_segment(data = dynamical_pct.df %>% dplyr::filter(type == "learned"), linetype = "solid", size = 0.5, alpha = 0.7,
							 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange2.v + ymin2.v, yend = pct * yrange2.v + ymin2.v), inherit.aes = FALSE, color = "#4DAF4A") +
	geom_segment(data = dynamical_pct.df %>% dplyr::filter(type == "true"), linetype = "solid", size = 0.5, alpha = 0.9,
							 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange2.v + ymin2.v, yend = pct * yrange2.v + ymin2.v), inherit.aes = FALSE, color = "gray") +
	scale_y_continuous(name = "NRMSE", 
										 sec.axis = sec_axis(~ (. - ymin2.v) / yrange2.v , name = "Percentage of recovered genes\nby dynamical model", 
										 										labels = function(b) { paste0(round(b * 100, 0), "%")})) +
	labs( x = "Noise level", 
				title = Tex(r'(NRMSE between estimated $v_g$ and true $v_g$ of each gene)')) +
	theme(legend.position = "none")




# 
# tmp2.l <- lapply(seq_len(10), function(noise_level) {
# 	out.l <- make_sim_eva(noise_level = noise_level, n_obs = 800L, n_vars = 300L)
# 	return(out.l)
# })
# 
# qsave(tmp2.l, file = here::here("data/tmp2.l.qs"))

# tmp2.l <- qread(here::here("data/tmp2.l.qs"))
# 
# dynamical_pct300.df <- do.call(rbind, lapply(tmp2.l, "[[", "dynamical_pct"))
# dynamical_pct300.df$noise <- as.numeric(factor(as.numeric(dynamical_pct300.df$noise)))
# 
# 
# pcc.df <- do.call(rbind, lapply(tmp2.l, "[[", "pcc"))[, -1]
# tmp.df <- pcc.df %>% pivot_longer(
# 	cols = 1:2,
# 	values_to = "pcc"
# )
# tmp.df$noise <- factor(as.numeric(tmp.df$noise))
# tmp.df$name <- factor(tmp.df$name, levels = c("steady", "dynamical"))
# ymin3.v <- min(tmp.df$pcc, na.rm = TRUE)
# ymax3.v <- max(tmp.df$pcc, na.rm = TRUE)
# yrange3.v <- ymax3.v - ymin3.v
# 
# pcc300.p <- ggplot(tmp.df, aes(x = noise, y = pcc, color = name,  dodge = name, fill = type)) +
# 	geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.7) +
# 	scale_color_manual(values = c("#377EB8", "#4DAF4A"), name = "Type", limits = c("steady", "dynamical"), labels = c("Steady v", "Dynamical V")) +
# 	scale_fill_manual(values = c("white", "grey80"), name = "kNN", labels = c("Learned", "True")) + 
# 	geom_segment(data = dynamical_pct300.df %>% dplyr::filter(type == "learned"), linetype = "solid", size = 0.5, alpha = 0.7,
# 							 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange3.v + ymin3.v, yend = pct * yrange3.v + ymin3.v), inherit.aes = FALSE, color = "#4DAF4A") +
# 	geom_segment(data = dynamical_pct300.df %>% dplyr::filter(type == "true"), linetype = "solid", size = 0.5, alpha = 0.9,
# 							 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange3.v + ymin3.v, yend = pct * yrange3.v + ymin3.v), inherit.aes = FALSE, color = "gray") +
# 	scale_y_continuous(name = "PCC", 
# 										 sec.axis = sec_axis(~ (. - ymin3.v) / yrange3.v , name = "Percentage of recovered genes\nby dynamical model", 
# 										 										labels = function(b) { paste0(round(b * 100, 0), "%")})) +
# 	labs( x = "Noise level", 
# 				title = "300 genes and 800 cells") +
# 	theme(legend.position = "none")
# 
# 
# 
# 
# 
# rmse.df <- do.call(rbind, lapply(tmp2.l, "[[", "rmse"))[, -1]
# tmp.df <- rmse.df %>% pivot_longer(
# 	cols = 1:2,
# 	values_to = "rmse"
# )
# tmp.df$noise <- factor(as.numeric(tmp.df$noise))
# tmp.df$name <- factor(tmp.df$name, levels = c("steady", "dynamical"))
# ymin4.v <- min(tmp.df$rmse, na.rm = TRUE)
# ymax4.v <- max(tmp.df$rmse, na.rm = TRUE)
# yrange4.v <- ymax4.v - ymin4.v
# 
# rmse300.p <- ggplot(tmp.df, aes(x = noise, y = rmse, color = name,  dodge = name, fill = type)) +
# 	geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.7) +
# 	scale_color_manual(values = c("#377EB8", "#4DAF4A"), name = "Type", limits = c("steady", "dynamical"), labels = c("Steady v", "Dynamical V")) +
# 	scale_fill_manual(values = c("white", "grey80"), name = "KNN", labels = c("Learned", "True")) + 
# 	geom_segment(data = dynamical_pct300.df %>% dplyr::filter(type == "learned"), linetype = "solid", size = 0.5, alpha = 0.7,
# 							 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange4.v + ymin4.v, yend = pct * yrange4.v + ymin4.v), inherit.aes = FALSE, color = "#4DAF4A") +
# 	geom_segment(data = dynamical_pct300.df %>% dplyr::filter(type == "true"), linetype = "solid", size = 0.5, alpha = 0.9,
# 							 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange4.v + ymin4.v, yend = pct * yrange4.v + ymin4.v), inherit.aes = FALSE, color = "gray") +
# 	scale_y_continuous(name = "NRMSE", 
# 										 sec.axis = sec_axis(~ (. - ymin4.v) / yrange4.v , name = "Percentage of recovered genes\nby dynamical model", 
# 										 										labels = function(b) { paste0(round(b * 100, 0), "%")})) +
# 	labs( x = "Noise level", 
# 				title = "300 genes and 800 cells") +
# 	theme(legend.position = "none")


### cos sim for each cell
cos.df <- do.call(rbind, lapply(tmp.l, "[[", "cos"))
tmp.df <- cos.df %>% pivot_longer(
	cols = 1:2,
	values_to = "cos"
)
tmp.df$noise <- factor(as.numeric(tmp.df$noise))
tmp.df$name <- factor(tmp.df$name, levels = c("steady", "dynamical"))


cos.p <- ggplot(tmp.df, aes(x = noise, y = cos, color = name,  dodge = name, fill = type)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.7) +
	scale_color_manual(values = c("#377EB8", "#4DAF4A"), name = "Type", limits = c("steady", "dynamical"), labels = c("Steady v", "Dynamical V")) +
	scale_fill_manual(values = c("white", "grey80"), name = "kNN", labels = c("Learned", "True")) + 
	# geom_segment(data = dynamical_pct.df %>% dplyr::filter(type == "learned"), linetype = "solid", size = 0.5, alpha = 0.7,
	# 						 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange.v + ymin.v, yend = pct * yrange.v + ymin.v), inherit.aes = FALSE, color = "#4DAF4A") +
	# geom_segment(data = dynamical_pct.df %>% dplyr::filter(type == "true"), linetype = "solid", size = 0.5, alpha = 0.9,
	# 						 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange.v + ymin.v, yend = pct * yrange.v + ymin.v), inherit.aes = FALSE, color = "gray") +
	# scale_y_continuous(name = "PCC", 
	# 									 sec.axis = sec_axis(~ (. - ymin.v) / yrange.v , name = "Percentage of recovered genes\nby dynamical model", 
	# 									 										labels = function(b) { paste0(round(b * 100, 0), "%")})) +
	labs( x = "Noise level", y = "Cosine similarity",
				title = Tex(r'(cos.sim between estimated $v_c$ and true $v_c$ of each cell)')) +
	theme(legend.position = "none")

### abs len diff for each cell
len_diff.df <- do.call(rbind, lapply(tmp.l, "[[", "len_diff"))
tmp.df <- len_diff.df %>% pivot_longer(
	cols = 1:2,
	values_to = "len_diff"
)
tmp.df$noise <- factor(as.numeric(tmp.df$noise))
tmp.df$name <- factor(tmp.df$name, levels = c("steady", "dynamical"))


len_diff.p <- ggplot(tmp.df, aes(x = noise, y = len_diff, color = name,  dodge = name, fill = type)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.7) +
	scale_color_manual(values = c("#377EB8", "#4DAF4A"), name = "Model", limits = c("steady", "dynamical"), labels = c("Steady V", "Dynamical V")) +
	scale_fill_manual(values = c("white", "grey80"), name = "k-NN", labels = c("Learned", "True")) + 
	# geom_segment(data = dynamical_pct.df %>% dplyr::filter(type == "learned"), linetype = "solid", size = 0.5, alpha = 0.7,
	# 						 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange.v + ymin.v, yend = pct * yrange.v + ymin.v), inherit.aes = FALSE, color = "#4DAF4A") +
	# geom_segment(data = dynamical_pct.df %>% dplyr::filter(type == "true"), linetype = "solid", size = 0.5, alpha = 0.9,
	# 						 aes(x = noise - 0.3, xend = noise + 0.3, y = pct * yrange.v + ymin.v, yend = pct * yrange.v + ymin.v), inherit.aes = FALSE, color = "gray") +
	# scale_y_continuous(name = "PCC", 
	# 									 sec.axis = sec_axis(~ (. - ymin.v) / yrange.v , name = "Percentage of recovered genes\nby dynamical model", 
	# 									 										labels = function(b) { paste0(round(b * 100, 0), "%")})) +
	labs( x = "Noise level", y = "Absolute length difference",
				title ="Absolute difference between estimated speed and true speed in high.dim for each cell") +
	theme(legend.position = "none")


### example showing 4 scatter plots of vector length
len_pcc.df <- do.call(rbind, lapply(tmp.l, "[[", "len_pcc"))
tmp.df <- len_pcc.df
tmp.df$noise <- factor(as.numeric(tmp.df$noise))
tmp.df$name <- factor(tmp.df$name, levels = c("steady", "dynamical"))


len_pcc.p <- ggplot(tmp.df, aes(x = noise, y = pcc, color = name,  shape = type)) +
	geom_point(size = 1.5, alpha = 0.6) +
	scale_color_manual(values = c("#377EB8", "#4DAF4A"), name = "Model", limits = c("steady", "dynamical"), labels = c("Steady", "Dynamical")) +
	scale_shape_discrete(name = "k-NN") +
	labs( x = "Noise level", y = "PCC",
				title = "PCC between estimated speed and true speed in high.dim for each cell") +
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0))




mp <- plot_grid(pcc.p, rmse.p, align = "hv", axis = "tblr",
														nrow = 2, ncol = 1, label_size = 10, labels = "auto")


save_plot(here::here("figs", "sfigs", "sfig.sim_pcc_rmse.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 3, device = cairo_pdf)


mp <- plot_grid(len_diff.p + theme(legend.position = c(0, 1),
																	 legend.justification = c(0, 1)),
								len_pcc.p, cos.p, align = "hv", axis = "tblr",
								nrow = 3, ncol = 1, label_size = 10, labels = "auto")


save_plot(here::here("figs", "sfigs", "sfig.sim_len_cos.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 3, ncol = 3, device = cairo_pdf)



