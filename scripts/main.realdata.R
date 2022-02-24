rm(list=ls())
source(here::here("scripts/utils.R"))

getR2 <- function(path, name, ntop = 300) {
	sce.o <- qread(path)
	rowData.df <- rowData(sce.o) %>% as.data.frame() %>%
		filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))
	sce.o <- sce.o[rownames(rowData.df)[seq_len(min(nrow(rowData.df), ntop))], ]
	r2_est.v <- sapply(seq_len(nrow(sce.o)), function(i) fit_t_loess(x = assay(sce.o, "fit_t")[i, ], y = assay(sce.o, "Ms")[i, ], span = 1)$rsquared)
	r2u_est.v <- sapply(seq_len(nrow(sce.o)), function(i) fit_t_loess(x = assay(sce.o, "fit_t")[i, ], y = assay(sce.o, "Mu")[i, ], span = 1)$rsquared)
	return(data.frame(r2s = r2_est.v, r2u = r2u_est.v, name))
}

# fucci.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))

# fucci.df <- getR2(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"), "fucci")
# pancreas.df <- getR2(here::here("data/Pancreas/dynamical.o.qs"), "pancreas")
# chromaffin.df <- getR2(here::here("data/chromaffin/dynamical.o.qs"), "chromaffin")
# bonemarrow.df <- getR2(here::here("data/bonemarrow/dynamical.o.qs"), "bonemarrow")
# dentategyrus_lamanno.df <- getR2(here::here("data/dentategyrus_lamanno/dynamical.o.qs"), "dentategyrus_lamanno")
# dentategyrus_hochgerner.df <- getR2(here::here("data/dentategyrus/dynamical.o.qs"), "dentategyrus_hochgerner")
# erythroid.df <- getR2(here::here("data/erythroid/dynamical.o.qs"), "erythroid")
# forebrain.df <- getR2(here::here("data/forebrain/dynamical.o.qs"), "forebrain")
# gastrulation_e75.df <- getR2(here::here("data/gastrulation_e75/dynamical.o.qs"), "gastrulation_e75")
# pbmc68k.df <- getR2(here::here("data/pbmc68k/dynamical.o.qs"), "pbmc68k")
# 
# tmp.df <- rbind(fucci.df, pancreas.df, chromaffin.df, bonemarrow.df, dentategyrus_lamanno.df,
# 								dentategyrus_hochgerner.df, erythroid.df, forebrain.df, gastrulation_e75.df, pbmc68k.df)
# qsave(tmp.df, here::here("data/realdata_r2.df.qs"))

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
	geom_boxplot(color = "#377EB8",outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.8) +
	scale_fill_manual(values = c("white", "grey80"), name = "S/U") + 
	labs( x = "", 
				y =  expression(italic(R)^2), 
				title = "Explained variations of real data") +
	theme(legend.position = "none", 
				axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 0.5, angle = 30),
				axis.title.x = element_blank())









