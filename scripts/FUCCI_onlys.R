
rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))

embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)

p1 <- plotEmbFucci(sce.o, "UMAP", title = "Vector fields", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)


### scale total counts
sce.o <- logNormCounts(sce.o, assay.type= "spliced")

### Do PCA on spliced counts using HVG
set.seed(200)
dec <- modelGeneVar(sce.o)
top.hvgs <- getTopHVGs(dec, n= 2000)

# library(org.Hs.eg.db)
# go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049",
# 															columns = "SYMBOL")[, "SYMBOL"]
# cc.df <- rowData(sce.o)[rownames(sce.o) %in% go.v, ] %>% as.data.frame() %>%
# 	filter(`velocity_genes`)

# tmp <- sin(  pi * rowRanks(assay(sce.o, "Ms")) / ncol(sce.o) ) 
# dimnames(tmp) <- dimnames(sce.o)
# assay(sce.o, "velocity") <- tmp
nan.idx <- which(!(rownames(sce.o) %in% top.hvgs)) #   which(!rowData(sce.o)$velocity_genes) # 
assay(sce.o, "velocity")[nan.idx, ] <- matrix(rep(NaN, length(nan.idx) * ncol(sce.o)), ncol = ncol(sce.o))


scv <- reticulate::import("scvelo")
adata <- zellkonverter::SCE2AnnData(sce.o)
do.call(scv$tl$velocity_graph, c(list(data=adata)))
new.o <- zellkonverter::AnnData2SCE(adata)


embedded2 <- embedVelocity(reducedDim(new.o, "UMAP"), new.o)
grid2.df <- gridVectors(reducedDim(new.o, "UMAP"), embedded2)

p2 <- plotEmbFucci(new.o, "UMAP", title = "Only using Ms", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid2.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)


mp <- wrap_plots(list(p1, p2), ncol = 2)

save_plot(here::here("figs", "FUCCI", str_c("only_using_Ms.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 1.5, ncol = 2 * 1.5, device = cairo_pdf)




