rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <-  as(readRDS(here::here('data/sourcedata/Mahdessian_v1.2/RNAData/a.rds')), "SingleCellExperiment")
colnames(sce.o) <- read_csv(here::here("data/sourcedata/Mahdessian_v1.2/RNAData/a.obs_names.csv"))$well_plate
rownames(sce.o) <- rowData(sce.o)$Gene

### keep spliced and unspliced assays
assays(sce.o) <- assays(sce.o)[c("spliced", "unspliced")]
assays(sce.o) <- lapply(assays(sce.o), function(x) as(x, "dgCMatrix"))

### filter genes with more than 20 counts
idx <- which((rowSums(assay(sce.o, "spliced")) > 20) & (rowSums(assay(sce.o, "unspliced")) > 20))
sce.o <- sce.o[idx, ]

### filter cells with more than 200 genes
idx <- which((colSums(assay(sce.o, "spliced") > 0) > 200) & (colSums(assay(sce.o, "unspliced") > 0) > 200))
sce.o <- sce.o[, idx]

### load in FUCCI data
new_fucci.df <- read_csv(here::here("data/sourcedata/Mahdessian_v1.2/fucci_coords.csv"))
idx <- match(colnames(sce.o), new_fucci.df$cell)
sum(is.na(idx))
new_fucci.df <- new_fucci.df[idx, ]

sce.o$Green530 <- new_fucci.df$green530_lognorm_rescale
sce.o$Red585 <- new_fucci.df$red585_lognorm_rescale
sce.o$fucci_time <- new_fucci.df$fucci_time_hrs
sce.o$fucci_time <- sce.o$fucci_time / max(sce.o$fucci_time, na.rm = TRUE)

### scale total counts
sce.o <- logNormCounts(sce.o, assay.type= "spliced")


### Do PCA on spliced counts using HVG
dec <- modelGeneVar(sce.o)
top.hvgs <- getTopHVGs(dec, n=2000)

set.seed(100)
sce.o <- runPCA(sce.o, subset_row = top.hvgs)
### do UMAP
library(umap)
umap.o <- umap::umap(reducedDim(sce.o, "PCA")[, 1:30], min_dist = 0.5)
umap.m <- umap.o$layout[, c(2, 1)]

reducedDim(sce.o, "UMAP") <- umap.m
plotUMAP(sce.o)

### Do PCA on GO cell cycle genes
library(org.Hs.eg.db)
go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
																						columns = "SYMBOL")[, "SYMBOL"]

sce.o <- runPCA(sce.o, subset_row = rownames(sce.o) %in% go.v, name = "PCA.cc")
sce.o <- runUMAP(sce.o, dimred="PCA.cc", name = "UMAP.cc")



### remove sample without FUCCI time
sce.o <- sce.o[, which(!is.na(sce.o$fucci_time))]


library(velociraptor)
steady.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA")
reducedDims(steady.o) <- reducedDims(sce.o)
rowData(steady.o) <- cbind(rowData(steady.o), rowData(sce.o))
colData(steady.o) <- cbind(colData(steady.o), colData(sce.o))

dynamical.o <- scvelo(sce.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA")
reducedDims(dynamical.o) <- reducedDims(sce.o)
rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(sce.o))
colData(dynamical.o) <- cbind(colData(dynamical.o), colData(sce.o))


qsave(sce.o, file = here::here("data/sourcedata/Mahdessian_v1.2/sce.o.qs"))
qsave(steady.o, file = here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))
qsave(dynamical.o, file = here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))
qsave(umap.o, file = here::here("data/sourcedata/Mahdessian_v1.2/umap.o.qs"))








plotUMAP(dynamical.o, colour_by="velocity_pseudotime")
plotUMAP(dynamical.o, colour_by="fucci_time")
plotUMAP(dynamical.o, colour_by="sizeFactor")

embedded <- embedVelocity(reducedDim(dynamical.o, "UMAP"), dynamical.o)
grid.df <- gridVectors(reducedDim(dynamical.o, "UMAP"), embedded)

plotUMAP(dynamical.o, colour_by="velocity_pseudotime") +
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, 
																				 xend=end.1, yend=end.2), arrow=arrow(length=unit(0.05, "inches")))

plotUMAP(dynamical.o, colour_by="fucci_time") +
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, 
																				 xend=end.1, yend=end.2), arrow=arrow(length=unit(0.05, "inches")))

plotVelocity(dynamical.o, "TOP2A", use.dimred = 2)
plot(dynamical.o$fucci_time, assay(dynamical.o, "velocity")["TOP2A", ])
plot(dynamical.o$fucci_time, assay(dynamical.o, "Ms")["TOP2A", ])
plot(dynamical.o$fucci_time, assay(dynamical.o, "fit_t")["TOP2A", ])

plot(assay(dynamical.o, "Ms")["TOP2A", ], assay(dynamical.o, "velocity")["TOP2A", ] )

plotVelocity(dynamical.o, "CDK1", use.dimred = 2)
plot(dynamical.o$fucci_time, assay(dynamical.o, "velocity")["CDK1", ])
plot(dynamical.o$fucci_time, assay(dynamical.o, "Ms")["CDK1", ])
plot(assay(dynamical.o, "Ms")["CDK1", ], assay(dynamical.o, "velocity")["CDK1", ] )

plot(assay(dynamical.o, "Ms")["CDK1", ], assay(dynamical.o, "Ms")["TOP2A", ])

plot(assay(dynamical.o, "Ms")["SMC2", ], assay(dynamical.o, "Ms")["TOP2A", ])
plotVelocity(dynamical.o, "SMC4", use.dimred = 2)
plotVelocity(dynamical.o, "SMC2", use.dimred = 2)
plot(dynamical.o$fucci_time, assay(dynamical.o, "Ms")["SMC2", ])
plot(dynamical.o$fucci_time, assay(dynamical.o, "velocity")["SMC2", ])
plot(dynamical.o$fucci_time, assay(dynamical.o, "Ms")["SMC4", ])
plot(dynamical.o$fucci_time, assay(dynamical.o, "velocity")["SMC4", ])
plot(assay(dynamical.o, "Ms")["SMC2", ], assay(dynamical.o, "Ms")["SMC4", ])











