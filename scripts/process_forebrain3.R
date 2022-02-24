rm(list=ls())
source(here::here("scripts/utils.R"))

library(velociraptor)

sce.o <- qread(here::here("data/forebrain/sce.o.qs"))
steady.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA", scvelo.params = list(moments = list(n_neighbors = 550L)))
reducedDims(steady.o) <- reducedDims(sce.o)
rowData(steady.o) <- cbind(rowData(steady.o), rowData(sce.o))
colData(steady.o) <- cbind(colData(steady.o), colData(sce.o))
metadata(steady.o) <- c(metadata(steady.o), metadata(sce.o))

dynamical.o <- scvelo(sce.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA", scvelo.params = list(moments = list(n_neighbors = 550L)))
reducedDims(dynamical.o) <- reducedDims(sce.o)
rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(sce.o))
colData(dynamical.o) <- cbind(colData(dynamical.o), colData(sce.o))
metadata(dynamical.o) <- c(metadata(dynamical.o), metadata(sce.o))


qsave(steady.o, file = here::here("data/forebrain/steady2.o.qs"))
qsave(dynamical.o, file = here::here("data/forebrain/dynamical2.o.qs"))


sce.o <- qread(here::here("data/forebrain/sce.o.qs"))
steady.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA", scvelo.params = list(moments = list(n_neighbors = 550L, method = 'sklearn')))
reducedDims(steady.o) <- reducedDims(sce.o)
rowData(steady.o) <- cbind(rowData(steady.o), rowData(sce.o))
colData(steady.o) <- cbind(colData(steady.o), colData(sce.o))
metadata(steady.o) <- c(metadata(steady.o), metadata(sce.o))
qsave(steady.o, file = here::here("data/forebrain/steady3.o.qs"))


