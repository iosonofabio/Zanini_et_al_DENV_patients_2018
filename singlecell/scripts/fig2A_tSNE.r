cpm = read.table(gzfile('data/fabio.tsv.gz'))
dim(cpm)  # 21863 8796

pbmc <- CreateSeuratObject(raw.data = cpm, min.cells = 10, min.genes = 250, project = "comb")  #  all cells are retained by this filter

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ]) / Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 1.575)  # this results in 500 variable genes

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)

ndims <- 5

pbmc <- RunTSNE(object = pbmc, dims.use = 1:ndims, do.fast = TRUE, perplexity=10)
TSNEPlot(object = pbmc)
