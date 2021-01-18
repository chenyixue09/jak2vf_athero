library(Seurat)

preprocess_data <- function(datadir, minCell = 10, minGene = 200, maxGene, maxUMI, pctMito = 0.1, id)
{
	sample.data <- Read10X(data.dir = datadir)
	sample.obj <- CreateSeuratObject(counts = sample.data, min.cells = minCell, min.features = minGene, project = id)
	
	# QC mitochondrial genes
	sample.obj[["percent.mt"]] <- PercentageFeatureSet(object = sample.obj, pattern = "^mt-")
	
	# Filter by maximum #genes, maximum #UMIs, and mito%
	idx <- sample.obj[["nFeature_RNA"]] >= minGene & sample.obj[["nFeature_RNA"]] <= maxGene & sample.obj[["nCount_RNA"]] <= maxUMI & sample.obj[["percent.mt"]] <= 100*pctMito
	sample.obj = subset(sample.obj, cells = colnames(sample.obj)[idx])
	
	return(sample.obj)
}
