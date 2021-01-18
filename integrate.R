library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
source("./preprocess.R")

dpi = 300

theme_set(theme_cowplot())

# Load and preprocess each dataset
Jak2VF_rep1.obj <- preprocess_data(datadir = "./data/Jak2VF_rep1", minCell = 10, minGene = 200, maxGene = 5000, maxUMI = 30000, pctMito = 0.1, id = "Jak2VF_rep1")
Jak2VF_GasDKO.obj <- preprocess_data(datadir = "./data/Jak2VF_GasDKO", minCell = 10, minGene = 200, maxGene = 4000, maxUMI = 25000, pctMito = 0.1, id = "Jak2VF_GasDKO")
Control.obj <- preprocess_data(datadir = "./data/Control", minCell = 10, minGene = 200, maxGene = 4000, maxUMI = 20000, pctMit = 0.1, id = "Control")
Jak2VF_rep2.obj <- preprocess_data(datadir = "./data/Jak2VF_rep2", minCell = 10, minGene = 200, maxGene = 4000, maxUMI = 20000, pctMito = 0.1, id = "Jak2VF_rep2")

# Set genotype
Jak2VF_rep1.obj[["genotype"]] <- "Jak2VF"
Jak2VF_GasDKO.obj[["genotype"]] <- "Jak2VFGasD-/-"
Control.obj[["genotype"]] <- "Control"
Jak2VF_rep2.obj[["genotype"]] <- "Jak2VF"

obj.list <- list(Jak2VF_rep1.obj, Jak2VF_GasDKO.obj, Control.obj, Jak2VF_rep2.obj)

# SCTransform for each sample
for (i in seq_along(obj.list))
{
	obj.list[[i]] <- SCTransform(obj.list[[i]], vars.to.regress = "percent.mt", variable.features.n = 500)
}

cell.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 1500)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = cell.features)

obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = cell.features, dims = 1:30)
obj.integrated <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", dims = 1:30)

# PCA
obj.integrated <- RunPCA(obj.integrated)

# UMAP
obj.integrated <- RunUMAP(obj.integrated, dims = 1:20)

# Cluster
obj.integrated <- FindNeighbors(obj.integrated, reduction = "pca", k.param = 20, dims = 1:20)
obj.integrated <- FindClusters(obj.integrated, resolution = 0.5)

png("UMAP_integrated.png", width = 12*dpi, height = 4*dpi, res = dpi)
DimPlot(obj.integrated, reduction = "umap", split.by = "genotype", label = TRUE, label.size = 4, ncol = 3)
dev.off()

# Differential expression analysis
DefaultAssay(obj.integrated) <- "RNA"

obj.integrated <- NormalizeData(obj.integrated, normalization.method = "LogNormalize", scale.factor = 10000)

DE.cluster <- FindAllMarkers(obj.integrated, only.pos = TRUE, test = "MAST", min.pct = 0.25, logfc.threshold = 0.4, min.diff.pct = -Inf)

# Top DE genes
top5 <- DE.cluster %>% group_by(cluster) %>% top_n(5, avg_logFC)

# Scale gene expression for heatmap
obj.integrated <- ScaleData(obj.integrated, features = top5$gene, vars.to.regress = "percent.mt")

png("Heatmap_top5_gene.png", width = 18*dpi, height = 11*dpi, res = dpi)
DoHeatmap(obj.integrated, features = top5$gene, angle = 0, size = 3.5) + scale_fill_gradient2(low = "dodgerblue2", high = "red", mid = "white", midpoint = 0, na.value = "white") + theme(axis.text.y = element_text(size = 9))
dev.off()

## Cell cycle scoring analysis

# Load cell cycle genes that have been mapped to mouse
cc.genes.mouse <- readRDS("./cell_cycle_gene/cc.genes.mouse.rds")
s.genes <- cc.genes.mouse$s.genes
g2m.genes <- cc.genes.mouse$g2m.genes

# Cell cycle scoring
obj.integrated <- CellCycleScoring(obj.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Plot predicted phase
png("cell_phase.png", width = 12*dpi, height = 4*dpi, res = dpi)
DimPlot(obj.integrated, reduction = "umap", group.by = "Phase", split.by = "genotype", ncol = 3)
dev.off()








