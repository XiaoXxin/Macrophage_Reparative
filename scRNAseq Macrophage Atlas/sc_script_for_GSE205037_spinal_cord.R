# global options
rm(list = ls())

# library
library(Seurat) # v4.3.0
library(scuttle) # v1.12.0
library(SummarizedExperiment) # v1.32.0
library(clusterProfiler) # v4.8.3
library(limma) # v3.58.1
library(org.Mm.eg.db) # v3.18.0
library(ALRA) # 0.0.0.9000
library(GSVA) # v1.44.2
library(tidyverse) # v2.0.0

# load data
d3 <- Read10X(data.dir = "./rawData/3d")
d7 <- Read10X(data.dir = "./rawData/7d")

# create Seurat object
seurat_d3 <- CreateSeuratObject(counts = d3, 
                                min.cells = 3, 
                                min.features = 200, 
                                names.field = 1, 
                                names.delim = "-", 
                                project = "day3")

seurat_d7 <- CreateSeuratObject(counts = d7, 
                                min.cells = 3, 
                                min.features = 200, 
                                names.field = 1, 
                                names.delim = "-", 
                                project = "day7")

seuratProject <- merge(seurat_d3, y = seurat_d7, add.cell.ids = c("d3", "d7"), project = "seuratProject")

# filter data
VlnPlot(object = seuratProject, features = c("nCount_RNA", "nFeature_RNA"))
seuratProject <- subset(x = seuratProject, subset = nCount_RNA < 40000 & nFeature_RNA < 6000)

# logNormalize
seuratProject <- NormalizeData(object = seuratProject, normalization.method = "LogNormalize", scale.factor = 10000)
# find variable features
seuratProject <- FindVariableFeatures(object = seuratProject, selection.method = "vst", nfeatures = 2000)
# scale Data
seuratProject <- ScaleData(object = seuratProject, features = rownames(seuratProject))
# PCA
seuratProject <- RunPCA(object = seuratProject, features = VariableFeatures(object = seuratProject))
# find clusters
seuratProject <- FindNeighbors(seuratProject, dims = 1:20)
seuratProject <- FindClusters(seuratProject, resolution = 0.5)
# UMAP
seuratProject <- RunUMAP(seuratProject, dims = 1:20)
DimPlot(seuratProject, reduction = "umap", label = T)

# save
saveRDS(seuratProject, file = "./RData/seuratProject.rds")
#seuratProject <- readRDS("./RData/seuratProject.rds")

# SingleR
data <- GetAssayData(seuratProject, layer="data")
clusters <- seuratProject@meta.data$seurat_clusters
refdata <- celldex::HumanPrimaryCellAtlasData()
cell_singleR <- SingleR(test = data, ref = refdata, labels = refdata$label.main,
                        clusters = clusters, 
                        assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cell_singleR ), 
                      celltype=cell_singleR$labels)

# find biomarkers
# Macro hi
VlnPlot(seuratProject, features = c("Ms4a6c", "Adgre1", "Ms4a7"))
VlnPlot(seuratProject, features = c("Cd68", "Cd14", "Csf1r", "Fcgr1"))
# Macro low
VlnPlot(seuratProject, features = c("P2ry12", "Tmem119"))

# extract clusters
mac = seuratProject[,seuratProject@meta.data$seurat_clusters %in% c(0,5,6,7,9,10,12)]
DimPlot(mac, reduction = "umap", label = T, split.by = "orig.ident")
dim(mac)
# 19823 macrophages


# re-analysis the selected clusters
mac <- NormalizeData(object = mac, normalization.method = "LogNormalize", scale.factor = 10000)
mac <- FindVariableFeatures(object = mac, selection.method = "vst", nfeatures = 2000)
mac <- ScaleData(object = mac, features = rownames(mac))
mac <- RunPCA(object = mac, features = VariableFeatures(object = mac))
mac <- FindNeighbors(mac, dims = 1:20)
mac <- FindClusters(mac, resolution = 0.5)
mac <- RunUMAP(mac, dims = 1:10)
DimPlot(mac, reduction = "umap", label = T, label.size = 5, label.box = FALSE, split.by = "orig.ident")


# save
mac <- JoinLayers(mac)
saveRDS(mac, file = "./RData/result_mac_spinal_cord.rds")
# mac <- readRDS("./RData/result_mac_spinal_cord.rds")

# imputation of missing values using ALRA
exprSet <- GetAssayData(mac, layer = "data") %>% as.matrix %>% t
exprSet_alra <- alra(exprSet)
exprSet_alra <- exprSet_alra[[3]]
rownames(exprSet_alra) <- rownames(exprSet)

# ssGSEA
buildGseaList <- function(species, selected) {
  require(clusterProfiler)
  require(tidyverse)
  
  keggList <- download_KEGG(species = species)
  
  keggRef <- keggList$KEGGPATHID2NAME %>% filter(to %in% selected)
  
  keggPath <- keggList$KEGGPATHID2EXTID %>% filter(from %in% keggRef$from)
  
  keggPath$from <- keggRef$to[match(keggPath$from, keggRef$from)]
  
  keggPath <- split(keggPath$to, keggPath$from)
  
  keggPath
  
}

geneList <- buildGseaList(species = "mmu", 
                          selected = c("Fatty acid biosynthesis - Mus musculus (house mouse)",
                                       "Fatty acid degradation - Mus musculus (house mouse)",
                                       "Glycolysis / Gluconeogenesis - Mus musculus (house mouse)",
                                       "Citrate cycle (TCA cycle) - Mus musculus (house mouse)",
                                       "Pentose phosphate pathway - Mus musculus (house mouse)",
                                       "Oxidative phosphorylation - Mus musculus (house mouse)"))

names(geneList) <- c("TCA", "FAS", "FAO", "Glycolysis", "OXPHOS", "PPP")

exprSet_alra <- t(exprSet_alra)
IDs <- bitr(rownames(exprSet_alra), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
entrezMatrix <- exprSet_alra[IDs$SYMBOL,]
rownames(entrezMatrix) <- IDs$ENTREZID


metaScore <- gsva(entrezMatrix, 
                  geneList,
                  method = "ssgsea",
                  kcdf = "Gaussian",
                  abs.ranking = T)  %>% t() %>% scale() %>% as.data.frame()


save(metaScore, file = "./RData/result_mac_metabolism.RData")
