# global options
rm(list = ls())
options(stringsAsFactors = FALSE)

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
files <- list.files(path = "./rawData", full.names = T)

d4_1 <- Read10X_h5(files[1])
d4_2 <- Read10X_h5(files[2])
d4_3 <- Read10X_h5(files[3])
d4_4 <- Read10X_h5(files[4])
d4_5 <- Read10X_h5(files[5])
d4_6 <- Read10X_h5(files[6])
d4_7 <- Read10X_h5(files[7])
d4_8 <- Read10X_h5(files[8])

d7_1 <- Read10X_h5(files[9])
d7_2 <- Read10X_h5(files[10])
d7_3 <- Read10X_h5(files[11])
d7_4 <- Read10X_h5(files[12])
d7_5 <- Read10X_h5(files[13])
d7_6 <- Read10X_h5(files[14])
d7_7 <- Read10X_h5(files[15])
d7_8 <- Read10X_h5(files[16])

# create Seurat object
fileList <- ls(pattern = "d.*")
seuratList <- list(d4_1,d4_2,d4_3,d4_4,d4_5,d4_6,d4_7,d4_8,
                   d7_1,d7_2,d7_3,d7_4,d7_5,d7_6,d7_7,d7_8)


for (i in 1:16) {
  assign(paste0("seurat_",fileList[i]),
         CreateSeuratObject(counts = seuratList[[i]], min.cells = 3, min.features = 200, names.field = 1, names.delim = "-", project = fileList[i]))
}


# merge data
seuratProject <- list(seurat_d4_1,seurat_d4_3,seurat_d4_4,seurat_d4_5,seurat_d4_6,seurat_d4_7,seurat_d4_8,
                      seurat_d7_1,seurat_d7_2,seurat_d7_3,seurat_d7_4,seurat_d7_5,seurat_d7_6,seurat_d7_7,seurat_d7_8)


seuratProject <- merge(seuratProject[[1]], seuratProject[-1], project = "seuratProject")

# filter data
VlnPlot(object = seuratProject, features = c("nCount_RNA", "nFeature_RNA"))
seuratProject <- subset(x = seuratProject, subset = nCount_RNA < 20000 & nFeature_RNA < 4000)
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
seuratProject <- FindClusters(seuratProject, resolution = 0.3)
# UMAP
seuratProject <- RunUMAP(seuratProject, dims = 1:20)
DimPlot(seuratProject, reduction = "umap", label = T)

# save
seuratProject <- JoinLayers(seuratProject)
saveRDS(seuratProject, file = "./RData/result_seuratProject_muscle.rds")
# seuratProject <- readRDS("./RData/result_seuratProject_muscle.rds")

# SingleR
data <- GetAssayData(seuratProject, layer="data")
clusters <- seuratProject@meta.data$seurat_clusters
refdata <- celldex::HumanPrimaryCellAtlasData()
cell_singleR <- SingleR(test = data, ref = refdata, labels = refdata$label.main,
                        clusters = clusters, 
                        assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cell_singleR ), 
                      celltype=cell_singleR$labels)

# macrophage markers
VlnPlot(seuratProject, features = c("Fcgr1", "Itgam", "Cd68", "Mertk", "Ccr2", "Apoe", "Lyz1", "Ccl4"))

# extract clusters
mac <- seuratProject[,seuratProject@meta.data$seurat_clusters %in% c(0,3,6,8,11)]
DimPlot(mac, reduction = "umap", label = T)
dim(mac)
# 13315 macrophages

# re-analysis the selected clusters
mac <- NormalizeData(object = mac, normalization.method = "LogNormalize", scale.factor = 10000)
mac <- FindVariableFeatures(object = mac, selection.method = "vst", nfeatures = 2000)
mac <- ScaleData(object = mac, features = rownames(mac))
mac <- RunPCA(object = mac, features = VariableFeatures(object = mac))
mac <- FindNeighbors(mac, dims = 1:10)
mac <- FindClusters(mac, resolution = 0.5)
mac <- RunUMAP(mac, dims = 1:10)
DimPlot(mac, reduction = "umap", label = T, label.size = 5, label.box = FALSE, split.by = "orig.ident")

phenoData <- mac@meta.data
phenoData$orig.ident <- gsub("d4.*", "day4", phenoData$orig.ident)
phenoData$orig.ident <- gsub("d7.*", "day7", phenoData$orig.ident)
mac <- AddMetaData(mac, phenoData)

# save
saveRDS(mac, file = "./RData/result_mac_muscle.rds")
# mac <- readRDS("./RData/result_mac_muscle.rds")

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

