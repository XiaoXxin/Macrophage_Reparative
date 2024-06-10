# global options
rm(list = ls())

# library
library(Seurat)
library(scuttle)
library(SummarizedExperiment)
library(clusterProfiler)
library(limma)
library(org.Mm.eg.db)
library(ALRA)
library(GSVA)
library(tidyverse)


# load data
meta <- read.table(file = "./rawData/Single_cell_Meta_data.txt", sep = ",", header = T)
meta$Cell_barcode2 <- gsub("_.*", "", meta$Cell_barcode)
meta_48 <- meta[meta$time_point==48,]
meta_96 <- meta[meta$time_point==96,]
meta_168 <- meta[meta$time_point==168,]
# remove duplicated cells
meta_48 <- meta_48[!meta_48$Cell_barcode2 %in% meta_48$Cell_barcode2[duplicated(meta_48$Cell_barcode2)],]
meta_96 <- meta_96[!meta_96$Cell_barcode2 %in% meta_96$Cell_barcode2[duplicated(meta_96$Cell_barcode2)],]
meta_168 <- meta_168[!meta_168$Cell_barcode2 %in% meta_168$Cell_barcode2[duplicated(meta_168$Cell_barcode2)],]

# read 10x
d48a <- Read10X_h5("./rawData/48a_raw_feature_bc_matrix.h5")
d48a <- d48a[[1]]
colnames(d48a) <- gsub("-.*", "", colnames(d48a))
d48a <- d48a[,colnames(d48a) %in% meta_48$Cell_barcode2]
colnames(d48a) <- paste0(colnames(d48a), "_48h")

d96w <- Read10X_h5("./rawData/96_1w_raw_feature_bc_matrix.h5")
d96w <- d96w[[1]]
colnames(d96w) <- gsub("-.*", "", colnames(d96w))
d96 <- d96w[,colnames(d96w) %in% meta_96$Cell_barcode2]
d1w <- d96w[,colnames(d96w) %in% meta_168$Cell_barcode2]
colnames(d96) <- paste0(colnames(d96), "_96h")
colnames(d1w) <- paste0(colnames(d1w), "_168h")

meta_48 <- meta_48 %>% 
  mutate(Cell_barcode2 = paste0(Cell_barcode2, "_48h")) %>%
  remove_rownames() %>%
  column_to_rownames("Cell_barcode2")
meta_48 <- meta_48[colnames(d48a),]
meta_96 <- meta_96 %>% 
  mutate(Cell_barcode2 = paste0(Cell_barcode2, "_96h")) %>%
  remove_rownames() %>%
  column_to_rownames("Cell_barcode2")
meta_96 <- meta_96[colnames(d96),]
meta_168 <- meta_168 %>% 
  mutate(Cell_barcode2 = paste0(Cell_barcode2, "_168h")) %>%
  remove_rownames() %>%
  column_to_rownames("Cell_barcode2")
meta_168 <- meta_168[colnames(d1w),]

# create Seurat object
seurat_48 <- CreateSeuratObject(counts = d48a, 
                                min.cells = 3, 
                                min.features = 200, 
                                names.field = 1, 
                                names.delim = "_", 
                                meta.data = meta_48,
                                project = "h48")
seurat_96 <- CreateSeuratObject(counts = d96, 
                                min.cells = 3, 
                                min.features = 200, 
                                names.field = 1, 
                                names.delim = "_", 
                                meta.data = meta_96,
                                project = "h96")
seurat_168 <- CreateSeuratObject(counts = d1w, 
                                min.cells = 3, 
                                min.features = 200, 
                                names.field = 1, 
                                names.delim = "_", 
                                meta.data = meta_168,
                                project = "h168")

seuratProject <- list(seurat_48, seurat_96, seurat_168)
seuratProject <- merge(seuratProject[[1]], y = seuratProject[-1], project = "seuratProject")

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
seuratProject <- FindClusters(seuratProject, resolution = 0.3)
# UMAP
seuratProject <- RunUMAP(seuratProject, dims = 1:20)
DimPlot(seuratProject, reduction = "umap", label = T)

# save
seuratProject <- JoinLayers(seuratProject)
saveRDS(seuratProject, file = "./RData/result_seuratProject_liver.rds")
# seuratProject <- readRDS("./RData/result_seuratProject_liver.rds")

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
DimPlot(seuratProject, reduction = "umap", label = T)
DimPlot(seuratProject, group.by = "cell_type")
VlnPlot(seuratProject, features = c("Cd68", "Cd14", "Csf1r", "Fcgr1"))

# extract clusters
mac <- subset(x = seuratProject, cell_type == "Macrophages")
dim(mac)
# 1262 macrophages

# re-analysis the selected clusters
mac <- NormalizeData(object = mac, normalization.method = "LogNormalize", scale.factor = 10000)
mac <- FindVariableFeatures(object = mac, selection.method = "vst", nfeatures = 2000)
mac <- ScaleData(object = mac, features = rownames(mac))
mac <- RunPCA(object = mac, features = VariableFeatures(object = mac))
mac <- FindNeighbors(mac, dims = 1:10)
mac <- FindClusters(mac, resolution = 0.5)
mac <- RunUMAP(mac, dims = 1:10)
DimPlot(mac, reduction = "umap", label = T, label.size = 5, label.box = FALSE)

# save
saveRDS(mac, file = "./RData/result_mac_liver.rds")
# mac <- readRDS("./RData/result_mac_liver.rds")

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
