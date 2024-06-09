# global options
rm(list = ls())

# library
library(Seurat)
library(SingleR)
library(ALRA)
library(GSVA)
library(monocle)
library(tidyverse)

# load raw data
w1 <- Read10X(data.dir = "./rawData/w1")
w5 <- Read10X(data.dir = "./rawData/w5")

# merge data
seurat_w1 <- CreateSeuratObject(counts = w1, 
                                    min.cells = 3, 
                                    min.features = 200, 
                                    names.field = 1, 
                                    names.delim = "-", 
                                    project = "day1")
seurat_w5 <- CreateSeuratObject(counts = w5, 
                                min.cells = 3, 
                                min.features = 200, 
                                names.field = 1, 
                                names.delim = "-", 
                                project = "day5")
seuratProject <- merge(seurat_w1, y = seurat_w5, add.cell.ids = c("w1", "w5"), project = "seuratProject")

# filter data
VlnPlot(object = seuratProject, features = c("nCount_RNA", "nFeature_RNA"))
seuratProject <- subset(x = seuratProject, subset = nCount_RNA < 20000 & nFeature_RNA < 6000)
# logNormalize
seuratProject <- NormalizeData(object = seuratProject, normalization.method = "LogNormalize", scale.factor = 10000)
# find variable features
seuratProject <- FindVariableFeatures(object = seuratProject, selection.method = "vst", nfeatures = 2000)
# scale Data
seuratProject <- ScaleData(object = seuratProject, features = rownames(seuratProject))
# PCA
seuratProject <- RunPCA(object = seuratProject, features = VariableFeatures(object = seuratProject))
# find clusters
seuratProject <- FindNeighbors(seuratProject, dims = 1:10)
seuratProject <- FindClusters(seuratProject, resolution = 0.5)
# UMAP
seuratProject <- RunUMAP(seuratProject, dims = 1:10)
DimPlot(seuratProject, reduction = "umap", label = T)

# save
seuratProject <- JoinLayers(seuratProject)
saveRDS(seuratProject, file = "./RData/result_seuratProject_skin.rds")
# seuratProject <- readRDS("./RData/result_seuratProject_skin.rds")

# SingleR
data <- GetAssayData(seuratProject, layer="data")
clusters <- seuratProject@meta.data$seurat_clusters
refdata <- celldex::HumanPrimaryCellAtlasData()
cell_singleR <- SingleR(test = data, ref = refdata, labels = refdata$label.main,
                        clusters = clusters, 
                        assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cell_singleR ), 
                      celltype=cell_singleR$labels)

#Macrophage
VlnPlot(seuratProject, features = c("Csf1r", "Fcgr1", "Adgre1", "Apoe"))

# extract clusters
mac = seuratProject[,seuratProject@meta.data$seurat_clusters %in% c(1,4)]
DimPlot(mac, reduction = "umap", label = T)
dim(mac)
# 2953 macrophages

# re-analysis the selected clusters
mac <- NormalizeData(object = mac, normalization.method = "LogNormalize", scale.factor = 10000)
mac <- FindVariableFeatures(object = mac, selection.method = "vst", nfeatures = 2000)
mac <- ScaleData(object = mac, features = rownames(mac))
mac <- RunPCA(object = mac, features = VariableFeatures(object = mac))
mac <- FindNeighbors(mac, dims = 1:20)
mac <- FindClusters(mac, resolution = 0.5)
mac <- RunUMAP(mac, dims = 1:10)
DimPlot(mac, reduction = "umap", label = T, label.size = 5, label.box = FALSE)

# save
saveRDS(mac, file = "./RData/result_mac_skin.rds")
# mac <- readRDS("./RData/result_mac_skin.rds")

# imputation of missing values using ALRA
exprSet <- GetAssayData(mac, slot = "data") %>% as.matrix %>% t
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

# transform Seurat object to monocle
rm(list = ls())
mac <- readRDS("./RData/result_mac_skin.rds")

data <- GetAssayData(mac, slot = "counts") %>% as.matrix

pData <- mac@meta.data %>% 
  subset(select = c(orig.ident)) %>%
  add_column("cluster" = mtype_GSE186986$subtype)
pd <- new("AnnotatedDataFrame", data = pData)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

# create cell data set
cds <- newCellDataSet(data, 
                      phenoData = pd, 
                      featureData = fd,
                      expressionFamily = negbinomial.size(),
                      lowerDetectionLimit = 0)

# estimate size factor and dispersion
cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds)

# filter low-quality samples
cds <- detectGenes(cds, min_expr = 0.1)

# choose genes that define a sample's progress
markers <- VariableFeatures(mac)
#markers <- FindMarkers(mac, ident.1 = "day1", ident.2 = "day5", group.by = "orig.ident") %>% rownames()
cds <- setOrderingFilter(cds, markers)
dim(cds)

# reduce data dimensionality
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

# order samples
cds <- orderCells(cds)

save(cds, file = "./RData/result_monocle_skin.RData")
# load("./RData/result_monocle.RData")

