# library
library(DESeq2) # v1.40.1
library(clusterProfiler) # v4.8.3
library(limma) # v3.58.1
library(fgsea) # v1.22.0
library(ggpubr) # v0.6.0
library(FactoMineR) # v2.11
library(factoextra) # v1.0.7
library(sva) # v3.50.0
library(org.Mm.eg.db) # v3.18.0
library(GSVA) # v1.44.2
library(GO.db) # v3.18.0
library(ggplotify) # v0.1.2
library(ComplexHeatmap) # v2.18.0
library(circlize) # v0.4.16
library(paletteer) # v1.6.0
library(tidyverse) # v2.0.0

# global options
rm(list = ls())

# create dir
dir.create("./rawData")
dir.create("./RData")

# download GSE252899 rawdata from GEO to ./rawData dir

# load data
count <- read.table(file = "./rawData/GSE252899_raw_counts.txt.gz", sep = "\t", header = T, row.names = 1)

# remove duplicated genes
IDs <- bitr(rownames(count), fromType = "REFSEQ", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")

count <- count[rownames(count) %in% IDs$REFSEQ,]
IDs <- IDs[match(rownames(count), IDs$REFSEQ),]

counts <- avereps(count, ID = IDs$SYMBOL) %>% round()

# sample annotation informations
phenoData <- data.frame(sample = colnames(counts),
                        group = factor(c("TA", "TA", "TA", "WT", "WT", "WT"), levels = c("WT", "TA")))
phenoData$batch <- factor(c("B1", "B1", "B2", "B1", "B1", "B2"))

# build DESeqDataSet
dds <- DESeqDataSetFromMatrix(counts, 
                              colData = phenoData, 
                              design = ~group+batch)

# remove low expressed genes
dim(dds)
cps <- fpm(dds)
keep1 <- rowSums(cps >= 1) >= ceiling(ncol(dds)*3/4)
dds <- dds[keep1, ]
dim(dds)

# call DESeq function
dds <- DESeq(dds)

# get log2 scaled matrix
exprSet <- rlogTransformation(dds, blind = FALSE) %>% assay()

IDs <- IDs[match(rownames(exprSet), IDs$SYMBOL),]

res <- results(dds, contrast = c("group","TA","WT")) 
dif <- as.data.frame(res) %>% na.omit()
dif$symbol <- rownames(dif)

save(counts, exprSet, dds, phenoData,IDs,dif,res, file = "./RData/result_DESeq2.RData")

#### Fig. S5b ####
rm(list = ls())
load("./RData/result_DESeq2.RData")

geneList <- dif$log2FoldChange;head(geneList)
StoE <- bitr(dif$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
names(geneList) <- as.character(dif$symbol);head(geneList)
geneList <- geneList[names(geneList) %in% StoE$SYMBOL]
names(geneList) <- as.character(StoE$ENTREZID[match(names(geneList), StoE$SYMBOL)]);head(geneList)
geneList <- sort(geneList,decreasing = T);head(geneList)

goList <- toTable(org.Mm.egGO2ALLEGS)
goList <- goList[goList$Ontology == "BP",]
goAnn <- toTable(GOTERM)[,2:3]
goList$term <- goAnn$Term[match(goList$go_id, goAnn$go_id)]

goList <- split(goList$gene_id, goList$term)

res.fgsea.go <- fgsea(pathways = goList,
                   stats = geneList,
                   minSize = 5,
                   maxSize = 500)
res.fgsea.go <- res.fgsea.go[res.fgsea.go$pval <=0.05,]
res.fgsea.go$pathway <- gsub("interleukin-", "IL", res.fgsea.go$pathway)
res.fgsea.go$pathway <- gsub("T-helper ", "Th", res.fgsea.go$pathway)
names(goList) <- gsub("interleukin-", "IL", names(goList))
names(goList) <- gsub("T-helper ", "Th", names(goList))

pathways <- c("negative regulation of innate immune response", "positive regulation of type 2 immune response",
              "inflammatory response to wounding",
              "negative regulation of IL10 production", 
              "acute inflammatory response",
              "cellular response to IL1"
              )

p.gsea <- plotGseaTable(goList[pathways], geneList, res.fgsea.go, colwidths = c(3, 1.5, 0.6, 0.8, 0), gseaParam = 0.5, render = F)
p.gsea <- as.ggplot(p.gsea)

pdf(file = "./fgsea.pdf", width = 8, height = 3)
p.gsea
dev.off()

#### Fig. S5a ####
rm(list = ls())
load("./RData/result_DESeq2.RData")

dat <- as.data.frame(t(exprSet))
pca <- PCA(dat, graph = FALSE)

fviz_pca_ind(pca, 
             label = "all",
             axes = c(1,2),
             title = "BMDM",
             pointsize = 2.5,
             pointshape = 19,
             addEllipses = F,
             mean.point = FALSE,
             habillage = as.factor(phenoData$group))


mod <- model.matrix(~as.factor(group), data = phenoData)
exprSet_batch <- ComBat(dat = exprSet, batch = phenoData$batch, mod = mod)

dat <- as.data.frame(t(exprSet_batch))
pca <- PCA(dat, graph = FALSE)

p.pca <- fviz_pca_ind(pca, 
             label = "none",
             axes = c(1,2),
             title = "BMDM",
             pointsize = 7,
             pointshape = 19,
             addEllipses = F,
             mean.point = FALSE,
             habillage = as.factor(phenoData$group))+
  scale_color_manual(values = c("#A0A0A4", "#E64B35"), breaks=c("WT", "TA"), labels=c("+/+", "A/A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(-200, 200)+ylim(-200,200)

pdf(file = "./pca.pdf", width = 5, height = 5)
p.pca
dev.off()


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

IDss <- bitr(rownames(exprSet_batch), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
entrezMatrix <- exprSet_batch[IDss$SYMBOL,]
rownames(entrezMatrix) <- IDss$ENTREZID

# ssGSEA
metaScore_BMDM <- gsva(entrezMatrix, 
                  geneList,
                  method = "ssgsea",
                  kcdf = "Gaussian",
                  abs.ranking = T)  %>% t() %>% scale() %>% as.data.frame()

save(metaScore_BMDM, file = "./RData/result_BMDM_metabolism.RData")
