# library
library(DESeq2)
library(clusterProfiler)
library(limma)
library(fgsea)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(sva)
library(org.Mm.eg.db)
library(GSVA)
library(GO.db)
library(ggplotify)
library(ComplexHeatmap)
library(circlize)
library(paletteer)
library(tidyverse)

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

#### Fig. XXX ####
rm(list = ls())
load("./RData/result_DESeq2.RData")

geneList <- dif$log2FoldChange;head(geneList)
StoE <- bitr(dif$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
names(geneList) <- as.character(dif$symbol);head(geneList)
geneList <- geneList[names(geneList) %in% StoE$SYMBOL]
names(geneList) <- as.character(StoE$ENTREZID[match(names(geneList), StoE$SYMBOL)]);head(geneList)
geneList <- sort(geneList,decreasing = T);head(geneList)
buildGseaList<- function(species) {
  require(clusterProfiler)
  require(tidyverse)
  R.utils::setOption("clusterProfiler.download.method",'auto')
  keggList <- download_KEGG(species = species, keggType = "KEGG", keyType = "kegg")
  keggRef <- keggList$KEGGPATHID2NAME
  keggPath <- keggList$KEGGPATHID2EXTID
  keggPath$from <- keggRef$to[match(keggPath$from, keggRef$from)]
  keggPath <- split(keggPath$to, keggPath$from)
  keggPath
}
keggList <- buildGseaList(species = "mmu")
names(keggList) <- gsub(" - Mus .*", "", names(keggList))

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

#### Fig. XXX ####
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
                                       "Oxidative phosphorylation - Mus musculus (house mouse)",
                                       "Glycolysis / Gluconeogenesis - Mus musculus (house mouse)",
                                       "Pentose phosphate pathway - Mus musculus (house mouse)"))

IDss <- bitr(rownames(exprSet_batch), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
entrezMatrix <- exprSet_batch[IDss$SYMBOL,]
rownames(entrezMatrix) <- IDss$ENTREZID

# ssGSEA
metaScore <- gsva(entrezMatrix, 
                  geneList,
                  method = "ssgsea",
                  kcdf = "Gaussian",
                  abs.ranking = T)  %>% t() %>% scale() %>% as.data.frame()

colnames(metaScore) <- gsub(" ", "_", colnames(metaScore))
colnames(metaScore) <- gsub("-", "_", colnames(metaScore))


metaScore$Group <- phenoData$group
metaScore <- gather(metaScore, metabolism, score, -Group)
metaScore$metabolism <- gsub("___Mus_musculus.*", "", metaScore$metabolism)

metaScore$metabolism <- gsub("Fatty_acid_biosynthesis", "FAS", metaScore$metabolism)
metaScore$metabolism <- gsub("Fatty_acid_degradation", "FAO", metaScore$metabolism)
metaScore$metabolism <- gsub("Glycolysis_/_Gluconeogenesis", "Glycolysis", metaScore$metabolism)
metaScore$metabolism <- gsub("Oxidative_phosphorylation", "OXPHOS", metaScore$metabolism)
metaScore$metabolism <- gsub("Pentose_phosphate_pathway", "PPP", metaScore$metabolism)

metaScore$metabolism <- factor(metaScore$metabolism, levels = c("Glycolysis", "OXPHOS", "PPP", "FAO", "FAS"))

# Get the name and the y position of each label
data <- metaScore %>% group_by(metabolism, Group) %>% summarize(mean = round(mean(score),3),
                                                                sd = sd(score),
                                                                se= (sd(score) / sqrt(length(score))))

data$hline <- -0.05
data$hline2 <- 1.2

label_data <- data.frame(name = levels(data$metabolism),
                         y = -0.2)

res <- list()
for (i in c("Glycolysis", "OXPHOS", "PPP", "FAO", "FAS")) {
  res[[i]] <- compare_means(
  score~Group,
  metaScore[metaScore$metabolism == i,],
  method = "t.test",
  paired = FALSE,
  symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")),
  p.adjust.method = "none")
}

label_data$sig <- lapply(res, function(x) x$p) %>% as.numeric() %>% round(5) %>% paste0("p=", .)

p.meta<- ggplot(data, aes(x = metabolism, y = mean, fill = Group))+
  geom_col(position="dodge")+
  scale_fill_manual(values = c("#A0A0A4", "#E64B35"), breaks=c("WT", "TA"), labels=c("+/+", "A/A"))+
  geom_errorbar(aes(y=hline, ymax=hline, ymin=hline))+
  geom_errorbar(aes(y=hline2, ymax=hline2, ymin=hline2), width = 0.5)+
  geom_errorbar(aes(y=mean, ymax=mean, ymin=mean-se), position = position_dodge(0.9), width = 0.5)+
  geom_linerange(aes(x=metabolism, ymin=mean+0.16, ymax=hline2), position = position_dodge(1))+
  xlab(NULL)+
  ylab(NULL)+
  geom_text(data = label_data, aes(x = name, y = y, label = name),hjust=c(0.5,0.5,0.4,0.5,0.45), angle = c(-36,-108,0,108,36), colour = "black", alpha=0.8, size=4.5, fontface="bold", inherit.aes = FALSE)+
  geom_text(data = label_data, aes(x = name, y = 1.28, label = sig), angle = c(-36,-108,0,108,36), colour = "black", alpha=0.8, size=5, fontface="bold", inherit.aes = FALSE)+
  geom_text(data = data, aes(x = metabolism, y = mean+0.08, label = mean), colour = "black", angle = c(-18,-54,-90,-126, -162, -198, 126, 90, 54, 18), position = position_dodge(1),alpha=0.8, size = 5)+
  ylim(-0.7,1.3)+
  theme_minimal() +
  theme(
    legend.text=element_text(size=12),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_polar()

pdf(file = "./plot/meta.pdf", width = 9, height = 7)
p.meta
dev.off()



