
# library
library(limma) # v3.58.1
library(clusterProfiler) # v4.8.3
library(org.Mm.eg.db) # v3.18.0
library(ggpubr) # v0.6.0
library(GSVA) # v1.44.2
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # v3.10.0
library(pheatmap) # v1.0.12
library(tidyverse) # v2.0.0

# global options
rm(list = ls())

# create dir
dir.create("./rawData")
dir.create("./rawData/GSE131364")
dir.create("./RData")

# download GSE131364 rawdata from GEO to ./rawData/GSE131364 dir

# load data
file <- list.files(path = "./rawData/GSE131364", pattern = "\\.txt$", full.names = T)
name <- str_split(file,"_",simplify = T)[,5] %>% gsub(".txt", "", .)
counts <- list()
for (i in 1:length(name)) {
  count <- read.table(file = file[i], sep = "\t", header = T)
  count <- count[,c("tracking_id", "condition", "replicate","raw_frags")]
  count <- count[grep("WT", count$condition),]
  count$sample <- paste(count$condition,name[i],count$replicate, sep = "_")
  count <- count[,c("tracking_id", "raw_frags", "sample")]
  count <- spread(count, sample, raw_frags)
  counts[[i]] <- count
}

counts <- Reduce(function(x,y) merge(x,y, by = "tracking_id"), counts)

counts <- column_to_rownames(counts, "tracking_id")
counts <- as.matrix(counts)
counts <- counts[,!colnames(counts) %in% paste("WT", name[2:5], rep(c(0,1,2),each = 4), sep = "_")]


rownames(counts) <- gsub("[.].*", "", rownames(counts))
IDs <- bitr(rownames(counts), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
dim(counts)
counts <- counts[rownames(counts) %in% IDs$ENSEMBL,]
dim(counts)
IDs <- IDs[match(rownames(counts),IDs$ENSEMBL), ]
counts <- avereps(counts, ID = IDs$SYMBOL)
dim(counts)
counts <- round(counts)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_len <- lapply(exons_gene, function(x) {sum(width(reduce(x)))})
gene_length <- sapply(exons_gene_len, function(x) {x})
id_length <- as.data.frame(gene_length)
id_length$ENTREZID <- row.names(id_length)
IDs <- bitr(id_length$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL",OrgDb = "org.Mm.eg.db")
id_length <- merge(IDs, id_length, by = "ENTREZID")


counts <- counts[rownames(counts) %in% id_length$SYMBOL,]

id_length <- id_length[match(rownames(counts), id_length$SYMBOL),] %>% as.data.frame()

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}


tpms <- apply(counts,2,function(x) countToTpm(x, effLen = id_length$gene_length))
tpms[1:3,]


pdata <- data.frame(sample = colnames(counts), 
                    group = str_split(colnames(counts), "_",simplify = T)[,2])
pdata$group <- gsub("Control", "M0", pdata$group)
pdata$group <- gsub("IFG", "M(IFNg)", pdata$group)
pdata$group <- gsub("IL10", "M(IL10)", pdata$group)
pdata$group <- gsub("IL4", "M(IL4)", pdata$group)
pdata$group <- gsub("LPS", "M(LPS)", pdata$group)

phenoData <- pdata

# filter low-expressed genes
tpm <- tpms
dim(tpm)
keep <- rowSums(tpm[,1:ncol(tpm)] >=1 ) >= 3
tpm <- tpm[keep, ]
dim(tpm)

ltpm <- log2(tpm+1)

save(ltpm, tpm, tpms, phenoData,IDs, file = "./RData/GSE131364.RData")
# load("./RData/GSE131364.RData")

#### Fig. 1G ####
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

# rename gene names
phenoData <- phenoData[!grepl("M[(]IFNg[)]", phenoData$group),]
ltpm <- ltpm[,phenoData$sample]
IDs <- bitr(rownames(ltpm), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
entrezMatrix <- ltpm[IDs$SYMBOL,]
rownames(entrezMatrix) <- IDs$ENTREZID

metaScore <- gsva(entrezMatrix, 
                  geneList,
                  method = "ssgsea",
                  kcdf = "Gaussian",
                  abs.ranking = T)  %>% t() %>% scale() %>% as.data.frame()

annotation_col <- data.frame(row.names = phenoData$sample, CellType = phenoData$group)
ann_colors = list(
  CellType = c('M0' = "#E64B35",
               "M(LPS)" = "#00A087",
               "M(IL4)" = "#3C5488",
               "M(IL10)" = "#F39B7F")
)


selected <- metaScore %>% t()


pdf(file = "./plot/mouse_ssGSEA_heatmap.pdf", width = 9, height = 4)
pheatmap(selected, 
         fontsize_row = 9,
         show_colnames = F,
         scale = "row", 
         cluster_cols = F,
         border_color = "black", 
         color = colorRampPalette(c("#5e90b8", "white","#f1939c"))(50),
         cellwidth = 23, 
         cellheight = 12,
         gaps_col = c(3,6,9),
         annotation_names_col =F,
         annotation_col = annotation_col,
         annotation_colors = ann_colors)
dev.off()


