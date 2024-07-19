
# library
library(limma) # v3.58.1
library(FactoMineR) # v2.11
library(factoextra) # v1.0.7
library(clusterProfiler) # v4.8.3
library(org.Hs.eg.db) # v3.18.0
library(GSVA) # v1.44.2
library(sva) # v3.50.0
library(ggpubr) # v0.6.0
library(paletteer) # v1.6.0
library(scattermore) # v1.2
library(ggrepel) # v0.9.5
library(pheatmap) # v1.0.12
library(tidyverse) # v2.0.0

# global options
rm(list = ls())

# create dir
dir.create("./rawData/GSE162698")

# download GSE162698 rawdata from GEO to ./rawData/GSE162698 dir

# load data
tpm <- read.table(file = "./rawData/GSE162698/GSE162698_tpm_summary_invitro_Macs.txt", sep = "\t", header = T)

tpm <- avereps(tpm[,3:17],ID = tpm$symbol)

# phenoData
phenoData <- data.frame(sample = colnames(tpm), 
                        ID = str_split(colnames(tpm), "_", simplify = T)[,1], 
                        group = str_split(colnames(tpm), "_", simplify = T)[,2])

phenoData <- phenoData[phenoData$group %in% c("M0", "M2IL10"),]
phenoData <- arrange(phenoData, group)

tpm <- tpm[,phenoData$sample]

# filter low-expressed genes
dim(tpm)
keep <- rowSums(tpm[,1:ncol(tpm)] >=1 ) >= 3
tpm <- tpm[keep, ]
dim(tpm)

ltpm <- log2(tpm+1)
rownames(ltpm) <- gsub("[.].*", "", rownames(ltpm))

# batch effects
dat <- as.data.frame(t(ltpm))
pca <- PCA(dat, graph = FALSE)

fviz_pca_ind(pca, 
             label = "none", 
             title = " ",
             pointsize = 2.5,
             pointshape = 19,
             addEllipses = F, 
             mean.point = FALSE,
             habillage = as.factor(phenoData$ID))

mod = model.matrix(~as.factor(group), data=phenoData)
ltpm_batch <- ComBat(dat = ltpm, batch = phenoData$ID, mod = mod)

dat <- as.data.frame(t(ltpm_batch))
pca <- PCA(dat, graph = FALSE)

fviz_pca_ind(pca, 
             label = "none", 
             title = " ",
             pointsize = 2.5,
             pointshape = 19,
             addEllipses = F, 
             mean.point = FALSE,
             habillage = as.factor(phenoData$group))

# DEGs analysis
group <- factor(phenoData[,"group"])

design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))

contrast.matrix <- makeContrasts(contrasts = c("M2IL10-M0"), levels = design)

fit <- lmFit(ltpm_batch, design) %>% contrasts.fit(contrast.matrix) %>% eBayes()

dif <- topTable(fit, coef = c("M2IL10-M0"), n = nrow(fit))

dif$symbol <- rownames(dif)


#### Fig. S1E ####
test <- dif[,c(1, 4, 7)]
test$change = ifelse(test$P.Value < 0.05 & abs(test$logFC) >= log2(1.5), ifelse(test$logFC> log2(1.5) ,'Up in M2IL10','Up in M1'), 'Stable')
test[,2] <- -log10(test$P.Value)
colnames(test)[2] <- "v"

test$name <- test$symbol

test$name[!test$name%in% c("FASN", "SREBF1", "ACSL1", "ALOX5", "MBOAT7", "CYP2S1")] = ""

p.v <- ggplot(test, aes(x = logFC, y = v, color = change))+
  geom_scattermore(pointsize = 8, pixels = c(2048, 2048))+
  scale_color_manual(values = c("#A0A0A499","#5e90b899",  "#f1939c99"))+
  xlim(-6,6)+
  labs(x = "log2FC (M(IL10)/M0)", y ="-log10(p value)")+
  geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed")+
  geom_vline(aes(xintercept= log2(1.5)), linetype="dashed")+
  geom_vline(aes(xintercept= -log2(1.5)), linetype="dashed")+
  geom_label_repel(aes(label=name), 
                   size = 3,
                   fontface="bold", color="black", 
                   box.padding=unit(0.35, "lines"), 
                   point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50",
                   max.overlaps = Inf,
                   force_pull = 3)+
  ggtitle("Up in M0                                Up in M(IL10)")+
  theme_classic()+
  theme(panel.grid=element_blank(),
        legend.title = element_blank(),
        plot.title=element_text(hjust=0.5), 
        axis.title = element_text(size = 15),
        legend.position = "none")


pdf(file = "./human_v_IL10_M0.pdf", width = 6.5, height = 4)
p.v
dev.off()


#### Fig. S1D ####
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

geneList <- buildGseaList(species = "hsa", 
                          selected = c("Fatty acid biosynthesis",
                                       "Fatty acid degradation",
                                       "Glycolysis / Gluconeogenesis",
                                       "Citrate cycle (TCA cycle)",
                                       "Pentose phosphate pathway",
                                       "Oxidative phosphorylation"))

names(geneList) <- c("TCA", "FAS", "FAO", "Glycolysis", "OXPHOS", "PPP")

# rename gene names
IDs <- bitr(rownames(ltpm_batch), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrezMatrix <- ltpm_batch[IDs$SYMBOL,]
rownames(entrezMatrix) <- IDs$ENTREZID

# ssGSEA
metaScore <- gsva(entrezMatrix, 
                  geneList,
                  method = "ssgsea",
                  kcdf = "Gaussian",
                  abs.ranking = T)  %>% t() %>% scale() %>% as.data.frame()


annotation_col <- data.frame(row.names = phenoData$sample, Group = phenoData$group)
annotation_col$Group <- gsub("M2IL10", "M(IL10)", annotation_col$Group)
ann_colors = list(
  Group = c('M0' = "#cad0d3",
            "M(IL10)" = "#5c8987")
)


selected <- metaScore %>% t()


pdf(file = "./human_ssGSEA_heatmap.pdf", width = 5, height = 5)
pheatmap(selected, 
         fontsize_row = 9,
         show_colnames = F,
         scale = "none", 
         border_color = "black", 
         cluster_cols = F,
         color = colorRampPalette(c("#5e90b8", "white","#f1939c"))(50),
         cellwidth = 23, 
         cellheight = 12,
         annotation_col = annotation_col,
         annotation_colors = ann_colors)
dev.off()




