
# library
library(Seurat)
library(SCENIC)
library(doParallel)
library(SCopeLoomR)
library(AUCell)
library(GSVA)
library(GSEABase)
library(ggsc)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(tidyverse)



# create dir
dir.create("./RData")
dir.create("./RData/cis")

# load data (SKIN)
# download cis data ("mm9-500bp-upstream-7species.mc9nr.feather", "mm9-tss-centered-10kb-7species.mc9nr.feather") from https://resources.aertslab.org
# copy "result_mac_skin.rds" file from ~/scRNAseq Macrophage Atlas/GSE186986/RData/
rm(list = ls())

mac <- readRDS("./RData/result_mac_skin.rds")
exprMat  <-  as.matrix(mac@assays$RNA@data)
cellInfo <-  mac@meta.data[,c(1,2,3)]
colnames(cellInfo)=c('CellType', 'nUMI', 'nGene')

### Initialize settings
dbname <- list.files("../RData/cis")
dbname <- dbname[c(3,5)]
names(dbname) <- c("500bp","10kb")

data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9

scenicOptions <- initializeScenic(org="mgi", dbDir="./RData/cis", nCores=20, dbs = defaultDbNames[["mgi"]]) 

saveRDS(cellInfo, file = "int/cellInfo.Rds")
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

rm(list = ls()) 
scenicOptions=readRDS(file="int/scenicOptions.Rds")

### Exploring output 
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))

###########################################################
#### Fig. 3A & Fig. S4B ####
rm(list = ls()) 

scenicOptions=readRDS(file="int/scenicOptions.Rds")
mac <- readRDS("./RData/result_mac_skin.rds")

scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)


a <- embeddings$`SCENIC: tSNE_AUC_50pcs_30perpl`
colnames(a) <- c("x","y")
b <- regulonsAUC@assays@data@listData$AUC
b <- t(b) %>% as.data.frame()

data <- cbind(a,b)
data$group <- gsub("_.*", "", rownames(data))

p2 <- ggplot(data, aes(x = group, y = `Pparg_extended (71g)`, fill = group))+
  geom_violin(scale = "width")+
  xlab(NULL)+ylab("PPARg activity")+
  stat_compare_means(method = "t.test", aes(label = paste0("p ", after_stat(p.format))))+
  scale_x_discrete(labels = c("w1" = "Day1", "w5" = "Day5"))+
  scale_fill_manual(values = c("black","#b84b48"))+
  theme_classic()+
  theme(legend.position = "none") 

# for revising
ggplot(data, aes(x = group, y = `Stat3_extended (379g)`, fill = group))+
  geom_violin(scale = "width")+
  geom_boxplot(color = "white", width = 0.1, outlier.shape = NA)+
  xlab(NULL)+ylab("STAT3 activity")+
  stat_compare_means(method = "t.test", aes(label = paste0("p ", after_stat(p.format))))+
  scale_x_discrete(labels = c("w1" = "Day1", "w5" = "Day5"))+
  scale_fill_manual(values = c("black","#b84b48"))+
  theme_classic()+
  theme(legend.position = "none")


dat <- data[,c("x", "y", "Pparg_extended (71g)", "Srebf2_extended (2418g)", "Cebpd_extended (417g)", "Stat3_extended (379g)")]
colnames(dat) <- c("SCENIC_x", "SCENIC_y", "PPARG", "SREBF2", "CEBPD", "STAT3")

mac <- AddMetaData(mac, metadata = dat)

p1 <- sc_feature(mac, features = "PPARG",pointsize = 15, pixels = c(2048, 2048))+
  ggtitle("PPARg activity")+
  theme(panel.grid=element_blank(),
        legend.title = element_blank(),
        plot.title=element_text(hjust=0.5), 
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15))+
  scale_color_gradientn(colors = c("#cad0d322","#b84b48", "#b84b48"),breaks = c(0,0.07), labels = c("Low", "High"))

p1.2 <- sc_feature(mac, features = "SREBF2",pointsize = 15, pixels = c(2048, 2048))+
  ggtitle("SREBF2 activity")+
  theme(panel.grid=element_blank(),
        legend.title = element_blank(),
        plot.title=element_text(hjust=0.5), 
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15))+
  scale_color_gradientn(colors = c("#cad0d322","#b84b48", "#b84b48"),breaks = c(0.04,0.3), labels = c("Low", "High"))

p1.3 <- sc_feature(mac, features = "CEBPD",pointsize = 15, pixels = c(2048, 2048))+
  ggtitle("CEBPD activity")+
  theme(panel.grid=element_blank(),
        legend.title = element_blank(),
        plot.title=element_text(hjust=0.5), 
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15))+
  scale_color_gradientn(colors = c("#cad0d322","#b84b48", "#b84b48"),breaks = c(0,0.09), labels = c("Low", "High"))


pdf(file = "./plot/scRNA_PPARG_acti.pdf", width = 5, height = 4)
p1
dev.off()

pdf(file = "./plot/scRNA_PPARG_acti_2.pdf", width = 3, height = 4)
p2
dev.off()


pdf(file = "./plot/scRNA_SREBF2_acti.pdf", width = 5, height = 4)
p1.2
dev.off()
pdf(file = "./plot/scRNA_CEBPD_acti.pdf", width = 5, height = 4)
p1.3
dev.off()



#### Fig. S4A ####
rm(list = ls())
scenicOptions=readRDS(file="int/scenicOptions.Rds")
cellInfo <- readRDS("int/cellInfo.Rds")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellInfo[colnames(regulonAUC), "CellType"])
rss <- rss[grep("extended",rownames(rss)),]
rssPlot <- plotRSS(rss, zThreshold = 0.6)
rssDif <- rssPlot$df

auc <- getAUC(regulonAUC)
dim(auc)

# parameters for column annotations
annotation_col <- data.frame(row.names = colnames(auc), CellType = gsub("_.*", "", colnames(auc)))
annotation_col$CellType <- gsub("w1", "Phase I", annotation_col$CellType)
annotation_col$CellType <- gsub("w5", "Phase II", annotation_col$CellType)

colMeta <- annotation_col

colAnn <- HeatmapAnnotation(
  df = colMeta,
  col = list(
    CellType = c('Phase I' = "#A0A0A4",
                 'Phase II' = "#5e90b8")
  ),
  which = 'column',
  
  # parameters for labels
  annotation_label = c("Source"),
  annotation_name_gp = gpar(fontsize = 10),
  
  # parameters for legend
  annotation_legend_param = list(
    CellType = list(
      title = 'Source',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 10)))
)

# parameters for row annotations
rowOrder <- gsub("_extended.*", "", rssPlot$rowOrder)

annotation_row <- data.frame(row.names = rowOrder, 
                             CTSR = c(rep("Phase II", 64), rep("Phase I", 61)))


rowMeta <- annotation_row

rowAnn <- HeatmapAnnotation(
  df = rowMeta,
  col = list(
    CTSR = c('Phase I' = "#dc943b", 'Phase II' = "#67b4b0")
  ),
  which = 'row',
  
  # parameters for labels
  annotation_label = c("CTSR"),
  annotation_name_gp = gpar(fontsize = 10),
  
  # parameters for legend
  show_legend = c(TRUE),
  annotation_legend_param = list(
    CTSR = list(
      title = 'CTSR',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 10))
  )
  )

# parameters for cell colors
bk <- seq(-0.8, 0.8, length.out = 100)
color <- colorRampPalette(c("#5e90b8", "white", "#f1939c"))(100)

# matrix
selected <- auc[rssPlot$rowOrder,]
rownames(selected) <- gsub("_extended.*", "", rownames(selected))
selected <- t(selected) %>% scale() %>% t()

matrix <- selected


matrix <- by(t(matrix), colMeta$CellType, colMeans) %>% Reduce(cbind,.) %>% magrittr::set_colnames(c("Phase I", "Phase II"))


# labels
genes <- c("Pparg", "Srebf2", "Cebpb", "Cebpd")

genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = match(genes,rownames(matrix)),
    labels = genes,
    labels_gp = gpar(fontsize = 12, fontface = 'bold'),
    padding = 0.75),
  width = unit(2.0, 'cm') +
    
    max_text_width(
      genes,
      gp = gpar(fontsize = 12,  fontface = 'bold')))


p.h <- Heatmap(
  matrix = matrix,
  
  col = colorRamp2(bk, color),
  
  border_gp = gpar(col = "black"),
  
  # legend parameters
  name = "AUC",
  heatmap_legend_param = list(
    title_position = 'topleft',
    legend_direction = 'vertical',
    legend_width = unit(30, 'cm'),
    title_gp = gpar(fontsize = 12, fontface = 'bold'),
    labels_gp = gpar(fontsize = 10)), 
  
  # row (gene) parameters
  cluster_rows = T,
  show_row_names = F,
  row_names_gp = gpar(fontsize = 10),
  cluster_row_slices = F,
  
  #row_title = c(rep("",8)),
  
  # column (sample) parameters
  cluster_columns = F,
  show_column_names = T,
  
  
  # specify top and bottom annotations
  #top_annotation = colAnn,
  left_annotation = rowAnn,
  
  width = unit(1.5, "inches")
)

pdf(file = "./plot/AUCheatmap.pdf", width = 5, height = 7)
draw(p.h + genelabels,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right',
     row_sub_title_side = 'right')

dev.off()

#### Fig. 3B ####
rm(list = ls())

mac <- readRDS("./RData/result_mac_skin.rds")

countexp <- mac@assays$RNA@counts %>% as.matrix()
countexp <- log2(countexp+1) %>% t()

# imputation
result.completed <- alra(countexp)
countexp2 <- result.completed[[3]]
rownames(countexp2) <- rownames(countexp)
countexp2 <- t(countexp2)

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
                                       "Fatty acid degradation - Mus musculus (house mouse)"))

IDss <- bitr(rownames(countexp2), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
countexp3 <- countexp2[IDss$SYMBOL,]
rownames(countexp3) <- IDss$ENTREZID


gsva_es <- gsva(as.matrix(countexp3), 
                geneList, 
                method = c("ssgsea"), 
                kcdf = c("Gaussian"))

signature_exp <- gsva_es %>% t() %>% magrittr::set_colnames(c("FAS", "FAO")) %>% data.frame()


scenicOptions=readRDS(file="int/scenicOptions.Rds")
cellInfo <- readRDS("int/cellInfo.Rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellInfo[colnames(regulonAUC), "CellType"])
rss <- rss[grep("extended",rownames(rss)),]
rssPlot <- plotRSS(rss, zThreshold = 0.6)
auc <- getAUC(regulonAUC)

selected <- auc[rssPlot$rowOrder, rownames(signature_exp)]
rownames(selected) <- gsub("_extended.*", "", rownames(selected))

selected <- selected[c("Pparg", "Srebf2", "Cebpd"),] %>% t() %>% as.data.frame()
selected[selected == 0] = NA
selected <- na.omit(selected)

rownames(signature_exp) <- gsub("[.]", "-", rownames(signature_exp))

dat <- data.frame(FAS = signature_exp[rownames(selected),]$FAS,
                  FAO = signature_exp[rownames(selected),]$FAO,
                  Pparg = selected[,"Pparg"],
                  Srebf2 = selected[,"Srebf2"],
                  Cebpd = selected[,"Cebpd"])


dat <- gather(dat, "TF", "act", -FAS, -FAO)

p <- ggplot(dat, aes(x = FAS, y = act))+
  geom_hex(bins = 15)+
  scale_fill_gradient(low = "#6163A366", high = "#6163A3")+
  theme_bw()+
  geom_smooth(method = 'lm', se = T, color = "#DB802D")+
  stat_cor(method = "spearman")+
  labs(x = "FAS", y = "TF activity")+
  theme(legend.position = "none",
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5))+
  facet_wrap(~TF, scales = "free", nrow = 1)

pdf(file = "./plot/FAS_cor.pdf", width = 7, height = 2.5)
p
dev.off()


