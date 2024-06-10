
# library
library(ggpubr)
library(patchwork)
library(ggsci)
library(ggsc)
library(gg.gap)
library(scattermore)
library(ggplotify)
library(ggthemes)
library(Seurat)
library(monocle)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)
library(tidyverse)


##### Fig. 1b #####
rm(list = ls())
mac <- readRDS("./GSE186986/RData/result_mac_skin.rds")

p.dim <- sc_feature(mac, features = "orig.ident",pointsize = 20, pixels = c(2048, 2048))+
  ggtitle("SKIN")+
  theme(panel.grid=element_blank(),
        legend.title = element_blank(),
        plot.title=element_text(hjust=0.5), 
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15))+
  scale_color_manual(values = c("#F39B7F55", "#8491B455"), labels = c("Phase I","Phase II"))

pdf(file = "./scRNA_dim.pdf", width = 5, height = 4)
p.dim
dev.off()

##### Fig. 1c #####
rm(list = ls())
mac <- readRDS("./GSE186986/RData/result_mac_skin.rds")

# find DEGs
markers <- FindMarkers(mac, ident.1 = "day1", ident.2 = "day5", group.by = "orig.ident")

fc <- markers$avg_log2FC
names(fc) <- rownames(markers)

# select genes
geneList_day1 = rownames(markers)[markers$avg_log2FC > 0]
geneList_day5 = rownames(markers)[markers$avg_log2FC < 0]

geneList <- list("Phase I" = geneList_day1, "Phase II" = geneList_day5)

# GO enrichment
x <- compareCluster(geneList,
                    fun="enrichGO", 
                    OrgDb = "org.Mm.eg.db",
                    pvalueCutoff=0.05,
                    ont = 'BP',
                    pAdjustMethod = 'BH',
                    keyType = 'SYMBOL')
save(x, file = "./GSE186986/RData/result_compareCluster.RData")
# load("./GSE186986/RData/result_compareCluster.RData")

categorys <- c("unsaturated fatty acid biosynthetic process",
               "fatty acid biosynthetic process",
               "lipid catabolic process",
               "glucose metabolic process",
               "pyruvate metabolic process",
               "ATP metabolic process"
)

pdf(file = "./scRNA_go.pdf", width = 5, height = 5)
cnetplot(x,
         label_format = 9999, 
         showCategory = categorys, 
         layout = "kk",
         node_label = "category",
         pie.params = list(pie = "Count", legend_n = 2),
         cex.params = list(category_label = 1.5))+scale_fill_manual(values =c("#4DBBD5FF", "#F39B7FFF","#4DBBD5FF", "#F39B7FFF"))
dev.off()

##### Fig. S1a-b & Fig. 1e  #####
rm(list = ls())

readscdata_FF <- function(GSE, tissue, label, groups){
  mac <- readRDS(paste0("./", GSE, "/RData/result_mac_", tissue, ".rds"))
  load(paste0("./", GSE, "/RData/result_mac_metabolism.RData"))
  
  phenoData <- mac@meta.data
  phenoData <- cbind(phenoData, metaScore)
  
  phenoData$group <- ifelse(phenoData[,groups[1]] == groups[2], "Phase I", "Phase II")
  
  dat <- phenoData[,c("group", "FAS", "FAO")]
  dat$tissue <- label
  
  dat
}

dat_6035873 <- readscdata_FF(GSE = "6035873", tissue = "liver", label = "LIVER",groups = c("orig.ident", "h48"))
dat_GSE141259 <- readscdata_FF(GSE = "GSE141259", tissue = "lung", label = "LUNG",groups = c("grouping", "d10", "d28"))
dat_GSE152501 <- readscdata_FF(GSE = "GSE152501", tissue = "air", label = "AWI", groups = c("orig.ident", "day1", "day7"))
dat_GSE163465 <- readscdata_FF(GSE = "GSE163465", tissue = "heart", label = "MCIN", groups = c("orig.ident", "day3"))
dat_GSE180420 <- readscdata_FF(GSE = "GSE180420", tissue = "kidney", label = "KDI", groups = c("pheno", "IRI_short_3", "IRI_short_14"))
dat_GSE186986 <- readscdata_FF(GSE = "GSE186986", tissue = "skin", label = "SKIN", groups = c("orig.ident", "day1", "day5"))
dat_GSE200843 <- readscdata_FF(GSE = "GSE200843", tissue = "joints", label = "ACLI", groups = c("orig.ident", "day1", "day7"))
dat_GSE205037 <- readscdata_FF(GSE = "GSE205037", tissue = "spinal_cord", label = "SCI", groups = c("orig.ident", "day3", "day7"))
dat_GSE205690 <- readscdata_FF(GSE = "GSE205690", tissue = "muscle", label = "SMI", groups = c("orig.ident", "day4", "day7"))

datall <- Reduce(rbind, list(dat_6035873,
                             dat_GSE141259,
                             dat_GSE152501,
                             dat_GSE163465,
                             dat_GSE180420,
                             dat_GSE186986,
                             dat_GSE200843,
                             dat_GSE205037,
                             dat_GSE205690))
datall3 <- datall %>% 
  group_by(tissue) %>%
  summarise(n=n())

pdf(file = "./stat1.pdf", width = 5, height = 3)
ggplot(datall3, aes(x=tissue, y=n))+
  geom_bar(stat = "identity")+
  geom_text(aes(x = tissue, y = n+1000, label = n),size=3)+
  theme_classic()+
  xlab(NULL)+
  ylab("Number of cells")+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()

datall2 <- datall %>% 
  group_by(tissue, group) %>%
  summarise(n=n(),
            mean_FAS = mean(FAS),
            mean_FAO = mean(FAO))

pdf(file = "./stat2.pdf", width = 7, height = 3)
ggplot(datall2, aes(fill=group, y=n, x=tissue)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = c("#f1939c","#5e90b8"))+
  xlab(NULL)+
  ylab("Proportion")+
  labs(fill = "Phase")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()

datall2 <- gather(datall2, metabolism, score, -group, -tissue, -n)
datall2$group <- factor(datall2$group, levels = c("Phase II", "Phase I"))
datall2$metabolism <- gsub("mean_", "", datall2$metabolism)


p <- ggplot(datall2, aes(x=metabolism, y=group, fill= score)) + 
  geom_tile(colour = "white", linewidth = 3)+
  scale_fill_gradientn(colours = c("#4DBBD5","#4DBBD5DD","#4DBBD599","white", "#E64B3599","#E64B35DD","#E64B35"), 
                       limits= c(-1.2,1.2), 
                       name=c("Metabolic score"))+
  theme_bw()+
  theme(axis.title = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(~tissue, nrow = 2, scales = "free_x")

pdf(file = "./FASFAO.pdf", width = 7, height = 3)
p
dev.off()

##### Fig. S1c  #####
rm(list = ls())

readscdata_Other <- function(GSE, tissue, label, groups){
  mac <- readRDS(paste0("./", GSE, "/RData/result_mac_", tissue, ".rds"))
  load(paste0("./", GSE, "/RData/result_mac_metabolism.RData"))
  
  phenoData <- mac@meta.data
  phenoData <- cbind(phenoData, metaScore)
  
  phenoData$group <- ifelse(phenoData[,groups[1]] == groups[2], "Phase I", "Phase II")
  
  dat <- phenoData[,c("group", "TCA", "Glycolysis", "OXPHOS", "PPP")]
  dat$tissue <- label
  
  dat
}

dat_6035873 <- readscdata_Other(GSE = "6035873", tissue = "liver", label = "LIVER",groups = c("orig.ident", "h48"))
dat_GSE141259 <- readscdata_Other(GSE = "GSE141259", tissue = "lung", label = "LUNG",groups = c("grouping", "d10", "d28"))
dat_GSE152501 <- readscdata_Other(GSE = "GSE152501", tissue = "air", label = "AWI", groups = c("orig.ident", "day1", "day7"))
dat_GSE163465 <- readscdata_Other(GSE = "GSE163465", tissue = "heart", label = "MCIN", groups = c("orig.ident", "day3"))
dat_GSE180420 <- readscdata_Other(GSE = "GSE180420", tissue = "kidney", label = "KDI", groups = c("pheno", "IRI_short_3", "IRI_short_14"))
dat_GSE186986 <- readscdata_Other(GSE = "GSE186986", tissue = "skin", label = "SKIN", groups = c("orig.ident", "day1", "day5"))
dat_GSE200843 <- readscdata_Other(GSE = "GSE200843", tissue = "joints", label = "ACLI", groups = c("orig.ident", "day1", "day7"))
dat_GSE205037 <- readscdata_Other(GSE = "GSE205037", tissue = "spinal_cord", label = "SCI", groups = c("orig.ident", "day3", "day7"))
dat_GSE205690 <- readscdata_Other(GSE = "GSE205690", tissue = "muscle", label = "SMI", groups = c("orig.ident", "day4", "day7"))

datall <- Reduce(rbind, list(dat_6035873,
                             dat_GSE141259,
                             dat_GSE152501,
                             dat_GSE163465,
                             dat_GSE180420,
                             dat_GSE186986,
                             dat_GSE200843,
                             dat_GSE205037,
                             dat_GSE205690))


datall2 <- datall %>% 
  group_by(tissue, group) %>%
  summarise(n=n(),
            mean_TCA = mean(TCA),
            mean_Glycolysis = mean(Glycolysis),
            mean_OXPHOS = mean(OXPHOS),
            mean_PPP = mean(PPP))

datall2 <- gather(datall2, metabolism, score, -group, -tissue, -n)
datall2$group <- factor(datall2$group, levels = c("Phase II", "Phase I"))
datall2$metabolism <- gsub("mean_", "", datall2$metabolism)


p <- ggplot(datall2, aes(x=metabolism, y=group, fill= score)) + 
  geom_tile(colour = "white", linewidth = 3)+
  scale_fill_gradientn(colours = c("#4DBBD5","#4DBBD5DD","#4DBBD599","white", "#E64B3599","#E64B35DD","#E64B35"), 
                       limits= c(-1.2,1.2), 
                       name=c("Metabolic score"))+
  theme_bw()+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
  facet_wrap(~tissue, nrow = 2, scales = "free_x")

pdf(file = "./metaOther.pdf", width = 12, height = 4)
p
dev.off()

##### Fig. 1F  #####
rm(list = ls())

# load data
autoMAC <- function(GSE, tissue, label, groups){
  require(tidyverse)
  require(Seurat)
  require(SingleR)
  require(SummarizedExperiment)
  require(scuttle)
  require(ALRA)
  
  
  mac <- readRDS(paste0("./", GSE, "/RData/result_mac_", tissue, ".rds"))
  m012ac_se <- readRDS("../scRNAseq SingleR/RData/result_ref_m012ac.rds")
  
  exprSet <- GetAssayData(mac, layer = "counts") %>% as.matrix %>% t
  exprSet_alra <- alra(exprSet)
  exprSet_alra <- exprSet_alra[[3]]
  rownames(exprSet_alra) <- rownames(exprSet)
  mac_counts <- t(exprSet_alra)
  
  # intersect
  common_gene <- intersect(rownames(mac_counts), rownames(m012ac_se))
  mac_counts <- mac_counts[common_gene,]
  m012ac_se <- m012ac_se[common_gene,]
  
  # create SummarizedExperiment object and log normalize
  mac_se <- SummarizedExperiment(assays = list(counts = mac_counts))
  mac_se <- logNormCounts(mac_se)
  
  # annotate
  mac_res <- SingleR(test = mac_se, ref = m012ac_se, labels = m012ac_se$ref_label)
  
  # export the result
  anno_df <- mac_res@listData$scores %>% 
    as.data.frame() %>% 
    add_column(subtype = mac_res$labels)%>% 
    add_column(index = rownames(mac_res))
  
  phenoData <- mac@meta.data %>% rownames_to_column("index")
  
  phenoData$group <- ifelse(phenoData[,groups[1]] == groups[2], "Phase I", "Phase II")
  
  res <- merge(anno_df, phenoData, by = "index")
  res <- res[,c("index", "subtype", "group")]
  res$tissue <- label
  
  res
  
}

mtype_GSE186986 <- autoMAC(GSE = "GSE186986", tissue = "skin", label = "SKIN", groups = c("orig.ident", "day1", "day5"))

# transform Seurat object to monocle
mac <- readRDS("./GSE186986/RData/result_mac_skin.rds")

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

save(cds, file = "./GSE186986/RData/result_monocle_skin.RData")
# load("./GSE186986/RData/result_monocle.RData")


dims <- cds@reducedDimS %>% t %>% as.data.frame %>%
  magrittr::set_colnames(c("Dim1", "Dim2")) %>%
  cbind(cds@phenoData@data)

p1 <- plot_cell_trajectory(cds, show_branch_points = F, cell_size = NA)+
  geom_scattermore(data = dims, aes(x= Dim1, y= Dim2, color = orig.ident), pointsize = 18, pixels = c(2048, 2048))+
  scale_color_manual(values = c("day1" = "#f1939c", "day5" = "#5e90b8"), labels = c("Phase I", "Phase II"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

p2 <- plot_cell_trajectory(cds, show_branch_points = F, cell_size = NA)+
  geom_scattermore(data = dims, aes(x= Dim1, y= Dim2, color = subtype), pointsize = 18, pixels = c(2048, 2048))+
  scale_color_manual(values = c("M0" = "#E64B35", "M1_IFNg" = "#4DBBD5", "M1_LPS" = "#00A087", "M2a_IL4" = "#3C5488", "M2c_IL10" = "#F39B7F"), 
                     labels = c("M0", "M(IFNg)", "M(LPS)", "M(IL4)", "M(IL10)"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

typall <- mtype_GSE186986

typall2 <- typall %>% 
  group_by(tissue, group, subtype) %>%
  summarise(n=n())


p3 <- ggplot(typall2, aes(fill=subtype, y=n, x=group)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = c("M0" = "#E64B35", "M1_IFNg" = "#4DBBD5", "M1_LPS" = "#00A087", "M2a_IL4" = "#3C5488", "M2c_IL10" = "#F39B7F"),guide='none')+
  xlab(NULL)+
  ylab("Proportion")+
  labs(fill = "Phase")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))



pdf(file = "./scRNA_mactype.pdf", width = 8, height = 3)
p1+p2+p3+plot_layout(nrow = 1, widths = c(3,3,1), guides = "collect")
dev.off()


##### Fig. S9D  #####
rm(list = ls())

autoCycle <- function(GSE, tissue, label, groups){

  require(scran)
  require(scuttle)
  require(org.Mm.eg.db)
  require(SummarizedExperiment)
  
  mac <- readRDS(paste0("./", GSE, "/RData/result_mac_", tissue, ".rds"))
  
  phenoData <- mac@meta.data
  phenoData$group <- ifelse(phenoData[,groups[1]] == groups[2], "Phase I", "Phase II")
  mac@meta.data$orig.ident <- phenoData$group
  
  mac_counts <- GetAssayData(mac, slot = "counts") %>% as.matrix
  
  mac_se <- SingleCellExperiment(list(counts = mac_counts))
  
  ensembl <- mapIds(org.Mm.eg.db, keys=rownames(mac_se), keytype="SYMBOL", column="ENSEMBL")
  
  mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  
  assigned <- cyclone(mac_se, pairs=mm.pairs, gene.names=ensembl)
  
  assigned 
  
  dat <- data.frame(phase = mac@meta.data$orig.ident,
                    cycle = assigned$phases)
  
  dat$tissue <- label
  dat
}

dat_6035873 <- autoCycle(GSE = "6035873", tissue = "liver", label = "LIVER",groups = c("orig.ident", "h48"))
dat_GSE141259 <- autoCycle(GSE = "GSE141259", tissue = "lung", label = "LUNG",groups = c("grouping", "d10", "d28"))
dat_GSE152501 <- autoCycle(GSE = "GSE152501", tissue = "air", label = "AWI", groups = c("orig.ident", "day1", "day7"))
dat_GSE163465 <- autoCycle(GSE = "GSE163465", tissue = "heart", label = "MCIN", groups = c("orig.ident", "day3"))
dat_GSE180420 <- autoCycle(GSE = "GSE180420", tissue = "kidney", label = "KDI", groups = c("pheno", "IRI_short_3", "IRI_short_14"))
dat_GSE186986 <- autoCycle(GSE = "GSE186986", tissue = "skin", label = "SKIN", groups = c("orig.ident", "day1", "day5"))
dat_GSE200843 <- autoCycle(GSE = "GSE200843", tissue = "joints", label = "ACLI", groups = c("orig.ident", "day1", "day7"))
dat_GSE205037 <- autoCycle(GSE = "GSE205037", tissue = "spinal_cord", label = "SCI", groups = c("orig.ident", "day3", "day7"))
dat_GSE205690 <- autoCycle(GSE = "GSE205690", tissue = "muscle", label = "SMI", groups = c("orig.ident", "day4", "day7"))


datall <- Reduce(rbind, list(dat_6035873,
                             dat_GSE141259,
                             dat_GSE152501,
                             dat_GSE163465,
                             dat_GSE180420,
                             dat_GSE186986,
                             dat_GSE200843,
                             dat_GSE205037,
                             dat_GSE205690))

save(datall, file = "result_cell_cyc.RData")


datall$num <- 1

res <- datall %>%
  group_by(tissue, phase) %>%
  mutate(sum = sum(num)) %>% 
  group_by(tissue, phase,cycle) %>% 
  mutate(sumc = sum(num))

res <- res[!duplicated(res),]
res$por <- res$sumc/res$sum

res <- na.omit(res)

res <- res[,c("phase", "cycle","tissue", "por")]

res <- spread(res, cycle, por)

res[is.na(res)] <- 0

res <- gather(res, cycle, num, -phase, -tissue)

res2 <- res %>%
  group_by(phase, cycle) %>% 
  summarise(mean = mean(num),
            n = 9,
            sd = sd(num),
            sem = sd/sqrt(n))
  


pdf(file = "./cycle.pdf", width = 7, height = 4)
ggplot(res2, aes(x = phase, y = mean, fill = cycle))+
  geom_bar(position = position_dodge(1), stat = "identity", color = "black")+
  geom_errorbar(aes(x = phase, ymin=mean-sem, ymax=mean+sem),
                position = position_dodge(1), width=0.4, colour="black", alpha=0.9, size=1)+
  geom_point(data = res, aes(y = num),size = 1.5, show.legend = FALSE, color = "black",
              position = position_dodge(1))+
  scale_fill_manual(values = c("#00A087", "#E64B35", "#4DBBD5"))+
  labs(x= "", y = "Proportion")+
  theme_classic()+
  theme(legend.title = element_blank())
dev.off()

