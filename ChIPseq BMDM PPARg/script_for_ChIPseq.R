# library
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(Gviz)
library(biomaRt)
library(rtracklayer)
library(GenomicFeatures)
library(tidyverse)

#### Fig. 5a ####
rm(list = ls())

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak_ta <- readPeakFile("./TA_summits.bed")
peak_wt <- readPeakFile("./WT_summits.bed")

# annotate peak
peakAnno_ta <- annotatePeak(peak_ta, 
                            tssRegion = c(-2500, 2500), 
                            TxDb = txdb,
                            annoDb = "org.Mm.eg.db")
peakAnno_wt <- annotatePeak(peak_wt, 
                            tssRegion = c(-2500, 2500), 
                            TxDb = txdb, 
                            annoDb = "org.Mm.eg.db")
res_ta <- as.data.frame(peakAnno_ta)
res_wt <- as.data.frame(peakAnno_wt)

res_ta <- res_ta[grep("Promoter", res_ta$annotation),]
res_wt <- res_wt[grep("Promoter", res_wt$annotation),]

goList <- list(WT = unique(res_wt$SYMBOL),
               TA = unique(res_ta$SYMBOL))

go <- compareCluster(goList, 
                     fun="enrichGO",
                     OrgDb = "org.Mm.eg.db", 
                     ont = "BP", 
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 0.05, 
                     qvalueCutoff = 0.99,
                     keyType = "SYMBOL")


term <- c("fatty acid beta-oxidation", 
          "lipid modification", 
          "lipid oxidation",
          "lipid localization",
          "regulation of lipid biosynthetic process",
          "regulation of steroid metabolic process",
          "phospholipid biosynthetic process",
          "positive regulation of lipid biosynthetic process")



pdf(file = "./ChIPseq_GO.pdf", width = 6, height = 4)
dotplot(go, showCategory = term, label_format = 9999, color = "qvalue")+
  scale_color_gradientn(colours = c("#4B74B2", "#90BEF4", "#E6F1F3", "#FFDF92", "#FC8C5A", "#DB3124"))
dev.off()


#### Fig. S9a ####

mm10 <- useMart(host='nov2020.archive.ensembl.org',
                biomart='ENSEMBL_MART_ENSEMBL', 
                dataset='mmusculus_gene_ensembl')


###################Dbi#########################
genome = "mm10"
chr_no <- "chr1" 
chr_start <- 120112760
chr_end <- 120123297
peak_start <- c(120120023, 120121227)

peak_start <- c(120120023)


track_bmt <- BiomartGeneRegionTrack(biomart = mm10,
                                    chromosome = chr_no,
                                    filter = list(with_refseq_mrna = TRUE),
                                    name=NULL,
                                    cex.title= 0.8,
                                    col.line = NULL,
                                    col = NULL,
                                    col.title="white",
                                    background.title= "white",
                                    stacking = "hide")

track_pos <- GenomeAxisTrack (lwd = 0.5, fontsize = 8)

track_ta <- DataTrack(range ="./TA_logLR.bigWig", 
                      type = "polygon",
                      name="T166A",
                      col.mountain ="#732e7e",
                      fill.mountain = c("#732e7e", "#732e7e"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 0.6), 
                      transformation = function(x) { x[x < 0] <- 0; x })

track_wt <- DataTrack(range ="./WT_logLR.bigWig", 
                      type = "polygon",
                      name="WT",
                      col.mountain = "#008e59",
                      fill.mountain = c("#008e59","#008e59"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 0.6),
                      transformation = function(x) { x[x < 0] <- 0; x })

ht <- HighlightTrack(trackList = list(track_ta, track_wt),
                     start = peak_start, width = c(800),
                     chromosome = chr_no)



pdf(file = "./ChIPseq_Dbi.pdf", width = 6, height = 2)
plotTracks(list(track_pos, ht, track_bmt),
           from = chr_start,
           to = chr_end, 
           chromosome = chr_no,
           transcriptAnnotation = "symbol",
           sizes = c(1,2,2,1),
           reverseStrand = TRUE)

dev.off()

###################Plin2#########################
genome = "mm10"
chr_no <- "chr4" 
chr_start <- 86650234
chr_end <- 86675943
peak_start <- 86671686


track_bmt <- BiomartGeneRegionTrack(biomart = mm10,
                                    chromosome = chr_no,
                                    filter = list(with_refseq_mrna = TRUE),
                                    name=NULL,
                                    cex.title= 0.8,
                                    col.line = NULL,
                                    col = NULL,
                                    col.title="white",
                                    background.title= "white",
                                    stacking = "hide")

track_pos <- GenomeAxisTrack (lwd = 0.5, fontsize = 8)

track_ta <- DataTrack(range ="./TA_logLR.bigWig", 
                      type = "polygon",
                      name="T166A",
                      col.mountain ="#732e7e",
                      fill.mountain = c("#732e7e", "#732e7e"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 1.5), 
                      transformation = function(x) { x[x < 0] <- 0; x })

track_wt <- DataTrack(range ="./WT_logLR.bigWig", 
                      type = "polygon",
                      name="WT",
                      col.mountain = "#008e59",
                      fill.mountain = c("#008e59","#008e59"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 1.5),
                      transformation = function(x) { x[x < 0] <- 0; x })

ht <- HighlightTrack(trackList = list(track_ta, track_wt),
                     start = peak_start, width = 1000,
                     chromosome = chr_no)



pdf(file = "./ChIPseq_Plin2.pdf", width = 6, height = 2)
plotTracks(list(track_pos, ht, track_bmt),
           from = chr_start,
           to = chr_end, 
           chromosome = chr_no,
           transcriptAnnotation = "symbol",
           sizes = c(1,2,2,1),
           reverseStrand = TRUE)
dev.off()


###################Cd163#########################
genome = "mm10"
chr_no <- "chr6" 
chr_start <- 124300366
chr_end <- 124335166


track_bmt <- BiomartGeneRegionTrack(biomart = mm10,
                                    chromosome = chr_no,
                                    filter = list(with_refseq_mrna = TRUE),
                                    name=NULL,
                                    cex.title= 0.8,
                                    col.line = NULL,
                                    col = NULL,
                                    col.title="white",
                                    background.title= "white",
                                    stacking = "hide")

track_pos <- GenomeAxisTrack (lwd = 0.5, fontsize = 8)

track_ta <- DataTrack(range ="./TA_logLR.bigWig", 
                      type = "polygon",
                      name="T166A",
                      col.mountain ="#732e7e",
                      fill.mountain = c("#732e7e", "#732e7e"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 1), 
                      transformation = function(x) { x[x < 0] <- 0; x })

track_wt <- DataTrack(range ="./WT_logLR.bigWig", 
                      type = "polygon",
                      name="WT",
                      col.mountain = "#008e59",
                      fill.mountain = c("#008e59","#008e59"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 1),
                      transformation = function(x) { x[x < 0] <- 0; x })


pdf(file = "./ChIPseq_Cd163.pdf", width = 6, height = 2)
plotTracks(list(track_pos, track_ta, track_wt, track_bmt),
           from = chr_start,
           to = chr_end, 
           chromosome = chr_no,
           transcriptAnnotation = "symbol",
           sizes = c(1,2,2,1),
           reverseStrand = F)
dev.off()

###################Mrc1#########################
genome = "mm10"
chr_no <- "chr2" 
chr_start <- 14220530
chr_end <- 14337731


track_bmt <- BiomartGeneRegionTrack(biomart = mm10,
                                    chromosome = chr_no,
                                    filter = list(with_refseq_mrna = TRUE),
                                    name=NULL,
                                    cex.title= 0.8,
                                    col.line = NULL,
                                    col = NULL,
                                    col.title="white",
                                    background.title= "white",
                                    stacking = "hide")

track_pos <- GenomeAxisTrack (lwd = 0.5, fontsize = 8)

track_ta <- DataTrack(range ="./TA_logLR.bigWig", 
                      type = "polygon",
                      name="T166A",
                      col.mountain ="#732e7e",
                      fill.mountain = c("#732e7e", "#732e7e"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 1), 
                      transformation = function(x) { x[x < 0] <- 0; x })

track_wt <- DataTrack(range ="./WT_logLR.bigWig", 
                      type = "polygon",
                      name="WT",
                      col.mountain = "#008e59",
                      fill.mountain = c("#008e59","#008e59"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 1),
                      transformation = function(x) { x[x < 0] <- 0; x })


pdf(file = "./ChIPseq_Mrc1.pdf", width = 6, height = 2)
plotTracks(list(track_pos, track_ta, track_wt, track_bmt),
           from = chr_start,
           to = chr_end, 
           chromosome = chr_no,
           transcriptAnnotation = "symbol",
           sizes = c(1,2,2,1),
           reverseStrand = F)
dev.off()

###################Arg1#########################
genome = "mm10"
chr_no <- "chr10" 
chr_start <- 24913811
chr_end <- 24929474


track_bmt <- BiomartGeneRegionTrack(biomart = mm10,
                                    chromosome = chr_no,
                                    filter = list(with_refseq_mrna = TRUE),
                                    name=NULL,
                                    cex.title= 0.8,
                                    col.line = NULL,
                                    col = NULL,
                                    col.title="white",
                                    background.title= "white",
                                    stacking = "hide")

track_pos <- GenomeAxisTrack (lwd = 0.5, fontsize = 8)

track_ta <- DataTrack(range ="./TA_logLR.bigWig", 
                      type = "polygon",
                      name="T166A",
                      col.mountain ="#732e7e",
                      fill.mountain = c("#732e7e", "#732e7e"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 1), 
                      transformation = function(x) { x[x < 0] <- 0; x })

track_wt <- DataTrack(range ="./WT_logLR.bigWig", 
                      type = "polygon",
                      name="WT",
                      col.mountain = "#008e59",
                      fill.mountain = c("#008e59","#008e59"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 1),
                      transformation = function(x) { x[x < 0] <- 0; x })


pdf(file = "./ChIPseq_Arg1.pdf", width = 6, height = 2)
plotTracks(list(track_pos, track_ta, track_wt, track_bmt),
           from = chr_start,
           to = chr_end, 
           chromosome = chr_no,
           transcriptAnnotation = "symbol",
           sizes = c(1,2,2,1),
           reverseStrand = T)
dev.off()
