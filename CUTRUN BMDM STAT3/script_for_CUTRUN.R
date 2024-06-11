# library
library(Gviz)
library(biomaRt)
library(rtracklayer)
library(GenomicFeatures)
library(tidyverse)

#### Fig. 6a ####

mm10 <- useMart(host='nov2020.archive.ensembl.org',
                biomart='ENSEMBL_MART_ENSEMBL', 
                dataset='mmusculus_gene_ensembl')


###################Pdgfb#########################
genome = "mm10"
chr_no <- "chr15" 
chr_start <- 79994354
chr_end <- 80027267
peak_start <- c(80011790)


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

track_stat3 <- DataTrack(range ="./BMDM_STAT3.bigWig", 
                      type = "polygon",
                      name="STAT3",
                      col.mountain ="#732e7e",
                      fill.mountain = c("#732e7e", "#732e7e"),
                      cex.title= 0.8,
                      cex.axis=0.5,
                      col.title="white",
                      background.title="lightgray",
                      showAxis=T,
                      ylim=c(0, 50), 
                      transformation = function(x) { x[x < 0] <- 0; x })

ht <- HighlightTrack(trackList = list(track_stat3),
                     start = peak_start, width = c(4000),
                     chromosome = chr_no)



pdf(file = "./CUTRUN_Pdgfb.pdf", width = 6, height = 1.5)
plotTracks(list(track_pos, ht, track_bmt),
           from = chr_start,
           to = chr_end, 
           chromosome = chr_no,
           transcriptAnnotation = "symbol",
           sizes = c(1,2,1),
           reverseStrand = TRUE)

dev.off()



###################Tgfb1#########################
genome = "mm10"
chr_no <- "chr7" 
chr_start <- 25685260
chr_end <- 25707335
peak_start <- c(25686361)


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

track_stat3 <- DataTrack(range ="./BMDM_STAT3.bigWig", 
                         type = "polygon",
                         name="STAT3",
                         col.mountain ="#732e7e",
                         fill.mountain = c("#732e7e", "#732e7e"),
                         cex.title= 0.8,
                         cex.axis=0.5,
                         col.title="white",
                         background.title="lightgray",
                         showAxis=T,
                         ylim=c(0, 50), 
                         transformation = function(x) { x[x < 0] <- 0; x })

ht <- HighlightTrack(trackList = list(track_stat3),
                     start = peak_start, width = c(2500),
                     chromosome = chr_no)



pdf(file = "./CUTRUN_Tgfb1.pdf", width = 6, height = 1.5)
plotTracks(list(track_pos, ht, track_bmt),
           from = chr_start,
           to = chr_end, 
           chromosome = chr_no,
           transcriptAnnotation = "symbol",
           sizes = c(1,2,1),
           reverseStrand = F)

dev.off()

###################Vegfb#########################
genome = "mm10"
chr_no <- "chr19" 
chr_start <- 6981650
chr_end <- 6989244
peak_start <- c(6987016)


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

track_stat3 <- DataTrack(range ="./BMDM_STAT3.bigWig", 
                         type = "polygon",
                         name="STAT3",
                         col.mountain ="#732e7e",
                         fill.mountain = c("#732e7e", "#732e7e"),
                         cex.title= 0.8,
                         cex.axis=0.5,
                         col.title="white",
                         background.title="lightgray",
                         showAxis=T,
                         ylim=c(0, 40), 
                         transformation = function(x) { x[x < 0] <- 0; x })

ht <- HighlightTrack(trackList = list(track_stat3),
                     start = peak_start, width = c(1500),
                     chromosome = chr_no)



pdf(file = "./CUTRUN_Vegfb.pdf", width = 6, height = 1.5)
plotTracks(list(track_pos, ht, track_bmt),
           from = chr_start,
           to = chr_end, 
           chromosome = chr_no,
           transcriptAnnotation = "symbol",
           sizes = c(1,2,1),
           reverseStrand = T)

dev.off()

