
# library
library(limma)
library(DMwR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(SummarizedExperiment)
library(scuttle)
library(tidyverse)

# global options
rm(list = ls())

# create dir
dir.create("./rawData")
dir.create("./rawData/GSE131364")
dir.create("./RData")

# download GSE131364 rawdata from GEO to ./rawData dir

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

table(duplicated(counts$tracking_id))

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

# create SummarizedExperiment object
pdata <- data.frame(sample = colnames(counts), 
                    group = str_split(colnames(counts), "_",simplify = T)[,2])
pdata$group <- gsub("Control", "M0", pdata$group)
pdata$group <- gsub("IFG", "M1_IFNg", pdata$group)
pdata$group <- gsub("IL10", "M2c_IL10", pdata$group)
pdata$group <- gsub("IL4", "M2a_IL4", pdata$group)
pdata$group <- gsub("LPS", "M1_LPS", pdata$group)

# SMOTE
ranPercent <- function(){
  per <- c()
  per[1] <- sample(1:98, size = 1, replace = F)
  per[2] <- sample(1:(99-per[1]), size = 1, replace = F)
  per[3] <- 100-per[1]-per[2]
  per/100
  
}

dat_new_list <- list()

for (i in unique(pdata$group)) {
  sample = pdata$sample[pdata$group == i]
  dat_select <- counts[,sample] %>% t()
  dat_new <- c()
  for (j in 1:17) {
    a <- ranPercent()
    dat_j <- dat_select*a
    dat_j <- colSums(dat_j)
    dat_new <- rbind(dat_new, dat_j)
  }
  dat_new <- as.data.frame(dat_new) %>% 
    remove_rownames() %>% 
    rbind(dat_select) %>% 
    magrittr::set_rownames(paste0(i,"_", "rep", 1:nrow(.)))
  
  dat_new_list[[i]] <- dat_new
}


dat_new_list <- Reduce(rbind, dat_new_list) %>% t() %>% round()

pdata <- data.frame(sample = colnames(dat_new_list), group= colnames(dat_new_list))
pdata$group <- str_split(pdata$group, "_rep", simplify = T)[,1]
rownames(pdata) <- pdata$sample
pdata$sample <- NULL
colnames(pdata) <- "ref_label"

m012ac <- SummarizedExperiment(assays = list(counts = dat_new_list),colData = pdata) 

# compute log-transformed normalized expression values
m012ac <- logNormCounts(m012ac)
saveRDS(m012ac, "./RData/result_ref_m012ac.rds")
