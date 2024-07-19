
# library
library(limma) # v3.58.1
library(clusterProfiler) # v4.10.1
library(org.Mm.eg.db) # v3.18.0
library(mixOmics) # v6.26.0
library(ggpubr) # v0.6.0
library(scattermore) # v1.2
library(ggrepel) # v0.9.5
library(tidyverse) # v2.0.0

# the input files "M2c_M0.all.txt" can be obtained from Supplementary Table 3 for the source article or from the corresponding author upon request

#### Fig. 2a ####
rm(list = ls())

all <- read.table(file = "M2c_M0.all.txt", header = T, sep = "\t", fill = T, quote = "")

test <- data.frame(logFC = log2(as.numeric(all$ratio)),
                   pvalue = all$t.test_p.value,
                   level = all$MS2classlevel1,
                   name = all$MS2Metabolite)

test$name[test$level != "Glycerolipids [GL]"] = ""
test$name[abs(test$logFC) < log2(1.5)] = ""
test$name[test$pvalue > 0.05] = ""

test$change = ifelse(test$pvalue < 0.05 & abs(test$logFC) > log2(1.5), ifelse(test$logFC> log2(1.5),'Up in M2c','Up in M0'), 'Stable')
test[,2] <- -log10(test$pvalue)
colnames(test)[2] <- "v"

p.v <- ggplot(test, aes(x = logFC, y = v, color = change))+
  geom_scattermore(pointsize = 8, pixels = c(2048, 2048))+
  scale_color_manual(values = c("#A0A0A499","#5e90b899",  "#f1939c99"))+
  xlim(-6,6)+
  labs(x = "log2FC (M(IL10)/M0)", y ="-log10(P value)")+
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



pdf(file = "./v_M2c_M0.pdf", width = 5, height = 4.2)
p.v
dev.off()

#### Fig. S3a ####
rm(list = ls())

all <- read.table(file = "M2c_M0.all.txt", header = T, sep = "\t", fill = T, quote = "")

X=all[,2:9] %>% t()
Y=c("M0", "M0", "M0", "M0", "M(IL10)", "M(IL10)", "M(IL10)", "M(IL10)")

plsda.datatm <-plsda(X, Y, ncomp = 3)
df <- unclass(plsda.datatm)

df1 = as.data.frame(df$variates$X)
df1$group = c("M0", "M0", "M0", "M0", "M(IL10)", "M(IL10)", "M(IL10)", "M(IL10)")
df1$samples = rownames(df1)

explain = df$prop_expl_var$X
x_lable <- round(explain[1],digits=3)
y_lable <- round(explain[2],digits=3)

col=c("#f1939c", "#5e90b8")
p1<-ggplot(df1,aes(x=comp1,y=comp2,
                   color=group,shape=group))+
  theme_bw()+
  geom_point(size=2)+
  theme(panel.grid = element_blank())+
  ggtitle("PLS-DA")+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("X-variate1 (",x_lable*100,"%)"),
       y=paste0("X-variate2 (",y_lable*100,"%)"))+
  stat_ellipse(data=df1,
               geom = "polygon",level = 0.95,
               linetype = 2,linewidth=0.5,
               aes(fill=group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = col,guide = "none") +
  scale_fill_manual(values = c("#f1939c", "#5e90b8"))+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        plot.title = element_text(hjust=0.5),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())

pdf(file = "./PLS-DA.pdf", width = 5, height = 4)
p1
dev.off()

#### Fig. 2b ####
rm(list = ls())

enrich <- read.table("enrichment.txt", header = T, sep = "\t")

enrich <- enrich[enrich$Level1 == "Metabolism",]
enrich$GeneRatio <- enrich$Diff.metbolites..NumCompound./27
enrich <- enrich[enrich$Pvalue <= 0.05,]

showCategory =12
font.size =12
p<-enrich %>% 
  arrange(Pvalue) %>% 
  slice(1:showCategory) %>% 
  ggplot(aes(x=forcats::fct_reorder(Pathway,Pvalue,.desc = T),y=GeneRatio,fill=Pvalue))+ 
  geom_bar(stat="identity")+
  coord_flip()+
  scale_fill_continuous(low="#f1939c",high="#5e90b8", guide=guide_colorbar(reverse=TRUE))+
  labs(x=NULL) +
  ggtitle("")+
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))



pdf(file = "GO.pdf", width = 6, height = 5)
p
dev.off()





