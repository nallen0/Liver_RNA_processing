
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gplots)
library(gtools)
library(clipr)
library(stringr)
library(VennDiagram)
if (!require('BiocManager', quietly=TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
BiocManager::install('SPIA')
BiocManager::install('pcaExplorer')
library(EnhancedVolcano)
library(SPIA)
library(org.Hs.eg.db)
library(sjPlot)
library(gt)
library(pcaExplorer)

R.version

glds242 <- read.csv('GLDS-242_rna_seq_differential_expression.csv', header = TRUE)

glds48 <- read.csv('GLDS-48_rna_seq_differential_expression.csv', header = TRUE)

glds47 <- read.csv('GLDS-47_rna_seq_differential_expression.csv', header = TRUE)

glds137 <- read.csv('GLDS-137_rna_seq_differential_expression.csv', header = TRUE)

glds202 <- read.csv('GLDS-202_rna_seq_differential_expression.csv', header = TRUE)

glds173 <- read.csv('GLDS-173_rna_seq_differential_expression.csv', header = TRUE)

glds25 <- read.csv('GLDS-25_array_differential_expression.csv', header = TRUE)

#glds173 is rnaseq data and glds25 is rna microarray data of the same samples

EnhancedVolcano(glds242,
                lab=glds242$SYMBOL, 
                x='Log2fc_FLT_C1vGC_C2', y= 'Adj.p.value_.FLT_C1.v.GC_C2.',
                pCutoff = 5e-2,
                FCcutoff = 2,
                xlim = c(-5.5,5.5),
                ylim = c(0.0,5.0),
                pointSize = 1.0,
                drawConnectors = TRUE, widthConnectors = 0.75,
                title = 'GLDS 242 Liver - FLT vs GC'
                )

EnhancedVolcano(glds48,
                lab=glds48$SYMBOL, 
                x='Log2fc_.FLT_I.v.GC_I.', y= 'Adj.p.value_.FLT_I.v.GC_I.',
                pCutoff = 5e-2,
                FCcutoff = 2,
                xlim = c(-10,10),
                ylim = c(0.0,15.0),
                pointSize = 1.0,
                drawConnectors = TRUE, widthConnectors = 0.75,
                title = 'GLDS 48 Liver - FLT_I vs GC_I'
)
EnhancedVolcano(glds48,
                lab=glds48$SYMBOL, 
                x='Log2fc_.FLT_C.v.GC_C.', y= 'Adj.p.value_.FLT_C.v.GC_C.',
                pCutoff = 5e-2,
                FCcutoff = 2,
                xlim = c(-10,10),
                ylim = c(0.0,15.0),
                pointSize = 1.0,
                drawConnectors = TRUE, widthConnectors = 0.75,
                title = 'GLDS 48 Liver - FLT_C vs GC_C'
)


EnhancedVolcano(glds47,
                lab=glds47$SYMBOL, 
                x='Log2fc_.FLT.v.GC.', y= 'Adj.p.value_.FLT.v.GC.',
                pCutoff = 5e-2,
                FCcutoff = 2,
                xlim = c(-10,10),
                ylim = c(0.0,15.0),
                pointSize = 1.0,
                #drawConnectors = TRUE, widthConnectors = 0.75,
                title = 'GLDS 47 Liver - FLT vs GC'
)

EnhancedVolcano(glds137,
                lab=glds137$SYMBOL, 
                x='Log2fc_.FLT.v.GC.', y= 'Adj.p.value_.FLT.v.GC.',
                pCutoff = 5e-2,
                FCcutoff = 2,
                xlim = c(-10,10),
                ylim = c(0.0,5.0),
                pointSize = 1.0,
                drawConnectors = TRUE, widthConnectors = 0.75,
                title = 'GLDS 137 Liver - FLT vs GC'
)

EnhancedVolcano(glds173,
                lab=glds173$SYMBOL, 
                #selectLab = c('Nr1d1', 'Col15a', 'Mapk15'),
                x='Log2fc_.FLT.v.GC.', y= 'Adj.p.value_.FLT.v.GC.',
                pCutoff = 5e-2,
                FCcutoff = 2,
               # xlim = c(-5.5,5.5),
                #ylim = c(0.0,5.0),
                pointSize = 1.0,
                drawConnectors = TRUE, widthConnectors = 0.75, maxoverlapsConnectors = 20, boxedLabels = FALSE,
                title = 'GLDS 173 Liver - FLT vs GC'
)

EnhancedVolcano(glds173,
                lab=glds173$SYMBOL, 
                #selectLab = c('Nr1d1', 'Col15a', 'Mapk15'),
                x='Log2fc_.GC.v.VIV.', y= 'Adj.p.value_.GC.v.VIV.',
                pCutoff = 5e-2,
                FCcutoff = 2,
                 xlim = c(-4,4),
                ylim = c(0.0,15.0),
                pointSize = 1.0,
                drawConnectors = TRUE, widthConnectors = 0.75, maxoverlapsConnectors = 20, boxedLabels = FALSE,
                title = 'GLDS 173 Liver - GC vs VIV'
)

EnhancedVolcano(glds202,
                lab=glds202$SYMBOL, 
                x='Log2fc_.HLU_IRC_1mon.v.HLLC_IRC_1mon.', y= 'Adj.p.value_.HLU_IRC_1mon.v.HLLC_IRC_1mon.',
                pCutoff = 5e-2,
                FCcutoff = 2,
                xlim = c(-5,5),
                ylim = c(0.0,5.0),
                pointSize = 1.0,
                drawConnectors = TRUE, widthConnectors = 0.75,
                title = 'GLDS 202 Brain - HU vs NL'
) 

EnhancedVolcano(glds25,
                lab=glds25$SYMBOL, 
                x='Log2fc_.Space.Flight.v.Ground.Control.', y= 'Adj.p.value_.Space.Flight.v.Ground.Control.',
                pCutoff = 5e-2,
                FCcutoff = 1.5,
                xlim = c(-5,5),
                ylim = c(0.0,5.0),
                pointSize = 1.0,
                drawConnectors = TRUE, widthConnectors = 0.75,
                title = 'GLDS 25 Liver -micro array - Flt vs GC'
) 



geneAkr <- upDEGglds48$Log2fc_.FLT_I.v.GC_I.[]
geneAkr2 <- upDEGglds48$Log2fc_.FLT_I.v.GC_I.[2]
a<- list(geneAkr,geneAkr2)
hist(geneAkr,geneAkr2)

print(upDEGglds173)

geneinterest <- glds242[ 'SYMBOL', Cdkn1a]
print (geneinterest)
print(upDEGglds173)

#venn diagram data prep
pthreshold <- 0.05
log2fcthreshold <- 1.5
print(log2fcthreshold)
log2fcneg <- -1.5
print(log2fcneg)

glds48DEGs <- filter(glds48, Adj.p.value_.FLT_I.v.GC_I. <=pthreshold & ( Log2fc_.FLT_I.v.GC_I. >= log2fcthreshold | Log2fc_.FLT_I.v.GC_I.<=log2fcneg) )
glds48DEGs_C <- filter(glds48, Adj.p.value_.FLT_C.v.GC_C. <=pthreshold & ( Log2fc_.FLT_C.v.GC_C. >= log2fcthreshold | Log2fc_.FLT_C.v.GC_C.<=log2fcneg) )
glds48DEGsymbol <- as.list(glds48DEGs$SYMBOL)
glds48DEGensembl <- as.list(glds48DEGs$ENSEMBL)
glds48DEGsymbol_C <- as.list(glds48DEGs_C$SYMBOL)
upDEGglds48 <-filter(glds48, Adj.p.value_.FLT_I.v.GC_I. <=pthreshold &  Log2fc_.FLT_I.v.GC_I. >= log2fcthreshold)
upDEGglds48 <- as.list(upDEGglds48$SYMBOL)
downDEGglds48 <- filter(glds48, Adj.p.value_.FLT_I.v.GC_I. <=pthreshold &  Log2fc_.FLT_I.v.GC_I. <= log2fcneg)
downDEGglds48 <- as.list(downDEGglds48$SYMBOL)

glds48_IvC <- venn.diagram(
  x=list(glds48DEGsymbol,glds48DEGsymbol_C),
  category.names = c("GLDS 48 I", "GLDS 48 C"),
  output = TRUE,
  filename = 'glds48_IvC.tiff',
  na = "remove"
)
print(glds48_IvC)
plot(glds48_IvC)
overlap <- calculate.overlap(
  x=list(glds48DEGsymbol,glds48DEGsymbol_C)
  
)
print(overlap)
plot(glds48_IvC)

glds242DEGs <- filter(glds242, Adj.p.value_.FLT_C1.v.GC_C2. <=pthreshold& ( Log2fc_FLT_C1vGC_C2 >= log2fcthreshold | Log2fc_FLT_C1vGC_C2 <= log2fcneg))
glds242DEGsymbol <- as.list(glds242DEGs$SYMBOL)
glds242DEGensembl <- as.list(glds242DEGs$ENSEMBL)

glds47DEGs <- filter(glds47, Adj.p.value_.FLT.v.GC. <=pthreshold & ( Log2fc_.FLT.v.GC. >= log2fcthreshold | Log2fc_.FLT.v.GC.<=log2fcneg))
glds47DEGsymbol <- as.list(glds47DEGs$SYMBOL)
glds47DEGensembl <- as.list(glds47DEGs$ENSEMBL)

glds137DEGs <- filter(glds137, Adj.p.value_.FLT.v.GC. <=pthreshold& (Log2fc_.FLT.v.GC. >= log2fcthreshold | Log2fc_.FLT.v.GC. <= log2fcneg))
glds137DEGsymbol <- as.list(glds137DEGs$SYMBOL)
glds137DEGensembl <- as.list(glds137DEGs$ENSEMBL)

glds173DEGs <- filter(glds173, Adj.p.value_.FLT.v.GC. <=pthreshold& (Log2fc_.FLT.v.GC. >= log2fcthreshold | 'Log2fc_.FLT.v.GC.' <= log2fcneg))
glds173DEGsymbol <- as.list(glds173DEGs$SYMBOL)
glds173DEGensembl <- as.list(glds173DEGs$ENSEMBL)
upDEGglds173 <- filter(glds173, Adj.p.value_.GC.v.FLT. <=pthreshold& Log2fc_.GC.v.FLT. >= log2fcthreshold)
upDEGglds173 <- as.list(upDEGglds173$SYMBOL)
downDEGglds173 <- filter(glds173, Adj.p.value_.GC.v.FLT. <=pthreshold& Log2fc_.GC.v.FLT. <= log2fcneg)
downDEGglds173 <- as.list(downDEGglds173$SYMBOL)

glds202DEGs <- filter(glds202, Adj.p.value_.HLU_IRC_1mon.v.HLLC_IRC_1mon. <=pthreshold& (Log2fc_.HLU_IRC_1mon.v.HLLC_IRC_1mon. >= log2fcthreshold | Log2fc_.HLU_IRC_1mon.v.HLLC_IRC_1mon. <= log2fcneg))
glds202DEGsymbol <- as.list(glds137DEGs$SYMBOL)
glds202DEGensembl <- as.list(glds202DEGs$ENSEMBL)


circglds173 <- glds173 %>%
  filter(str_detect(GOSLIM_IDS, "0043153"))
circglds173exp <- data.frame( SYMBOL=circglds173$SYMBOL,Log2FC173=circglds173$Log2fc_.FLT.v.GC.)

circglds48 <- glds48 %>%
  filter(str_detect(GOSLIM_IDS, "0007623"))
circglds48exp <- data.frame( SYMBOL=circglds48$SYMBOL,Log2Fc=circglds48$Log2fc_.FLT_I.v.GC_I.)

circglds242 <- glds242 %>%
  filter(str_detect(GOSLIM_IDS, "0007623"))
circglds242exp <- data.frame( SYMBOL=circglds242$SYMBOL,Log2Fc=circglds242$Log2fc_FLT_C1vGC_C2)

circglds47 <- glds47 %>%
  filter(str_detect(GOSLIM_IDS, "0007623"))
circglds47exp <- data.frame( SYMBOL=circglds47$SYMBOL,Log2Fc=circglds47$Log2fc_.FLT.v.GC.)

circglds137 <- glds137 %>%
  filter(str_detect(GOSLIM_IDS, "0007623"))
circglds137exp <- data.frame( SYMBOL=circglds137$SYMBOL,Log2Fc=circglds137$Log2fc_.FLT.v.GC.)

m0 <- rbind(circglds173exp, circglds48exp, all = FALSE)
m1 <- rbind(m0, circglds242exp,all = FALSE)
m2 <- rbind(m1, circglds47exp,all = FALSE)
m3 <- rbind(m2, circglds137exp,all = FALSE)

rownames(m3) <- m3[,1]
m3[,1] <- NULL
m3 <- as.matrix(m3)
m3

heatmap.2(m3,na.rm = TRUE, key = TRUE, trace = c("none"), scale = c("row"),
          #key 
          density.info = c("none"),
          key.xlab = "Log2fc"
        
          )

venn.diagram(
  x=list(glds242DEGensembl,glds48DEGensembl,glds47DEGensembl,glds137DEGensembl,glds173DEGensembl),
  category.names = c("GLDS 242", "GLDS 48", "GLDS 47","GLDS 137", "GLDS 173"),
  output = TRUE,
  filename = 'venn1.tiff',
  na = "remove"
)

overlap <- calculate.overlap(
  x=list(glds173DEGsymbol,glds48DEGsymbol)

)

#vtable <-
venn(list(glds242DEGsymbol,glds48DEGsymbol,glds137DEGsymbol,glds47DEGsymbol,glds173DEGsymbol))
itemslist <- as.list( venn(glds242DEGsymbol,glds48DEGsymbol,glds137DEGsymbol,glds47DEGsymbol), show.plot = FALSE)
print(vtable)

rankedUP173 <- upDEGglds173[with(upDEGglds173,order("Log2fc_.GC.v.FLT."))]
rankedDOWN173 <- downDEGglds173[with(downDEGglds173, order("Log2fc_.GC.V.FLT."))]
write_clip(rankedUP173)

vtable2 <-venn(list(upDEGglds173,upDEGglds48))
venn(list(downDEGglds173, downDEGglds48))
print(vtable2)

#Signaling Pathway Anlaysis 
genes <-c(glds48DEGs$Log2fc_.FLT_I.v.GC_I.)
genes2 <- c(glds48DEGs$ENTREZID)
genes3 <- c(genes,genes2)
genes$V2 <- glds48DEGs$ENTREZID
genes$V1 <- glds48DEGs$Log2fc_.FLT_I.v.GC_I.



sig_genes <- glds48DEGs
sig_genes <- sig_genes[!is.na(sig_genes$ENTREZID),]
sig_genes <- sig_genes[!duplicated(sig_genes$ENTREZID),]
tg1 <- sig_genes[sig_genes$Adj.p.value_.FLT_I.v.GC_I. <=.05]
de_genes <- tg1$Log2fc_.FLT_I.v.GC_I.
names(de_genes) <- as.vector(tg1$ENTREZID)
names(sig_genes) <- subset(glds48, Adj.p.value_.FLT_I.v.GC_I. <=pthreshold, select = ENTREZID)
sig_genes <- na.omit(sig_genes)
all_genes <- glds48$ENTREZID
spia_result <- spia(de=de_genes, all=all_genes, organism="mmu", plots = TRUE)
plotP(spia_result)
print(spia_result)

sig_genes2 <- glds173DEGs
sig_genes2 <- sig_genes2[!is.na(sig_genes2$ENTREZID),]
sig_genes2 <- sig_genes2[!duplicated(sig_genes2$ENTREZID),]
tg2 <- sig_genes2[(sig_genes2$Adj.p.value_.FLT.v.GC.<=.05),]
de_genes2 <- tg2$Log2fc_.FLT.v.GC.
names(de_genes2) <- as.vector(tg2$ENTREZID)
names(sig_genes2) <- subset(glds173, Adj.p.value_.FLT.v.GC. <=pthreshold, select = ENTREZID)
all_genes2 <- glds173$ENTREZID
spia_result2 <- spia(de=de_genes2, all=all_genes2, organism="mmu", plots = TRUE)
plotP(spia_result2)
print(spia_result2)
colnames(spia_result2) <- c("Name", "Pathway ID", "Pathway Size", "Number DE Genes", "Total Perturbation", "Probability of DE genes", "Probability of Perturbation", "Combined Probability p-value", "FDR", "Adj. p-value", "Status", "KEGGLink")
subset(spia_result2, n=20)

gt_tbl <- gt(spia_result2[1:5,1:10])
gt_tbl <- 
  gt_tbl%>%
  tab_header(title = md( "**Signaling Pathway Impact Analysis - GLDS 173**"), subtitle = "Top 10 pathways shown")%>%
  tab_source_note(source_note = "Source: R v4.2.0, SPIA Package v2.24.0")%>%
  tab_source_note(source_note = md("Ref. Romero R. et. al.,(2009) *Bioinformatics*, Oxford University Press."))%>%
  cols_align(align = c( "auto"), columns = everything())

gt_tbl




sig_genes3 <- glds242DEGs
sig_genes3 <- sig_genes3[!is.na(sig_genes3$ENTREZID),]
sig_genes3 <- sig_genes3[!duplicated(sig_genes3$ENTREZID),]
tg3 <- sig_genes3[(sig_genes3$Adj.p.value_.FLT_C1.v.GC_C2. <=.05),]
de_genes3 <- tg3$Log2fc_FLT_C1vGC_C2
names(de_genes3) <- as.vector(tg3$ENTREZID)
names(sig_genes3) <- subset(glds242, Adj.p.value_.FLT_C1.v.GC_C2. <=pthreshold, select = ENTREZID)
all_genes3 <- glds242$ENTREZID
spia_result3 <- spia(de=de_genes3, all=all_genes3, organism="mmu", plots = TRUE)
plotP(spia_result3)
print(spia_result3)
subset(spia_result3, ID=="04064")

sig_genes4 <- glds137DEGs
sig_genes4 <- sig_genes4[!is.na(sig_genes4$ENTREZID),]
sig_genes4 <- sig_genes4[!duplicated(sig_genes4$ENTREZID),]
tg4 <- sig_genes4[(sig_genes4$Adj.p.value_.FLT.v.GC. <=.05),]
de_genes4 <- tg4$Log2fc_.FLT.v.GC.
names(de_genes4) <- as.vector(tg4$ENTREZID)
names(sig_genes4) <- subset(glds137, Adj.p.value_.FLT.v.GC. <=pthreshold, select = ENTREZID)
all_genes4 <- glds137$ENTREZID
spia_result4 <- spia(de=de_genes4, all=all_genes4, organism="mmu", plots = TRUE)
plotP(spia_result4)
print(spia_result4)
subset(spia_result4, ID=="04064")

sig_genes5 <- glds202DEGs
sig_genes5 <- sig_genes5[!is.na(sig_genes5$ENTREZID),]
sig_genes5 <- sig_genes5[!duplicated(sig_genes5$ENTREZID),]
tg5 <- sig_genes5[(sig_genes5$Adj.p.value_.HLU_IRC_1mon.v.HLLC_IRC_1mon. <=.05),]
de_genes5 <- tg5$Log2fc_.HLU_IRC_1mon.v.HLLC_IRC_1mon.
names(de_genes5) <- as.vector(tg5$ENTREZID)
names(sig_genes5) <- subset(glds202, Adj.p.value_.HLU_IRC_1mon.v.HLLC_IRC_1mon. <=pthreshold, select = ENTREZID)
all_genes5 <- glds202$ENTREZID
spia_result5 <- spia(de=de_genes5, all=all_genes5, organism="mmu", plots = TRUE)
plotP(spia_result5)
print(spia_result5)

gt_tbl2 <- gt(spia_result5[1:10,])
gt_tbl2 <- 
  gt_tbl2%>%
  tab_header(title = md( "**Signaling Pathway Impact Analysis**"), subtitle = "Top 10 pathways shown")%>%
  tab_source_note(source_note = "Source: R v4.2.0, SPIA Package v2.24.0")%>%
  tab_source_note(source_note = md("Ref. Romero R. et. al.,(2009) *Bioinformatics*, Oxford University Press."))
gt_tbl2

head(spia_result2)

keg2 <- spia_result2[order(spia_result2$ID),]
keg1 <- spia_result[order(spia_result$ID),]
keg3 <- spia_result5[order(spia_result5$ID),]


keg22<- keg2[(keg2$tA >= 1), ]
keg22 <- rbind(keg2[(keg2$tA<=-1),])

keg11 <- subset(keg1, ID %in% keg22$ID)
keg22 <- subset(keg22, ID %in% keg11$ID)
keg33 <- subset(keg3, ID %in% keg22$ID)
keg22 <- subset(keg22, ID %in% keg33$ID)
keg11 <- subset(keg11, ID %in% keg22$ID)

keg_tAd <- cbind(keg11$tA, keg22$tA)
rownames(keg_tAd) <- keg22$Name
colnames(keg_tAd) <- c("GLDS 173", "GLDS 48", "GLDS 202 brain")


keg_pG <- cbind(keg2$pG,keg3$pG)

rownames(keg_tA) <- spia_result2$Name
colnames(keg_tA) <- c("GLDS 173", "GLDS 48")
rownames(keg_pG) <- spia_result2$Name
colnames(keg_pG) <- c("GLDS 173", "GLDS 48")
keg_tA <- as.matrix(keg_tA)
keg_pG <- as.matrix(keg_pG)

heatmap.2(keg_tA, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")

heatmap.2(keg_tAd, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")

heatmap.2(keg_pG, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")

write_clip(glds202DEGs)

countmatrix <- as.matrix(glds202[])
pcaExplorer(countmatrix = countmatrix, coldata = coldata)
pcaExplorer()
