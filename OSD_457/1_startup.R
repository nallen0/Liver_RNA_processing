# Packages to install if necessary

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require(devtools)) install.packages("devtools")

# Packages to install if necessary
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("SPIA")
BiocManager::install("pcaExplorer")
BiocManager::install("ggfortify")
BiocManager::install("ggplot")
BiocManager::install("VennDiagram")
BiocManager::install("ggVennDiagram")
BiocManager::install("GEOquery")
BiocManager::install("PCAtools")
BiocManager::install("bind.r")
BiocManager::install("DEGreport")
BiocManager::install("lasso2")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("venneuler")
BiocManager::install("pathview")
BiocManager::install("DOSE")
BiocManager::install("stringr")
BiocManager::install("enrichR")

devtools::install_github("gaospecial/ggVennDiagram")

remotes::install_github("YuLab-SMU/ggtree")
remotes::install_version("matrixStats", version="1.1.0") # restart your session and run previous scripts

install.packages("babelgene")
install.packages("tidyverse")
install.packages("progress")
install.packages("webshot2")
webshot2::install_phantomjs()

# Import libraries 
library(webshot2)
library(progress)
library(BiocGenerics)
library(tximport)
library(DESeq2)
library(tidyverse)
library(stringr)
library(DESeq2)
library(EnhancedVolcano)
library(SPIA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(sjPlot)
library(gt)
library(pcaExplorer)
library(ggfortify)
library(ggplot2)
library(genefilter)
library(vsn)
library(VennDiagram)
library(ggVennDiagram)
library(gt)
library(enrichR)
library(DEP)
library(dplyr)
library(pcaExplorer)
library(PCAtools)
library(matrixStats)
library(bind.r)
library(rowr)
library(purrr)
library(tibble)
library(DEGreport)
library(clusterProfiler)
library(ggtree)
library(glue)
library(UpSetR)
library(venneuler)
library(pathview)
library(biomaRt)
library(DOSE)
library(babelgene)
library(DEGreport)
library(ggplot2)
library(pheatmap)


# Key Functions
#cbind function
bind_cols_fill <- function(df_list) {
  
  max_rows <- map_int(df_list, nrow) %>% max()
  
  map(df_list, function(df) {
    if(nrow(df) == max_rows) return(df)
    first <- names(df)[1] %>% sym()
    df %>% add_row(!!first := rep(NA, max_rows - nrow(df)))
  }) %>% bind_cols()
}


### Define which organism is used in the study - this should be consistent with the name in the "name" column of the GL-DPPD-7110_annotations.csv file, which matches the abbreviations used in the Panther database for each organism ###

organism <- "MOUSE"

### Pull in the GeneLab annotation table (GL-DPPD-7110_annotations.csv) file ###

org_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv"

org_table <- read.table(org_table_link, sep = ",", header = TRUE)


### Define the link to the GeneLab annotation table for the organism of interest ###

annotations_link <- org_table[org_table$name == organism, "genelab_annots_link"]

annot <- read.table(annotations_link, sep = "\t", header = TRUE, quote = "", comment.char = "", row.names = 1)


### Pull all factors for each sample in the study from the runsheet created in Step 9a ###

compare_csv_from_runsheet <- function(runsheet_path) {
  df = read.csv(runsheet_path)
  # get only Factor Value columns
  factors = as.data.frame(df[,grep("Factor.Value", colnames(df), ignore.case=TRUE)])
  colnames(factors) = paste("factor",1:dim(factors)[2], sep= "_")
  result = data.frame(sample_id = df[,c("Original.Sample.Name")], factors)	
  return(result)
}


### print session info ###

print("Session Info below: ")
sessionInfo()
BiocManager::valid()
version