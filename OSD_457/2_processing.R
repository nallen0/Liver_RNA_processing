###processing BNL1###
##Starting from featureCounts instead of RSEM counts, preparing data for GL pipeline

## This is the runsheet created in GL RCP processing
runsheet_path = "./GLDS-457_rna_seq_bulkRNASeq_v1_runsheet_liver.csv"
counts_dir = "./GLDS-457_rna_seq_RSEM_Unnormalized_Counts_liver.csv"
norm_output = "./norm_output"
DGE_output = "./dge_output"
figure_output = "./figures"

#read in counts
cts <- read.csv(counts_dir, 
                header = TRUE, row.names = 'X')
head(cts)

str(cts)
summary(cts)
any(cts %% 1 != 0)  # TRUE if there are non-integer values
cts <- round(cts)

compare_csv_from_runsheet <- function(runsheet_path) {
  df = read.csv(runsheet_path)
  # get only Factor Value columns
  factors = as.data.frame(df[,grep("Factor.Value", colnames(df), ignore.case=TRUE)])
  colnames(factors) = paste("factor",1:dim(factors)[2], sep= "_")
  result = data.frame(sample_id = df[,c("Original.Sample.Name")], factors)	
  return(result)
}

### Load metadata from runsheet csv file ###

compare_csv <- compare_csv_from_runsheet(runsheet_path)
cdata <- compare_csv
rownames(cdata) <- cdata[,1]
colnames(cdata) <- c('sample_id', 'condition')

all(rownames(cdata) == colnames(cts))
cts <- cts[, rownames(cdata)]
all(rownames(cdata) == colnames(cts))

study <- as.data.frame(compare_csv[,2:dim(compare_csv)[2]])
colnames(study) <- colnames(compare_csv)[2:dim(compare_csv)[2]]
rownames(study) <- compare_csv[,1]
dim(study)

if (dim(study) >= 2){
  group<-apply(study,1,paste,collapse = " & ") ## concatenate multiple factors into one condition per sample
} else{
  group<-study[,1]
}

group_names <- paste0("(",group,")",sep = "") ## human readable group names
group <- sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", group))) # group naming compatible with R models, this maintains the default behaviour of make.names with the exception that 'X' is never prepended to group names
names(group) <- group_names
rm(group_names)

contrast.names <- combn(levels(factor(names(group))),2) ## generate matrix of pairwise group combinations for comparison
contrasts <- apply(contrast.names, MARGIN=2, function(col) sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2))))) # limited make.names call for each group (also removes leading parentheses)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) ## format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names) 
contrasts

### Make DESeqDataSet object ###
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = cdata,
                              design = ~condition)

sampleTable <- data.frame(rownames(study))
sampleTable$Condition <- (condition=factor(group))
colnames(sampleTable) <- c("id","condition")

### Filter out genes with counts less than 10 in any sample ###
keepGenes <- apply(counts(dds),1, function(x){all(x>10)})

dds_1 <- dds[keepGenes,]
summary(dds_1)

### Run DESeq analysis ###
dds_1 <- DESeq(dds_1)
res <- results(dds_1)

head(results(dds_1, tidy=TRUE))
plotMA(dds_1, ylim=c(-10,10))

###Generate metadata
metad <- sampleTable
rownames(metad) <- sampleTable$id

### Generate F statistic p-value (similar to ANOVA p-value) using DESeq2 likelihood ratio test (LRT) design ###
dds_1_lrt <- DESeq(dds_1, test = "LRT", reduced = ~1)
res_1_lrt <- results(dds_1_lrt)
plotMA(dds_1_lrt, ylim=c(-10,10))

### Create a data frame containing normalized counts ###
normCounts <- as.data.frame(counts(dds_1, normalized=TRUE))

### Add 1 to all normalized counts to avoid issues with downstream calculations ###
normCounts <- normCounts +1


### Start the DGE output table with the normalized counts for all samples ###
## reduced output table 1 will be used to generate human-readable DGE table
reduced_output_table_1 <- normCounts

## output tables 1 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots
output_table_1 <- normCounts

### Iterate through Wald Tests to generate pairwise comparisons of all groups ###
for (i in 1:dim(contrasts)[2]){
  res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
  res_1 <- as.data.frame(res_1@listData)[,c(2,4,5,6)]
  colnames(res_1) <- c(paste0("Log2fc_", colnames(contrasts)[i]), paste0("Stat_",colnames(contrasts)[i]), paste0("P.value_",colnames(contrasts)[i]), paste0("Adj.p.value_",colnames(contrasts)[i]))
  output_table_1 <- cbind(output_table_1,res_1)
  reduced_output_table_1 <- cbind(reduced_output_table_1,res_1)
  rm(res_1)
}

### Generate and add all sample mean column to the DGE table ###
output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)
reduced_output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)

### Generate and add all sample stdev column to the DGE table ###
output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)
reduced_output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)

### Add F statistic p-value (similar to ANOVA p-value) column to the DGE table ###
output_table_1$LRT.p.value <- res_1_lrt@listData$padj
reduced_output_table_1$LRT.p.value <- res_1_lrt@listData$padj

### Generate and add group mean and stdev columns to the DGE table ###
tcounts <- as.data.frame(t(normCounts))
tcounts$group <- names(group) # Used final table group name formatting (e.g. '( Space Flight & Blue Light )' )

group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean))) # Compute group name group-wise means
colnames(group_means) <- paste0("Group.Mean_", group_means['group',]) # assign group name as column names

group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd))) # Compute group name group-wise standard deviation
colnames(group_stdev) <- paste0("Group.Stdev_", group_stdev['group',]) # assign group name as column names

group_means <- group_means[-c(1),] # Drop group name row from data rows (now present as column names)
group_stdev <- group_stdev[-c(1),] # Drop group name row from data rows (now present as column names)
output_table_1 <- cbind(output_table_1,group_means, group_stdev) # Column bind the group-wise data
reduced_output_table_1 <- cbind(reduced_output_table_1,group_means, group_stdev) # Column bind the group-wise data

rm(group_stdev,group_means,tcounts)

### Add columns needed to generate GeneLab visualization plots to the DGE table ###
## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison ##
updown_table <- sign(output_table_1[,grep("Log2fc_",colnames(output_table_1))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,updown_table)
rm(updown_table)

## Add column to indicate contrast significance with p <= 0.1 ##
sig.1_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.1_table)
rm(sig.1_table)

## Add column to indicate contrast significance with p <= 0.05 ##
sig.05_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.05_table)
rm(sig.05_table)

## Add columns for the volcano plot with p-value and adjusted p-value ##
log_pval_table <- log2(output_table_1[,grep("P.value_",colnames(output_table_1))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_1 <- cbind(output_table_1,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_1[,grep("Adj.p.value_",colnames(output_table_1))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_1 <- cbind(output_table_1,log_adj_pval_table)
rm(log_adj_pval_table)

### Read in GeneLab annotation table for the organism of interest ###
annot <- read.table(annotations_link, sep = "\t", header = TRUE, quote = "", comment.char = "", row.names = 1)

### Combine annotations table and the DGE table ###
output_table_1 <- merge(annot, output_table_1, by='row.names', all.y=TRUE)
output_table_1 <- output_table_1 %>% 
  dplyr::rename(
    ENSEMBL = Row.names ## Change ENSEMBL to TAIR for plant studies ##
  )

reduced_output_table_1 <- merge(annot, reduced_output_table_1, by='row.names', all.y=TRUE)
reduced_output_table_1 <- reduced_output_table_1 %>% 
  dplyr::rename(
    ENSEMBL = Row.names ## Change ENSEMBL to TAIR for plant studies ##
  )

view(reduced_output_table_1%>%
       filter(SYMBOL=="Cdkn1a")%>%
       dplyr::select(contains("v(WT_NL)"))%>%
       dplyr::select(contains("_IR_"))
)

### Export unnormalized and normalized counts tables ###
normCounts_exp <- as.data.frame(counts(dds_1, normalized=TRUE))

write.csv(txi.rsem$counts,file.path(norm_output, "RSEM_Unnormalized_Counts.csv"))
write.csv(normCounts_exp,file.path(norm_output, "Normalized_Counts.csv"))

### Export sample grouping and contrasts tables ###
write.csv(sampleTable,file.path(DGE_output, "SampleTable.csv"))
write.csv(contrasts,file.path(DGE_output, "contrasts.csv"))

### Export human-readable DGE table ###
write.csv(reduced_output_table_1,file.path(DGE_output, "differential_expression.csv"), row.names = FALSE)

### Export computer-readable DGE and PCA tables used for GeneLab visualization ###
write.csv(output_table_1,file.path(DGE_output, "visualization_output_table.csv"), row.names = FALSE)
write.csv(PCA_raw$x,file.path(DGE_output, "visualization_PCA_table.csv"), row.names = TRUE)

### print session info ###
write.csv(org_table,file.path("org_table.csv"), row.names = FALSE)
write.csv(annot,file.path("annot.csv"), row.names = FALSE)

print("Session Info below: ")
sessionInfo()