## Prepare PCA table for GeneLab visualization plots ##
normCounts1 <- data.frame(normCounts+1)
exp_raw <- log2(normCounts1)
PCA_raw <- prcomp(t(exp_raw), scale = F)

autoplot(PCA_raw, data=metad, colour="condition", frame=T, frame.type="norm")

# Create a tibble for LRT results
res_1_lrt_tb <- res_1_lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset to return genes with padj < 0.05
pthreshold <-  0.05
sigLRT_genes <- res_1_lrt_tb %>% 
  dplyr::filter(padj < pthreshold)

# Check number of significant genes
ngenes <- nrow(sigLRT_genes)

# Reformat meta data
metad2 <- metad
metad
metad2 <- metad2 %>%
  separate(condition, into = c("genotype", "treatment"), sep = "_")
metad3 <- metad
metad3[["genotype"]] <- metad2$genotype
metad3[["treatment"]] <- metad2$treatment
metad3

#Pull data
cluster_norm <- normCounts[sigLRT_genes$gene, ]
#cluster_norm <- cluster_norm[complete.cases(cluster_norm), ] #for issue with unfilterd data
colnames(cluster_norm) <- metad$id

#Run clustering
clusters <- degPatterns(cluster_norm, metadata = metad3,time = 'genotype', col='treatment', minc = 10)
clusters$plot
clusters$normalized

#Enhancing cluster plots
degPlotCluster(clusters$normalized, time = 'treatment', col='genotype',
               min_genes = 10, smooth = F)+
  scale_color_manual(values=c("#8C8943","#3072AA"))+
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Significant LRT Clustering - Relative Gene Expression",
       subtitle = "Groups based on pairwise clustering of genes",
       caption = paste0("Data source: ", ngenes, " genes LRT adjusted p<0.05, clusters >10 genes"))
ggsave("osd457_LRT.pdf", width = 6, height = 6)


###DEG parameters###
pthreshold <- .05
log2fcthreshold <- 1.5
print(log2fcthreshold)
log2fcneg <- -1.5
print(log2fcneg)

###DEG list building###
contrasts
input1 <- "Wild Type_FLT"
input2 <- "Wild Type_GC"

###
ap_term <- glue("Adj.p.value_(",input1,")v(",input2,")")
log2fc_term <- glue("Log2fc_(",input1,")v(",input2,")")
ap_term <- as.symbol(ap_term)
log2fc_term <- as.symbol(log2fc_term)

degs <- output_table_1 %>%
  filter(!!ap_term <= pthreshold 
         & ( !!log2fc_term >= log2fcthreshold | !!log2fc_term <=log2fcneg))%>%
  dplyr::select(ENSEMBL, SYMBOL, all_of(ap_term), all_of(log2fc_term), ENTREZID)
degs_UP <-filter(degs, !!log2fc_term >= log2fcthreshold)
degs_UP <- list(degs_UP$ENSEMBL)
degs_DWN <-filter(degs, !!log2fc_term <= log2fcthreshold)
degs_DWN <- list(degs_DWN$ENSEMBL)
###
#check specific genes
degs%>%
  dplyr::filter(SYMBOL == "Apoa4")%>%
  dplyr::select(contains("Wild Type"))

###Ontology Gathering###

UP <- as.character(unlist(degs_UP))
DWN <- as.character(unlist(degs_DWN))

eGO_BP3_UP <- enrichGO(gene     = UP,
                       keyType = "ENSEMBL",
                       OrgDb    = org.Mm.eg.db,
                       ont      = "BP",
                       readable = TRUE)
eGO_BP3_UP
as.data.frame(eGO_BP3_UP)%>%
  dplyr::select(contains("ID")|contains  ("Description"))%>%
  head()

dotplot(eGO_BP3_UP)


eGO_BP3_DWN <- enrichGO(gene     = DWN,
                       keyType = "ENSEMBL",
                       OrgDb    = org.Mm.eg.db,
                       ont      = "BP",
                       readable = TRUE)
eGO_BP3_DWN
as.data.frame(eGO_BP3_DWN)%>%
  dplyr::select(contains("ID")|contains  ("Description"))%>%
  head()
dotplot(eGO_BP3_DWN)


###Generate single UP/DWN table of GO terms###    
eGO_table <- eGO_BP3_UP@result
eGO_table[["Direction"]] <- "+"
eGO_table <- head(arrange(eGO_table,+p.adjust),5)
eGO_table2 <- eGO_BP3_DWN@result
eGO_table2[["Direction"]] <- "-"
eGO_table2 <- head(arrange(eGO_table2,+p.adjust),5)
eGO_table <- rbind(eGO_table, eGO_table2)

###Prepare publication ready table of GO terms### Does not include gene names
n_vals <- paste0("Up regulated genes Log2FC >1.5. Down regulated genes Log2FC< -1.5.")
gt_tbl <- gt(eGO_table %>% 
               arrange(Direction))
gGO_BP3 <- 
  gt_tbl%>%
  tab_header(title = md( paste0("**Top GO Biological Processes - ",input1,"v",input2,"**")), subtitle = "Up and Down Regulated Processes: +/-")%>%
  tab_source_note(source_note = n_vals)%>%
  tab_source_note(source_note = "Source: R v4.1.2, clusterProfiler v4.2.3")%>%
  cols_align(align = c( "auto"), columns = everything())%>%
  fmt_scientific(columns = c('pvalue', 'p.adjust', 'qvalue'), rows = everything(),decimals = 2, exp_style = 'e')%>%
  cols_hide('geneID')

gGO_BP3



##Volcano Plotting
group1 <- "Wild Type_FLT"
group2 <- "Wild Type_GC"

x <- paste('Log2fc_(', group1,')v(',group2,')', sep = "")
y <- paste('Adj.p.value_(', group1,')v(',group2,')', sep = "")
title <- paste(group1, " v ", group2, sep = '')

ev <- EnhancedVolcano(output_table_1,
                      lab=output_table_1$SYMBOL, 
                      x=x, y=y,
                      pCutoff = 5e-2,
                      FCcutoff = 1.5,
                      xlim = c(-6,6),
                      ylim = c(0.0,10),
                      pointSize = 1.0,
                      drawConnectors = TRUE, widthConnectors = 0.75,
                      title = title
)
ev_out <- paste0(figure_output, group1, 'v', group2, '.png', sep='')
png(ev_out)
#dev.new(width = 20, height = 20)
plot(ev)
dev.off()

#Signaling Pathway Anlaysis 
genes <- degs
rownames(genes) <- degs$ENTREZID
head(genes)
de_genes <- genes[[log2fc_term]]
names(de_genes) <- rownames(genes)
all_genes <- output_table_1$ENTREZID

spia_result <- spia(de=de_genes, all=all_genes, organism="mmu", plots = TRUE)
plotP(spia_result)
print(spia_result)

colnames(spia_result) <- c("Name", "Pathway ID", "Pathway Size", "Number DE Genes", "Total Perturbation", "Probability of DE genes", "Probability of Perturbation", "Combined Probability p-value", "FDR", "Adj. p-value", "Status", "KEGGLink")
subset(spia_result, n=20)

gt_tbl <- gt(spia_result[1:5,1:10])
gt_tbl <- 
  gt_tbl%>%
  tab_header(title = md( "**Signaling Pathway Impact Analysis - BNL Liver Aim 2 HU_IR v NL**"), subtitle = "Top 10 pathways shown")%>%
  tab_source_note(source_note = "Source: R v4.2.0, SPIA Package v2.24.0")%>%
  tab_source_note(source_note = md("Ref. Romero R. et. al.,(2009) *Bioinformatics*, Oxford University Press."))%>%
  cols_align(align = c( "auto"), columns = everything())

gt_tbl

show(de_genes)

BiocManager::install('ggupset')
library(ggupset)
upsetplot(eGOhu)
gse_prep <- output_table_1%>%
              dplyr::select(ENSEMBL, all_of(log2fc_term))%>%
              arrange(desc({{log2fc_term}}))

gse_input <- setNames(gse_prep[[log2fc_term]], gse_prep$ENSEMBL)

degs
gse <- gseGO(gene     = gse_input,
             keyType = "ENSEMBL",
             OrgDb    = org.Mm.eg.db,
             ont      = "BP")

data.frame(gse)
g1 <- gseaplot(gse, geneSetID =1)
g2 <- gseaplot(gse, geneSetID =2)
cowplot::plot_grid(g1,g2)
gseaplot(gse, geneSetID =3)+gseaplot(gse, geneSetID =4)


gse.r <- setReadable(gse, 'org.Mm.eg.db', 'ENSEMBL')
cnetplot(gse.r, foldChange=gse_prep[[log2fc_term]], circular=T, colorEdge=T,
         showCategory = 2)
