run_analysis <- function(input1, input2, output_table_1, pthreshold, log2fcthreshold, log2fcneg, figure_output) {
  # Dynamic column names
  ap_term <- glue("Adj.p.value_({input1})v({input2})")
  log2fc_term <- glue("Log2fc_({input1})v({input2})")
  ap_term <- as.symbol(ap_term)
  log2fc_term <- as.symbol(log2fc_term)
  
  # Filter DEG list
  degs <- output_table_1 %>%
    filter(
      !!ap_term <= pthreshold &
        ( !!log2fc_term >= log2fcthreshold | !!log2fc_term <= log2fcneg)
    ) %>%
    dplyr::select(ENSEMBL, SYMBOL, all_of(ap_term), all_of(log2fc_term), ENTREZID)
  
  # Up and Down-regulated Genes
  degs_UP <- list(filter(degs, !!log2fc_term >= log2fcthreshold)$ENSEMBL)
  degs_DWN <- list(filter(degs, !!log2fc_term <= log2fcneg)$ENSEMBL)

  # Enrichment Analysis
  UP <- as.character(unlist(degs_UP))
  DWN <- as.character(unlist(degs_DWN))
  
  eGO_BP3_UP <- enrichGO(
    gene     = UP,
    keyType  = "ENSEMBL",
    OrgDb    = org.Mm.eg.db,
    ont      = "BP",
    readable = TRUE
  )
  
  eGO_BP3_DWN <- enrichGO(
    gene     = DWN,
    keyType  = "ENSEMBL",
    OrgDb    = org.Mm.eg.db,
    ont      = "BP",
    readable = TRUE
  )
  
  # Save Dotplots
  dotplot(eGO_BP3_UP)
  ggsave(filename = file.path(figure_output, paste0(input1, "_vs_", input2, "_UP_dotplot.png")))
  
  dotplot(eGO_BP3_DWN)
  ggsave(filename = file.path(figure_output, paste0(input1, "_vs_", input2, "_DWN_dotplot.png")))
  
  # Combine UP and DOWN GO results
  eGO_table <- eGO_BP3_UP@result
  eGO_table[["Direction"]] <- "+"
  eGO_table <- head(arrange(eGO_table, +p.adjust), 5)
  eGO_table2 <- eGO_BP3_DWN@result
  eGO_table2[["Direction"]] <- "-"
  eGO_table2 <- head(arrange(eGO_table2, +p.adjust), 5)
  eGO_table <- rbind(eGO_table, eGO_table2)
  
  # Prepare table for publication-ready format
  n_vals <- paste0("Up regulated genes Log2FC >1.5. Down regulated genes Log2FC< -1.5.")
  gt_tbl <- gt(eGO_table %>% arrange(Direction))
  gGO_BP3 <- gt_tbl %>%
    tab_header(
      title = md(paste0("**Top GO Biological Processes - ", input1, " v ", input2, "**")),
      subtitle = "Up and Down Regulated Processes: +/-"
    ) %>%
    tab_source_note(source_note = n_vals) %>%
    tab_source_note(source_note = "Source: R v4.1.2, clusterProfiler v4.2.3") %>%
    cols_align(align = c("auto"), columns = everything()) %>%
    fmt_scientific(columns = c('pvalue', 'p.adjust', 'qvalue'), rows = everything(), decimals = 2, exp_style = 'e') %>%
    cols_hide('geneID')
  
  # Save the GO table as PNG
  go_table_file <- file.path(figure_output, paste0(input1, "_vs_", input2, "_GO_table.png"))
  gtsave(gGO_BP3, filename = go_table_file, path = figure_output)
  
  # Volcano Plot
  x <- paste('Log2fc_(', input1, ')v(', input2, ')', sep = "")
  y <- paste('Adj.p.value_(', input1, ')v(', input2, ')', sep = "")
  title <- paste(input1, " v ", input2, sep = '')
  
  ev <- EnhancedVolcano(
    output_table_1,
    lab = output_table_1$SYMBOL,
    x = x, y = y,
    pCutoff = 5e-2,
    FCcutoff = 1.5,
    xlim = c(-6, 6),
    ylim = c(0.0, 10),
    pointSize = 1.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    title = title
  )
  
  ev_out <- file.path(figure_output, paste0(input1, "_vs_", input2, ".png"))
  png(ev_out)
  plot(ev)
  dev.off()
  
  # SPIA Analysis with Error Handling
#  de_genes <- degs[[log2fc_term]]
#  names(de_genes) <- degs$ENTREZID
#  all_genes <- unique(output_table_1$ENTREZID)
#  spia_result <- spia(de = de_genes, all = all_genes, organism = "mmu", plots = FALSE)
#  
#  # Attempt to Plot SPIA Results
#  spia_plot_path <- file.path(figure_output, paste0(input1, "_vs_", input2, "_SPIA.png"))
#  tryCatch({
#    png(spia_plot_path)
#    plotP(spia_result)
#    dev.off()
#  }, error = function(e) {
#    message(glue("SPIA plotP() failed for {input1} vs {input2}: {e$message}"))
#  })
  
  return(list(
    degs = degs,
    eGO_BP3_UP = eGO_BP3_UP,
    eGO_BP3_DWN = eGO_BP3_DWN,
    spia_result = spia_result
  ))
}

comparisons <- list(
  c("Wild Type_FLT", "Wild Type_GC"),
  c("Nrf2KO_FLT", "Nrf2KO_GC"),
  c("Nrf2KO_GC", "Wild Type_GC"),
  c("Nrf2KO_FLT", "Wild Type_FLT")
  # Add more comparisons as needed
)

# Create a progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta [:current/:total]",
  total = length(comparisons),
  clear = FALSE,
  width = 60
)

output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)

results <- lapply(comparisons, function(groups) {
  input1 <- groups[1]
  input2 <- groups[2]
  
  # Update progress bar
  pb$tick()
  
  run_analysis(
    input1 = input1,
    input2 = input2,
    output_table_1 = output_table_1,
    pthreshold = 0.05,
    log2fcthreshold = 1.5,
    log2fcneg = -1.5,
    figure_output = output_dir
  )
})
