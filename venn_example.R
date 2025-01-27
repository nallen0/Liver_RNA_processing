##Monocyte Venn diagram
shiftAPP_degs <- cv1 %>%
  filter(abs(Monocytes_avg_log2FC) > 1.2 & Monocytes_p_val_adj < 0.05) %>%  # Adjust column names accordingly
  pull(Gene)  # Count the number of rows (DEGs)
dayAPP_degs <- cv2 %>%
  filter(abs(Monocytes_avg_log2FC) > 1.2 & Monocytes_p_val_adj < 0.05) %>%  # Adjust column names accordingly
  pull(Gene)  # Count the number of rows (DEGs)
shiftLM_degs <- cv3 %>%
  filter(abs(Monocytes_avg_log2FC) > 1.2 & Monocytes_p_val_adj < 0.05) %>%  # Adjust column names accordingly
  pull(Gene)  # Count the number of rows (DEGs)

bind_cols_fill <- function(df_list) {
  
  max_rows <- purrr::map_int(df_list, nrow) %>% max()
  
  purrr::map(df_list, function(df) {
    if(nrow(df) == max_rows) return(df)
    first <- names(df)[1] %>% sym()
    df %>% add_row(!!first := rep(NA, max_rows - nrow(df)))
  }) %>% bind_cols()
}
monocyte_degs_list <- list("shiftAPP_degs", "dayAPP_degs", "shiftLM_degs")
monocyte_degs <- bind_cols_fill(list(shiftAPP = as.data.frame(shiftAPP_degs),
                                     dayAPP = as.data.frame(dayAPP_degs),
                                     shiftLM = as.data.frame(shiftLM_degs)))
library(ggVennDiagram)
ggVennDiagram(monocyte_degs,set_size = 3, label = "count")+
  scale_x_continuous(expand = expansion(mult = .2)) + scale_fill_distiller(palette = "Blues", direction = 1)
ggsave("monocyte_venn.pdf")