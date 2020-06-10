# Understanding Gender Differential Susceptibility in Cancer


library(WGCNA)
library(tidyverse)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)



# load R workspace
load("./r_workspaces/tcga_gtex_thyroid_wgcna_networks.RData")



network <- tibble(gene = names(females_tumour_network$colors), module = unname(females_tumour_network$colors), color = labels2colors(females_tumour_network$colors))
network <- network %>%
  filter(color != "grey") %>%
  arrange(module)

cor_matrix <- females_tumour[, network$gene]
cor_matrix <- cor(cor_matrix)
cor_matrix <- abs(cor_matrix)

pdf(file="./plots/wgcna_networks/wgcna_thyroid_tumours_females_heatmap.pdf", w=6, h=5)

colors = structure(network$color, names = network$module)

rows <- rowAnnotation(
  modules = network$module,
  col = list(modules = colors),
  show_annotation_name = FALSE,
  show_legend = c("modules" = FALSE),
  border = TRUE)

columns <- columnAnnotation(
  modules = network$module,
  col = list(modules = colors),
  show_annotation_name = FALSE,
  show_legend = c("modules" = FALSE),
  border = TRUE)

ht <- Heatmap(
  matrix = cor_matrix,
  show_heatmap_legend = TRUE,
  name = "Absolute\nPearson's r",
  col = circlize::colorRamp2(c(0, 1), c("white", "red")),
  column_title = "Modules",
  row_title = "Modules",
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = rows,
  top_annotation = columns,
  border = T,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title_gp = gpar(fontsize = 12),
  row_title_gp = gpar(fontsize = 12))

print(ht)

dev.off()
