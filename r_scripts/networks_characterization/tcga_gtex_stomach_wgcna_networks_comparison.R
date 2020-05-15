# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data




# -- Stomach

# load R workspace
load("./r_workspaces/tcga_gtex_stomach_wgcna_networks.RData")

source("./r_scripts/utils.R")



library(WGCNA)
library(tidyverse)
library(gplots)
library(ComplexHeatmap)





#Males tumour vs females tumour
males_females_tumour_overlapTable <- overlap_table_wgcna(
    males_tumour_network,
    females_tumour_network,
    "Stomach Tumours (TCGA) - Men",
    "Stomach Tumours (TCGA) - Women",
    "./plots/wgcna_networks/tcga_stomach_tumours_wgcna_networks_comparison_malesVSfemales.pdf", 20, 20)

males_females_tumour_comparison <- networks_comparison(
    males_tumour,
    females_tumour,
    males_tumour_network,
    females_tumour_network,
    males_females_tumour_overlapTable)
write.table(males_females_tumour_comparison$corTable1, "./files/stomach_males_females_tumour_comparison_corTable1.txt", quote=F, sep="\t", row.names=T)
write.table(males_females_tumour_comparison$corTable2, "./files/stomach_males_females_tumour_comparison_corTable2.txt", quote=F, sep="\t", row.names=T)


pdf(file="./plots/wgcna_networks/tcga_stomach_tumours_wgcna_networks_comparison_malesVSfemales2.pdf", w=6, h=6)

mat <- as.matrix(males_females_tumour_comparison$corTable1)
mat2 <- mat
mat2[which(mat2 == 1)] = 0
mat2[which(mat2 == 2)] = 1
mat2[which(mat2 == 3)] = 2
mat2[which(mat2 == 4)] = 3

colors = structure(c("#bdbdbd", "#fee8c8", "#fdbb84", "#e34a33"), names = c("0", "1", "2", "3"))
column_labels = structure(
  table(labels2colors(females_tumour_network$colors))[names(table(labels2colors(females_tumour_network$colors))) != "grey"],
  names = colnames(mat2))
rown_labels = structure(
  table(labels2colors(males_tumour_network$colors))[names(table(labels2colors(males_tumour_network$colors))) != "grey"],
  names = rownames(mat2))

code <- c("0" = "gender-specific", "1" = "lowly", "2" = "moderately", "3" = "highly")
ha1 <- as.character(code[as.character(apply(mat2,1,max))])
ha1 <- rowAnnotation(
  Preservation = ha1,
  col = list(Preservation = c("gender-specific" = "#31a354", "lowly" = "#9ecae1", "moderately" = "#4292c6", "highly" = "#08306b")),
  show_annotation_name = FALSE,
  show_legend = c("Preservation" = FALSE),
  border = TRUE)

ha2 <- as.character(code[as.character(apply(mat2,2,max))])
ha2 <- columnAnnotation(
  Preservation = ha2,
  col = list(Preservation = c("gender-specific" = "#31a354", "lowly" = "#9ecae1", "moderately" = "#4292c6", "highly" = "#08306b")),
  show_annotation_name = FALSE,
  show_legend = c("Preservation" = FALSE),
  border = TRUE)

lgd1 <- Legend(at = 0:3, title = "Overlap", legend_gp = gpar(fill = c("#bdbdbd", "#fee8c8", "#fdbb84", "#e34a33")))
lgd2 <- Legend(labels = c("gender-specific", "lowly", "moderately", "highly"), title = "Preservation", legend_gp = gpar(fill = c("#31a354", "#9ecae1", "#4292c6", "#08306b")))
pd <- packLegend(lgd1, lgd2)

ht <- Heatmap(
  matrix = mat2,
  show_heatmap_legend = FALSE,
  name = "Overlap",
  column_title = "Female modules",
  row_title = "Male modules",
  column_labels = column_labels[colnames(mat2)],
  row_labels = rown_labels[rownames(mat2)],
  left_annotation = ha1,
  top_annotation = ha2,
  col = colors,
  border = T,
  cluster_rows = T,
  cluster_columns = T,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 12),
  row_title_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
  if(mat[i, j] == 1)
    grid.text("*", x, y, gp = gpar(fontsize = 10))})
draw(ht,annotation_legend_list = pd)

dev.off()

male_code <- c("1" = "male-specific", "2" = "lowly-conserved", "3" = "moderately-conserved", "4" = "highly-conserved")
female_code <- c("1" = "female-specific", "2" = "lowly-conserved", "3" = "moderately-conserved", "4" = "highly-conserved")

tumour_male_modules <- tibble(network = "males",
  module = names(rowSums(males_females_tumour_comparison$corTable1)),
  state = apply(males_females_tumour_comparison$corTable1, 1, max)) %>%
  mutate(state = male_code[state])

tumour_female_modules <- tibble(network = "females",
  module = names(colSums(males_females_tumour_comparison$corTable1)),
  state = apply(males_females_tumour_comparison$corTable1, 2, max)) %>%
  mutate(state = female_code[state])

stomach_tumour_gender_diff_networks <- bind_rows(tumour_male_modules, tumour_female_modules)
write.table(stomach_tumour_gender_diff_networks, "./files/stomach_tumour_gender_diff_networks.txt", quote=F, sep="\t", row.names=F)


# topology analysis
get_values <- function(ind1, ind2, mat, flag){
  values = c()
  if(flag == 1){
    for(i in 1:length(ind2)){
      values = c(values, mat[ind2[i], ind1[i]])
    }
  }
  else if(flag == 2){
    for(i in 1:length(ind2)){
      values = c(values, mat[ind1[i], ind2[i]])
    }
  }
  else(stop("unknown flag"))
  return(values)
}

tumour_male_modules2 <- tumour_male_modules %>%
  mutate(indexes1 = apply(males_females_tumour_comparison$corTable1, 1, function(x) which(x == max(x)))) %>%
  unnest() %>%
  mutate(indexes2 = as.numeric(as.factor(module))) %>%
  mutate(cor.kME = get_values(indexes1, indexes2, males_females_tumour_comparison$corTable2, 1)) %>%
  group_by(network, module, state) %>%
  summarise(cor.kME = max(cor.kME)) %>%
  ungroup() %>%
  mutate(topology = if_else(cor.kME > 0.7, "preserved", "not_preserved")) %>%
  mutate(topology = replace(topology, (topology == "not_preserved" & state == "male-specific"), NA))


tumour_female_modules2 <- tumour_female_modules %>%
  mutate(indexes1 = apply(males_females_tumour_comparison$corTable1, 2, function(x) which(x == max(x)))) %>%
  unnest() %>%
  mutate(indexes2 = as.numeric(as.factor(module))) %>%
  mutate(cor.kME = get_values(indexes1, indexes2, males_females_tumour_comparison$corTable2, 2)) %>%
  group_by(network, module, state) %>%
  summarise(cor.kME = max(cor.kME)) %>%
  ungroup() %>%
  mutate(topology = if_else(cor.kME > 0.7, "preserved", "not_preserved")) %>%
  mutate(topology = replace(topology, (topology == "not_preserved" & state == "female-specific"), NA))

stomach_tumour_gender_diff_networks_topology <- bind_rows(tumour_male_modules2, tumour_female_modules2)
write.table(stomach_tumour_gender_diff_networks_topology, "./files/stomach_tumour_gender_diff_networks_topology.txt", quote=F, sep="\t", row.names=F)





#Males normal vs females normal
males_females_normal_overlapTable <- overlap_table_wgcna(
    males_normal_network,
    females_normal_network,
    "Stomach Normal (GTEx) - Men",
    "Stomach Normal (GTEx) - Women",
    "./plots/wgcna_networks/gtex_stomach_wgcna_networks_comparison_malesVSfemales.pdf", 12, 12)

males_females_normal_comparison <- networks_comparison(
    males_normal,
    females_normal,
    males_normal_network,
    females_normal_network,
    males_females_normal_overlapTable)
write.table(males_females_normal_comparison$corTable1, "./files/stomach_males_females_normal_comparison_corTable1.txt", quote=F, sep="\t", row.names=T)
write.table(males_females_normal_comparison$corTable2, "./files/stomach_males_females_normal_comparison_corTable2.txt", quote=F, sep="\t", row.names=T)


pdf(file="./plots/wgcna_networks/gtex_stomach_wgcna_networks_comparison_malesVSfemales2.pdf", w=6, h=6)

mat <- as.matrix(males_females_normal_comparison$corTable1)
mat2 <- mat
mat2[which(mat2 == 1)] = 0
mat2[which(mat2 == 2)] = 1
mat2[which(mat2 == 3)] = 2
mat2[which(mat2 == 4)] = 3

colors = structure(c("#bdbdbd", "#fee8c8", "#fdbb84", "#e34a33"), names = c("0", "1", "2", "3"))
column_labels = structure(
  table(labels2colors(females_normal_network$colors))[names(table(labels2colors(females_normal_network$colors))) != "grey"],
  names = colnames(mat2))
rown_labels = structure(
  table(labels2colors(males_normal_network$colors))[names(table(labels2colors(males_normal_network$colors))) != "grey"],
  names = rownames(mat2))

code <- c("0" = "gender-specific", "1" = "lowly", "2" = "moderately", "3" = "highly")
ha1 <- as.character(code[as.character(apply(mat2,1,max))])
ha1 <- rowAnnotation(
  Preservation = ha1,
  col = list(Preservation = c("gender-specific" = "#31a354", "lowly" = "#9ecae1", "moderately" = "#4292c6", "highly" = "#08306b")),
  show_annotation_name = FALSE,
  show_legend = c("Preservation" = FALSE),
  border = TRUE)

ha2 <- as.character(code[as.character(apply(mat2,2,max))])
ha2 <- columnAnnotation(
  Preservation = ha2,
  col = list(Preservation = c("gender-specific" = "#31a354", "lowly" = "#9ecae1", "moderately" = "#4292c6", "highly" = "#08306b")),
  show_annotation_name = FALSE,
  show_legend = c("Preservation" = FALSE),
  border = TRUE)

lgd1 <- Legend(at = 0:3, title = "Overlap", legend_gp = gpar(fill = c("#bdbdbd", "#fee8c8", "#fdbb84", "#e34a33")))
lgd2 <- Legend(labels = c("gender-specific", "lowly", "moderately", "highly"), title = "Preservation", legend_gp = gpar(fill = c("#31a354", "#9ecae1", "#4292c6", "#08306b")))
pd <- packLegend(lgd1, lgd2)

ht <- Heatmap(
  matrix = mat2,
  show_heatmap_legend = FALSE,
  name = "Overlap",
  column_title = "Female modules",
  row_title = "Male modules",
  column_labels = column_labels[colnames(mat2)],
  row_labels = rown_labels[rownames(mat2)],
  left_annotation = ha1,
  top_annotation = ha2,
  col = colors,
  border = T,
  cluster_rows = T,
  cluster_columns = T,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  column_title_gp = gpar(fontsize = 12),
  row_title_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
  if(mat[i, j] == 1)
    grid.text("*", x, y, gp = gpar(fontsize = 10))})
draw(ht,annotation_legend_list = pd)

dev.off()

normal_male_modules <- tibble(network = "males",
  module = names(rowSums(males_females_normal_comparison$corTable1)),
  state = apply(males_females_normal_comparison$corTable1, 1, max)) %>%
  mutate(state = male_code[state])

normal_female_modules <- tibble(network = "females",
  module = names(colSums(males_females_normal_comparison$corTable1)),
  state = apply(males_females_normal_comparison$corTable1, 2, max)) %>%
  mutate(state = female_code[state])

stomach_normal_gender_diff_networks <- bind_rows(normal_male_modules, normal_female_modules)
write.table(stomach_normal_gender_diff_networks, "./files/stomach_normal_gender_diff_networks.txt", quote=F, sep="\t", row.names=F)



# topology analysis
normal_male_modules2 <- normal_male_modules %>%
  mutate(indexes1 = apply(males_females_normal_comparison$corTable1, 1, function(x) which(x == max(x)))) %>%
  unnest() %>%
  mutate(indexes2 = as.numeric(as.factor(module))) %>%
  mutate(cor.kME = get_values(indexes1, indexes2, males_females_normal_comparison$corTable2, 1)) %>%
  group_by(network, module, state) %>%
  summarise(cor.kME = max(cor.kME)) %>%
  ungroup() %>%
  mutate(topology = if_else(cor.kME > 0.7, "preserved", "not_preserved")) %>%
  mutate(topology = replace(topology, (topology == "not_preserved" & state == "male-specific"), NA))


normal_female_modules2 <- normal_female_modules %>%
  mutate(indexes1 = apply(males_females_normal_comparison$corTable1, 2, function(x) which(x == max(x)))) %>%
  unnest() %>%
  mutate(indexes2 = as.numeric(as.factor(module))) %>%
  mutate(cor.kME = get_values(indexes1, indexes2, males_females_normal_comparison$corTable2, 2)) %>%
  group_by(network, module, state) %>%
  summarise(cor.kME = max(cor.kME)) %>%
  ungroup() %>%
  mutate(topology = if_else(cor.kME > 0.7, "preserved", "not_preserved")) %>%
  mutate(topology = replace(topology, (topology == "not_preserved" & state == "female-specific"), NA))

stomach_normal_gender_diff_networks_topology <- bind_rows(normal_male_modules2, normal_female_modules2)
write.table(stomach_normal_gender_diff_networks_topology, "./files/stomach_normal_gender_diff_networks_topology.txt", quote=F, sep="\t", row.names=F)
