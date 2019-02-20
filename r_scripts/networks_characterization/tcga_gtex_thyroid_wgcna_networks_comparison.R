# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data




# -- Thyroid

# load R workspace
load("./r_workspaces/tcga_gtex_thyroid_wgcna_networks.RData")

source("./r_scripts/utils.R")



library(WGCNA)
library(tidyverse)
library(gplots)





#Males tumour vs females tumour
males_females_tumour_overlapTable <- overlap_table_wgcna(
    males_tumour_network,
    females_tumour_network,
    "Thyroid Tumours (TCGA) - Men",
    "Thyroid Tumours (TCGA) - Women",
    "./plots/wgcna_networks/tcga_thyroid_tumours_wgcna_networks_comparison_malesVSfemales.pdf", 20, 10)

males_females_tumour_comparison <- networks_comparison(
    males_tumour,
    females_tumour,
    males_tumour_network,
    females_tumour_network,
    males_females_tumour_overlapTable)
write.table(males_females_tumour_comparison$corTable1, "./files/thyroid_males_females_tumour_comparison_corTable1.txt", quote=F, sep="\t", row.names=T)
write.table(males_females_tumour_comparison$corTable2, "./files/thyroid_males_females_tumour_comparison_corTable2.txt", quote=F, sep="\t", row.names=T)


pdf(file="./plots/wgcna_networks/tcga_thyroid_tumours_wgcna_networks_comparison_malesVSfemales2.pdf", w=6, h=5)
heatmap.2( x=as.matrix(males_females_tumour_comparison$corTable1),
    Rowv=T,
    Colv=T,
    dendrogram="both",
    scale = "none",
    trace = "none",
    key = FALSE,
    margins=c(8,8),
    col = c("grey", "#1c9099", "#deebf7", "#9ecae1", "#3182bd"),
    main = "",
    xlab = "Women modules",
    ylab = "Men modules",
    labRow = paste(rownames(as.matrix(males_females_tumour_comparison$corTable1)), table(labels2colors(males_tumour_network$colors))[names(table(labels2colors(males_tumour_network$colors))) != "grey"]),
    labCol = paste(colnames(as.matrix(males_females_tumour_comparison$corTable1)), table(labels2colors(females_tumour_network$colors))[names(table(labels2colors(females_tumour_network$colors))) != "grey"]))
legend("topleft", fill = c("#1c9099", "#deebf7", "#9ecae1", "#3182bd"),
  legend = c("Network-specific module", "Module pair lowly-conserved", "Module pair moderately-conserved", "Module pair highly-conserved"), y.intersp = 1.5, cex=0.6 )
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

thyroid_tumour_gender_diff_networks <- bind_rows(tumour_male_modules, tumour_female_modules)
write.table(thyroid_tumour_gender_diff_networks, "./files/thyroid_tumour_gender_diff_networks.txt", quote=F, sep="\t", row.names=F)





#Males normal vs females normal
males_females_normal_overlapTable <- overlap_table_wgcna(
    males_normal_network,
    females_normal_network,
    "Thyroid Normal (GTEx) - Men",
    "Thyroid Normal (GTEx) - Women",
    "./plots/wgcna_networks/gtex_thyroid_wgcna_networks_comparison_malesVSfemales.pdf", 20, 20)

males_females_normal_comparison <- networks_comparison(
    males_normal,
    females_normal,
    males_normal_network,
    females_normal_network,
    males_females_normal_overlapTable)
write.table(males_females_normal_comparison$corTable1, "./files/thyroid_males_females_normal_comparison_corTable1.txt", quote=F, sep="\t", row.names=T)
write.table(males_females_normal_comparison$corTable2, "./files/thyroid_males_females_normal_comparison_corTable2.txt", quote=F, sep="\t", row.names=T)


pdf(file="./plots/wgcna_networks/gtex_thyroid_wgcna_networks_comparison_malesVSfemales2.pdf", w=8, h=8)
heatmap.2( x=as.matrix(males_females_normal_comparison$corTable1),
    Rowv=T,
    Colv=T,
    dendrogram="both",
    scale = "none",
    trace = "none",
    key = FALSE,
    margins=c(8,8),
    col = c("grey", "#1c9099", "#deebf7", "#9ecae1", "#3182bd"),
    main = "",
    xlab = "Women modules",
    ylab = "Men modules",
    labRow = paste(rownames(as.matrix(males_females_normal_comparison$corTable1)), table(labels2colors(males_normal_network$colors))[names(table(labels2colors(males_normal_network$colors))) != "grey"]),
    labCol = paste(colnames(as.matrix(males_females_normal_comparison$corTable1)), table(labels2colors(females_normal_network$colors))[names(table(labels2colors(females_normal_network$colors))) != "grey"]))
legend("topleft", fill = c("#1c9099", "#deebf7", "#9ecae1", "#3182bd"),
  legend = c("Network-specific module", "Module pair lowly-conserved", "Module pair moderately-conserved", "Module pair highly-conserved"), y.intersp = 1.5, cex=0.6 )
dev.off()

normal_male_modules <- tibble(network = "males",
  module = names(rowSums(males_females_normal_comparison$corTable1)),
  state = apply(males_females_normal_comparison$corTable1, 1, max)) %>%
  mutate(state = male_code[state])

normal_female_modules <- tibble(network = "females",
  module = names(colSums(males_females_normal_comparison$corTable1)),
  state = apply(males_females_normal_comparison$corTable1, 2, max)) %>%
  mutate(state = female_code[state])

thyroid_normal_gender_diff_networks <- bind_rows(normal_male_modules, normal_female_modules)
write.table(thyroid_normal_gender_diff_networks, "./files/thyroid_normal_gender_diff_networks.txt", quote=F, sep="\t", row.names=F)
