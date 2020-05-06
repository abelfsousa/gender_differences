# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Correlation between MEs and sample phenotypes


# -- Thyroid


library(tidyverse)
library(WGCNA)
library(data.table)
library(gplots)
library(ComplexHeatmap)
library(circlize)



# load R workspace
load("./r_workspaces/tcga_gtex_thyroid_wgcna_networks.RData")


# load networks comparison Male vs Female in TCGA tumours
color_code <- c("male-specific" = "orange", "female-specific" = "red", "lowly-conserved" = "#deebf7", "moderately-conserved" = "#9ecae1", "highly-conserved" = "#3182bd")
thyroid_tumour_gender_diff_networks <- read_tsv("./files/thyroid_tumour_gender_diff_networks.txt") %>%
  mutate(colors = unname(color_code[state]))


# load TCGA thyroid metadata
thyroid_tcga_clinical1 <- read.delim("./files/tcga_thca_meta.txt", row.names = c(1), stringsAsFactors=F) %>%
  rownames_to_column(var = "sample") %>%
  as.tibble() %>%
  filter(sample_type_id == 1) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  dplyr::select(sample, gender, histological_type, tumor_stage)



# load additional TCGA clinical data
thyroid_tcga_clinical2 <- fread("./data/tcga/tcga_clinical_data.txt") %>%
  as.tibble() %>%
  mutate(bcr_patient_barcode = str_replace_all(bcr_patient_barcode, "-", ".")) %>%
  filter(type == "THCA") %>%
  dplyr::select(sample = bcr_patient_barcode, gender, OS.time, histological_type, ajcc_pathologic_tumor_stage) %>%
  filter(histological_type != "Other, specify" & ajcc_pathologic_tumor_stage != "[Not Available]")



thyroid_tcga_meta <- inner_join(thyroid_tcga_clinical1, thyroid_tcga_clinical2, by = "sample")






tumourME_males <- moduleEigengenes(males_tumour, labels2colors(males_tumour_network$colors))$eigengenes
tumourME_males <- cbind(sample = rownames(males_tumour), tumourME_males) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  dplyr::select(-MEgrey) %>%
  inner_join(thyroid_tcga_clinical2, by = "sample") %>%
  dplyr::rename(ajcc_tumor_stage = ajcc_pathologic_tumor_stage) %>%
  mutate(ajcc_tumor_stage2 = as.numeric(as.factor(ajcc_tumor_stage)), histological_type2 = as.numeric(as.factor(histological_type)))
write.table(tumourME_males, "./files/thyroid_tumourME_traits_males.txt", sep="\t", quote=F, row.names = F)


#male modules (all highly/moderately-conserved in female modules)
tumourME_males_cor <- cor(
  tumourME_males %>% dplyr::select(OS.time, histological_type2, ajcc_tumor_stage2),
  tumourME_males %>% dplyr::select(starts_with("ME")),
  method = "pearson") %>%
  abs()


pdf(file="./plots/wgcna_networks_traits/thyroid_tumourME_males_cor.pdf", height=6, width=6)
#par(oma = c(3,0,0,5))
heatmap.2(
  x = tumourME_males_cor,
  trace="none",
  Rowv=TRUE,
  Colv=TRUE,
  dendrogram="both",
  cexRow=1.2,
  cexCol = 1.2,
  key=TRUE,
  keysize=2,
  symkey=TRUE,
  key.title="Pearson's r",
  key.xlab = NA,
  key.ylab = NA,
  density.info = "none",
  col=bluered(100),
  ColSideColors = thyroid_tumour_gender_diff_networks %>% filter(network == "males") %>% pull(colors),
  labRow = c("Overall survival", "Histological type", "AJCC tumor stage"))
dev.off()



tumourME_females <- moduleEigengenes(females_tumour, labels2colors(females_tumour_network$colors))$eigengenes
tumourME_females <- cbind(sample = rownames(females_tumour), tumourME_females) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  dplyr::select(-MEgrey) %>%
  inner_join(thyroid_tcga_clinical2, by = "sample") %>%
  dplyr::rename(ajcc_tumor_stage = ajcc_pathologic_tumor_stage) %>%
  mutate(ajcc_tumor_stage2 = as.numeric(as.factor(ajcc_tumor_stage)), histological_type2 = as.numeric(as.factor(histological_type)))
write.table(tumourME_females, "./files/thyroid_tumourME_traits_females.txt", sep="\t", quote=F, row.names = F)


#female-specific tumour modules
tumourME_females_cor <- cor(
  tumourME_females %>% dplyr::select(OS.time, histological_type2, ajcc_tumor_stage2),
  tumourME_females %>% dplyr::select(thyroid_tumour_gender_diff_networks %>% filter(state == "lowly-conserved") %>% pull(module) %>% paste("ME", ., sep="")),
  method = "pearson") %>%
  abs()

tumourME_females_y <- tumourME_females %>%
  dplyr::select(sample, MEcyan, MEgreenyellow) %>%
  pivot_longer(-sample, names_to = "ME", values_to = "ME_value")

tumourME_females_x <- tumourME_females %>%
  dplyr::select(sample, OS.time, histological_type2, ajcc_tumor_stage2) %>%
  pivot_longer(-sample, names_to = "feature", values_to = "feature_value")

linear_model <- function(mat){
  model <- lm(ME_value ~ feature_value, data = mat)

  rsq <- summary(model)$r.squared
  anovaP <- broom::tidy(aov(model)) %>%
    filter(term == "feature_value") %>%
    pull(p.value)

  res <- tibble(rsquared = rsq, pvalue = anovaP)

}

tumourME_females_mod <- inner_join(tumourME_females_y, tumourME_females_x, by = "sample")
tumourME_females_mod <- tumourME_females_mod %>%
  group_by(ME, feature) %>%
  nest() %>%
  ungroup() %>%
  mutate(data = map2(feature, data, ~ if(.x != "OS.time"){mutate_at(.y, 3, as.character)}else{.y})) %>%
  mutate(model = map(.x = data, .f = linear_model)) %>%
  select(-data) %>%
  unnest(model)

tumourME_females_pval <- tumourME_females_mod %>%
  select(-rsquared) %>%
  mutate(pvalue = -log10(pvalue)) %>%
  pivot_wider(names_from = "feature", values_from = "pvalue") %>%
  column_to_rownames(var = "ME") %>%
  t()

tumourME_females_rsq <- tumourME_females_mod %>%
  select(-pvalue) %>%
  mutate(rsquared = round(rsquared, 2)) %>%
  pivot_wider(names_from = "feature", values_from = "rsquared") %>%
  column_to_rownames(var = "ME") %>%
  t()


pdf(file="./plots/wgcna_networks_traits/thyroid_tumourME_females_cor.pdf", height=1.5, width=3.5)
mat <- tumourME_females_pval
mat2 <- tumourME_females_rsq
column_labels = structure(c("M1", "M2"), names = colnames(mat))
row_labels = structure(c("Survival (days)", "Histological type", "Tumour stage"), names = rownames(mat))

#col_fun = colorRamp2(c(0, 10, 18), brewer.pal(3, "Reds"))
col_fun = colorRamp2(c(0, 18), c("white", "red"))
print(Heatmap(
  matrix = mat,
  name = "P-value\n(-log10)",
  column_title = "",
  row_title = "",
  column_labels = column_labels[colnames(mat)],
  row_labels = row_labels[rownames(mat)],
  col = col_fun,
  heatmap_legend_param = list(at = c(0, 4, 8, 12, 18)),
  border = T,
  cluster_rows = T,
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 0,
  column_names_centered = TRUE,
  rect_gp = gpar(col = "#d9d9d9", lwd = 2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(mat2[i, j], x, y, gp = gpar(fontsize = 10))}))
dev.off()








# melted datasets
males_tumour_network_MEs <- cbind(sample = rownames(males_tumour), males_tumour_network$MEs) %>%
  as.tibble() %>%
  gather(key = "moduleN", value = "ME", -sample) %>%
  mutate(moduleN = as.numeric(str_replace_all(moduleN, "[ME]", "")), sample = str_sub(sample, 1, 12)) %>%
  mutate(moduleL = labels2colors(moduleN)) %>%
  filter(moduleL != "grey") %>%
  inner_join(thyroid_tumour_gender_diff_networks %>% filter(network == "males"), by=c("moduleL" = "module")) %>%
  arrange(moduleL) %>%
  dplyr::select(network, moduleL, moduleN, state, sample, ME) %>%
  inner_join(thyroid_tcga_clinical2, by=c("sample")) %>%
  mutate(ajcc_pathologic_tumor_stage = as.numeric(as.factor(ajcc_pathologic_tumor_stage)), histological_type = as.numeric(as.factor(histological_type)))
write.table(males_tumour_network_MEs, "./files/thyroid_males_tumour_network_MEs_traits.txt", quote=F, sep="\t", row.names=F)


females_tumour_network_MEs <- cbind(sample = rownames(females_tumour), females_tumour_network$MEs) %>%
  as.tibble() %>%
  gather(key = "moduleN", value = "ME", -sample) %>%
  mutate(moduleN = as.numeric(str_replace_all(moduleN, "[ME]", "")), sample = str_sub(sample, 1, 12)) %>%
  mutate(moduleL = labels2colors(moduleN)) %>%
  filter(moduleL != "grey") %>%
  inner_join(thyroid_tumour_gender_diff_networks %>% filter(network == "females"), by=c("moduleL" = "module")) %>%
  arrange(moduleL) %>%
  dplyr::select(network, moduleL, moduleN, state, sample, ME) %>%
  inner_join(thyroid_tcga_clinical2, by=c("sample")) %>%
  mutate(ajcc_pathologic_tumor_stage = as.numeric(as.factor(ajcc_pathologic_tumor_stage)), histological_type = as.numeric(as.factor(histological_type)))
write.table(females_tumour_network_MEs, "./files/thyroid_females_tumour_network_MEs_traits.txt", quote=F, sep="\t", row.names=F)


save(list=ls(), file="r_workspaces/thyroid_networks_MEs_traits.RData")
