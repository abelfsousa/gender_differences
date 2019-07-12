# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Correlation between MEs and sample phenotypes


# -- Stomach


library(tidyverse)
library(WGCNA)
library(data.table)
library(gplots)



# load R workspace
load("./r_workspaces/tcga_gtex_stomach_wgcna_networks.RData")


# load networks comparison Male vs Female in TCGA tumours
color_code <- c("male-specific" = "orange", "female-specific" = "red", "lowly-conserved" = "#deebf7", "moderately-conserved" = "#9ecae1", "highly-conserved" = "#3182bd")
stomach_tumour_gender_diff_networks <- read_tsv("./files/stomach_tumour_gender_diff_networks.txt") %>%
  mutate(colors = unname(color_code[state]))


# load TCGA stomach metadata
stomach_tcga_clinical1 <- read.delim("./files/tcga_stad_meta.txt", row.names = c(1), stringsAsFactors=F) %>%
  rownames_to_column(var = "sample") %>%
  as.tibble() %>%
  filter(sample_type_id == 1) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  dplyr::select(sample, gender, Lauren_Class, tumor_stage)



# load additional TCGA clinical data
stomach_tcga_clinical2 <- fread("./data/tcga/tcga_clinical_data.txt") %>%
  as.tibble() %>%
  mutate(bcr_patient_barcode = str_replace_all(bcr_patient_barcode, "-", ".")) %>%
  filter(type == "STAD") %>%
  dplyr::select(sample = bcr_patient_barcode, gender, OS.time, histological_type, ajcc_pathologic_tumor_stage, histological_grade) %>%
  filter(!(histological_type %in% c("[Discrepancy]", "[Not Available]", "Stomach, Adenocarcinoma, Not Otherwise Specified (NOS)")) & !(ajcc_pathologic_tumor_stage %in% c("[Discrepancy]", "[Not Available]"))) %>%
  filter(!is.na(OS.time))



stomach_tcga_meta <- inner_join(stomach_tcga_clinical1, stomach_tcga_clinical2, by = "sample")






tumourME_males <- moduleEigengenes(males_tumour, labels2colors(males_tumour_network$colors))$eigengenes
tumourME_males <- cbind(sample = rownames(males_tumour), tumourME_males) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  dplyr::select(-MEgrey) %>%
  inner_join(stomach_tcga_clinical2, by = "sample") %>%
  dplyr::rename(ajcc_tumor_stage = ajcc_pathologic_tumor_stage) %>%
  mutate(ajcc_tumor_stage2 = as.numeric(as.factor(ajcc_tumor_stage)), histological_type2 = as.numeric(as.factor(histological_type)), histological_grade2 = as.numeric(as.factor(histological_grade)))
write.table(tumourME_males, "./files/stomach_tumourME_traits_males.txt", sep="\t", quote=F, row.names = F)


#male-specific tumour modules
tumourME_males_cor <- cor(
  tumourME_males %>% dplyr::select(OS.time, histological_type2, ajcc_tumor_stage2, histological_grade2),
  tumourME_males %>% dplyr::select(stomach_tumour_gender_diff_networks %>% filter(state == "male-specific") %>% pull(module) %>% paste("ME", ., sep="")),
  method = "pearson") %>%
  abs()


pdf(file="./plots/wgcna_networks_traits/stomach_tumourME_males_cor.pdf", height=3, width=5)
par(oma = c(0,0,0,5))
heatmap.2(
  x = tumourME_males_cor,
  trace="none",
  Rowv=TRUE,
  Colv=TRUE,
  dendrogram="both",
  cexRow=1.3,
  cexCol = 1.6,
  key=TRUE,
  keysize=2,
  symkey=FALSE,
  key.title="Pearson's r",
  key.xlab = NA,
  key.ylab = NA,
  key.par = list(cex.axis = 1),
  #srtCol=0,
  density.info = "none",
  #col=bluered(100),
  col = colorpanel(100, "white", "orange", "red"),
  #ColSideColors = stomach_tumour_gender_diff_networks %>% filter(network == "males") %>% pull(colors),
  #labCol = tumourME_males_cor %>% colnames %>% str_replace("ME", ""),
  labCol = c("M1", "M3", "M2"),
  labRow = c("Survival (days)", "Histological type", "Tumour stage", "Histological grade"))
dev.off()



tumourME_females <- moduleEigengenes(females_tumour, labels2colors(females_tumour_network$colors))$eigengenes
tumourME_females <- cbind(sample = rownames(females_tumour), tumourME_females) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  dplyr::select(-MEgrey) %>%
  inner_join(stomach_tcga_clinical2, by = "sample") %>%
  dplyr::rename(ajcc_tumor_stage = ajcc_pathologic_tumor_stage) %>%
  mutate(ajcc_tumor_stage2 = as.numeric(as.factor(ajcc_tumor_stage)), histological_type2 = as.numeric(as.factor(histological_type)), histological_grade2 = as.numeric(as.factor(histological_grade)))
write.table(tumourME_females, "./files/stomach_tumourME_traits_females.txt", sep="\t", quote=F, row.names = F)


#female modules (all conserved in male modules)
tumourME_females_cor <- cor(
  tumourME_females %>% dplyr::select(OS.time, histological_type2, ajcc_tumor_stage2, histological_grade2),
  tumourME_females %>% dplyr::select(starts_with("ME")),
  method = "pearson") %>%
  abs()


pdf(file="./plots/wgcna_networks_traits/stomach_tumourME_females_cor.pdf", height=4, width=8)
par(oma = c(4,0,0,5))
heatmap.2( x = tumourME_females_cor,
  trace="none",
  Rowv=TRUE,
  Colv=TRUE,
  dendrogram="both",
  cexRow=1.2,
  cexCol = 1.2,
  key=TRUE,
  keysize=1,
  symkey=TRUE,
  key.title="Pearson's r",
  key.xlab = NA,
  key.ylab = NA,
  density.info = "none",
  col=bluered(100),
  ColSideColors = stomach_tumour_gender_diff_networks %>% filter(network == "females") %>% pull(colors),
  labRow = c("Overall survival", "Histological type", "AJCC tumor stage", "Histological grade"))
dev.off()








# melted datasets
males_tumour_network_MEs <- cbind(sample = rownames(males_tumour), males_tumour_network$MEs) %>%
  as.tibble() %>%
  gather(key = "moduleN", value = "ME", -sample) %>%
  mutate(moduleN = as.numeric(str_replace_all(moduleN, "[ME]", "")), sample = str_sub(sample, 1, 12)) %>%
  mutate(moduleL = labels2colors(moduleN)) %>%
  filter(moduleL != "grey") %>%
  inner_join(stomach_tumour_gender_diff_networks %>% filter(network == "males"), by=c("moduleL" = "module")) %>%
  arrange(moduleL) %>%
  dplyr::select(network, moduleL, moduleN, state, sample, ME) %>%
  inner_join(stomach_tcga_clinical2, by=c("sample")) %>%
  mutate(ajcc_pathologic_tumor_stage = as.numeric(as.factor(ajcc_pathologic_tumor_stage)), histological_type = as.numeric(as.factor(histological_type)), histological_grade = as.numeric(as.factor(histological_grade)))
write.table(males_tumour_network_MEs, "./files/stomach_males_tumour_network_MEs_traits.txt", quote=F, sep="\t", row.names=F)


females_tumour_network_MEs <- cbind(sample = rownames(females_tumour), females_tumour_network$MEs) %>%
  as.tibble() %>%
  gather(key = "moduleN", value = "ME", -sample) %>%
  mutate(moduleN = as.numeric(str_replace_all(moduleN, "[ME]", "")), sample = str_sub(sample, 1, 12)) %>%
  mutate(moduleL = labels2colors(moduleN)) %>%
  filter(moduleL != "grey") %>%
  inner_join(stomach_tumour_gender_diff_networks %>% filter(network == "females"), by=c("moduleL" = "module")) %>%
  arrange(moduleL) %>%
  dplyr::select(network, moduleL, moduleN, state, sample, ME) %>%
  inner_join(stomach_tcga_clinical2, by=c("sample")) %>%
  mutate(ajcc_pathologic_tumor_stage = as.numeric(as.factor(ajcc_pathologic_tumor_stage)), histological_type = as.numeric(as.factor(histological_type)), histological_grade = as.numeric(as.factor(histological_grade)))
write.table(females_tumour_network_MEs, "./files/stomach_females_tumour_network_MEs_traits.txt", quote=F, sep="\t", row.names=F)