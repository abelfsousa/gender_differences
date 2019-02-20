# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# MEs and sample traits


# -- Thyroid


library(tidyverse)
library(WGCNA)
library(data.table)
library(gplots)



# load R workspace
load("./r_workspaces/tcga_gtex_thyroid_wgcna_networks.RData")


# load networks comparison Male vs Female in TCGA tumours
thyroid_tumour_gender_diff_networks <- read_tsv("./files/thyroid_tumour_gender_diff_networks.txt")


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
  dplyr::select(-MEgrey)


tumourME_males <- inner_join(tumourME_males, thyroid_tcga_clinical2, by = "sample") %>%
  dplyr::select(-c(sample, gender)) %>%
  mutate(ajcc_pathologic_tumor_stage = as.numeric(as.factor(ajcc_pathologic_tumor_stage)), histological_type = as.numeric(as.factor(histological_type)))


tumourME_males_cor <- cor(tumourME_males %>% dplyr::select(OS.time, histological_type, ajcc_pathologic_tumor_stage), tumourME_males %>% dplyr::select(-c(OS.time, histological_type, ajcc_pathologic_tumor_stage)), method = "pearson")


pdf(file="./plots/wgcna_networks_traits/thyroid_tumourME_males_cor.pdf", height=4, width=6)
par(oma = c(2,0,0,4))
heatmap.2( x = tumourME_males_cor, trace="none", Rowv=TRUE, Colv=TRUE, dendrogram="both", cexRow=1.2, cexCol = 1.2, key=TRUE, keysize=1, symkey=TRUE, key.title="Pearson's r", key.xlab = NA, key.ylab = NA, density.info = "none", col=bluered(100) )
dev.off()



tumourME_females <- moduleEigengenes(females_tumour, labels2colors(females_tumour_network$colors))$eigengenes
tumourME_females <- cbind(sample = rownames(females_tumour), tumourME_females) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  dplyr::select(-MEgrey)


tumourME_females <- inner_join(tumourME_females, thyroid_tcga_clinical2, by = "sample") %>%
  dplyr::select(-c(sample, gender)) %>%
  mutate(ajcc_pathologic_tumor_stage = as.numeric(as.factor(ajcc_pathologic_tumor_stage)), histological_type = as.numeric(as.factor(histological_type)))


tumourME_females_cor <- cor(tumourME_females %>% dplyr::select(OS.time, histological_type, ajcc_pathologic_tumor_stage), tumourME_females %>% dplyr::select(-c(OS.time, histological_type, ajcc_pathologic_tumor_stage)), method = "pearson")


pdf(file="./plots/wgcna_networks_traits/thyroid_tumourME_females_cor.pdf", height=4, width=6)
par(oma = c(2,0,0,4))
heatmap.2( x = tumourME_females_cor, trace="none", Rowv=TRUE, Colv=TRUE, dendrogram="both", cexRow=1.2, cexCol = 1.2, key=TRUE, keysize=1, symkey=TRUE, key.title="Pearson's r", key.xlab = NA, key.ylab = NA, density.info = "none", col=bluered(100) )
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
