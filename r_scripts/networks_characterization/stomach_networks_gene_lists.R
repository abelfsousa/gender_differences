# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data




# -- Stomach

# load R workspace
load("./r_workspaces/tcga_gtex_stomach_wgcna_networks.RData")


library(WGCNA)
library(tidyverse)


stomach_tumour_gender_diff_networks <- read_tsv("./files/stomach_tumour_gender_diff_networks.txt")
stomach_normal_gender_diff_networks <- read_tsv("./files/stomach_normal_gender_diff_networks.txt")



# males tumour network
males_tumor_kme <- signedKME(males_tumour, males_tumour_network$MEs)
males_tumor_kme <- cbind(gene = rownames(males_tumor_kme), males_tumor_kme)
males_tumor_kme <- males_tumor_kme %>%
  as.tibble() %>%
  gather(key = "M", value = "kme", -gene) %>%
  mutate(M = str_replace_all(M, "[kME]", ""))

males_tumour_modules <- tibble(genes_ens = colnames(males_tumour)) %>%
  inner_join(tcga.geneIDs.annot, by = c("genes_ens" = "geneID")) %>%
  mutate(moduleN = as.character(males_tumour_network$colors), moduleL = labels2colors(males_tumour_network$colors)) %>%
  inner_join(males_tumor_kme, by=c("genes_ens" = "gene", "moduleN" = "M")) %>%
  filter(moduleL != "grey") %>%
  inner_join(stomach_tumour_gender_diff_networks %>% filter(network == "males") %>% dplyr::select(-network), by=c("moduleL" = "module")) %>%
  arrange(moduleL)
write.table(males_tumour_modules, "./files/stomach_males_tumour_modules.txt", quote=F, sep="\t", row.names=F)


# females tumour network
females_tumor_kme <- signedKME(females_tumour, females_tumour_network$MEs)
females_tumor_kme <- cbind(gene = rownames(females_tumor_kme), females_tumor_kme)
females_tumor_kme <- females_tumor_kme %>%
  as.tibble() %>%
  gather(key = "M", value = "kme", -gene) %>%
  mutate(M = str_replace_all(M, "[kME]", ""))

females_tumour_modules <- tibble(genes_ens = colnames(females_tumour)) %>%
  inner_join(tcga.geneIDs.annot, by = c("genes_ens" = "geneID")) %>%
  mutate(moduleN = as.character(females_tumour_network$colors), moduleL = labels2colors(females_tumour_network$colors)) %>%
  inner_join(females_tumor_kme, by=c("genes_ens" = "gene", "moduleN" = "M")) %>%
  filter(moduleL != "grey") %>%
  inner_join(stomach_tumour_gender_diff_networks %>% filter(network == "females") %>% dplyr::select(-network), by=c("moduleL" = "module")) %>%
  arrange(moduleL)
write.table(females_tumour_modules, "./files/stomach_females_tumour_modules.txt", quote=F, sep="\t", row.names=F)


# males normal network
males_normal_kme <- signedKME(males_normal, males_normal_network$MEs)
males_normal_kme <- cbind(gene = rownames(males_normal_kme), males_normal_kme)
males_normal_kme <- males_normal_kme %>%
  as.tibble() %>%
  gather(key = "M", value = "kme", -gene) %>%
  mutate(M = str_replace_all(M, "[kME]", ""))

males_normal_modules <- tibble(genes_ens = colnames(males_normal)) %>%
  inner_join(tcga.geneIDs.annot, by = c("genes_ens" = "geneID")) %>%
  mutate(moduleN = as.character(males_normal_network$colors), moduleL = labels2colors(males_normal_network$colors)) %>%
  inner_join(males_normal_kme, by=c("genes_ens" = "gene", "moduleN" = "M")) %>%
  filter(moduleL != "grey") %>%
  inner_join(stomach_normal_gender_diff_networks %>% filter(network == "males") %>% dplyr::select(-network), by=c("moduleL" = "module")) %>%
  arrange(moduleL)
write.table(males_normal_modules, "./files/stomach_males_normal_modules.txt", quote=F, sep="\t", row.names=F)


# females normal network
females_normal_kme <- signedKME(females_normal, females_normal_network$MEs)
females_normal_kme <- cbind(gene = rownames(females_normal_kme), females_normal_kme)
females_normal_kme <- females_normal_kme %>%
  as.tibble() %>%
  gather(key = "M", value = "kme", -gene) %>%
  mutate(M = str_replace_all(M, "[kME]", ""))

females_normal_modules <- tibble(genes_ens = colnames(females_normal)) %>%
  inner_join(tcga.geneIDs.annot, by = c("genes_ens" = "geneID")) %>%
  mutate(moduleN = as.character(females_normal_network$colors), moduleL = labels2colors(females_normal_network$colors)) %>%
  inner_join(females_normal_kme, by=c("genes_ens" = "gene", "moduleN" = "M")) %>%
  filter(moduleL != "grey") %>%
  inner_join(stomach_normal_gender_diff_networks %>% filter(network == "females") %>% dplyr::select(-network), by=c("moduleL" = "module")) %>%
  arrange(moduleL)
write.table(females_normal_modules, "./files/stomach_females_normal_modules.txt", quote=F, sep="\t", row.names=F)
