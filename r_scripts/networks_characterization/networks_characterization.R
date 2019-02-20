# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Characterization



library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(viridis)




# load gene modules


# -- Stomach

sto_males_tumour_modules <- read_tsv("./files/stomach_males_tumour_modules.txt") %>%
  dplyr::select(genes_ens, geneName, geneType, chrom, moduleN, moduleL) %>%
  mutate(tissue = "Stomach", type = "Tumour", gender = "Male")

sto_females_tumour_modules <- read_tsv("./files/stomach_females_tumour_modules.txt") %>%
  dplyr::select(genes_ens, geneName, geneType, chrom, moduleN, moduleL) %>%
  mutate(tissue = "Stomach", type = "Tumour", gender = "Female")

sto_males_normal_modules <- read_tsv("./files/stomach_males_normal_modules.txt") %>%
  dplyr::select(genes_ens, geneName, geneType, chrom, moduleN, moduleL) %>%
  mutate(tissue = "Stomach", type = "Normal", gender = "Male")

sto_females_normal_modules <- read_tsv("./files/stomach_females_normal_modules.txt") %>%
  dplyr::select(genes_ens, geneName, geneType, chrom, moduleN, moduleL) %>%
  mutate(tissue = "Stomach", type = "Normal", gender = "Female")


# -- Thyroid

thy_males_tumour_modules <- read_tsv("./files/thyroid_males_tumour_modules.txt") %>%
  dplyr::select(genes_ens, geneName, geneType, chrom, moduleN, moduleL) %>%
  mutate(tissue = "Thyroid", type = "Tumour", gender = "Male")

thy_females_tumour_modules <- read_tsv("./files/thyroid_females_tumour_modules.txt") %>%
  dplyr::select(genes_ens, geneName, geneType, chrom, moduleN, moduleL) %>%
  mutate(tissue = "Thyroid", type = "Tumour", gender = "Female")

thy_males_normal_modules <- read_tsv("./files/thyroid_males_normal_modules.txt") %>%
  dplyr::select(genes_ens, geneName, geneType, chrom, moduleN, moduleL) %>%
  mutate(tissue = "Thyroid", type = "Normal", gender = "Male")

thy_females_normal_modules <- read_tsv("./files/thyroid_females_normal_modules.txt") %>%
  dplyr::select(genes_ens, geneName, geneType, chrom, moduleN, moduleL) %>%
  mutate(tissue = "Thyroid", type = "Normal", gender = "Female")



networks <- bind_rows(
  sto_males_tumour_modules,
  sto_females_tumour_modules,
  sto_males_normal_modules,
  sto_females_normal_modules,
  thy_males_tumour_modules,
  thy_females_tumour_modules,
  thy_males_normal_modules,
  thy_females_normal_modules)



# number of modules by tissue and gender
wgcna_modules <- networks %>%
  group_by(tissue, type, gender) %>%
  summarise(mod = length(unique(moduleL))) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = gender, y = mod, fill = gender )) +
    geom_bar(stat = "identity") +
    facet_grid(tissue ~ type) +
    scale_fill_manual(values=c("#fbb4b9", "#74a9cf"), guide=F) +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    labs(x = "Gender", y = "Number of modules")
ggsave(filename="wgcna_modules.png", plot = wgcna_modules, path = "./plots/wgcna_networks", width=4, height=6)
unlink("wgcna_modules.png")




# number of genes by module
wgcna_modules_genes <- networks %>%
  group_by(tissue, type, gender, moduleL, moduleN) %>%
  summarise( n_genes = n() ) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(moduleN = as.character(moduleN)) %>%
  mutate(moduleN = fct_reorder(moduleN, n_genes)) %>%
  ggplot() +
    geom_bar(mapping = aes(x = gender, y = n_genes, fill = moduleN), stat = "identity", position = "dodge") +
    #geom_text(mapping = aes(x = gender, y = n_genes, group = moduleN, label = moduleN), position = position_dodge(width = 1), vjust = -0.5, size = 2, color="black") +
    scale_fill_viridis(discrete=T, option = "C", guide=F) +
    facet_grid(tissue ~ type, scales = "fixed") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    labs(x = "Gender", y = "Number of genes per module")
ggsave(filename="wgcna_modules_genes.png", plot = wgcna_modules_genes, path = "./plots/wgcna_networks", width=10, height=6)
unlink("wgcna_modules_genes.png")








# load enrichment tables


# -- Stomach

stomach_hyp <- read_tsv("./files/stomach_modules_hypergeo_enr.txt") %>%
  mutate(source = "Stomach")


stomach_gse <- read_tsv("./files/stomach_modules_gsea_enr.txt") %>%
  mutate(source = "Stomach")



# -- Thyroid

thyroid_hyp <- read_tsv("./files/thyroid_modules_hypergeo_enr.txt") %>%
  mutate(source = "Thyroid")


thyroid_gse <- read_tsv("./files/thyroid_modules_gsea_enr.txt") %>%
  mutate(source = "Thyroid")




thyroid_hyp_bp <- thyroid_hyp %>%
  filter(p.adjust < 0.05) %>%
  filter(state2 == "female-specific" | state2 == "male-specific") %>%
  ggplot(mapping = aes(x = state2, y = ..count.., fill = Description)) +
    geom_bar(position = "dodge") +
    facet_wrap(~ tissue) +
    scale_fill_discrete(name = "Term") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=13),
      axis.text.x = element_text(colour="black", size=10),
      axis.text.y = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    scale_x_discrete(labels = c("Female\nspecific", "Male\nspecific")) +
    labs(y = "Enriched terms (FDR < 5%)", x = "Sex")
ggsave(filename="thyroid_networks_number_terms_hyp.png", plot = thyroid_hyp_bp, path = "./plots/wgcna_modules_characterization", width=5, height=3)
unlink("thyroid_networks_number_terms_hyp.png")



stomach_hyp_bp <- stomach_hyp %>%
  filter(p.adjust < 0.05) %>%
  filter(state2 == "female-specific" | state2 == "male-specific") %>%
  ggplot(mapping = aes(x = state2, y = ..count.., fill = Description)) +
    geom_bar(position = "dodge") +
    facet_wrap(~ tissue) +
    scale_fill_discrete(name = "Term") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=13),
      axis.text.x = element_text(colour="black", size=10),
      axis.text.y = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    scale_x_discrete(labels = c("Male-specific")) +
    labs(y = "Enriched terms (FDR < 5%)", x = "Sex")
ggsave(filename="stomach_networks_number_terms_hyp.png", plot = stomach_hyp_bp, path = "./plots/wgcna_modules_characterization", width=5, height=3)
unlink("stomach_networks_number_terms_hyp.png")
