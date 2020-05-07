# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between tumour and normal samples
# Normal samples from GTEx

# Thyroid and Stomach




library(tidyverse)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)


# load TCGA annotation
tcga_annotation <- read_tsv("./data/annotation/geneAnnot.gencode.v22.txt") %>%
  mutate(geneID = str_replace(geneID, "\\.[0-9]+", "")) %>%
  dplyr::rename(ensemblGeneID = geneID)

chrom_genes <- tcga_annotation %>%
  filter(geneType == "protein_coding" | geneType == "lincRNA") %>%
  group_by(chrom) %>%
  summarise(total_genes = n())



# load degs tables
stomach_degs <- read_tsv("./files/stomach_tumour_normal_signf_degs.txt") %>%
  dplyr::select(-c(FDR_tumour, FDR_normal)) %>%
  mutate(tissue = "Stomach")


thyroid_degs <- read_tsv("./files/thyroid_tumour_normal_signf_degs.txt") %>%
  dplyr::select(-c(FDR_tumour, FDR_normal)) %>%
  mutate(tissue = "Thyroid")


degs <- bind_rows(stomach_degs, thyroid_degs) %>%
  inner_join(tcga_annotation %>% filter(geneType %in% c("protein_coding", "lincRNA")) %>% dplyr::select(geneName, geneType) %>% distinct(), by = "geneName")


degs2 <- degs %>%
  dplyr::select(tissue, genes, geneName, chrom, state, log2FC_tumour, log2FC_normal) %>%
  gather(key = "fc_type", "value", -c(tissue, genes, geneName, chrom, state)) %>%
  na.exclude()


# load enrichment tables
stomach_enr <- read_tsv("./files/stomach_tumour_normal_degs_enr.txt") %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(-c(GeneRatio, BgRatio)) %>%
  mutate(tissue = "Stomach")

thyroid_enr <- read_tsv("./files/thyroid_tumour_normal_degs_enr.txt") %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(-c(GeneRatio, BgRatio)) %>%
  mutate(tissue = "Thyroid")


enrch <- bind_rows(stomach_enr, thyroid_enr)




# barplot of number of up-regulated DEGs in tumour/normal tissue by DEG type
logFC_barpl <- degs %>%
  dplyr::select(tissue, state, log2FC_tumour, log2FC_normal) %>%
  gather(-c(tissue, state), key = "fc_type", value = "log2FC") %>%
  na.exclude() %>%
  mutate(signal = if_else(log2FC > 0, "up_male", "up_female")) %>%
  group_by(tissue, state, fc_type, signal) %>%
  count() %>%
  ungroup() %>%
  ggplot(mapping = aes(x = state, y = n, fill = signal, color = fc_type)) +
    geom_bar(stat = "identity", position = "dodge", lwd=1) +
    facet_wrap(~ tissue) +
    scale_color_manual(values=c("green", "#d95f02"), name="Tissue", labels = c("Normal", "Tumour")) +
    scale_fill_manual(values=c("#fbb4b9", "#74a9cf"), name="DEG state", labels = c("Up-regulated female", "Up-regulated male")) +
    theme_classic()
logFC_barpl <- logFC_barpl +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=13),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15),
      legend.position = "bottom") +
    scale_x_discrete(labels = c("Common", "Normal\nspecific", "Tumour\nspecifc")) +
    labs(y = "Number of DEGs", x = "DEGs group") +
    guides(fill=guide_legend(nrow=2), color=guide_legend(nrow=2))

ggsave(filename="tumour_normal_degs_MvsF_logFC_barpl.png", plot = logFC_barpl, path = "./plots/diff_expression_maleVSfemale_gtex_normal", width=7, height=5)
unlink("tumour_normal_degs_MvsF_logFC_barpl.png")




# density plot of log2FC distribution of normal-specific GO-BP/CM enriched categories

normal_specific_enrch <- enrch %>%
  filter(state == "normal_specific") %>%
  group_by(tissue, Description) %>%
  top_n(-5, p.adjust) %>%
  ungroup() %>%
  dplyr::select(tissue, state, Description, ID, geneID) %>%
  mutate(geneName = str_split(geneID, "/")) %>%
  unnest() %>%
  dplyr::select(-geneID) %>%
  inner_join(degs2, by = c("tissue", "state", "geneName"))



normal_specific_enrch_thyroid <- normal_specific_enrch %>%
  filter(tissue == "Thyroid") %>%
  ggplot(mapping = aes(x = value, color = ID, fill = ID)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
    facet_wrap(
      ~ Description,
      ncol = 2,
      nrow = 2,
      labeller=labeller(Description = c("GO_BP" = "GO BP", "KEGG" = "KEGG"))) +
    scale_color_viridis(discrete = TRUE, option = "C", name = "Term") +
    scale_fill_viridis(discrete = TRUE, option = "C", name = "Term") +
    #scale_colour_brewer(palette = "Set3", name = "Term") +
    #scale_fill_brewer(palette = "Set3", name = "Term") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=13),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=16),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15),
      plot.title=element_text(colour="black", size=18, hjust = 0.5)) +
    labs(x = "log2FC", y = "Density", title = "Thyroid")

ggsave(filename="normal_specific_enrch_thyroid.png", plot = normal_specific_enrch_thyroid, path = "./plots/diff_expression_maleVSfemale_gtex_normal", width=10, height=4)
ggsave(filename="normal_specific_enrch_thyroid.pdf", plot = normal_specific_enrch_thyroid, path = "./plots/diff_expression_maleVSfemale_gtex_normal", width=10, height=4)
unlink("normal_specific_enrch_thyroid.png")
unlink("normal_specific_enrch_thyroid.pdf")




normal_specific_enrch_stomach <- normal_specific_enrch %>%
  filter(tissue == "Stomach", Description != "POS") %>%
  ggplot(mapping = aes(x = value, color = ID, fill = ID)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
    facet_wrap(
      ~ Description,
      labeller=labeller(Description = c("GO_BP" = "GO BP", "KEGG" = "KEGG", "POS" = "Positional"))) +
    #scale_color_viridis(discrete = TRUE, option = "C", name = "Term") +
    #scale_fill_viridis(discrete = TRUE, option = "C", name = "Term") +
    scale_colour_brewer(palette = "Set1", name = "Term") +
    scale_fill_brewer(palette = "Set1", name = "Term") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=12),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=16),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15),
      plot.title=element_text(colour="black", size=18, hjust = 0.5),
      panel.spacing = unit(0.5, "lines")) +
    labs(x = "log2FC", y = "Density", title = "Stomach")

ggsave(filename="normal_specific_enrch_stomach.png", plot = normal_specific_enrch_stomach, path = "./plots/diff_expression_maleVSfemale_gtex_normal", width=10, height=4)
ggsave(filename="normal_specific_enrch_stomach.pdf", plot = normal_specific_enrch_stomach, path = "./plots/diff_expression_maleVSfemale_gtex_normal", width=10, height=4)
unlink("normal_specific_enrch_stomach.png")
unlink("normal_specific_enrch_stomach.pdf")



# barplot of number of DEGs by chromosome
degs_chr_bp_stomach <- degs %>%
  filter(tissue == "Stomach") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(chrom = fct_relevel(chrom, c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))) %>%
  ggplot(mapping = aes(x=chrom, y = ..count.., fill = state)) +
  geom_bar() +
  scale_fill_manual(values = c("#bf812d", "#a1d76a", "#ca0020"), labels = c("Common", "Normal", "Tumour"), name = "SBG type") +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=18),
    axis.text.y=element_text(colour="black", size=16),
    axis.text.x=element_text(colour="black", size=13),
    plot.title = element_blank(),
    strip.text.y = element_text(colour="black", size=10),
    strip.text.x = element_text(colour="black", size=14),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=16),
    legend.title = element_text(colour="black", size=18)) +
  scale_y_continuous(name = "Number of genes", limits = c(NA, 40)) +
  scale_x_discrete(name = "Chromosome")

ggsave(filename="stomach_degs_chr_bp2.png", plot=degs_chr_bp_stomach, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 13, height = 3)
ggsave(filename="stomach_degs_chr_bp2.pdf", plot=degs_chr_bp_stomach, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 13, height = 3)
unlink("stomach_degs_chr_bp2.png")
unlink("stomach_degs_chr_bp2.pdf")


degs_chr_bp_stomach <- degs %>%
  filter(tissue == "Stomach") %>%
  mutate(chrom_type = map_chr(.x = chrom, .f = ~ if(!.x %in% c("chrX", "chrY")){if(.x == "chrM"){"chrM"}else{"Autosomes"}}else{.x})) %>%
  mutate_if(is.character, as.factor) %>%
  ggplot(mapping = aes(x=chrom_type, fill = state)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("#bf812d", "#a1d76a", "#ca0020"), labels = c("Common", "Normal", "Tumour"), name = "SBG type") +
  #theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=10),
    axis.text.y=element_text(colour="black", size=8),
    axis.text.x=element_text(colour="black", size=8),
    plot.title = element_blank(),
    strip.text.y = element_text(colour="black", size=8),
    strip.text.x = element_text(colour="black", size=8),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=8),
    legend.title = element_text(colour="black", size=10),
    legend.key.size = unit(0.2, "cm")) +
  scale_y_continuous(name = "Count", limits = c(NA, 40)) +
  scale_x_discrete(name = "Chromosome")

ggsave(filename="stomach_degs_chr_bp3.png", plot=degs_chr_bp_stomach, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 3, height = 2)
ggsave(filename="stomach_degs_chr_bp3.pdf", plot=degs_chr_bp_stomach, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 3, height = 2)
unlink("stomach_degs_chr_bp3.png")
unlink("stomach_degs_chr_bp3.pdf")

write_rds(degs_chr_bp_stomach, "./r_objects/plots/figure2/stomach_degs_chr_barplot.rds")


degs_chr_bp_thyroid <- degs %>%
  filter(tissue == "Thyroid") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(chrom = fct_relevel(chrom, c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))) %>%
  ggplot(mapping = aes(x=chrom, y = ..count.., fill = state)) +
  geom_bar() +
  scale_fill_manual(values = c("#bf812d", "#a1d76a", "#ca0020"), labels = c("Common", "Normal", "Tumour"), name = "SBG type") +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=18),
    axis.text.y=element_text(colour="black", size=16),
    axis.text.x=element_text(colour="black", size=13),
    plot.title = element_blank(),
    strip.text.y = element_text(colour="black", size=10),
    strip.text.x = element_text(colour="black", size=14),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=16),
    legend.title = element_text(colour="black", size=18)) +
  scale_y_continuous(name = "Number of genes", limits = c(NA, 100), breaks = seq(0, 120, by = 20)) +
  scale_x_discrete(name = "Chromosome")

ggsave(filename="thyroid_degs_chr_bp2.png", plot=degs_chr_bp_thyroid, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 15, height = 3)
ggsave(filename="thyroid_degs_chr_bp2.pdf", plot=degs_chr_bp_thyroid, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 15, height = 3)
unlink("thyroid_degs_chr_bp2.png")
unlink("thyroid_degs_chr_bp2.pdf")


degs_chr_bp_thyroid <- degs %>%
  filter(tissue == "Thyroid") %>%
  mutate(chrom_type = map_chr(.x = chrom, .f = ~ if(!.x %in% c("chrX", "chrY")){if(.x == "chrM"){"chrM"}else{"Autosomes"}}else{.x})) %>%
  mutate_if(is.character, as.factor) %>%
  filter(chrom_type != "chrM") %>%
  ggplot(mapping = aes(x=chrom_type, fill = state)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("#bf812d", "#a1d76a", "#ca0020"), labels = c("Common", "Normal", "Tumour"), name = "SBG type") +
  #theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=10),
    axis.text.y=element_text(colour="black", size=8),
    axis.text.x=element_text(colour="black", size=8),
    plot.title = element_blank(),
    strip.text.y = element_text(colour="black", size=8),
    strip.text.x = element_text(colour="black", size=8),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=8),
    legend.title = element_text(colour="black", size=10),
    legend.key.size = unit(0.2, "cm")) +
  scale_y_continuous(name = "Count") +
  scale_x_discrete(name = "Chromosome")

ggsave(filename="thyroid_degs_chr_bp3.png", plot=degs_chr_bp_thyroid, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 3, height = 2)
ggsave(filename="thyroid_degs_chr_bp3.pdf", plot=degs_chr_bp_thyroid, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 3, height = 2)
unlink("thyroid_degs_chr_bp3.png")
unlink("thyroid_degs_chr_bp3.pdf")

write_rds(degs_chr_bp_thyroid, "./r_objects/plots/figure2/thyroid_degs_chr_barplot.rds")


# density plot of log2FC distribution of normal/tumour-specific DEGs in the sexual chromosomes
stomach_normal_spec_degs_sex_chrom <- degs2 %>%
  filter(tissue == "Stomach", state != "common" & chrom == "chrX") %>%
  dplyr::select(geneName, state, value) %>%
  ggplot(mapping = aes(x = value, color = state, fill = state)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=13),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    labs(x = "log2FC", y = "Density")

ggsave(filename="stomach_normal_spec_degs_sex_chrom.png", plot = stomach_normal_spec_degs_sex_chrom, path = "./plots/diff_expression_maleVSfemale_gtex_normal", width=5, height=3)
unlink("stomach_normal_spec_degs_sex_chrom.png")



thyroid_normal_spec_degs_sex_chrom <- degs2 %>%
  filter(tissue == "Thyroid", state != "common" & chrom == "chrX") %>%
  dplyr::select(geneName, state, value) %>%
  ggplot(mapping = aes(x = value, color = state, fill = state)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=13),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    labs(x = "log2FC", y = "Density")

ggsave(filename="thyroid_normal_spec_degs_sex_chrom.png", plot = thyroid_normal_spec_degs_sex_chrom, path = "./plots/diff_expression_maleVSfemale_gtex_normal", width=5, height=3)
unlink("thyroid_normal_spec_degs_sex_chrom.png")



# cancer drivers gene list
driver_genes <- data.table::fread("./data/tcga/cancer_driver_genes.txt") %>%
  as_tibble()

driver_genes_summary <- driver_genes %>%
  group_by(Gene, Decision) %>%
  summarise(n = n()) %>%
  ungroup()


thyroid_enr %>%
  filter(state == "normal_specific") %>%
  select(tissue, state, gene=geneID) %>%
  mutate(gene = str_split(gene, "/")) %>%
  unnest(cols = gene) %>%
  distinct() %>%
  inner_join(
    thyroid_degs %>% filter(state == "normal_specific") %>% select(gene=geneName, log2fc=log2FC_normal),
    by = "gene") %>%
  filter(log2fc < 0) %>%
  filter(gene %in% driver_genes_summary$Gene)


stomach_enr %>%
  filter(state == "normal_specific") %>%
  select(tissue, state, gene=geneID) %>%
  mutate(gene = str_split(gene, "/")) %>%
  unnest(cols = gene) %>%
  distinct() %>%
  inner_join(
    stomach_degs %>% filter(state == "normal_specific") %>% select(gene=geneName, log2fc=log2FC_normal),
    by = "gene") %>%
  filter(log2fc < 0) %>%
  filter(gene %in% driver_genes_summary$Gene)





save(list=ls(), file="r_workspaces/thyroid_stomach_MvsF_tcga_gtex_characterization.RData")
