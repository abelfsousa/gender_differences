# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between tumour and normal samples
# Normal samples from TCGA

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
stomach_degs <- read_tsv("./files/stomach_tumour_normal_signf_degs2.txt") %>%
  dplyr::select(-c(FDR_tumour, FDR_normal)) %>%
  mutate(tissue = "Stomach")


thyroid_degs <- read_tsv("./files/thyroid_tumour_normal_signf_degs2.txt") %>%
  dplyr::select(-c(FDR_tumour, FDR_normal)) %>%
  mutate(tissue = "Thyroid")


degs <- bind_rows(stomach_degs, thyroid_degs) %>%
  inner_join(tcga_annotation %>% filter(geneType %in% c("protein_coding", "lincRNA")) %>% dplyr::select(geneName, geneType) %>% distinct(), by = "geneName")

degs2 <- degs %>%
  dplyr::select(tissue, genes, geneName, chrom, state, log2FC_tumour, log2FC_normal) %>%
  gather(key = "fc_type", "value", -c(tissue, genes, geneName, chrom, state)) %>%
  na.exclude()



# load enrichment tables
stomach_enr <- read_tsv("./files/stomach_tumour_normal_degs_enr2.txt") %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(-c(GeneRatio, BgRatio)) %>%
  mutate(tissue = "Stomach")

thyroid_enr <- read_tsv("./files/thyroid_tumour_normal_degs_enr2.txt") %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(-c(GeneRatio, BgRatio)) %>%
  mutate(tissue = "Thyroid")


enrch <- bind_rows(stomach_enr, thyroid_enr)



# load normal differential expression results M vs F
thyroid_normal_table <- read_tsv("./files/diff_expr_thyroid_normal_tcga_maleVSfemale_edgeR.txt") %>%
  mutate(tissue = "Thyroid")

stomach_normal_table <- read_tsv("./files/diff_expr_stomach_normal_tcga_maleVSfemale_edgeR.txt") %>%
  mutate(tissue = "Stomach")

normal_diff_expr <- bind_rows(thyroid_normal_table, stomach_normal_table) %>%
  dplyr::select(tissue, geneName, logFC, FDR)




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
    theme_classic() +
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
ggsave(filename="tumour_normal_degs_MvsF_logFC_barpl.png", plot = logFC_barpl, path = "./plots/diff_expression_maleVSfemale_tcga_normal", width=7, height=5)
unlink("tumour_normal_degs_MvsF_logFC_barpl.png")



# density plot of log2FC distribution of the tumour-specific DEGs enriched in X chromosome.
# the log2FC distribution in the normal samples is also shown

tumour_specific_chrx_tumour <- enrch %>%
  dplyr::select(tissue, state, ID, geneID) %>%
  filter(state == "tumour_specific") %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest() %>%
  inner_join(degs2[, c("tissue", "state", "geneName", "fc_type", "value")], by = c("tissue", "state", "geneID" = "geneName"))


tumour_specific_chrx_normal <- enrch %>%
  dplyr::select(tissue, state, ID, geneID) %>%
  filter(state == "tumour_specific") %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest() %>%
  inner_join(normal_diff_expr[, c("tissue", "geneName", "logFC")], by = c("tissue", "geneID" = "geneName")) %>%
  mutate(fc_type = "log2FC_normal") %>%
  dplyr::rename(value = logFC)


tumour_specific_chrx <- bind_rows(tumour_specific_chrx_tumour, tumour_specific_chrx_normal) %>%
  ggplot(mapping = aes(x = value, color = fc_type, fill = fc_type)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
    facet_grid(tissue ~ ID) +
    #scale_color_viridis(discrete = TRUE, option = "E") +
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
ggsave(filename="tumour_specific_chrx_log2FC.png", plot = tumour_specific_chrx, path = "./plots/diff_expression_maleVSfemale_tcga_normal", width=10, height=5)
unlink("tumour_specific_chrx_log2FC.png")





save(list=ls(), file="r_workspaces/thyroid_stomach_MvsF_tcga_characterization.RData")
