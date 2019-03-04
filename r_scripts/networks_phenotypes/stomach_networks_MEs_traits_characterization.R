# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Differential hub genes expression between histological types


# -- Stomach


library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(gplots)
library(mclust)
library(survival)
library(data.table)
library(ranger)
library(ggfortify)

set.seed(123)


# load male tumour modules with kmE
males_tumour_modules <- read_tsv("./files/stomach_males_tumour_modules.txt")


# load MEs and traits for female tumours
tumourME_males <- read_tsv("./files/stomach_tumourME_traits_males.txt")


# load fpkm matrix
stomach_males_fpkm <- read.table("./files/stomach_males_tumour_fpkm_wgcna.txt", sep="\t", h=T)
stomach_males_fpkm <- t(stomach_males_fpkm)
colnames(stomach_males_fpkm) <- substr(colnames(stomach_males_fpkm), 1, 12)



# load TCGA clinical data
stomach_tcga_clinical <- fread("./data/tcga/tcga_clinical_data.txt") %>%
  as.tibble() %>%
  mutate(bcr_patient_barcode = str_replace_all(bcr_patient_barcode, "-", ".")) %>%
  filter(type == "STAD") %>%
  dplyr::select(sample = bcr_patient_barcode, gender, OS, OS.time, histological_type, ajcc_pathologic_tumor_stage, histological_grade) %>%
  filter(!(histological_type %in% c("[Discrepancy]", "[Not Available]")) & !(ajcc_pathologic_tumor_stage %in% c("[Discrepancy]", "[Not Available]"))) %>%
  filter(!is.na(OS.time)) %>%
  filter(sample %in% tumourME_males$sample)




# skyblue
# male-specific tumour module correlated with histological type
# select hub genes (kme >= 0.8)
skyblue <- males_tumour_modules %>%
  filter(moduleL == "skyblue", kme >= 0.8)


# select hub genes and samples used for correlation between eigenenes and histological type
# test for expression differences between histological types using a kruskal wallis test
skyblue_genes <- stomach_males_fpkm[skyblue$genes_ens, tumourME_males$sample]
skyblue_genes <- apply(skyblue_genes, 1, function(x) kruskal.test(as.numeric(x), as.factor(tumourME_males$histological_type))$p.value) %>%
  as.data.frame() %>%
  setNames(c("kruskal_pval")) %>%
  rownames_to_column(var = "gene") %>%
  as.tibble() %>%
  arrange(kruskal_pval) %>%
  inner_join(skyblue[, c("genes_ens", "geneName", "geneType", "chrom")], by=c("gene" = "genes_ens"))


skyblue_genes_expr <- stomach_males_fpkm %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "sample", value = "fpkm", -gene) %>%
  as.tibble() %>%
  inner_join(tumourME_males[, c("sample", "histological_type")], by = "sample") %>%
  inner_join(skyblue_genes[, c("gene", "geneName", "kruskal_pval")], by = "gene") %>%
  arrange(kruskal_pval)





# plot boxplot for top 5 genes
stomach_males_fpkm_plot <- skyblue_genes_expr %>%
  filter(geneName %in% c("HSP90AB1", "XPO5", "POLR1C", "MAD2L1BP", "MED20")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(geneName = fct_reorder(geneName, kruskal_pval)) %>%
  ggplot(mapping = aes(x = histological_type, y = fpkm, fill = histological_type)) +
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap( ~ geneName) +
  scale_fill_brewer(type = "qual", palette = "Set2", name = "Histological type") +
  theme_classic() +
  theme(
  axis.title = element_text(colour="black", size=14),
  axis.text.y = element_text(colour="black", size=14),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  legend.text=element_text(colour="black", size=13),
  legend.title=element_text(colour="black", size=14),
  strip.background = element_blank(),
  strip.text = element_text(colour="black", size=14),
  legend.position = "bottom") +
  labs(x = "Histological type", y = "FPKM (log2)", title = "") +
  guides(fill=guide_legend(nrow=6))
ggsave(filename="stomach_tumour_skyblue_hist_types.png", plot=stomach_males_fpkm_plot, path = "./plots/wgcna_networks_traits/", width=9, height=10)
unlink("stomach_tumour_skyblue_hist_types.png")



#
stomach_males_fpkm_plot2 <- skyblue_genes_expr %>%
  group_by(histological_type, gene, geneName) %>%
  summarise(fpkm = median(fpkm)) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = histological_type, y = fpkm, fill = histological_type)) +
  geom_boxplot() +
  #stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
  scale_fill_brewer(type = "qual", palette = "Set2", name = "Histological type") +
  theme_classic() +
  theme(
  axis.title = element_text(colour="black", size=14),
  axis.text.y = element_text(colour="black", size=14),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  legend.text=element_text(colour="black", size=13),
  legend.title=element_text(colour="black", size=14),
  strip.background = element_blank(),
  strip.text = element_text(colour="black", size=14),
  legend.position = "bottom") +
  labs(x = "Histological type", y = "Median log2 FPKM", title = "") +
  guides(fill=guide_legend(nrow=6))
ggsave(filename="stomach_tumour_skyblue_hist_types2.png", plot=stomach_males_fpkm_plot2, path = "./plots/wgcna_networks_traits/", width=8, height=6)
unlink("stomach_tumour_skyblue_hist_types2.png")



# plot boxplots for all hub genes
stomach_males_fpkm_plot3 <- skyblue_genes_expr %>%
  mutate_if(is.character, as.factor) %>%
  mutate(geneName = fct_reorder(geneName, kruskal_pval)) %>%
  ggplot(mapping = aes(x = histological_type, y = fpkm, fill = histological_type)) +
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap( ~ geneName) +
  scale_fill_brewer(type = "qual", palette = "Set2", name = "Histological type") +
  theme_classic() +
  theme(
  axis.title = element_text(colour="black", size=14),
  axis.text.y = element_text(colour="black", size=14),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  legend.text=element_text(colour="black", size=13),
  legend.title=element_text(colour="black", size=14),
  strip.background = element_blank(),
  strip.text = element_text(colour="black", size=14),
  legend.position = "bottom") +
  labs(x = "Histological type", y = "FPKM (log2)", title = "") +
  guides(fill=guide_legend(nrow=6))
ggsave(filename="stomach_tumour_skyblue_hist_types3.png", plot=stomach_males_fpkm_plot3, path = "./plots/wgcna_networks_traits/", width=18, height=20)
unlink("stomach_tumour_skyblue_hist_types3.png")



# hub genes heatmap
# samples colored by histological type
skyblue_genes_heatmap <- stomach_males_fpkm[skyblue$genes_ens, tumourME_males$sample]
hm_color_code <- c("Stomach Adenocarcinoma, Signet Ring Type" = "#66c2a5", "Stomach, Adenocarcinoma, Diffuse Type" = "#fc8d62", "Stomach, Intestinal Adenocarcinoma, Mucinous Type" = "#8da0cb", "Stomach, Intestinal Adenocarcinoma, Not Otherwise Specified (NOS)" = "#e78ac3", "Stomach, Intestinal Adenocarcinoma, Papillary Type" = "#a6d854", "Stomach, Intestinal Adenocarcinoma, Tubular Type" = "#ffd92f")
skyblue_samples_col <- as.character(hm_color_code[tumourME_males$histological_type])
pdf(file="./plots/wgcna_networks_traits/stomach_skyblue_genes_heatmap.pdf", w=10, h=6)
heatmap.2( x = as.matrix(skyblue_genes_heatmap),
  scale = "none",
  trace="none",
  #hclustfun=function(x) hclust(x, method="average"),
  Rowv=TRUE,
  Colv=TRUE,
  dendrogram = "both",
  cexRow=0.4, key=TRUE,
  keysize=1.2,
  symkey=TRUE,
  key.title = "FPKM",
  col=greenred(100),
  labCol=F,
  labRow = F,
  ColSideColors = skyblue_samples_col )
dev.off()





# sample clustering in each gene using FPKMs
# use gaussian mixture model

skyblue_genes_expr_cluster <- skyblue_genes_expr %>%
  group_by(gene, geneName) %>%
  mutate(class = as.character(mclust::Mclust(fpkm, G=2)$classification)) %>%
  ungroup()



# density plot for top 5 genes
skyblue_genes_expr_cluster_density <- skyblue_genes_expr_cluster %>%
  filter(geneName %in% c("HSP90AB1", "XPO5", "POLR1C", "MAD2L1BP", "MED20")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(geneName = fct_reorder(geneName, kruskal_pval)) %>%
  ggplot(mapping = aes(x = fpkm, fill = class, color = class)) +
  geom_density(alpha = 0.4) +
  facet_wrap( ~ geneName) +
  theme_classic() +
  theme(
  axis.title = element_text(colour="black", size=14),
  axis.text = element_text(colour="black", size=14),
  legend.text=element_text(colour="black", size=13),
  legend.title=element_text(colour="black", size=14),
  strip.background = element_blank(),
  strip.text = element_text(colour="black", size=14),
  legend.position = "bottom") +
  labs(x = "FPKM (log2)", title = "")
ggsave(filename="skyblue_genes_expr_cluster_density.png", plot=skyblue_genes_expr_cluster_density, path = "./plots/wgcna_networks_traits/", width=8, height=8)
unlink("skyblue_genes_expr_cluster_density.png")



# density plot for all hub genes
skyblue_genes_expr_cluster_density2 <- skyblue_genes_expr_cluster %>%
  mutate_if(is.character, as.factor) %>%
  mutate(geneName = fct_reorder(geneName, kruskal_pval)) %>%
  ggplot(mapping = aes(x = fpkm, fill = class, color = class)) +
  geom_density(alpha = 0.4) +
  facet_wrap( ~ geneName) +
  theme_classic() +
  theme(
  axis.title = element_text(colour="black", size=14),
  axis.text = element_text(colour="black", size=14),
  legend.text=element_text(colour="black", size=13),
  legend.title=element_text(colour="black", size=14),
  strip.background = element_blank(),
  strip.text = element_text(colour="black", size=14),
  legend.position = "bottom") +
  labs(x = "FPKM (log2)", title = "")
ggsave(filename="skyblue_genes_expr_cluster_density2.png", plot=skyblue_genes_expr_cluster_density2, path = "./plots/wgcna_networks_traits/", width=10, height=8)
unlink("skyblue_genes_expr_cluster_density2.png")




# Kaplan Meier Survival Curves

# HSP90AB1 gene

meta <- stomach_tcga_clinical %>%
  inner_join(skyblue_genes_expr_cluster %>% filter(geneName == "HSP90AB1") %>% dplyr::select(sample, geneName, class), by = c("sample"))

km_fit <- survfit(Surv(OS.time, OS) ~ class, data=meta)
km_pval <- survdiff(Surv(OS.time, OS) ~ class, data=meta)
autoplot(km_fit)
