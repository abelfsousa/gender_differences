# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Differential hub genes expression between histological types


# -- Thyroid


library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(gplots)
library(mclust)
library(survival)
library(data.table)
library(ranger)
library(ggfortify)
library(survminer)
library(mygene)

set.seed(123)


# load female tumour modules with kmE
females_tumour_modules <- read_tsv("./files/thyroid_females_tumour_modules.txt")


# load MEs and traits for female tumours
tumourME_females <- read_tsv("./files/thyroid_tumourME_traits_females.txt")


# load fpkm matrix
thyroid_females_fpkm <- read.table("./files/thyroid_females_tumour_fpkm_wgcna.txt", sep="\t", h=T)
thyroid_females_fpkm <- t(thyroid_females_fpkm)
colnames(thyroid_females_fpkm) <- substr(colnames(thyroid_females_fpkm), 1, 12)

thyroid_males_fpkm <- read.table("./files/thyroid_males_tumour_fpkm_wgcna.txt", sep="\t", h=T)
thyroid_males_fpkm <- t(thyroid_males_fpkm)
colnames(thyroid_males_fpkm) <- substr(colnames(thyroid_males_fpkm), 1, 12)



# load TCGA clinical data
thyroid_tcga_clinical <- fread("./data/tcga/tcga_clinical_data.txt") %>%
  as.tibble() %>%
  mutate(bcr_patient_barcode = str_replace_all(bcr_patient_barcode, "-", ".")) %>%
  filter(type == "THCA") %>%
  dplyr::select(sample = bcr_patient_barcode, gender, OS, OS.time, histological_type, ajcc_pathologic_tumor_stage) %>%
  filter(histological_type != "Other, specify" & ajcc_pathologic_tumor_stage != "[Not Available]") %>%
  filter(sample %in% tumourME_females$sample)

thyroid_tcga_clinical_males <- fread("./data/tcga/tcga_clinical_data.txt") %>%
  as.tibble() %>%
  mutate(bcr_patient_barcode = str_replace_all(bcr_patient_barcode, "-", ".")) %>%
  filter(type == "THCA") %>%
  dplyr::select(sample = bcr_patient_barcode, gender, OS, OS.time, histological_type, ajcc_pathologic_tumor_stage) %>%
  filter(histological_type != "Other, specify" & ajcc_pathologic_tumor_stage != "[Not Available]") %>%
  filter(sample %in% colnames(thyroid_males_fpkm))




# greenyellow
# female-specific tumour module correlated with histological type
# select hub genes (kme >= 0.8)
greenyel <- females_tumour_modules %>%
  filter(moduleL == "greenyellow", abs(kme) >= 0.8)


# select hub genes and samples used for correlation between eigenenes and histological type
# test for expression differences between histological types using a kruskal wallis test

#females
greenyel_genes <- thyroid_females_fpkm[greenyel$genes_ens, tumourME_females$sample]
greenyel_genes <- apply(greenyel_genes, 1, function(x) kruskal.test(as.numeric(x), as.factor(tumourME_females$histological_type))$p.value) %>%
  as.data.frame() %>%
  setNames(c("kruskal_pval")) %>%
  rownames_to_column(var = "gene") %>%
  as.tibble() %>%
  arrange(kruskal_pval) %>%
  inner_join(greenyel[, c("genes_ens", "geneName", "geneType", "chrom")], by=c("gene" = "genes_ens")) %>%
  dplyr::select(gene, geneName, geneType, chrom, kruskal_pval)


greenyel_genes_expr <- thyroid_females_fpkm %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "sample", value = "fpkm", -gene) %>%
  as.tibble() %>%
  inner_join(tumourME_females[, c("sample", "histological_type")], by = "sample") %>%
  inner_join(greenyel_genes[, c("gene", "geneName", "kruskal_pval")], by = "gene") %>%
  arrange(kruskal_pval)


greenyel_genes_expr_c <- greenyel_genes_expr %>%
  filter(geneName %in% c("SUCLG1", "GBAS", "MPC1", "GOT2", "SLC25A4")) %>%
  group_by(gene, geneName, kruskal_pval, histological_type) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(geneName = fct_reorder(geneName, kruskal_pval))


greenyel_genes_export <- greenyel_genes_expr %>%
  group_by(gene, geneName, histological_type, kruskal_pval) %>%
  summarise(n = n(), median_fpkm = median(fpkm)) %>%
  ungroup() %>%
  arrange(kruskal_pval)

greenyel_genes_annot <- queryMany(unique(greenyel_genes_export$gene), scopes='ensembl.gene', fields=c("symbol", "name", "summary"), species='human', return.as = "DataFrame") %>%
  as_tibble() %>%
  dplyr::select(query, name, summary)

greenyel_genes_export <- greenyel_genes_export %>%
  inner_join(greenyel_genes_annot, by = c("gene" = "query")) %>%
  dplyr::select(gene, geneName, name, summary, histological_type, kruskal_pval, n, median_fpkm)

write_tsv(greenyel_genes_export, "./files/thyroid_females_greenyellow_module_hub_genes_kruskal.txt")


#males
greenyel_genes_males <- thyroid_males_fpkm[greenyel$genes_ens, thyroid_tcga_clinical_males$sample]
greenyel_genes_males <- apply(greenyel_genes_males, 1, function(x) kruskal.test(as.numeric(x), as.factor(thyroid_tcga_clinical_males$histological_type))$p.value) %>%
  as.data.frame() %>%
  setNames(c("kruskal_pval")) %>%
  rownames_to_column(var = "gene") %>%
  as.tibble() %>%
  arrange(kruskal_pval) %>%
  inner_join(greenyel[, c("genes_ens", "geneName", "geneType", "chrom")], by=c("gene" = "genes_ens"))


greenyel_genes_males_expr <- thyroid_males_fpkm %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "sample", value = "fpkm", -gene) %>%
  as.tibble() %>%
  inner_join(thyroid_tcga_clinical_males[, c("sample", "histological_type")], by = "sample") %>%
  inner_join(greenyel_genes_males[, c("gene", "geneName", "kruskal_pval")], by = "gene") %>%
  arrange(kruskal_pval)


greenyel_genes_males_expr_c <- greenyel_genes_males_expr %>%
  filter(geneName %in% c("SUCLG1", "GBAS", "MPC1", "GOT2", "SLC25A4")) %>%
  group_by(gene, geneName, kruskal_pval, histological_type) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(geneName = fct_reorder(geneName, kruskal_pval))



# compare p-value distribution between genders
p_dis_genders <- bind_rows(
  greenyel_genes %>% dplyr::select(geneName, kruskal_pval) %>% mutate(gender = "females"),
  greenyel_genes_males %>% dplyr::select(geneName, kruskal_pval) %>% mutate(gender = "males")) %>%
  ggplot(mapping = aes(x = gender, y = -log10(kruskal_pval), fill = gender)) +
  geom_boxplot() +
  stat_compare_means(size = 3, label.y = 16) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(colour="black", size=14),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black", size=12),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.text=element_text(colour="black", size=10),
    legend.title=element_text(colour="black", size=14),
    legend.position = "bottom") +
  scale_fill_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  labs(x = "", y = "P-value (-log10)", title = "") +
  guides(fill=guide_legend(nrow=2))
ggsave(filename="thyroid_p_dist_genders.png", plot=p_dis_genders, path = "./plots/wgcna_networks_traits/", width=2, height=5)
ggsave(filename="thyroid_p_dist_genders.pdf", plot=p_dis_genders, path = "./plots/wgcna_networks_traits/", width=2, height=5)
unlink("thyroid_p_dist_genders.png")
unlink("thyroid_p_dist_genders.pdf")




# plot boxplots for top 5 genes

#females
thyroid_females_fpkm_plot <- greenyel_genes_expr %>%
  filter(geneName %in% c("SUCLG1", "GBAS", "MPC1", "GOT2", "SLC25A4")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(geneName = fct_reorder(geneName, kruskal_pval)) %>%
  ggplot(mapping = aes(x = histological_type, y = fpkm, fill = histological_type)) +
  geom_boxplot() +
  geom_text(data=greenyel_genes_expr_c, mapping=aes(x = histological_type, y = 1, label = counts)) +
  stat_compare_means(label.y = 8) +
  facet_wrap( ~ geneName, ncol = 5, nrow = 1) +
  scale_fill_brewer(type = "qual", palette = "Set2", name = "Histological type") +
  theme_classic() +
  theme(
    axis.title = element_text(colour="black", size=14),
    axis.text.y = element_text(colour="black", size=14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text=element_text(colour="black", size=11),
    legend.title=element_text(colour="black", size=14),
    strip.background = element_blank(),
    strip.text = element_text(colour="black", size=14),
    legend.position = "bottom") +
  labs(x = "Histological type", y = "FPKM (log2)", title = "") +
  guides(fill=guide_legend(ncol=2))
ggsave(filename="thyroid_tumour_greenyellow_hist_types.png", plot=thyroid_females_fpkm_plot, path = "./plots/wgcna_networks_traits/", width=15, height=5)
ggsave(filename="thyroid_tumour_greenyellow_hist_types.pdf", plot=thyroid_females_fpkm_plot, path = "./plots/wgcna_networks_traits/", width=15, height=5)
unlink("thyroid_tumour_greenyellow_hist_types.png")
unlink("thyroid_tumour_greenyellow_hist_types.pdf")


#males
thyroid_males_fpkm_plot <- greenyel_genes_males_expr %>%
  filter(geneName %in% c("SUCLG1", "GBAS", "MPC1", "GOT2", "SLC25A4")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(geneName = fct_reorder(geneName, kruskal_pval)) %>%
  ggplot(mapping = aes(x = histological_type, y = fpkm, fill = histological_type)) +
  geom_boxplot() +
  geom_text(data=greenyel_genes_males_expr_c, mapping=aes(x = histological_type, y = 1, label = counts)) +
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
  guides(fill=guide_legend(nrow=3))
ggsave(filename="thyroid_tumour_greenyellow_hist_types_males.png", plot=thyroid_males_fpkm_plot, path = "./plots/wgcna_networks_traits/", width=8, height=8)
unlink("thyroid_tumour_greenyellow_hist_types_males.png")



# plot median fpkm by histological subtype

#females
thyroid_females_fpkm_plot2 <- greenyel_genes_expr %>%
  group_by(histological_type, gene, geneName) %>%
  summarise(fpkm = median(fpkm)) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = histological_type, y = fpkm, fill = histological_type)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
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
  guides(fill=guide_legend(nrow=3))
ggsave(filename="thyroid_tumour_greenyellow_hist_types2.png", plot=thyroid_females_fpkm_plot2, path = "./plots/wgcna_networks_traits/", width=6, height=6)
unlink("thyroid_tumour_greenyellow_hist_types2.png")



# plot boxplots for all hub genes

#females
thyroid_females_fpkm_plot3 <- greenyel_genes_expr %>%
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
  guides(fill=guide_legend(nrow=3))
ggsave(filename="thyroid_tumour_greenyellow_hist_types3.png", plot=thyroid_females_fpkm_plot3, path = "./plots/wgcna_networks_traits/", width=20, height=30)
unlink("thyroid_tumour_greenyellow_hist_types3.png")



# hub genes heatmap
# samples colored by histological type
greenyel_genes_heatmap <- thyroid_females_fpkm[greenyel$genes_ens, tumourME_females$sample]
greenyel_samples_col <- if_else(tumourME_females$histological_type == "Thyroid Papillary Carcinoma - Classical/usual", "#66c2a5", if_else(tumourME_females$histological_type == "Thyroid Papillary Carcinoma - Tall Cell (>= 50% tall cell features)", "#8da0cb", "#fc8d62"))
pdf(file="./plots/wgcna_networks_traits/thyroid_greenyel_genes_heatmap.pdf", w=10, h=6)
heatmap.2( x = as.matrix(greenyel_genes_heatmap),
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
  ColSideColors = greenyel_samples_col )
dev.off()





# sample clustering in each gene using FPKMs
# use gaussian mixture model

#females
greenyel_genes_expr_cluster_females <- greenyel_genes_expr %>%
  group_by(gene, geneName) %>%
  nest() %>%
  ungroup() %>%
  mutate(class = purrr::map(.x = data, .f = ~ mclust::Mclust(.x$fpkm, G=2)$classification )) %>%
  unnest()

#males
greenyel_genes_expr_cluster_males <- greenyel_genes_males_expr %>%
  group_by(gene, geneName) %>%
  nest() %>%
  ungroup() %>%
  mutate(class = purrr::map(.x = data, .f = ~ mclust::Mclust(.x$fpkm, G=2)$classification )) %>%
  unnest()



# density plot for top 5 genes
greenyel_genes_expr_cluster_density <- greenyel_genes_expr_cluster_females %>%
  filter(geneName %in% c("SUCLG1", "GBAS", "MPC1", "GOT2", "SLC25A4")) %>%
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
ggsave(filename="greenyel_genes_expr_cluster_density.png", plot=greenyel_genes_expr_cluster_density, path = "./plots/wgcna_networks_traits/", width=8, height=8)
unlink("greenyel_genes_expr_cluster_density.png")



# density plot for all hub genes
greenyel_genes_expr_cluster_density2 <- greenyel_genes_expr_cluster_females %>%
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
ggsave(filename="greenyel_genes_expr_cluster_density2.png", plot=greenyel_genes_expr_cluster_density2, path = "./plots/wgcna_networks_traits/", width=30, height=20)
unlink("greenyel_genes_expr_cluster_density2.png")




# Kaplan Meier Survival Curves

# GBAS gene

meta <- thyroid_tcga_clinical %>%
  inner_join(greenyel_genes_expr_cluster_females %>% filter(geneName == "GBAS") %>% dplyr::select(sample, geneName, class), by = c("sample"))

km_fit <- survfit(Surv(OS.time, OS) ~ class, data=meta)
km_pval <- survdiff(Surv(OS.time, OS) ~ class, data=meta)
autoplot(km_fit)
ggsurvplot(km_fit, data = meta, risk.table = F, conf.int = TRUE,
  pval = TRUE) + labs(title = "Women")

meta <- thyroid_tcga_clinical_males %>%
  inner_join(greenyel_genes_expr_cluster_males %>% filter(geneName == "GBAS") %>% dplyr::select(sample, geneName, class), by = c("sample"))

km_fit <- survfit(Surv(OS.time, OS) ~ class, data=meta)
km_pval <- survdiff(Surv(OS.time, OS) ~ class, data=meta)
autoplot(km_fit)
ggsurvplot(km_fit, data = meta, risk.table = F, conf.int = TRUE,
  pval = TRUE) + labs(title = "Men")



meta <- thyroid_tcga_clinical %>%
  inner_join(greenyel_genes_expr_cluster_females %>% filter(geneName == "GOT2") %>% dplyr::select(sample, geneName, class), by = c("sample"))

km_fit <- survfit(Surv(OS.time, OS) ~ class, data=meta)
km_pval <- survdiff(Surv(OS.time, OS) ~ class, data=meta)
autoplot(km_fit)


meta <- thyroid_tcga_clinical_males %>%
  inner_join(greenyel_genes_expr_cluster_males %>% filter(geneName == "GOT2") %>% dplyr::select(sample, geneName, class), by = c("sample"))

km_fit <- survfit(Surv(OS.time, OS) ~ class, data=meta)
km_pval <- survdiff(Surv(OS.time, OS) ~ class, data=meta)
autoplot(km_fit)
