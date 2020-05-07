# Understanding Gender Differential Susceptibility in Cancer


# Comparison of the methylation status of differentially expressed genes between sample groups


library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(cowplot)



# load datasets

# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))


# metadata
thyroid_tcga_meta <- read_tsv("./files/tcga_thca_meta.txt") %>%
  dplyr::rename(sample = submitter_id) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  filter(sample_type_id != "06")

stomach_tcga_meta <- read_tsv("./files/tcga_stad_meta.txt") %>%
  dplyr::rename(sample = submitter_id) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  filter(sample_type_id != "06")


# differentially expressed genes
thca_TvsN <- read_tsv("./files/thyroid_males_females_signf_degs.txt")
thca_MvsF <- read_tsv("./files/thyroid_tumour_normal_signf_degs.txt")

stad_TvsN <- read_tsv("./files/stomach_males_females_signf_degs.txt")
stad_MvsF <- read_tsv("./files/stomach_tumour_normal_signf_degs.txt")


# differential expression tables

# Thyroid - Tumour vs Normal

thca_diff_expr_TvsN_males <- read_tsv("./files/diff_expr_thyroid_males_tcga_tumourVSnormal_edgeR.txt")
thca_diff_expr_TvsN_females <- read_tsv("./files/diff_expr_thyroid_females_tcga_tumourVSnormal_edgeR.txt")

thca_TvsN_males_signf_degs <- thca_diff_expr_TvsN_males %>%
    filter(FDR <= 0.05, abs(logFC) > 1)

thca_TvsN_females_signf_degs <- thca_diff_expr_TvsN_females %>%
    filter(FDR <= 0.05, abs(logFC) > 1)


# Thyroid - Male vs Female in tumours
thca_diff_expr_MvsF_tumours <- read_tsv("./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_edgeR.txt")

thca_MvsF_tumours_signf_degs <- thca_diff_expr_MvsF_tumours %>%
    filter(FDR <= 0.05)


# Stomach - Male vs Female in tumours
stad_diff_expr_MvsF_tumours <- read_tsv("./files/diff_expr_stomach_tumour_tcga_maleVSfemale_edgeR.txt")

stad_MvsF_tumours_signf_degs <- stad_diff_expr_MvsF_tumours %>%
    filter(FDR <= 0.05)



# methylation data

# thyroid
thca_meth <- read_tsv("./data/tcga/firebrowser/gdac.broadinstitute.org_THCA.Methylation_Preprocess.Level_3.2016012800.0.0/THCA_methylation450_mean_beta_promoters.txt.gz") %>%
  rename_all(.funs = ~ str_sub(.x, 1, 15)) %>%
  pivot_longer(-Gene_Symbol, names_to = "sample", values_to = "beta") %>%
  rename(gene = Gene_Symbol) %>%
  filter(sample %in% thyroid_tcga_meta$sample) %>%
  filter(str_detect(gene, ";", negate = T))


# stomach
stad_meth <- read_tsv("./data/tcga/firebrowser/gdac.broadinstitute.org_STAD.Methylation_Preprocess.Level_3.2016012800.0.0/STAD_methylation450_mean_beta_promoters.txt.gz") %>%
  rename_all(.funs = ~ str_sub(.x, 1, 15)) %>%
  pivot_longer(-Gene_Symbol, names_to = "sample", values_to = "beta") %>%
  rename(gene = Gene_Symbol) %>%
  filter(sample %in% stomach_tcga_meta$sample) %>%
  filter(str_detect(gene, ";", negate = T))




# compare fold-change/FDR distributions between DEGs / SBGs (T vs N and M vs F)

wilcoxon_diff_methyl <- function(group1, group2, expression.data){

  #differential methylation analysis using Wilcoxon rank sum test (Mann-Whitney)
  #compares methylation differences between two experimental groups

  p.val <- c()
  log2FC <- c()

  for(gene in rownames(expression.data)){
    group1.data <- as.numeric( expression.data[ gene, colnames(expression.data) %in% group1 ] )
    group2.data <- as.numeric( expression.data[ gene, colnames(expression.data) %in% group2 ] )

    pvalue <- wilcox.test(group1.data, group2.data)$p.value
    fc <- log2( (median(group1.data, na.rm = T)+0.05) / (median(group2.data, na.rm = T)+0.05) )

    p.val <- c(p.val, pvalue)
    log2FC <- c(log2FC, fc)
  }

  FDR <- p.adjust(p.val, method="BH")
  gene.table <- cbind.data.frame(rownames(expression.data), log2FC, p.val, FDR)
  colnames(gene.table)[1] <- "genes"

  return(gene.table)
}

thca_meth_mat <- thca_meth %>%
  group_by(gene) %>%
  filter(sum(is.na(beta)) == 0) %>%
  ungroup() %>%
  pivot_wider(names_from="sample", values_from="beta") %>%
  column_to_rownames(var = "gene")

stad_meth_mat <- stad_meth %>%
  group_by(gene) %>%
  filter(sum(is.na(beta)) == 0) %>%
  ungroup() %>%
  pivot_wider(names_from="sample", values_from="beta") %>%
  column_to_rownames(var = "gene")


# Thyroid T vs N females
females <- thyroid_tcga_meta %>%
  filter(gender == "female")

thca_TvsN_females <- wilcoxon_diff_methyl(
  females[females$sample_type_id == "01", "sample"]$sample,
  females[females$sample_type_id == "11", "sample"]$sample,
  thca_meth_mat)

# Thyroid T vs N males
males <- thyroid_tcga_meta %>%
  filter(gender == "male")

thca_TvsN_males <- wilcoxon_diff_methyl(
  males[males$sample_type_id == "01", "sample"]$sample,
  males[males$sample_type_id == "11", "sample"]$sample,
  thca_meth_mat)


# Thyroid M vs F tumour
tumour <- thyroid_tcga_meta %>%
  filter(sample_type_id == "01")

thca_MvsF_tumour <- wilcoxon_diff_methyl(
  tumour[tumour$gender == "male", "sample"]$sample,
  tumour[tumour$gender == "female", "sample"]$sample,
  thca_meth_mat)

# Stomach M vs F tumour
tumour <- stomach_tcga_meta %>%
  filter(sample_type_id == "01")

stad_MvsF_tumour <- wilcoxon_diff_methyl(
  tumour[tumour$gender == "male", "sample"]$sample,
  tumour[tumour$gender == "female", "sample"]$sample,
  stad_meth_mat)




# volcano plots

# Thyroid - T vs N in females

volcano <- thca_TvsN_females %>%
  as_tibble() %>%
  semi_join(thca_TvsN_females_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = log2FC, y = -log10(FDR), color = diff_methyl)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", lwd = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", lwd = 0.3) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 14, hjust = 0.5)) +
  scale_colour_manual(values = c("grey60", "red"), labels = c("No", "Yes"), name = "Differentially\nmethylated") +
  scale_x_continuous(limits = c(-2,2)) +
  labs(x = "Fold-change (log2)", y = "FDR (-log10)", title = "Thyroid females - Tumour Vs Normal\nDifferential methylation of DEGs")

ggsave(filename="01_thyroid_TvsN_females_degs_methylation_state_volcano.png", plot=volcano, path = "./plots/methylation_analysis_differential/", width=5, height=5)
unlink("01_thyroid_TvsN_females_degs_methylation_state_volcano.png")


barplot1 <- thca_TvsN_females %>%
  as_tibble() %>%
  semi_join(thca_TvsN_females_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = factor(0), fill = diff_methyl)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("grey60", "red"), guide=F) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color="black")) +
  labs(y = "Percentage")

ggsave(filename="02_thyroid_TvsN_females_degs_methylation_state_barplotP.png", plot=barplot1, path = "./plots/methylation_analysis_differential/", width=3, height=3)
unlink("02_thyroid_TvsN_females_degs_methylation_state_barplotP.png")


barplot2 <- thca_TvsN_females %>%
  as_tibble() %>%
  semi_join(thca_TvsN_females_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = factor(0), fill = diff_methyl)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("grey60", "red"), guide = F) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color="black")) +
  labs(y = "Count")

ggsave(filename="03_thyroid_TvsN_females_degs_methylation_state_barplotC.png", plot=barplot2, path = "./plots/methylation_analysis_differential/", width=3, height=3)
unlink("03_thyroid_TvsN_females_degs_methylation_state_barplotC.png")


volcano2 <- ggdraw(volcano) +
  draw_plot(barplot1, .55, .65, .25, .25) +
  draw_plot(barplot2, .55, .35, .25, .25)

ggsave(filename="thyroid_TvsN_females_degs_methylation_state_volcano.pdf", plot=volcano2, path = "./plots/methylation_analysis_differential/", width=5, height=5)
unlink("thyroid_TvsN_females_degs_methylation_state_volcano.pdf")



# Thyroid - T vs N in males

volcano <- thca_TvsN_males %>%
  as_tibble() %>%
  semi_join(thca_TvsN_males_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = log2FC, y = -log10(FDR), color = diff_methyl)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", lwd = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", lwd = 0.3) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 14, hjust = 0.5)) +
  scale_colour_manual(values = c("grey60", "red"), labels = c("No", "Yes"), name = "Differentially\nmethylated") +
  scale_x_continuous(limits = c(-2,2)) +
  labs(x = "Fold-change (log2)", y = "FDR (-log10)", title = "Thyroid males - Tumour Vs Normal\nDifferential methylation of DEGs")

ggsave(filename="04_thyroid_TvsN_males_degs_methylation_state_volcano.png", plot=volcano, path = "./plots/methylation_analysis_differential/", width=5, height=5)
unlink("04_thyroid_TvsN_males_degs_methylation_state_volcano.png")


barplot1 <- thca_TvsN_males %>%
  as_tibble() %>%
  semi_join(thca_TvsN_males_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = factor(0), fill = diff_methyl)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("grey60", "red"), guide=F) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color="black")) +
  labs(y = "Percentage")

ggsave(filename="05_thyroid_TvsN_males_degs_methylation_state_barplotP.png", plot=barplot1, path = "./plots/methylation_analysis_differential/", width=3, height=3)
unlink("05_thyroid_TvsN_males_degs_methylation_state_barplotP.png")


barplot2 <- thca_TvsN_males %>%
  as_tibble() %>%
  semi_join(thca_TvsN_males_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = factor(0), fill = diff_methyl)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("grey60", "red"), guide = F) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color="black")) +
  labs(y = "Count")

ggsave(filename="06_thyroid_TvsN_males_degs_methylation_state_barplotC.png", plot=barplot2, path = "./plots/methylation_analysis_differential/", width=3, height=3)
unlink("06_thyroid_TvsN_males_degs_methylation_state_barplotC.png")


volcano2 <- ggdraw(volcano) +
  draw_plot(barplot1, .55, .65, .25, .25) +
  draw_plot(barplot2, .55, .35, .25, .25)

ggsave(filename="thyroid_TvsN_males_degs_methylation_state_volcano.pdf", plot=volcano2, path = "./plots/methylation_analysis_differential/", width=5, height=5)
unlink("thyroid_TvsN_males_degs_methylation_state_volcano.pdf")



# Thyroid - M vs F in tumours

volcano <- thca_MvsF_tumour %>%
  as_tibble() %>%
  semi_join(thca_MvsF_tumours_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = log2FC, y = -log10(FDR), color = diff_methyl)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", lwd = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", lwd = 0.3) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 14, hjust = 0.5)) +
  scale_colour_manual(values = c("grey60", "red"), labels = c("No", "Yes"), name = "Differentially\nmethylated") +
  scale_x_continuous(limits = c(-2,2)) +
  labs(x = "Fold-change (log2)", y = "FDR (-log10)", title = "Thyroid tumours - Male Vs Female\nDifferential methylation of DEGs")

ggsave(filename="07_thyroid_MvsF_tumours_degs_methylation_state_volcano.png", plot=volcano, path = "./plots/methylation_analysis_differential/", width=5, height=5)
unlink("07_thyroid_MvsF_tumours_degs_methylation_state_volcano.png")


barplot1 <- thca_MvsF_tumour %>%
  as_tibble() %>%
  semi_join(thca_MvsF_tumours_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = factor(0), fill = diff_methyl)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("grey60", "red"), guide=F) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color="black")) +
  labs(y = "Percentage")

ggsave(filename="08_thyroid_MvsF_tumours_degs_methylation_state_barplotP.png", plot=barplot1, path = "./plots/methylation_analysis_differential/", width=3, height=3)
unlink("08_thyroid_MvsF_tumours_degs_methylation_state_barplotP.png")


barplot2 <- thca_MvsF_tumour %>%
  as_tibble() %>%
  semi_join(thca_MvsF_tumours_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = factor(0), fill = diff_methyl)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("grey60", "red"), guide = F) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color="black")) +
  labs(y = "Count")

ggsave(filename="09_thyroid_MvsF_tumours_degs_methylation_state_barplotC.png", plot=barplot2, path = "./plots/methylation_analysis_differential/", width=3, height=3)
unlink("09_thyroid_MvsF_tumours_degs_methylation_state_barplotC.png")


volcano2 <- ggdraw(volcano) +
  draw_plot(barplot1, .55, .65, .25, .25) +
  draw_plot(barplot2, .55, .35, .25, .25)

ggsave(filename="thyroid_MvsF_tumours_degs_methylation_state_volcano.pdf", plot=volcano2, path = "./plots/methylation_analysis_differential/", width=5, height=5)
unlink("thyroid_MvsF_tumours_degs_methylation_state_volcano.pdf")




# Stomach - M vs F in tumours

volcano <- stad_MvsF_tumour %>%
  as_tibble() %>%
  semi_join(stad_MvsF_tumours_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = log2FC, y = -log10(FDR), color = diff_methyl)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", lwd = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", lwd = 0.3) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 14, hjust = 0.5)) +
  scale_colour_manual(values = c("grey60", "red"), labels = c("No", "Yes"), name = "Differentially\nmethylated") +
  scale_x_continuous(limits = c(-2,2)) +
  labs(x = "Fold-change (log2)", y = "FDR (-log10)", title = "Stomach tumours - Male Vs Female\nDifferential methylation of DEGs")

ggsave(filename="10_stomach_MvsF_tumours_degs_methylation_state_volcano.png", plot=volcano, path = "./plots/methylation_analysis_differential/", width=5, height=5)
unlink("10_stomach_MvsF_tumours_degs_methylation_state_volcano.png")


barplot1 <- stad_MvsF_tumour %>%
  as_tibble() %>%
  semi_join(stad_MvsF_tumours_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = factor(0), fill = diff_methyl)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("grey60", "red"), guide=F) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color="black")) +
    labs(y = "Percentage")

ggsave(filename="11_stomach_MvsF_tumours_degs_methylation_state_barplotP.png", plot=barplot1, path = "./plots/methylation_analysis_differential/", width=3, height=3)
unlink("11_stomach_MvsF_tumours_degs_methylation_state_barplotP.png")


barplot2 <- stad_MvsF_tumour %>%
  as_tibble() %>%
  semi_join(stad_MvsF_tumours_signf_degs, by = c("genes" = "geneName")) %>%
  mutate(diff_methyl = if_else(FDR < 0.05, "1", "0")) %>%
  ggplot(mapping = aes(x = factor(0), fill = diff_methyl)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("grey60", "red"), guide = F) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color="black")) +
  labs(y = "Count")

ggsave(filename="12_stomach_MvsF_tumours_degs_methylation_state_barplotC.png", plot=barplot2, path = "./plots/methylation_analysis_differential/", width=3, height=3)
unlink("12_stomach_MvsF_tumours_degs_methylation_state_barplotC.png")


volcano2 <- ggdraw(volcano) +
  draw_plot(barplot1, .55, .65, .25, .25) +
  draw_plot(barplot2, .55, .35, .25, .25)

ggsave(filename="stomach_MvsF_tumours_degs_methylation_state_volcano.pdf", plot=volcano2, path = "./plots/methylation_analysis_differential/", width=5, height=5)
unlink("stomach_MvsF_tumours_degs_methylation_state_volcano.pdf")






# comparison <- thca_TvsN_females %>%
#   as_tibble() %>%
#   left_join(thca_TvsN[, c("geneName", "state")], by = c("genes" = "geneName")) %>%
#   mutate(state = replace_na(state, "background")) %>%
#   filter(state != "male_specific") %>%
#   ggplot(mapping = aes(x = state, y = -log10(FDR), fill = state)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.8) +
#   stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3))) +
#   scale_fill_brewer(type = "div", palette = "BrBG") +
#   theme(
#     axis.text = element_text(size = 8, color = "black")) +
#   labs(title = "Thyroid - Females")
# ggsave(filename="thca_TvsN_females_fdr_distribution.png", plot=comparison, height=6, width=4, path = "./plots/methylation_analysis2/")
# unlink("thca_TvsN_females_fdr_distribution.png")
#
#
# comparison <- thca_TvsN_males %>%
#   as_tibble() %>%
#   left_join(thca_TvsN[, c("geneName", "state")], by = c("genes" = "geneName")) %>%
#   mutate(state = replace_na(state, "background")) %>%
#   filter(state != "female_specific") %>%
#   ggplot(mapping = aes(x = state, y = -log10(FDR), fill = state)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.8) +
#   stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3))) +
#   scale_fill_brewer(type = "div", palette = "BrBG") +
#   theme(
#     axis.text = element_text(size = 8, color = "black")) +
#   labs(title = "Thyroid - Males")
# ggsave(filename="thca_TvsN_males_fdr_distribution.png", plot=comparison, height=6, width=4, path = "./plots/methylation_analysis2/")
# unlink("thca_TvsN_males_fdr_distribution.png")
#
#
# comparison <- thca_MvsF_tumour %>%
#   as_tibble() %>%
#   left_join(thca_MvsF[, c("geneName", "state")], by = c("genes" = "geneName")) %>%
#   mutate(state = replace_na(state, "background")) %>%
#   filter(state != "normal_specific") %>%
#   ggplot(mapping = aes(x = state, y = -log10(FDR), fill = state)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.8) +
#   stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3))) +
#   scale_fill_brewer(type = "div", palette = "BrBG") +
#   theme(
#       axis.text = element_text(size = 8, color = "black")) +
#   labs(title = "Thyroid - Tumour")
# ggsave(filename="thca_MvsF_tumour_fdr_distribution.png", plot=comparison, height=6, width=4, path = "./plots/methylation_analysis2/")
# unlink("thca_MvsF_tumour_fdr_distribution.png")
#
#
# comparison <- stad_MvsF_tumour %>%
#   as_tibble() %>%
#   left_join(stad_MvsF[, c("geneName", "state")], by = c("genes" = "geneName")) %>%
#   mutate(state = replace_na(state, "background")) %>%
#   filter(state != "normal_specific") %>%
#   ggplot(mapping = aes(x = state, y = -log10(FDR), fill = state)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.8) +
#   stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3))) +
#   scale_fill_brewer(type = "div", palette = "BrBG") +
#   theme(
#     axis.text = element_text(size = 8, color = "black")) +
#   labs(title = "Stomach - Tumour")
# ggsave(filename="stad_MvsF_tumour_fdr_distribution.png", plot=comparison, height=6, width=4, path = "./plots/methylation_analysis2/")
# unlink("stad_MvsF_tumour_fdr_distribution.png")



save(list=ls(), file="r_workspaces/diff_methylation_analysis_promoters.RData")
