# Understanding Gender Differential Susceptibility in Cancer


# Comparison of DEGs between the normal tissue of males and females from TCGA and GTEx



library(tidyverse)
library(RColorBrewer)



# -- Thyroid



# load DEGs dataset
thyroid_tcga_degs_gender <- read_tsv("./files/diff_expr_thyroid_normal_gtex_maleVSfemale_edgeR.txt") %>%
  dplyr::rename(adj.P.Val = FDR)
thyroid_gtex_degs_gender <- read_tsv("./files/diff_expr_thyroid_normal_tcga_maleVSfemale_edgeR.txt") %>%
  dplyr::rename(adj.P.Val = FDR)



thyroid_tcga_gtex <- inner_join(thyroid_tcga_degs_gender[, c("geneName", "logFC", "adj.P.Val")], thyroid_gtex_degs_gender[, c("geneName", "logFC", "adj.P.Val")], by="geneName") %>%
	dplyr::rename(logFC_tcga = logFC.x, fdr_tcga = adj.P.Val.x, logFC_gtex = logFC.y, fdr_gtex = adj.P.Val.y) %>%
  mutate(fdr_tcga = -log10(fdr_tcga), fdr_gtex = -log10(fdr_gtex))


thyroid_tcga_gtex <- thyroid_tcga_gtex %>%
    filter(abs(logFC_tcga) < 2.5 & abs(logFC_gtex) < 2.5) %>%
    filter(fdr_gtex < 1, fdr_tcga < 10)


# scatterplot
thyroid_tcga_gtex_plot1 <- ggplot(data=thyroid_tcga_gtex, mapping=aes(x=logFC_gtex, y=logFC_tcga)) +
    geom_point() +
    geom_smooth(method = "lm", se=T, colour="blue") +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=15),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    labs(x = "log2FC GTEx", y = "log2FC TCGA", title=paste("Thyroid log2FC correlation:", round(cor(thyroid_tcga_gtex$logFC_tcga, thyroid_tcga_gtex$logFC_gtex),4)))
ggsave(filename="thyroid_tcga_gtex_degs_normal_cor_logFC.png", plot=thyroid_tcga_gtex_plot1, path = "./plots/maleVSfemale_tcga_gtex/")
unlink("thyroid_tcga_gtex_degs_normal_cor_logFC.png")

thyroid_tcga_gtex_plot2 <- ggplot(data=thyroid_tcga_gtex, mapping=aes(x=fdr_gtex, y=fdr_tcga)) +
    geom_point() +
    geom_smooth(method = "lm", se=T, colour="blue") +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=15),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    labs(x = "-log10 FDR GTEx", y = "-log10 FDR TCGA", title=paste("Thyroid -log10 FDR correlation:", round(cor(thyroid_tcga_gtex$fdr_tcga, thyroid_tcga_gtex$fdr_gtex),4)))
ggsave(filename="thyroid_tcga_gtex_degs_normal_cor_log10fdr.png", plot=thyroid_tcga_gtex_plot2, path = "./plots/maleVSfemale_tcga_gtex/")
unlink("thyroid_tcga_gtex_degs_normal_cor_log10fdr.png")






# -- Stomach



# load DEGs dataset
stomach_tcga_degs_gender <- read_tsv("./files/diff_expr_stomach_normal_gtex_maleVSfemale_edgeR.txt") %>%
    dplyr::rename(adj.P.Val = FDR)
stomach_gtex_degs_gender <- read_tsv("./files/diff_expr_stomach_normal_tcga_maleVSfemale_edgeR.txt") %>%
    dplyr::rename(adj.P.Val = FDR)



stomach_tcga_gtex <- inner_join(stomach_tcga_degs_gender[, c("geneName", "logFC", "adj.P.Val")], stomach_gtex_degs_gender[, c("geneName", "logFC", "adj.P.Val")], by="geneName") %>%
	dplyr::rename(logFC_tcga = logFC.x, fdr_tcga = adj.P.Val.x, logFC_gtex = logFC.y, fdr_gtex = adj.P.Val.y) %>%
  mutate(fdr_tcga = -log10(fdr_tcga), fdr_gtex = -log10(fdr_gtex))


stomach_tcga_gtex <- stomach_tcga_gtex %>%
    filter(abs(logFC_gtex) < 5 & abs(logFC_tcga) < 5) %>%
    filter(fdr_gtex < 10, fdr_tcga < 10)


# scatterplot
stomach_tcga_gtex_plot1 <- ggplot(data=stomach_tcga_gtex, mapping=aes(x=logFC_gtex, y=logFC_tcga)) +
    geom_point() +
    geom_smooth(method = "lm", se=T, colour="blue") +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=15),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    labs(x = "log2FC GTEx", y = "log2FC TCGA", title=paste("stomach log2FC correlation:", round(cor(stomach_tcga_gtex$logFC_tcga, stomach_tcga_gtex$logFC_gtex),4)))
ggsave(filename="stomach_tcga_gtex_degs_normal_cor_logFC.png", plot=stomach_tcga_gtex_plot1, path = "./plots/maleVSfemale_tcga_gtex/")
unlink("stomach_tcga_gtex_degs_normal_cor_logFC.png")

stomach_tcga_gtex_plot2 <- ggplot(data=stomach_tcga_gtex, mapping=aes(x=fdr_gtex, y=fdr_tcga)) +
    geom_point() +
    geom_smooth(method = "lm", se=T, colour="blue") +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=15),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    labs(x = "-log10 FDR GTEx", y = "-log10 FDR TCGA", title=paste("stomach -log10 FDR correlation:", round(cor(stomach_tcga_gtex$fdr_tcga, stomach_tcga_gtex$fdr_gtex),4)))
ggsave(filename="stomach_tcga_gtex_degs_normal_cor_log10fdr.png", plot=stomach_tcga_gtex_plot2, path = "./plots/maleVSfemale_tcga_gtex/")
unlink("stomach_tcga_gtex_degs_normal_cor_log10fdr.png")
