# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between tumour and normal samples


library(tidyverse)
library(RColorBrewer)
library(VennDiagram)





# -- Thyroid

# load datasets
diff_expr_tumour_table <- read_tsv("./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_limma.txt") %>%
    #dplyr::rename(adj.P.Val = FDR) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))

diff_expr_normal_table <- read_tsv("./files/diff_expr_thyroid_normal_gtex_maleVSfemale_limma.txt") %>%
    #dplyr::rename(adj.P.Val = FDR) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))



# significant DEGs between males and females
tumour_signf_degs <- diff_expr_tumour_table %>%
    filter(adj.P.Val <= 0.05, abs(logFC) > 1)

normal_signf_degs <- diff_expr_normal_table %>%
    filter(adj.P.Val <= 0.05, abs(logFC) > 1)

tumour_signf_degs <- tumour_signf_degs %>%
    mutate(state = if_else(genes %in% normal_signf_degs$genes, "common", "tumour_specific")) %>%
    dplyr::select(logFC, adj.P.Val, genes, geneName, geneType, chrom, fdr, state)

normal_signf_degs <- normal_signf_degs %>%
    mutate(state = if_else(genes %in% tumour_signf_degs$genes, "common", "normal_specific")) %>%
    dplyr::select(logFC, adj.P.Val, genes, geneName, geneType, chrom, fdr, state)


# all significant DEGs between tumour and normal
tumour_normal_signf_degs <- full_join(
    tumour_signf_degs[, c("genes", "geneName", "chrom", "state", "logFC", "adj.P.Val")],
    normal_signf_degs[, c("genes", "geneName", "chrom", "state", "logFC", "adj.P.Val")],
    by=c("genes","geneName","chrom","state")) %>%
    dplyr::rename(log2FC_tumour = logFC.x, log2FC_normal = logFC.y, FDR_tumour = adj.P.Val.x, FDR_normal = adj.P.Val.y) %>%
    mutate(log2FC_tumour = round(log2FC_tumour, 2), log2FC_normal = round(log2FC_normal, 2), FDR_tumour = round(FDR_tumour, 2), FDR_normal = round(FDR_normal, 2))
write.table(tumour_normal_signf_degs, "./files/thyroid_tumour_normal_signf_degs.txt", sep="\t", quote=F, row.names=F)


# Venn Diagram
males_females_signf_degs_venn <- draw.pairwise.venn(
    sum(table(tumour_normal_signf_degs$state)[c("common", "normal_specific")]),
    sum(table(tumour_normal_signf_degs$state)[c("common", "tumour_specific")]),
    sum(table(tumour_normal_signf_degs$state)[c("common")]),
    category = c("", ""),
    lty = rep("blank", 2),
    fill = c("green", "red"),
    alpha = rep(0.5, 2),
    cex = rep(1.4, 3),
    cat.cex = rep(1.5, 2))
ggsave(filename="thyroid_tumour_normal_signf_degs_venn.png", plot=males_females_signf_degs_venn, path = "./plots/diff_expression_maleVSfemale/", width=3, height=3)
unlink("thyroid_tumour_normal_signf_degs_venn.png")





save(list=ls(), file="r_workspaces/tcga_gtex_thyroid_diffExpr_malesVSfemales_characterization.RData")

