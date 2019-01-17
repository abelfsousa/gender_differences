# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between tumour and normal samples


library(tidyverse)
library(RColorBrewer)
library(VennDiagram)





# -- Stomach

# load datasets
diff_expr_tumour_table <- read_tsv("./files/diff_expr_stomach_tumour_tcga_maleVSfemale_limma.txt")
    #dplyr::rename(adj.P.Val = FDR) %>%
    #mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))

diff_expr_normal_table <- read_tsv("./files/diff_expr_stomach_normal_tcga_maleVSfemale_limma.txt")
    #dplyr::rename(adj.P.Val = FDR) %>%
    #mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))



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
write.table(tumour_normal_signf_degs, "./files/stomach_tumour_normal_signf_degs2.txt", sep="\t", quote=F, row.names=F)



# Venn Diagram
males_females_signf_degs_venn <- draw.pairwise.venn(
    sum(table(tumour_normal_signf_degs$state)[c("common")]),
    sum(table(tumour_normal_signf_degs$state)[c("common", "tumour_specific")]),
    sum(table(tumour_normal_signf_degs$state)[c("common")]),
    category = c("", ""),
    lty = rep("blank", 2),
    fill = c("green", "red"),
    alpha = rep(0.5, 2),
    cex = rep(1.5, 3),
    cat.cex = rep(1.5, 2))
ggsave(filename="stomach_tumour_normal_signf_degs_venn2.png", plot=males_females_signf_degs_venn, path = "./plots/diff_expression_maleVSfemale/", width=3, height=3)
unlink("stomach_tumour_normal_signf_degs_venn2.png")



# append DEGs state
diff_expr_tumour_table <- diff_expr_tumour_table %>%
  left_join(tumour_signf_degs[,c("genes", "state")], by="genes")

diff_expr_normal_table <- diff_expr_normal_table %>%
  left_join(normal_signf_degs[,c("genes", "state")], by="genes")



# volcano plots
diff_expr_tumour_table_vp <- ggplot( data = diff_expr_tumour_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state, fill=fdr) ) +
  geom_point(pch = 21) +
  scale_fill_manual(values=rev(brewer.pal(name="OrRd", n=4)[2:4]), name = "Significance") +
  scale_colour_manual(values=c("#3182bd", "green"), na.value="#bdbdbd", labels=c("Common", "Tumour-specific", "Not DEGs"), name = "DEG type") +
  geom_line(aes(x=0), color="black", linetype=2, size = 0.3) +
  geom_line(aes(x=1), color="black", linetype=2, size = 0.3) +
  geom_line(aes(x=-1), color="black", linetype=2, size = 0.3) +
  geom_line(aes(y=-log10(0.05)), color="black", linetype=2, size = 0.3) +
  theme_classic() +
  theme(axis.title = element_text(colour="black", size=18),
    axis.text = element_text(colour="black", size=16),
    legend.text=element_text(colour="black", size=16),
    legend.title=element_text(colour="black", size=18),
    plot.title=element_text(colour="black", size=20)) +
  labs(x = "logFC", y = "FDR (-log10)", title="Tumour")
ggsave(filename="diff_expr_stomach_tumour_tcga_maleVSfemale.png", plot=diff_expr_tumour_table_vp, path = "./plots/diff_expression_maleVSfemale/", width=6, height=7)
unlink("diff_expr_stomach_tumour_tcga_maleVSfemale.png")

diff_expr_normal_table_vp <- ggplot( data = diff_expr_normal_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state, fill=fdr) ) +
    geom_point(pch = 21) +
    scale_fill_manual(values=rev(brewer.pal(name="OrRd", n=4)[2:4]), name = "Significance") +
    scale_colour_manual(values=c("#3182bd"), na.value="#bdbdbd", labels=c("Common", "Not DEGs"), name = "DEG type") +
    geom_line(aes(x=0), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=-1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(y=-log10(0.05)), color="black", linetype=2, size = 0.3) +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=18),
      axis.text = element_text(colour="black", size=16),
      legend.text=element_text(colour="black", size=16),
      legend.title=element_text(colour="black", size=18),
      plot.title=element_text(colour="black", size=20)) +
    labs(x = "logFC", y = "FDR (-log10)", title="Normal")
ggsave(filename="diff_expr_stomach_normal_tcga_maleVSfemale.png", plot=diff_expr_normal_table_vp, path = "./plots/diff_expression_maleVSfemale/", width=6, height=7)
unlink("diff_expr_stomach_normal_tcga_maleVSfemale.png")





save(list=ls(), file="r_workspaces/tcga_stomach_diffExpr_malesVSfemales_characterization.RData")
