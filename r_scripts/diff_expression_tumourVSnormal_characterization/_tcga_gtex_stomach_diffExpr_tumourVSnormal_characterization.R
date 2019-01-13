# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between tumour and normal samples


library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)




# -- Stomach

# load datasets
diff_expr_males_table <- read_tsv("./files/diff_expr_stomach_males_tcga_gtex_tumourVSnormal.txt")
diff_expr_females_table <- read_tsv("./files/diff_expr_stomach_females_tcga_gtex_tumourVSnormal.txt")



# significant DEGs between males and females
males_signf_degs <- diff_expr_males_table %>%
    filter(adj.P.Val <= 0.05, abs(logFC) > 1)

females_signf_degs <- diff_expr_females_table %>%
    filter(adj.P.Val <= 0.05, abs(logFC) > 1)

males_signf_degs <- males_signf_degs %>%
    mutate(state = if_else(genes %in% females_signf_degs$genes, "common", "male_specific")) %>%
    dplyr::select(logFC, adj.P.Val, genes, geneName, geneType, chrom, fdr, state)

females_signf_degs <- females_signf_degs %>%
    mutate(state = if_else(genes %in% males_signf_degs$genes, "common", "female_specific")) %>%
    dplyr::select(logFC, adj.P.Val, genes, geneName, geneType, chrom, fdr, state)



# all significant DEGs between males and females
males_females_signf_degs <- full_join(
    males_signf_degs[, c("genes", "geneName", "chrom", "state", "logFC", "adj.P.Val")],
    females_signf_degs[, c("genes", "geneName", "chrom", "state", "logFC", "adj.P.Val")],
    by=c("genes","geneName","chrom","state"))

males_females_signf_degs <- males_females_signf_degs %>%
    dplyr::rename(logFC_males = logFC.x, logFC_females = logFC.y, adj.P.Val_males = adj.P.Val.x, adj.P.Val_females = adj.P.Val.y)
write.table(males_females_signf_degs, "./files/stomach_males_females_signf_degs.txt", sep="\t", quote=F, row.names=F)



# barplot of all significant DEGs between males and females: common, male-specific, female-specific
males_females_signf_degs_bp <- ggplot(data=males_females_signf_degs, mapping=aes(x=state, y=..count.., fill=state)) +
    geom_bar() +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=15),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    scale_x_discrete(labels = c("Common", "Female-specific", "Male-specific")) +
    labs(x = "", y = "Number of DEGs", title="stomach")
ggsave(filename="stomach_males_females_signf_degs_bp.png", plot=males_females_signf_degs_bp, path = "./plots/")
unlink("stomach_males_females_signf_degs_bp.png")



# append DEGs state
diff_expr_males_table <- diff_expr_males_table %>%
    left_join(males_signf_degs[,c("genes", "state")], by="genes")

diff_expr_females_table <- diff_expr_females_table %>%
    left_join(females_signf_degs[,c("genes", "state")], by="genes")



# volcano plots
diff_expr_males_vp <- ggplot( data = diff_expr_males_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state, fill=fdr) ) +
    geom_point(pch = 21) +
    scale_fill_manual(values=rev(brewer.pal(name="OrRd", n=4)[2:4])) +
    scale_colour_manual(values=c("#91bfdb", "green"), na.value="#FDCC8A", labels=c("Common", "Male-specific", "")) +
    geom_line(aes(x=0), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=-1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(y=-log10(0.05)), color="black", linetype=2, size = 0.3) +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text = element_text(colour="black", size=13),
        legend.text=element_text(colour="black", size=13),
        legend.title=element_text(colour="black", size=15),
        plot.title=element_text(colour="black", size=15)) +
    labs(x = "logFC", y = "FDR (-log10)", title="stomach - Males")
ggsave(filename="diff_expr_stomach_males_tcga_gtex_tumourVSnormal.png", plot=diff_expr_males_vp, path = "./plots/")
unlink("diff_expr_stomach_males_tcga_gtex_tumourVSnormal.png")

diff_expr_females_vp <- ggplot( data = diff_expr_females_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state, fill=fdr) ) +
    geom_point(pch = 21) +
    scale_fill_manual(values=rev(brewer.pal(name="OrRd", n=4)[2:4])) +
    scale_colour_manual(values=c("#91bfdb", "green"), na.value="#FDCC8A", labels=c("Common", "Female-specific", "")) +
    geom_line(aes(x=0), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=-1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(y=-log10(0.05)), color="black", linetype=2, size = 0.3) +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text = element_text(colour="black", size=13),
        legend.text=element_text(colour="black", size=13),
        legend.title=element_text(colour="black", size=15),
        plot.title=element_text(colour="black", size=15)) +
    labs(x = "logFC", y = "FDR (-log10)", title="stomach - Females")
ggsave(filename="diff_expr_stomach_females_tcga_gtex_tumourVSnormal.png", plot=diff_expr_females_vp, path = "./plots/")
unlink("diff_expr_stomach_females_tcga_gtex_tumourVSnormal.png")






# hypergeometric test of GO terms and KEGG pathways

# common DEGs: GO terms
common_degs_go <- enrichGO(
    gene = males_females_signf_degs[(males_females_signf_degs$state == "common"), "genes"]$genes,
    universe = unique(diff_expr_males_table$genes, diff_expr_females_table$genes),
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    pool=TRUE)


common_degs_go_dot <- dotplot(common_degs_go, x = "GeneRatio", color = "p.adjust", showCategory = 10, split = "ONTOLOGY", font.size = 12, title = "Common DEGs (Tumour vs Normal)\nGO terms: BP, CC, MF")
ggsave(filename="stomach_common_degs_go_dot.png", plot=common_degs_go_dot, path = "./plots/", width=8)
unlink("stomach_common_degs_go_dot.png")



# common DEGs: KEGG pathways
common_degs_keggPath <- enrichKEGG(
    gene=bitr(males_females_signf_degs[(males_females_signf_degs$state == "common"), "genes"]$genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
    organism = "hsa",
    keyType = "ncbi-geneid",
    universe = bitr(unique(diff_expr_males_table$genes, diff_expr_females_table$genes), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    use_internal_data = FALSE)


common_degs_keggPath_dot <- dotplot(common_degs_keggPath, x = "GeneRatio", color = "p.adjust", showCategory = 10, font.size = 12, title = "Common DEGs (Tumour vs Normal)\nKEGG pathways")
ggsave(filename="stomach_common_degs_keggPath_dot.png", plot=common_degs_keggPath_dot, path = "./plots/", width=8)
unlink("stomach_common_degs_keggPath_dot.png")





# female-specific DEGs: GO terms
femaleSpecific_degs_go <- enrichGO(
    gene = males_females_signf_degs[(males_females_signf_degs$state == "female_specific"), "genes"]$genes,
    universe = unique(diff_expr_males_table$genes, diff_expr_females_table$genes),
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    pool=TRUE)


femaleSpecific_degs_go_dot <- dotplot(femaleSpecific_degs_go, x = "GeneRatio", color = "p.adjust", showCategory = 10, split = "ONTOLOGY", font.size = 12, title = "Female - specific DEGs (Tumour vs Normal)\nGO terms: BP, CC, MF")
ggsave(filename="stomach_femaleSpecific_degs_go_dot.png", plot=femaleSpecific_degs_go_dot, path = "./plots/", width=11)
unlink("stomach_femaleSpecific_degs_go_dot.png")



# female-specific DEGs: KEGG pathways
femaleSpecific_degs_keggPath <- enrichKEGG(
    gene=bitr(males_females_signf_degs[(males_females_signf_degs$state == "female_specific"), "genes"]$genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
    organism = "hsa",
    keyType = "ncbi-geneid",
    universe = bitr(unique(diff_expr_males_table$genes, diff_expr_females_table$genes), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.2,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    use_internal_data = FALSE)


femaleSpecific_degs_keggPath_dot <- dotplot(femaleSpecific_degs_keggPath, x = "GeneRatio", color = "p.adjust", showCategory = 10, font.size = 12, title = "Female - specific DEGs (Tumour vs Normal)\nKEGG pathways")
ggsave(filename="stomach_femaleSpecific_degs_keggPath_dot.png", plot=femaleSpecific_degs_keggPath_dot, path = "./plots/", width=7)
unlink("stomach_femaleSpecific_degs_keggPath_dot.png")




# male-specific DEGs: GO terms
maleSpecific_degs_go <- enrichGO(
    gene = males_females_signf_degs[(males_females_signf_degs$state == "male_specific"), "genes"]$genes,
    universe = unique(diff_expr_males_table$genes, diff_expr_females_table$genes),
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    pool=TRUE)


maleSpecific_degs_go_dot <- dotplot(maleSpecific_degs_go, x = "GeneRatio", color = "p.adjust", showCategory = 10, split = "ONTOLOGY", font.size = 12, title = "Male - specific DEGs (Tumour vs Normal)\nGO terms: BP, CC, MF")
ggsave(filename="stomach_maleSpecific_degs_go_dot.png", plot=maleSpecific_degs_go_dot, path = "./plots/", width=11)
unlink("stomach_maleSpecific_degs_go_dot.png")



# female-specific DEGs: KEGG pathways
maleSpecific_degs_keggPath <- enrichKEGG(
    gene=bitr(males_females_signf_degs[(males_females_signf_degs$state == "male_specific"), "genes"]$genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
    organism = "hsa",
    keyType = "ncbi-geneid",
    universe = bitr(unique(diff_expr_males_table$genes, diff_expr_females_table$genes), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    use_internal_data = FALSE)


maleSpecific_degs_keggPath_dot <- dotplot(maleSpecific_degs_keggPath, x = "GeneRatio", color = "p.adjust", showCategory = 10, font.size = 12, title = "Male - specific DEGs (Tumour vs Normal)\nKEGG pathways")
ggsave(filename="stomach_maleSpecific_degs_keggPath_dot.png", plot=maleSpecific_degs_keggPath_dot, path = "./plots/", width=7)
unlink("stomach_maleSpecific_degs_keggPath_dot.png")






save(list=ls(), file="r_workspaces/tcga_gtex_stomach_diffExpr_tumourVSnormal_characterization.RData")




