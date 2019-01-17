# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between tumour and normal samples


library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)
library(viridis)




# -- Thyroid

# load datasets
diff_expr_males_table <- read_tsv("./files/diff_expr_thyroid_males_tcga_tumourVSnormal_edgeR.txt") %>%
  dplyr::rename(adj.P.Val = FDR)

diff_expr_females_table <- read_tsv("./files/diff_expr_thyroid_females_tcga_tumourVSnormal_edgeR.txt") %>%
  dplyr::rename(adj.P.Val = FDR)



# significant DEGs between males and females
males_signf_degs <- diff_expr_males_table %>%
    filter(adj.P.Val <= 0.05, abs(logFC) > 1)

females_signf_degs <- diff_expr_females_table %>%
    filter(adj.P.Val <= 0.05, abs(logFC) > 1)

males_signf_degs <- males_signf_degs %>%
    mutate(state = if_else(genes %in% females_signf_degs$genes, "common", "male_specific")) %>%
    dplyr::select(logFC, adj.P.Val, genes, geneName, geneType, chrom, state)

females_signf_degs <- females_signf_degs %>%
    mutate(state = if_else(genes %in% males_signf_degs$genes, "common", "female_specific")) %>%
    dplyr::select(logFC, adj.P.Val, genes, geneName, geneType, chrom, state)



# all significant DEGs between males and females
males_females_signf_degs <- full_join(
    males_signf_degs[, c("genes", "geneName", "chrom", "state", "logFC", "adj.P.Val")],
    females_signf_degs[, c("genes", "geneName", "chrom", "state", "logFC", "adj.P.Val")],
    by=c("genes","geneName","chrom","state")) %>%
    dplyr::rename(logFC_males = logFC.x, logFC_females = logFC.y, adj.P.Val_males = adj.P.Val.x, adj.P.Val_females = adj.P.Val.y)




#comparecluster
thca_compCluster <- list(up.tumour.male = NULL, up.normal.male = NULL, up.tumour.female = NULL, up.normal.female = NULL)
thca_compCluster$up.tumour.male <- males_signf_degs %>% filter(state == "male_specific", logFC > 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID)
thca_compCluster$up.normal.male <- males_signf_degs %>% filter(state == "male_specific", logFC < 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID)
thca_compCluster$up.tumour.female <- females_signf_degs %>% filter(state == "female_specific", logFC > 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID)
thca_compCluster$up.normal.female <- females_signf_degs %>% filter(state == "female_specific", logFC < 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID)

universe_enr <- diff_expr_males_table$geneName

thca_compCluster_enrichGO <- compareCluster(
  thca_compCluster,
  fun = "enrichGO",
  universe = bitr(universe_enr, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
  ont="BP",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE)

thca_compCluster_enrichGO_dot <- dotplot(thca_compCluster_enrichGO, font.size = 20, title = "Comparison of GO terms across gender-specific DEGs\nGO terms: BP")
ggsave(filename="thyroid_compareCluster_degs_GO_dot_edger.png", plot=thca_compCluster_enrichGO_dot, path = "./plots/diff_expression_tumourVSnormal/", width=17, height=8)
unlink("thyroid_compareCluster_degs_GO_dot_edger.png")
