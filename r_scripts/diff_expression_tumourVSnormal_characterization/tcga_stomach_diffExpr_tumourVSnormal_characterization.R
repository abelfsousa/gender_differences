# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between tumour and normal samples


library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)
library(viridis)




# -- Stomach

# load datasets
diff_expr_males_table <- read_tsv("./files/diff_expr_stomach_males_tcga_tumourVSnormal_limma.txt") %>%
    #dplyr::rename(adj.P.Val = FDR) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))

diff_expr_females_table <- read_tsv("./files/diff_expr_stomach_females_tcga_tumourVSnormal_limma.txt") %>%
    #dplyr::rename(adj.P.Val = FDR) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))



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
    by=c("genes","geneName","chrom","state")) %>%
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
ggsave(filename="stomach_males_females_signf_degs_bp.png", plot=males_females_signf_degs_bp, path = "./plots/diff_expression_tumourVSnormal/")
unlink("stomach_males_females_signf_degs_bp.png")




# Venn Diagram
males_females_signf_degs_venn <- draw.pairwise.venn(
    sum(table(males_females_signf_degs$state)[c("common", "male_specific")]),
    sum(table(males_females_signf_degs$state)[c("common", "female_specific")]),
    sum(table(males_females_signf_degs$state)[c("common")]),
    category = c("", ""),
    lty = rep("blank", 2),
    fill = c("light blue", "pink"),
    alpha = rep(0.5, 2),
    cex = rep(1.4, 3),
    cat.cex = rep(1.5, 2))
ggsave(filename="stomach_males_females_signf_degs_venn.png", plot=males_females_signf_degs_venn, path = "./plots/diff_expression_tumourVSnormal/", width=3, height=3)
unlink("stomach_males_females_signf_degs_venn.png")











# append DEGs state
diff_expr_males_table <- diff_expr_males_table %>%
    left_join(males_females_signf_degs[,c("genes", "state")], by="genes")

diff_expr_females_table <- diff_expr_females_table %>%
    left_join(males_females_signf_degs[,c("genes", "state")], by="genes")



# volcano plots
diff_expr_males_vp <- ggplot( data = diff_expr_males_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state, fill=fdr) ) +
    geom_point(pch = 21) +
    scale_fill_manual(values=rev(brewer.pal(name="OrRd", n=4)[2:4]), name = "Significance") +
    scale_colour_manual(values=c("#3182bd", "green", "yellow"), na.value="#bdbdbd", labels=c("Common", "Female-specific", "Male-specific", "Not DEGs"), name = "DEG type") +
    geom_line(aes(x=0), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=2), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=-1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=-2), color="black", linetype=2, size = 0.3) +
    geom_line(aes(y=-log10(0.05)), color="black", linetype=2, size = 0.3) +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=15),
      axis.text = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      plot.title=element_text(colour="black", size=15)) +
    labs(x = "logFC", y = "FDR (-log10)", title="Stomach - Males")
ggsave(filename="diff_expr_stomach_males_tcga_tumourVSnormal.png", plot=diff_expr_males_vp, path = "./plots/diff_expression_tumourVSnormal/", width=6, height=6)
unlink("diff_expr_stomach_males_tcga_tumourVSnormal.png")

diff_expr_females_vp <- ggplot( data = diff_expr_females_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state, fill=fdr) ) +
    geom_point(pch = 21) +
    scale_fill_manual(values=rev(brewer.pal(name="OrRd", n=4)[2:4]), name = "Significance") +
    scale_colour_manual(values=c("#3182bd", "green", "yellow"), na.value="#bdbdbd", labels=c("Common", "Female-specific", "Male-specific", "Not DEGs"), name = "DEG type") +
    geom_line(aes(x=0), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=2), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=-1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=-2), color="black", linetype=2, size = 0.3) +
    geom_line(aes(y=-log10(0.05)), color="black", linetype=2, size = 0.3) +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=15),
      axis.text = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      plot.title=element_text(colour="black", size=15)) +
    labs(x = "logFC", y = "FDR (-log10)", title="Stomach - Females")
ggsave(filename="diff_expr_stomach_females_tcga_tumourVSnormal.png", plot=diff_expr_females_vp, path = "./plots/diff_expression_tumourVSnormal/", width=6, height=6)
unlink("diff_expr_stomach_females_tcga_tumourVSnormal.png")






# hypergeometric test of GO terms and KEGG pathways

# enrichment function
enrich_test <- function(gene_set, universe, terms1, terms2, p_adj, q_value){
  enr <- enricher(
    gene = gene_set,
    universe = universe,
    TERM2GENE = terms1,
    TERM2NAME = terms2,
    pvalueCutoff = p_adj,
    qvalueCutoff = q_value,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500)
}

# load gene lists
kegg <- read.gmt("./data/gene_lists/c2.cp.kegg.v6.2.symbols.gmt") %>% mutate(ont = str_replace(ont, "KEGG_", ""))
kegg2 <- data.frame(ont = kegg$ont, name = "KEGG")
go_bp <- read.gmt("./data/gene_lists/c5.bp.v6.2.symbols.gmt") %>% mutate(ont = str_replace(ont, "GO_", ""))
go_bp2 <- data.frame(ont = go_bp$ont, name = "GO_BP")
go_mf <- read.gmt("./data/gene_lists/c5.mf.v6.2.symbols.gmt") %>% mutate(ont = str_replace(ont, "GO_", ""))
go_mf2 <- data.frame(ont = go_mf$ont, name = "GO_MF")
go_cc <- read.gmt("./data/gene_lists/c5.cc.v6.2.symbols.gmt") %>% mutate(ont = str_replace(ont, "GO_", ""))
go_cc2 <- data.frame(ont = go_cc$ont, name = "GO_CC")
onco <- read.gmt("./data/gene_lists/c6.all.v6.2.symbols.gmt")
onco2 <- data.frame(ont = onco$ont, name = "ONCO")
immuno <- read.gmt("./data/gene_lists/c7.all.v6.2.symbols.gmt")
immuno2 <- data.frame(ont = immuno$ont, name = "IMMUNO")

all_terms <- bind_rows(kegg, go_bp, go_mf, go_cc, onco, immuno) %>% mutate(ont = str_replace_all(ont, "_", " "))
all_terms2 <- bind_rows(kegg2, go_bp2, go_mf2, go_cc2, onco2, immuno2) %>% mutate(ont = str_replace_all(ont, "_", " "))


universe_enr <- diff_expr_males_table$geneName





# common DEGs
common_degs <- males_females_signf_degs %>%
  filter(state == "common") %>%
  pull(geneName)

common_degs_enr <- enrich_test(common_degs, universe_enr, all_terms, all_terms2, 1, 1)
common_degs_enr_tab <- common_degs_enr@result


common_degs_enr_barplot <- common_degs_enr_tab %>%
  as.tibble() %>%
  dplyr::select(ID, Description, p_adj=p.adjust) %>%
  filter(p_adj < 0.05) %>%
  mutate(p_adj = -log10(p_adj)) %>%
  group_by(Description) %>%
  top_n(n=10, wt=p_adj) %>%
  ungroup() %>%
  filter(!str_detect(Description, "GO_MF")) %>%
  filter(!str_detect(Description, "GO_CC")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, p_adj)) %>%
  ggplot(mapping = aes(x=ID, y = p_adj, fill = p_adj)) +
    geom_bar(stat="identity") +
    theme_classic() +
    facet_wrap( ~ Description, scales = "free") +
    theme(
      axis.title.x=element_text(colour="black", size=15),
      axis.title.y=element_blank(),
      axis.text.y=element_text(colour="black", size=10),
      axis.text.x=element_text(colour="black", size=12),
      plot.title = element_blank(),
      strip.text.x = element_text(size=12),
      strip.background = element_blank()) +
    coord_flip() +
    scale_fill_viridis(option="D", name="") +
    scale_y_continuous(name = "Adjusted P-value (-log10)")
ggsave(filename="stomach_common_degs_enr_barplot.png", plot=common_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 16, height = 8)
ggsave(filename="stomach_common_degs_enr_barplot.pdf", plot=common_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 16, height = 8)
unlink("stomach_common_degs_enr_barplot.png")
unlink("stomach_common_degs_enr_barplot.pdf")






# female-specific DEGs
female_specific_degs <- males_females_signf_degs %>%
  filter(state == "female_specific") %>%
  pull(geneName)


female_specific_degs_enr <- enrich_test(female_specific_degs, universe_enr, all_terms, all_terms2, 1, 1)
female_specific_degs_enr_tab <- female_specific_degs_enr@result


female_specific_degs_enr_barplot <- female_specific_degs_enr_tab %>%
  as.tibble() %>%
  dplyr::select(ID, Description, p_adj=p.adjust) %>%
  filter(p_adj < 0.05) %>%
  mutate(p_adj = -log10(p_adj)) %>%
  group_by(Description) %>%
  top_n(n=10, wt=p_adj) %>%
  ungroup() %>%
  filter(!str_detect(Description, "GO_MF")) %>%
  filter(!str_detect(Description, "GO_CC")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, p_adj)) %>%
  ggplot(mapping = aes(x=ID, y = p_adj, fill = p_adj)) +
    geom_bar(stat="identity") +
    theme_classic() +
    facet_wrap( ~ Description, scales = "free", nrow=3) +
    theme(
      axis.title.x=element_text(colour="black", size=15),
      axis.title.y=element_blank(),
      axis.text.y=element_text(colour="black", size=10),
      axis.text.x=element_text(colour="black", size=12),
      plot.title = element_blank(),
      strip.text.x = element_text(size=12),
      strip.background = element_blank()) +
  coord_flip() +
  scale_fill_viridis(option="D", name="") +
  scale_y_continuous(name = "Adjusted P-value (-log10)")
ggsave(filename="stomach_female_specific_degs_enr_barplot.png", plot=female_specific_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 10, height = 10)
ggsave(filename="stomach_female_specific_degs_enr_barplot.pdf", plot=female_specific_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 10, height = 10)
unlink("stomach_female_specific_degs_enr_barplot.png")
unlink("stomach_female_specific_degs_enr_barplot.pdf")







# male-specific DEGs: GO terms
male_specific_degs <- males_females_signf_degs %>%
  filter(state == "male_specific") %>%
  pull(geneName)


male_specific_degs_enr <- enrich_test(male_specific_degs, universe_enr, all_terms, all_terms2, 1, 1)
male_specific_degs_enr_tab <- male_specific_degs_enr@result


male_specific_degs_enr_barplot <- male_specific_degs_enr_tab %>%
  as.tibble() %>%
  dplyr::select(ID, Description, p_adj=p.adjust) %>%
  filter(p_adj < 0.05) %>%
  mutate(p_adj = -log10(p_adj)) %>%
  group_by(Description) %>%
  top_n(n=10, wt=p_adj) %>%
  ungroup() %>%
  filter(!str_detect(Description, "GO_MF")) %>%
  filter(!str_detect(Description, "GO_CC")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, p_adj)) %>%
  ggplot(mapping = aes(x=ID, y = p_adj, fill = p_adj)) +
    geom_bar(stat="identity") +
    theme_classic() +
    facet_wrap( ~ Description, scales = "free", nrow=3) +
    theme(
      axis.title.x=element_text(colour="black", size=15),
      axis.title.y=element_blank(),
      axis.text.y=element_text(colour="black", size=10),
      axis.text.x=element_text(colour="black", size=12),
      plot.title = element_blank(),
      strip.text.x = element_text(size=12),
      strip.background = element_blank()) +
  coord_flip() +
  scale_fill_viridis(option="D", name="") +
  scale_y_continuous(name = "Adjusted P-value (-log10)")
ggsave(filename="stomach_male_specific_degs_enr_barplot.png", plot=male_specific_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 10, height = 10)
ggsave(filename="stomach_male_specific_degs_enr_barplot.pdf", plot=male_specific_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 10, height = 10)
unlink("stomach_male_specific_degs_enr_barplot.png")
unlink("stomach_male_specific_degs_enr_barplot.pdf")



#comparecluster
stad_compCluster <- list(up.tumour.male = NULL, up.normal.male = NULL, up.tumour.female = NULL, up.normal.female = NULL)
stad_compCluster$up.tumour.male <- males_signf_degs %>% filter(state == "male_specific", logFC > 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID)
stad_compCluster$up.normal.male <- males_signf_degs %>% filter(state == "male_specific", logFC < 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID)
stad_compCluster$up.tumour.female <- females_signf_degs %>% filter(state == "female_specific", logFC > 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID)
stad_compCluster$up.normal.female <- females_signf_degs %>% filter(state == "female_specific", logFC < 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID)



stad_compCluster_enrichGO <- compareCluster(
	stad_compCluster,
	fun = "enrichGO",
    universe = bitr(universe_enr, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
	ont="BP",
	OrgDb = "org.Hs.eg.db",
	pvalueCutoff = 0.05,
	qvalueCutoff = 0.05,
	pAdjustMethod = "BH",
	readable = TRUE)

stad_compCluster_enrichGO_dot <- dotplot(stad_compCluster_enrichGO, font.size = 20, title = "Comparison of GO terms across gender-specific DEGs\nGO terms: BP")
ggsave(filename="stomach_compareCluster_degs_GO_dot.png", plot=stad_compCluster_enrichGO_dot, path = "./plots/diff_expression_tumourVSnormal/", width=19, height=8)
unlink("stomach_compareCluster_degs_GO_dot.png")


stad_compCluster_enrichKEGG <- compareCluster(
	stad_compCluster,
	fun = "enrichKEGG",
    universe = bitr(universe_enr, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
	organism = "hsa",
	pvalueCutoff = 0.05,
	qvalueCutoff = 0.05,
	pAdjustMethod = "BH",
	use_internal_data = FALSE)

stad_compCluster_enrichKEGG_dot <- dotplot(stad_compCluster_enrichKEGG, font.size = 20, title = "Comparison of KEGG pathways across gender-specific DEGs")
ggsave(filename="stomach_compareCluster_degs_keggPath_dot.png", plot=stad_compCluster_enrichKEGG_dot, path = "./plots/diff_expression_tumourVSnormal/", width=15, height=8)
unlink("stomach_compareCluster_degs_keggPath_dot.png")




save(list=ls(), file="r_workspaces/tcga_stomach_diffExpr_tumourVSnormal_characterization.RData")
