# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between genders


library(tidyverse)
library(RColorBrewer)
library(VennDiagram)
library(clusterProfiler)
library(viridis)





# -- Stomach

# load datasets
diff_expr_tumour_table <- read_tsv("./files/diff_expr_stomach_tumour_tcga_maleVSfemale_edgeR.txt") %>%
  dplyr::rename(adj.P.Val = FDR) %>%
  mutate(fdr = if_else(adj.P.Val < 0.05, "FDR < 0.05", "FDR >= 0.05"))

diff_expr_normal_table <- read_tsv("./files/diff_expr_stomach_normal_tcga_maleVSfemale_edgeR.txt") %>%
  dplyr::rename(adj.P.Val = FDR) %>%
  mutate(fdr = if_else(adj.P.Val < 0.05, "FDR < 0.05", "FDR >= 0.05"))



# significant DEGs between males and females
tumour_signf_degs <- diff_expr_tumour_table %>%
    filter(adj.P.Val <= 0.05)

normal_signf_degs <- diff_expr_normal_table %>%
    filter(adj.P.Val <= 0.05)

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
    sum(table(tumour_normal_signf_degs$state)[c("common", "normal_specific")]),
    sum(table(tumour_normal_signf_degs$state)[c("common", "tumour_specific")]),
    sum(table(tumour_normal_signf_degs$state)[c("common")]),
    category = c("", ""),
    lty = rep("blank", 2),
    fill = c("green", "red"),
    alpha = rep(0.5, 2),
    cex = rep(2.5, 3),
    cat.cex = rep(1.5, 2))
ggsave(filename="stomach_tumour_normal_signf_degs_venn.png", plot=males_females_signf_degs_venn, path = "./plots/diff_expression_maleVSfemale_tcga_normal/", width=3, height=3)
unlink("stomach_tumour_normal_signf_degs_venn.png")



# append DEGs state
diff_expr_tumour_table <- diff_expr_tumour_table %>%
  left_join(tumour_signf_degs[,c("genes", "state")], by="genes")

diff_expr_normal_table <- diff_expr_normal_table %>%
  left_join(normal_signf_degs[,c("genes", "state")], by="genes")




diff_expr_table <- diff_expr_tumour_table %>%
  dplyr::select(genes, logFC, adj.P.Val, state, fdr) %>%
  mutate(tissue = "Tumour") %>%
  bind_rows(diff_expr_normal_table %>% dplyr::select(genes, logFC, adj.P.Val, state, fdr) %>% mutate(tissue = "Normal"))

diff_expr_vp <- ggplot( data = diff_expr_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state) ) +
    geom_point() +
    #scale_fill_manual(values=c("#D7301F", "#FDCC8A"), name = "Significance") +
    scale_colour_manual(values=c("#fdbb84", "green", "#d95f02"), na.value="#bdbdbd", labels=c("Common", "Normal-specific", "Tumour-specific", "Not DEGs"), name = "DEG type") +
    facet_wrap( ~ tissue, scales = "free") +
    geom_line(aes(x=0), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(x=-1), color="black", linetype=2, size = 0.3) +
    geom_line(aes(y=-log10(0.05)), color="black", linetype=2, size = 0.3) +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=18),
      axis.text = element_text(colour="black", size=16),
      legend.text=element_text(colour="black", size=16),
      legend.title=element_text(colour="black", size=18),
      strip.background = element_blank(),
      strip.text.x = element_text(colour="black", size=18)) +
    labs(x = "logFC", y = "FDR (-log10)")
ggsave(filename="diff_expr_stomach_all_tcga_maleVSfemale.png", plot=diff_expr_vp, path = "./plots/diff_expression_maleVSfemale_tcga_normal/", width=10, height=5)
unlink("diff_expr_stomach_all_tcga_maleVSfemale.png")







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

  enr <- enr@result %>% dplyr::select(Description, ID, Count, p.adjust, GeneRatio, BgRatio, geneID)
  return(list(enr))

}




# load gene lists
kegg <- read.gmt("./data/gene_lists/c2.cp.kegg.v6.2.symbols.gmt") %>% mutate(ont = str_replace(ont, "KEGG_", ""))
kegg2 <- data.frame(ont = kegg$ont, name = "KEGG") %>% dplyr::distinct()
go_bp <- read.gmt("./data/gene_lists/c5.bp.v6.2.symbols.gmt") %>% mutate(ont = str_replace(ont, "GO_", ""))
go_bp2 <- data.frame(ont = go_bp$ont, name = "GO_BP") %>% dplyr::distinct()
go_mf <- read.gmt("./data/gene_lists/c5.mf.v6.2.symbols.gmt") %>% mutate(ont = str_replace(ont, "GO_", ""))
go_mf2 <- data.frame(ont = go_mf$ont, name = "GO_MF") %>% dplyr::distinct()
go_cc <- read.gmt("./data/gene_lists/c5.cc.v6.2.symbols.gmt") %>% mutate(ont = str_replace(ont, "GO_", ""))
go_cc2 <- data.frame(ont = go_cc$ont, name = "GO_CC") %>% dplyr::distinct()
onco <- read.gmt("./data/gene_lists/c6.all.v6.2.symbols.gmt")
onco2 <- data.frame(ont = onco$ont, name = "ONCO") %>% dplyr::distinct()
immuno <- read.gmt("./data/gene_lists/c7.all.v6.2.symbols.gmt")
immuno2 <- data.frame(ont = immuno$ont, name = "IMMUNO") %>% dplyr::distinct()
pos <- read.gmt("./data/gene_lists/c1.all.v6.2.symbols.gmt")
pos2 <- data.frame(ont = pos$ont, name = "POS") %>% dplyr::distinct()
cm <- read.gmt("./data/gene_lists/c4.cm.v6.2.symbols.gmt")
cm2 <- data.frame(ont = cm$ont, name = "CM") %>% dplyr::distinct()

all_terms <- bind_rows(kegg, go_bp, onco, immuno, pos, cm) %>% mutate(ont = str_replace_all(ont, "_", " "))
all_terms2 <- bind_rows(kegg2, go_bp2, onco2, immuno2, pos2, cm2) %>% mutate(ont = str_replace_all(ont, "_", " "))


universe_enr <- diff_expr_tumour_table$geneName





# enrichment analysis
all_diff_genes <- tumour_normal_signf_degs %>%
  filter(state == "tumour_specific") %>%
  dplyr::select(state, geneName) %>%
  group_by(state) %>%
  mutate(geneName = list(geneName)) %>%
  unique() %>%
  rowwise() %>%
  mutate(enr = enrich_test(geneName, universe_enr, all_terms, all_terms2, 1, 1)) %>%
  dplyr::select(-geneName) %>%
  unnest()
write.table(all_diff_genes, "./files/stomach_tumour_normal_degs_enr2.txt", sep="\t", quote=F, row.names=F)





tn_enr_bp <- all_diff_genes %>%
  filter(p.adjust < 0.05) %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(state, Description) %>%
  top_n(5, log10_p) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count), Description = fct_infreq(Description)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_grid(Description ~ state,
    space = "free_y",
    scales = "free",
    labeller=labeller(
      state = c("tumour_specific" = "Tumour-specific"),
      Description = c("GO_BP" = "GO biological processes", "KEGG" = "KEGG", "ONCO" = "Onco", "IMMUNO" = "Immunogenic", "POS" = "Positional", "CM" = "Cancer\nmodules"))) +
  theme(
    axis.title.x=element_text(colour="black", size=15),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_text(colour="black", size=13),
    plot.title = element_blank(),
    strip.text.y = element_text(colour="black", size=13),
    strip.text.x = element_text(colour="black", size=15),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=13),
    legend.title = element_text(colour="black", size=13)) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
  scale_y_continuous(name = "Number of genes")
ggsave(filename="stomach_tumour_normal_enr_bp.png", plot=tn_enr_bp, path="./plots/diff_expression_maleVSfemale_tcga_normal/", width = 5, height = 2)
ggsave(filename="stomach_tumour_normal_enr_bp.pdf", plot=tn_enr_bp, path="./plots/diff_expression_maleVSfemale_tcga_normal/", width = 5, height = 2)
unlink("stomach_tumour_normal_enr_bp.png")
unlink("stomach_tumour_normal_enr_bp.pdf")



# barplot of number of DEGs by chromosome
degs_chr_bp <- tumour_normal_signf_degs %>%
  group_by(state, chrom) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(chrom = fct_relevel(chrom, c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))) %>%
  mutate(chrom = fct_rev(chrom)) %>%
  ggplot(mapping = aes(x=chrom, y = counts)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  coord_flip() +
  theme_classic() +
  facet_grid(. ~ state,
    scales = "fixed",
    labeller=labeller(
      state = c("common" = "Common", "normal_specific" = "Normal\nspecific", "tumour_specific" = "Tumour\nspecific"))) +
  theme(
    axis.title.x=element_text(colour="black", size=15),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=14),
    axis.text.x=element_text(colour="black", size=13),
    plot.title = element_blank(),
    strip.text.y = element_text(colour="black", size=10),
    strip.text.x = element_text(colour="black", size=14),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=13),
    legend.title = element_text(colour="black", size=15)) +
  scale_y_continuous(name = "Number of genes")
ggsave(filename="stomach_degs_chr_bp.png", plot=degs_chr_bp, path="./plots/diff_expression_maleVSfemale_tcga_normal/", width = 4, height = 5)
ggsave(filename="stomach_degs_chr_bp.pdf", plot=degs_chr_bp, path="./plots/diff_expression_maleVSfemale_tcga_normal/", width = 4, height = 5)
unlink("stomach_degs_chr_bp.png")
unlink("stomach_degs_chr_bp.pdf")





save(list=ls(), file="r_workspaces/tcga_stomach_diffExpr_malesVSfemales_characterization.RData")
