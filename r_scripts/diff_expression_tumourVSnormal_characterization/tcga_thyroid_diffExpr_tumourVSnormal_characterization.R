# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between tumour and normal samples


library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)
library(viridis)
library(mygene)




# -- Thyroid

# load datasets
diff_expr_males_table <- read_tsv("./files/diff_expr_thyroid_males_tcga_tumourVSnormal_edgeR.txt") %>%
  dplyr::rename(adj.P.Val = FDR) %>%
  mutate(fdr = if_else(adj.P.Val < 0.05, "FDR < 0.05", "FDR >= 0.05"))

diff_expr_females_table <- read_tsv("./files/diff_expr_thyroid_females_tcga_tumourVSnormal_edgeR.txt") %>%
  dplyr::rename(adj.P.Val = FDR) %>%
  mutate(fdr = if_else(adj.P.Val < 0.05, "FDR < 0.05", "FDR >= 0.05"))



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
write.table(males_females_signf_degs, "./files/thyroid_males_females_signf_degs.txt", sep="\t", quote=F, row.names=F)


#gene annotation using MyGene.Info services
#http://mygene.info/
#http://mygene.info/metadata/fields


males_females_signf_degs_annot <- queryMany(males_females_signf_degs$genes, scopes='ensembl.gene', fields=c("symbol", "name", "summary"), species='human', return.as = "DataFrame") %>%
  as_tibble() %>%
  dplyr::select(query, name, summary) %>%
  inner_join(males_females_signf_degs, by = c("query" = "genes")) %>%
  dplyr::select(ensembl_ID = query, gene_symbol = geneName, chrom, state, logFC_males, FDR_males = adj.P.Val_males, logFC_females, FDR_females = adj.P.Val_females, name, summary)
write.table(males_females_signf_degs_annot, "./files/thyroid_males_females_signf_degs_annot.txt", sep="\t", quote=F, row.names=F)



# barplot of all significant DEGs between males and females: common, male-specific, female-specific
males_females_signf_degs_bp <- ggplot(data=males_females_signf_degs, mapping=aes(x=state, y=..count.., fill=state)) +
    geom_bar() +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=15),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    scale_x_discrete(labels = c("Common", "Female-specific", "Male-specific")) +
    labs(x = "", y = "Number of DEGs", title="Thyroid")
ggsave(filename="thyroid_males_females_signf_degs_bp.png", plot=males_females_signf_degs_bp, path = "./plots/diff_expression_tumourVSnormal/")
unlink("thyroid_males_females_signf_degs_bp.png")




# Venn Diagram
males_females_signf_degs_venn <- draw.pairwise.venn(
  sum(table(males_females_signf_degs$state)[c("common", "male_specific")]),
  sum(table(males_females_signf_degs$state)[c("common", "female_specific")]),
  sum(table(males_females_signf_degs$state)[c("common")]),
  scaled = T,
	category = c("", ""),
	lty = rep("blank", 2),
	fill = c("#3182bd", "#dd1c77"),
	alpha = rep(0.5, 2),
	cex = rep(2.5, 3),
	cat.cex = rep(1.5, 2))
ggsave(filename="thyroid_males_females_signf_degs_venn.png", plot=males_females_signf_degs_venn, path = "./plots/diff_expression_tumourVSnormal/", width=4, height=4)
ggsave(filename="thyroid_males_females_signf_degs_venn.pdf", plot=males_females_signf_degs_venn, path = "./plots/diff_expression_tumourVSnormal/", width=4, height=4)
unlink("thyroid_males_females_signf_degs_venn.png")
unlink("thyroid_males_females_signf_degs_venn.pdf")











# append DEGs state
diff_expr_males_table <- diff_expr_males_table %>%
    left_join(males_signf_degs[,c("genes", "state")], by="genes")

diff_expr_females_table <- diff_expr_females_table %>%
    left_join(females_signf_degs[,c("genes", "state")], by="genes")



# volcano plots
diff_expr_males_vp <- ggplot( data = diff_expr_males_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state, fill=fdr) ) +
    geom_point(pch = 21) +
    scale_fill_manual(values=c("#D7301F", "#FDCC8A"), name = "Significance") +
    scale_colour_manual(values=c("#3182bd", "green"), na.value="#bdbdbd", labels=c("Common", "Male-specific", "Not DEGs"), name = "DEG type") +
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
    labs(x = "logFC", y = "FDR (-log10)", title="Males")
ggsave(filename="diff_expr_thyroid_males_tcga_tumourVSnormal.png", plot=diff_expr_males_vp, path = "./plots/diff_expression_tumourVSnormal/", width=6, height=6)
unlink("diff_expr_thyroid_males_tcga_tumourVSnormal.png")

diff_expr_females_vp <- ggplot( data = diff_expr_females_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state, fill=fdr) ) +
    geom_point(pch = 21) +
    scale_fill_manual(values=c("#D7301F", "#FDCC8A"), name = "Significance") +
    scale_colour_manual(values=c("#3182bd", "green"), na.value="#bdbdbd", labels=c("Common", "Female-specific", "Not DEGs"), name = "DEG type") +
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
    labs(x = "logFC", y = "FDR (-log10)", title="Females")
ggsave(filename="diff_expr_thyroid_females_tcga_tumourVSnormal.png", plot=diff_expr_females_vp, path = "./plots/diff_expression_tumourVSnormal/", width=6, height=6)
unlink("diff_expr_thyroid_females_tcga_tumourVSnormal.png")


diff_expr_table <- diff_expr_males_table %>%
  dplyr::select(genes, logFC, adj.P.Val, state, fdr) %>%
  mutate(sex = "Males") %>%
  bind_rows(diff_expr_females_table %>% dplyr::select(genes, logFC, adj.P.Val, state, fdr) %>% mutate(sex = "Females"))

diff_expr_vp <- ggplot( data = diff_expr_table, mapping = aes(x=logFC, y=-log10(adj.P.Val), colour=state) ) +
    geom_point(pch = 20) +
    #scale_fill_manual(values=c("#D7301F", "#FDCC8A"), name = "Significance") +
    #scale_colour_manual(values=c("#3182bd", "green", "yellow"), na.value="#bdbdbd", labels=c("Common", "Female-specific", "Male-specific", "Not DEGs"), name = "DEG type") +
    scale_colour_manual(values=c("#fdbb84", "#dd1c77", "#3182bd"), na.value="#bdbdbd", labels=c("Common", "Female-specific", "Male-specific", "Not DEG"), name = "DEG type") +
    facet_wrap( ~ sex, scales = "free") +
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
    labs(x = "log2FC", y = "FDR (-log10)") +
    guides(color = guide_legend(override.aes = list(size=4)))
ggsave(filename="diff_expr_thyroid_all_tcga_tumourVSnormal.png", plot=diff_expr_vp, path = "./plots/diff_expression_tumourVSnormal/", width=10, height=5)
unlink("diff_expr_thyroid_all_tcga_tumourVSnormal.png")



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


all_terms <- bind_rows(kegg, go_bp, onco, cm, immuno, pos) %>% mutate(ont = str_replace_all(ont, "_", " "))
all_terms2 <- bind_rows(kegg2, go_bp2, onco2, cm2, immuno2, pos2) %>% mutate(ont = str_replace_all(ont, "_", " "))


universe_enr <- diff_expr_males_table$geneName





# enrichment analysis
all_diff_genes <- males_females_signf_degs %>%
  dplyr::select(state, geneName) %>%
  group_by(state) %>%
  mutate(geneName = list(geneName)) %>%
  unique() %>%
  rowwise() %>%
  mutate(enr = enrich_test(geneName, universe_enr, all_terms, all_terms2, 1, 1)) %>%
  dplyr::select(-geneName) %>%
  unnest()
write.table(all_diff_genes, "./files/thyroid_male_female_degs_enr.txt", sep="\t", quote=F, row.names=F)


mf_enr_bp <- all_diff_genes %>%
  filter(p.adjust < 0.05 & !(Description %in% c("POS", "IMMUNO", "ONCO"))) %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(state, Description) %>%
  top_n(5, log10_p) %>%
  ungroup() %>%
  mutate(ID = str_replace_all(ID,
    pattern=c(
      "MODULE 88" = "METABOLIC AND XENOBIOTIC RESPONSE GENES (MODULE 88)",
      "MODULE 55" = "EXTRACELLULAR SPACE (MODULE 55)",
      "\\bMODULE 6\\b" = "TRACHEA GENES (MODULE 6)",
      "MODULE 12" = "SPINAL CORD GENES (MODULE 12)",
      "\\bMODULE 2\\b" = "DORSAL ROOT GANGLIA GENES (MODULE 2)",
      "MODULE 45" = "WHOLE BLOOD GENES (MODULE 45)",
      "MODULE 84" = "INFLAMMATORY RESPONSE (MODULE 84)",
      "\\bMODULE 27\\b" = "RECEPTOR ACTIVITY (MODULE 27)",
      "\\bMODULE 64\\b" = "MEMBRANE RECEPTORS (MODULE 64)",
      "MODULE 75" = "IMMUNE RESPONSE (MODULE 75)",
      "MODULE 46" = "IMMUNE RESPONSE (MODULE 46)",
      "MODULE 44" = "THYMUS (MODULE 44)") )) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_grid(Description ~ state,
    scales = "free",
    space = "free_y",
    labeller=labeller(
      state = c("common" = "Common", "female_specific" = "Female-specific", "male_specific" = "Male-specific"),
      Description = c(
      "GO_BP" = "GO biological\nprocesses",
      "KEGG" = "KEGG\npathways",
      "CM" = "Cancer modules"))) +
  theme(
    axis.title.x=element_text(colour="black", size=20),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=15),
    axis.text.x=element_text(colour="black", size=16),
    plot.title = element_blank(),
    strip.text = element_text(colour="black", size=18),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=15),
    legend.title = element_text(colour="black", size=20),
    panel.spacing = unit(1.5, "lines")) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
  scale_y_continuous(name = "Number of genes")
ggsave(filename="thyroid_male_female_enr_bp.png", plot=mf_enr_bp, path="./plots/diff_expression_tumourVSnormal/", width = 15, height = 10)
ggsave(filename="thyroid_male_female_enr_bp.pdf", plot=mf_enr_bp, path="./plots/diff_expression_tumourVSnormal/", width = 15, height = 10)
unlink("thyroid_male_female_enr_bp.png")
unlink("thyroid_male_female_enr_bp.pdf")




# common DEGs
common_degs_enr_barplot <- all_diff_genes %>%
  filter(p.adjust < 0.05, state == "common") %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(state, Description) %>%
  top_n(5, log10_p) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_grid(Description ~ state,
    scales = "free", space = "free_y",
    labeller=labeller(state = c("common" = "Common DEGs"),
      Description = c(
        "GO_BP" = "GO-BP",
        "KEGG" = "KEGG pathways",
        "CM" = "Cancer modules",
        "ONCO" = "Oncogenic",
        "IMMUNO" = "Immunogenic",
        "POS" = "Positional"
      ))) +
  theme(
    axis.title.x=element_text(colour="black", size=16),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=14),
    axis.text.x=element_text(colour="black", size=14),
    plot.title = element_blank(),
    strip.text = element_text(colour="black", size=15),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=13),
    legend.title = element_text(colour="black", size=15)) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
  scale_y_continuous(name = "Number of genes")
ggsave(filename="thyroid_common_degs_enr_barplot.png", plot=common_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 10, height = 11)
ggsave(filename="thyroid_common_degs_enr_barplot.pdf", plot=common_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 10, height = 11)
unlink("thyroid_common_degs_enr_barplot.png")
unlink("thyroid_common_degs_enr_barplot.pdf")






# female-specific DEGs
female_specific_degs_enr_barplot <- all_diff_genes %>%
  filter(p.adjust < 0.05, state == "female_specific") %>%
  filter(!str_detect(Description, "IMMUNO")) %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(state, Description) %>%
  top_n(10, log10_p) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_grid(Description ~ state,
    scales = "free", space = "free_y",
    labeller=labeller(state = c("female_specific" = "Female-specific DEGs"),
      Description = c(
        "GO_BP" = paste("GO BP", all_diff_genes %>% filter(p.adjust < 0.05, state == "female_specific", Description == "GO_BP") %>% nrow),
        "KEGG" = paste("KEGG pathways", all_diff_genes %>% filter(p.adjust < 0.05, state == "female_specific", Description == "KEGG") %>% nrow),
        "CM" = paste("Cancer modules", all_diff_genes %>% filter(p.adjust < 0.05, state == "female_specific", Description == "CM") %>% nrow),
        "ONCO" = paste("ONCOGENIC sets", all_diff_genes %>% filter(p.adjust < 0.05, state == "female_specific", Description == "ONCO") %>% nrow)))) +
  theme(
    axis.title.x=element_text(colour="black", size=15),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_text(colour="black", size=13),
    plot.title = element_blank(),
    strip.text = element_text(colour="black", size=14),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=13),
    legend.title = element_text(colour="black", size=15)) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
  scale_y_continuous(name = "Number of genes")
ggsave(filename="thyroid_female_specific_degs_enr_barplot.png", plot=female_specific_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 8, height = 11)
ggsave(filename="thyroid_female_specific_degs_enr_barplot.pdf", plot=female_specific_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 8, height = 11)
unlink("thyroid_female_specific_degs_enr_barplot.png")
unlink("thyroid_female_specific_degs_enr_barplot.pdf")







# male-specific DEGs: GO terms
male_specific_degs_enr_barplot <- all_diff_genes %>%
  filter(p.adjust < 0.05, state == "male_specific") %>%
  filter(!str_detect(Description, "IMMUNO")) %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(state, Description) %>%
  top_n(10, log10_p) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_grid(Description ~ state,
    scales = "free", space = "free_y",
    labeller=labeller(state = c("male_specific" = "Male-specific DEGs"),
      Description = c(
        "GO_BP" = paste("GO BP", all_diff_genes %>% filter(p.adjust < 0.05, state == "male_specific", Description == "GO_BP") %>% nrow),
        "KEGG" = paste("KEGG pathways", all_diff_genes %>% filter(p.adjust < 0.05, state == "male_specific", Description == "KEGG") %>% nrow),
        "CM" = paste("Cancer modules", all_diff_genes %>% filter(p.adjust < 0.05, state == "male_specific", Description == "CM") %>% nrow),
        "ONCO" = paste("ONCOGENIC sets", all_diff_genes %>% filter(p.adjust < 0.05, state == "male_specific", Description == "ONCO") %>% nrow)))) +
  theme(
    axis.title.x=element_text(colour="black", size=15),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_text(colour="black", size=13),
    plot.title = element_blank(),
    strip.text = element_text(colour="black", size=14),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=13),
    legend.title = element_text(colour="black", size=15)) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
  scale_y_continuous(name = "Number of genes")
ggsave(filename="thyroid_male_specific_degs_enr_barplot.png", plot=male_specific_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 8, height = 11)
ggsave(filename="thyroid_male_specific_degs_enr_barplot.pdf", plot=male_specific_degs_enr_barplot, path="./plots/diff_expression_tumourVSnormal/", width = 8, height = 11)
unlink("thyroid_male_specific_degs_enr_barplot.png")
unlink("thyroid_male_specific_degs_enr_barplot.pdf")



#comparecluster
thca_compCluster <- list(
  `up tumour\nmale` = males_signf_degs %>% filter(state == "male_specific", logFC > 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID),
  `up normal\nmale` = males_signf_degs %>% filter(state == "male_specific", logFC < 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID),
  `up tumour\nfemale` = females_signf_degs %>% filter(state == "female_specific", logFC > 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID),
  `up normal\nfemale` = females_signf_degs %>% filter(state == "female_specific", logFC < 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID))


#thca_compCluster <- thca_compCluster[c(1,4)]
code <- c("up tumour\nmale" = "up tumour male", "up normal\nmale" = "up normal male", "up tumour\nfemale" = "up tumour female", "up normal\nfemale" = "up normal female")

thca_compCluster_enrichGO <- compareCluster(
	thca_compCluster,
	fun = "enrichGO",
  universe = bitr(diff_expr_males_table$genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
	ont="BP",
	OrgDb = "org.Hs.eg.db",
	pvalueCutoff = 0.05,
	qvalueCutoff = 0.05,
	pAdjustMethod = "BH",
	readable = TRUE)
thca_compCluster_enrichGO2 <- thca_compCluster_enrichGO@compareClusterResult
thca_compCluster_enrichGO2$Cluster <- unname(code[thca_compCluster_enrichGO2$Cluster])
write.table(thca_compCluster_enrichGO2, "./files/thyroid_male_female_compCluster_GO.txt", sep="\t", quote=F, row.names=F)

thca_compCluster_enrichGO_dot <- dotplot(thca_compCluster_enrichGO, font.size = 20, title = "")
thca_compCluster_enrichGO_dot$labels$colour = "Adjusted\nP-value"
thca_compCluster_enrichGO_dot$labels$size = "Gene ratio"
thca_compCluster_enrichGO_dot$theme$axis.text.y$size = 22
thca_compCluster_enrichGO_dot$theme$axis.text.y$size = 22
thca_compCluster_enrichGO_dot$theme$legend.text$size = rel(1.5)
thca_compCluster_enrichGO_dot$theme$legend.title$size = rel(2)
ggsave(filename="thyroid_compareCluster_degs_GO_dot.png", plot=thca_compCluster_enrichGO_dot, path = "./plots/diff_expression_tumourVSnormal/", width=16, height=7)
unlink("thyroid_compareCluster_degs_GO_dot.png")



thca_compCluster_enrichKEGG <- compareCluster(
	thca_compCluster,
	fun = "enrichKEGG",
  universe = bitr(diff_expr_males_table$genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
	organism = "hsa",
	pvalueCutoff = 0.05,
	qvalueCutoff = 0.05,
	pAdjustMethod = "BH",
	use_internal_data = FALSE)
thca_compCluster_enrichKEGG2 <- thca_compCluster_enrichKEGG@compareClusterResult
thca_compCluster_enrichKEGG2$Cluster <- unname(code[thca_compCluster_enrichKEGG2$Cluster])
write.table(thca_compCluster_enrichKEGG2, "./files/thyroid_male_female_compCluster_KEGG.txt", sep="\t", quote=F, row.names=F)


thca_compCluster_enrichKEGG_dot <- dotplot(thca_compCluster_enrichKEGG, font.size = 20, title = "KEGG pathways")
ggsave(filename="thyroid_compareCluster_degs_keggPath_dot.png", plot=thca_compCluster_enrichKEGG_dot, path = "./plots/diff_expression_tumourVSnormal/", width=12, height=6)
unlink("thyroid_compareCluster_degs_keggPath_dot.png")




# cancer genes list
cancer_genes <- read_tsv("./data/Census_allMon_May_13_17_05_42_2019.tsv")

thca_degs_TvsN_cancer_genes <- males_females_signf_degs %>%
 dplyr::select(geneName, state, logFC_males, logFC_females) %>%
 inner_join(cancer_genes %>% dplyr::select(`Gene Symbol`, `Genome Location`, `Tumour Types(Somatic)`, `Tumour Types(Germline)`, `Cancer Syndrome`, `Role in Cancer`), by = c("geneName" = "Gene Symbol"))



# cancer drivers gene list
driver_genes <- data.table::fread("./data/tcga/cancer_driver_genes.txt") %>%
 as_tibble()

driver_genes_summary <- driver_genes %>%
 group_by(Gene, Decision) %>%
 summarise(n = n()) %>%
 ungroup()

thca_degs_TvsN_drivers <- males_females_signf_degs %>%
 dplyr::select(geneName, state, logFC_males, logFC_females) %>%
 inner_join(driver_genes_summary, by = c("geneName" = "Gene"))




save(list=ls(), file="r_workspaces/tcga_thyroid_diffExpr_tumourVSnormal_characterization.RData")
