# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between genders in tumour and normal samples


library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)
library(viridis)
library(mygene)
library(data.table)
library(ggVennDiagram)


source("./r_scripts/utils.R")


# -- Thyroid

# load datasets
diff_expr_tumour_table <- read_tsv("./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_edgeR.txt") %>%
  dplyr::rename(adj.P.Val = FDR) %>%
  mutate(fdr = if_else(adj.P.Val < 0.05, "FDR < 0.05", "FDR >= 0.05"))

diff_expr_normal_table <- read_tsv("./files/diff_expr_thyroid_normal_gtex_maleVSfemale_edgeR.txt") %>%
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
    dplyr::rename(log2FC_tumour = logFC.x, log2FC_normal = logFC.y, FDR_tumour = adj.P.Val.x, FDR_normal = adj.P.Val.y)
    #mutate(log2FC_tumour = round(log2FC_tumour, 2), log2FC_normal = round(log2FC_normal, 2), FDR_tumour = round(FDR_tumour, 2), FDR_normal = round(FDR_normal, 2))
write.table(tumour_normal_signf_degs, "./files/thyroid_tumour_normal_signf_degs.txt", sep="\t", quote=F, row.names=F)



#gene annotation using MyGene.Info services
#http://mygene.info/
#http://mygene.info/metadata/fields


tumour_normal_signf_degs_annot <- queryMany(tumour_normal_signf_degs$genes, scopes='ensembl.gene', fields=c("symbol", "name", "summary"), species='human', return.as = "DataFrame") %>%
  as_tibble() %>%
  dplyr::select(query, name, summary) %>%
  inner_join(tumour_normal_signf_degs, by = c("query" = "genes")) %>%
  dplyr::select(ensembl_ID = query, gene_symbol = geneName, chrom, state, log2FC_tumour, FDR_tumour, log2FC_normal, FDR_normal, name, summary)
write.table(tumour_normal_signf_degs_annot, "./files/thyroid_tumour_normal_signf_degs_annot.txt", sep="\t", quote=F, row.names=F)



# Venn Diagram
males_females_signf_degs_venn <- draw.pairwise.venn(
  sum(table(tumour_normal_signf_degs$state)[c("common", "normal_specific")]),
  sum(table(tumour_normal_signf_degs$state)[c("common", "tumour_specific")]),
  sum(table(tumour_normal_signf_degs$state)[c("common")]),
  scaled = T,
  category = c("", ""),
  lty = rep("blank", 2),
  fill = c("#a1d76a", "#ca0020"),
  alpha = rep(0.5, 2),
  cex = rep(2.5, 3),
  cat.cex = rep(1.5, 2))
ggsave(filename="thyroid_tumour_normal_signf_degs_venn.png", plot=males_females_signf_degs_venn, path = "./plots/diff_expression_maleVSfemale_gtex_normal/", width=4, height=4)
ggsave(filename="thyroid_tumour_normal_signf_degs_venn.pdf", plot=males_females_signf_degs_venn, path = "./plots/diff_expression_maleVSfemale_gtex_normal/", width=4, height=4)
unlink("thyroid_tumour_normal_signf_degs_venn.png")
unlink("thyroid_tumour_normal_signf_degs_venn.pdf")

ggVenn <- bind_rows(
  normal_signf_degs %>% dplyr::select(geneName) %>% mutate(category = "normal\nspecific"),
  tumour_signf_degs %>% dplyr::select(geneName) %>% mutate(category = "tumour\nspecific")) %>%
  group_by(category) %>%
  summarise(geneName = list(geneName)) %>%
  ungroup() %>%
  mutate(geneName = set_names(geneName, category)) %>%
  pull(geneName) %>%
  ggVennDiagram(color = "black", category.names = NA) +
  scale_fill_gradientn(colors=c("#bf812d", "#ca0020", "#a1d76a"), values = c(0,0.05,1), guide = F)
ggVenn
ggsave(filename="thyroid_tumour_normal_signf_degs_venn_gg.png", plot=ggVenn, path = "./plots/diff_expression_maleVSfemale_gtex_normal/", width=4, height=4)
ggsave(filename="thyroid_tumour_normal_signf_degs_venn_gg.pdf", plot=ggVenn, path = "./plots/diff_expression_maleVSfemale_gtex_normal/", width=4, height=4)
unlink("thyroid_tumour_normal_signf_degs_venn_gg.png")
unlink("thyroid_tumour_normal_signf_degs_venn_gg.pdf")

write_rds(ggVenn, "./r_objects/plots/figure2/thyroid_tumour_normal_signf_degs_maleVsfemale_venn_ggplot.rds")



# append DEGs state
diff_expr_tumour_table <- diff_expr_tumour_table %>%
  left_join(tumour_signf_degs[,c("genes", "state")], by="genes")

diff_expr_normal_table <- diff_expr_normal_table %>%
  left_join(normal_signf_degs[,c("genes", "state")], by="genes")


# volcano plots


diff_expr_table <- diff_expr_tumour_table %>%
  dplyr::select(genes, geneName, chrom, logFC, adj.P.Val, state, fdr) %>%
  mutate(tissue = "Tumour") %>%
  bind_rows(diff_expr_normal_table %>% dplyr::select(genes, geneName, chrom, logFC, adj.P.Val, state, fdr) %>% mutate(tissue = "Normal")) %>%
  mutate(adj.P.Val = -log10(adj.P.Val+1e-100)) %>%
  mutate(adj.P.Val = pmap_dbl(.l = ., .f = ~ if(..3 == "chrY" | ..2 == "XIST"){NA}else{..5})) %>%
  mutate(logFC = pmap_dbl(.l = ., .f = ~ if(..3 == "chrY" | ..2 == "XIST"){NA}else{..4}))

diff_expr_vp <- ggplot( data = diff_expr_table, mapping = aes(x=logFC, y=adj.P.Val, colour=state) ) +
  geom_point(size = 1) +
  #scale_fill_manual(values=c("#D7301F", "#FDCC8A"), name = "Significance") +
  scale_colour_manual(values=c("#bf812d", "#a1d76a", "#ca0020"), na.value="#bdbdbd", labels=c("Common", "Normal-specific", "Tumour-specific", "Not SBG"), name = "SBG type") +
  facet_wrap( ~ tissue, scales = "fixed") +
  geom_line(aes(x=0), color="black", linetype=2, size = 0.3) +
  #geom_line(aes(x=1), color="black", linetype=2, size = 0.3) +
  #geom_line(aes(x=-1), color="black", linetype=2, size = 0.3) +
  geom_line(aes(y=-log10(0.05)), color="black", linetype=2, size = 0.3) +
  #theme_classic() +
  theme(axis.title = element_text(colour="black", size=16),
    axis.text = element_text(colour="black", size=12),
    legend.text=element_text(colour="black", size=12),
    legend.title=element_text(colour="black", size=16),
    plot.title = element_text(colour="black", size=18, hjust = 0.5),
    #strip.background = element_blank(),
    strip.text.x = element_text(colour="black", size=16),
    legend.position = "bottom") +
  labs(x = "Fold-change (log2)", y = "FDR (-log10)", title = "Thyroid\nMale vs Female") +
  guides(color=guide_legend(nrow=2))

ggsave(filename="diff_expr_thyroid_all_tcga_maleVSfemale.png", plot=diff_expr_vp, path = "./plots/diff_expression_maleVSfemale_gtex_normal/", width=6, height=5)
ggsave(filename="diff_expr_thyroid_all_tcga_maleVSfemale.pdf", plot=diff_expr_vp, path = "./plots/diff_expression_maleVSfemale_gtex_normal/", width=6, height=5)
unlink("diff_expr_thyroid_all_tcga_maleVSfemale.png")
unlink("diff_expr_thyroid_all_tcga_maleVSfemale.pdf")

write_rds(diff_expr_vp, "./r_objects/plots/figure2/diff_expr_thyroid_all_tcga_maleVSfemale_volcano_plot.rds")


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

all_terms <- bind_rows(kegg, go_bp, pos) %>% mutate(ont = str_replace_all(ont, "_", " "))
all_terms2 <- bind_rows(kegg2, go_bp2, pos2) %>% mutate(ont = str_replace_all(ont, "_", " "))


universe_enr <- diff_expr_tumour_table$geneName





# enrichment analysis
all_diff_genes <- tumour_normal_signf_degs %>%
  dplyr::select(state, geneName) %>%
  group_by(state) %>%
  summarise(geneName = list(geneName)) %>%
  mutate(enr = map(.x = geneName, .f = enrich_test, universe = universe_enr, terms1 = all_terms, terms2 = all_terms2, p_adj=1, q_value=1)) %>%
  dplyr::select(-geneName) %>%
  unnest(cols = c(enr))

write.table(all_diff_genes, "./files/thyroid_tumour_normal_degs_enr.txt", sep="\t", quote=F, row.names=F)





tn_enr_bp <- all_diff_genes %>%
  filter(p.adjust < 0.05 & state == "common") %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(state, Description) %>%
  top_n(5, log10_p) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count), Description = fct_infreq(Description)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_grid(Description ~ .,
    space = "free_y",
    scales = "free",
    labeller=labeller(
      state = c("common" = "Common", "normal_specific" = "Normal-specific"),
      Description = c("GO_BP" = "GO BP", "KEGG" = "KEGG", "ONCO" = "Onco", "IMMUNO" = "Immunogenic", "POS" = "Positional", "CM" = "Cancer\nmodules"))) +
  theme(
    axis.title.x=element_text(colour="black", size=15),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_text(colour="black", size=13),
    plot.title = element_blank(),
    strip.text.y = element_text(colour="black", size=14),
    strip.text.x = element_text(colour="black", size=14),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=13),
    legend.title = element_text(colour="black", size=15)) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adjusted P-value\n(-log10)") +
  scale_y_continuous(name = "Count") +
  labs(title = "Thyroid")

ggsave(filename="thyroid_tumour_normal_enr_bp_common.png", plot=tn_enr_bp, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 10, height = 5)
ggsave(filename="thyroid_tumour_normal_enr_bp_common.pdf", plot=tn_enr_bp, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 10, height = 5)
unlink("thyroid_tumour_normal_enr_bp_common.png")
unlink("thyroid_tumour_normal_enr_bp_common.pdf")


tn_enr_bp2 <- all_diff_genes %>%
  filter(p.adjust < 0.05 & state == "normal_specific" & Description != "POS") %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(state, Description) %>%
  top_n(5, log10_p) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count), Description = fct_infreq(Description)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_grid(Description ~ .,
    scales = "free",
    space = "free_y",
    labeller=labeller(Description = c("GO_BP" = "GO BP", "KEGG" = "KEGG", "POS" = "Pos"))) +
  theme(
    axis.title.x=element_text(colour="black", size=17),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=14),
    axis.text.x=element_text(colour="black", size=17),
    plot.title = element_blank(),
    strip.text.y = element_text(colour="black", size=17),
    strip.text.x = element_text(colour="black", size=20),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=14),
    legend.title = element_text(colour="black", size=15)) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adjusted P-value\n(-log10)")

ggsave(filename="thyroid_tumour_normal_enr_bp_normal_specific.png", plot=tn_enr_bp2, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 10, height = 4)
ggsave(filename="thyroid_tumour_normal_enr_bp_normal_specific.pdf", plot=tn_enr_bp2, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 10, height = 4)
unlink("thyroid_tumour_normal_enr_bp_normal_specific.png")
unlink("thyroid_tumour_normal_enr_bp_normal_specific.pdf")

write_rds(tn_enr_bp2, "./r_objects/plots/figure2/thyroid_tumour_normal_enrichment_barplot_normal_specific.rds")





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
      state = c("common" = "Common", "normal_specific" = "Normal-specific", "tumour_specific" = "Tumour-specific"))) +
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

ggsave(filename="thyroid_degs_chr_bp.png", plot=degs_chr_bp, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 6, height = 6)
ggsave(filename="thyroid_degs_chr_bp.pdf", plot=degs_chr_bp, path="./plots/diff_expression_maleVSfemale_gtex_normal/", width = 6, height = 6)
unlink("thyroid_degs_chr_bp.png")
unlink("thyroid_degs_chr_bp.pdf")




#comparecluster
compCluster <- list(
  `tumour-specific\nup-male` = tumour_normal_signf_degs %>% filter(state == "tumour_specific", log2FC_tumour > 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID),
  `tumour-specific\nup-female` = tumour_normal_signf_degs %>% filter(state == "tumour_specific", log2FC_tumour < 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID),
  `normal-specific\nup-male` = tumour_normal_signf_degs %>% filter(state == "normal_specific", log2FC_normal > 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID),
  `normal-specific\nup-female` = tumour_normal_signf_degs %>% filter(state == "normal_specific", log2FC_normal < 0) %>% pull(genes) %>% bitr(., "ENSEMBL", "ENTREZID", "org.Hs.eg.db") %>% pull(ENTREZID))


#code <- c("up tumour\nmale" = "up tumour male", "up normal\nmale" = "up normal male", "up tumour\nfemale" = "up tumour female", "up normal\nfemale" = "up normal female")
compCluster_enrichGO <- compareCluster(
	compCluster,
	fun = "enrichGO",
  universe = bitr(diff_expr_tumour_table$genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
	ont="BP",
	OrgDb = "org.Hs.eg.db",
	pvalueCutoff = 0.05,
	qvalueCutoff = 0.05,
	pAdjustMethod = "BH",
	readable = TRUE)

dotplot(compCluster_enrichGO, font.size = 20, title = "")



# cancer genes list
cancer_genes <- read_tsv("./data/Census_allMon_May_13_17_05_42_2019.tsv")

thca_degs_MvsF_cancer_genes <- tumour_normal_signf_degs %>%
  dplyr::select(geneName, state, log2FC_normal, log2FC_tumour) %>%
  inner_join(cancer_genes %>% dplyr::select(`Gene Symbol`, `Genome Location`, `Tumour Types(Somatic)`, `Tumour Types(Germline)`, `Cancer Syndrome`, `Role in Cancer`), by = c("geneName" = "Gene Symbol"))
write.table(thca_degs_MvsF_cancer_genes, "./files/thca_degs_MvsF_cancer_genes.txt", sep="\t", quote=F, row.names=F)


thca_degs_MvsF_normal_specific_enriched_cancer_genes <- all_diff_genes %>%
  filter(p.adjust < 0.05 & state == "normal_specific") %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(state, Description) %>%
  top_n(5, log10_p) %>%
  ungroup() %>%
  dplyr::select(state, Description, ID, geneID) %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest() %>%
  group_by(state, geneID) %>%
  summarise(ID = paste(ID, collapse="/"), Description = paste(unique(Description), collapse="/")) %>%
  ungroup() %>%
  inner_join(tumour_normal_signf_degs %>% dplyr::select(geneName, chrom, log2FC_normal), by = c("geneID" = "geneName")) %>%
  inner_join(cancer_genes %>% dplyr::select(`Gene Symbol`, `Genome Location`, `Tumour Types(Somatic)`, `Tumour Types(Germline)`, `Cancer Syndrome`, `Role in Cancer`), by = c("geneID" = "Gene Symbol"))
write.table(thca_degs_MvsF_normal_specific_enriched_cancer_genes, "./files/thca_degs_MvsF_normal_specific_enriched_cancer_genes.txt", sep="\t", quote=F, row.names=F)




# cancer drivers gene list
driver_genes <- fread("./data/tcga/cancer_driver_genes.txt") %>%
  as_tibble()

driver_genes_summary <- driver_genes %>%
  group_by(Gene, Decision) %>%
  summarise(n = n()) %>%
  ungroup()

thca_degs_MvsF_drivers <- tumour_normal_signf_degs %>%
  dplyr::select(geneName, state, log2FC_normal, log2FC_tumour) %>%
  inner_join(driver_genes_summary, by = c("geneName" = "Gene"))


diff_expr_normal_table %>%
  filter( !geneName %in% (tumour_normal_signf_degs %>% filter(state == "normal_specific" | state == "common") %>% pull(geneName)) ) %>%
  inner_join(driver_genes_summary, by = c("geneName" = "Gene"))




save(list=ls(), file="r_workspaces/tcga_gtex_thyroid_diffExpr_malesVSfemales_characterization.RData")
