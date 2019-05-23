# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially expressed genes between tumour and normal samples
# Thyroid and Stomach



library(tidyverse)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)


# load TCGA annotation
tcga_annotation <- read_tsv("./data/annotation/geneAnnot.gencode.v22.txt") %>%
  mutate(geneID = str_replace(geneID, "\\.[0-9]+", "")) %>%
  dplyr::rename(ensemblGeneID = geneID)

chrom_genes <- tcga_annotation %>%
  filter(geneType == "protein_coding" | geneType == "lincRNA") %>%
  group_by(chrom) %>%
  summarise(total_genes = n())


# load tumour-supressor genes table
tsgs <- read_tsv("./data/human_tsgs/Human_TSGs.txt")


# load onco-genes table
ocgs <- read_tsv("./data/human_oncogenes/ongene_human.txt")



# load degs tables
stomach_degs <- read_tsv("./files/stomach_males_females_signf_degs.txt") %>%
  dplyr::select(-c(adj.P.Val_males, adj.P.Val_females)) %>%
  mutate(tissue = "Stomach")


thyroid_degs <- read_tsv("./files/thyroid_males_females_signf_degs.txt") %>%
  dplyr::select(-c(adj.P.Val_males, adj.P.Val_females)) %>%
  mutate(tissue = "Thyroid")


degs <- bind_rows(stomach_degs, thyroid_degs) %>%
  inner_join(tcga_annotation %>% filter(geneType %in% c("protein_coding", "lincRNA")) %>% dplyr::select(ensemblGeneID, geneType), by = c("genes" = "ensemblGeneID"))



# load enrichment tables
stomach_enr <- read_tsv("./files/stomach_male_female_degs_enr.txt") %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(-c(GeneRatio, BgRatio, geneID)) %>%
  mutate(tissue = "Stomach")

thyroid_enr <- read_tsv("./files/thyroid_male_female_degs_enr.txt") %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(-c(GeneRatio, BgRatio, geneID)) %>%
  mutate(tissue = "Thyroid")


enrch <- bind_rows(stomach_enr, thyroid_enr)



# barplot of number of enriched terms by category
enrch_bp <- enrch %>%
  group_by(tissue, state, Description) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Description = fct_relevel(Description, c("GO_BP", "KEGG", "ONCO", "CM", "IMMUNO", "POS"))) %>%
  ggplot(mapping = aes(x = Description, y = counts, fill = state)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ tissue) +
    scale_fill_manual(values=c("#fdbb84", "#D7301F", "#3182bd"), labels=c("Common", "Female-specific", "Male-specific"), name = "DEGs group") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=10),
      axis.text.y = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    scale_x_discrete(labels = c("GO BP", "KEGG\npath", "Onco\nsets", "Cancer\nmod", "Immuno\nsets", "Position\nsets")) +
    scale_y_continuous(limits = c(NA, 600)) +
    labs(y = "Enriched terms (FDR < 5%)", x = "DEGs group")
ggsave(filename="thyroid_stomach_number_terms.png", plot = enrch_bp, path = "./plots/thyroid_stomach_tumourVSnormal", width=10, height=4)
unlink("thyroid_stomach_number_terms.png")



# barplot of number of genes by state and type
geneType_bp <- degs %>%
  group_by(tissue, state, geneType) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  ggplot(mapping = aes(x = state, y = counts, fill = geneType)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ tissue) +
    scale_fill_manual(values=c("#e41a1c", "#4daf4a"), labels=c("lincRNA", "Protein-coding"), name = "Gene biotype") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=10),
      axis.text.y = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    scale_x_discrete(labels = c("Common", "Female\nspecific", "Male\nspecifc")) +
    labs(y = "Number of genes", x = "DEGs group")
ggsave(filename="thyroid_stomach_number_genes_type.png", plot = geneType_bp, path = "./plots/thyroid_stomach_tumourVSnormal", width=7, height=4)
unlink("thyroid_stomach_number_genes_type.png")




# barplot of number of genes by chromosome and type
chrom_bp <- degs %>%
  group_by(tissue, chrom, state) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  inner_join(chrom_genes, by = "chrom") %>%
  mutate(fraction = counts/total_genes) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(chrom = fct_relevel(
    chrom,
    rev(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY", "chrM")))) %>%
  ggplot(mapping = aes(x = chrom, y = fraction, fill = state)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ tissue) +
    scale_fill_manual(values=c("#fdbb84", "#D7301F", "#3182bd"), labels=c("Common", "Female-specific", "Male-specific"), name = "DEG type") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=13),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    labs(y = "Fraction of genes", x = "Chromosome") +
    coord_flip()
ggsave(filename="thyroid_stomach_number_genes_chrom.png", plot = chrom_bp, path = "./plots/thyroid_stomach_tumourVSnormal", width=8, height=5)
unlink("thyroid_stomach_number_genes_chrom.png")




# boxplot of log2FC distribution by DEG type
logFC_dist <- degs %>%
  dplyr::select(tissue, state, logFC_males, logFC_females) %>%
  gather(-c(tissue, state), key = "gender", value = "logFC") %>%
  na.exclude() %>%
  mutate(signal = if_else(logFC > 0, "pos", "neg")) %>%
  ggplot(mapping = aes(x = state, y = logFC, fill = gender)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(color = "grey", alpha=0.2, width=0.1) +
    facet_wrap(~ tissue) +
    scale_fill_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
    scale_color_manual(values=c("#fbb4b9", "#74a9cf"), guide=F) +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=13),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15),
      legend.position = "bottom") +
    scale_x_discrete(labels = c("Common", "Female\nspecific", "Male\nspecifc")) +
    scale_y_continuous(limits = c(-8, NA)) +
    labs(y = "logFC distribution", x = "DEGs group")
ggsave(filename="thyroid_stomach_degs_logFC_dist.png", plot = logFC_dist, path = "./plots/thyroid_stomach_tumourVSnormal", width=7, height=5)
unlink("thyroid_stomach_degs_logFC_dist.png")




# barplot of number of up-regulated DEGs in tumour/normal tissue by DEG type
logFC_barpl <- degs %>%
  dplyr::select(tissue, state, logFC_males, logFC_females) %>%
  gather(-c(tissue, state), key = "gender", value = "logFC") %>%
  na.exclude() %>%
  mutate(signal = if_else(logFC > 0, "up_tumour", "up_normal")) %>%
  filter(state != "common") %>%
  group_by(tissue, state, gender, signal) %>%
  count() %>%
  ungroup() %>%
  ggplot(mapping = aes(x = state, y = n, fill = signal)) +
    geom_bar(stat = "identity", position = "dodge", lwd=1) +
    facet_wrap(~ tissue) +
    scale_fill_manual(values=c("#a1d76a", "#ca0020"), name="Expression", labels = c("Up-regulated normal", "Up-regulated tumour")) +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=13),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=18),
      legend.position = "bottom") +
    scale_x_discrete(labels = c("Female\nspecific", "Male\nspecifc")) +
    scale_y_continuous(limits = c(0, 500)) +
    labs(y = "Number of genes", x = "DEGs type", title = "All DEGs") +
    guides(fill=guide_legend(nrow=2), color=guide_legend(nrow=2))
ggsave(filename="thyroid_stomach_degs_logFC_barpl.png", plot = logFC_barpl, path = "./plots/thyroid_stomach_tumourVSnormal", width=5, height=5)
unlink("thyroid_stomach_degs_logFC_barpl.png")




# fraction of DEGs that are TSGs
degs_ts <- degs %>%
  inner_join(tsgs %>% dplyr::select(GeneSymbol), by = c("geneName" = "GeneSymbol")) %>%
  group_by(tissue, state) %>%
  count() %>%
  ungroup()

fraction_degs_tsg <- degs %>%
  group_by(tissue, state) %>%
  count() %>%
  ungroup() %>%
  mutate(tsgs = degs_ts$n) %>%
  mutate(fraction = tsgs/n)


# fisher-test
# null hypothesis:
# the number of TSGs by gene group (common, female-specific, male-specific) is similar

thyroid_contingency <- fraction_degs_tsg %>%
  filter(tissue == "Thyroid", state != "common") %>%
  mutate(not_tsgs = n-tsgs) %>%
  dplyr::select(state, tsgs, not_tsgs) %>%
  as.data.frame %>%
  column_to_rownames("state") %>%
  fisher.test()

stomach_contingency <- fraction_degs_tsg %>%
  filter(tissue == "Stomach", state != "common") %>%
  mutate(not_tsgs = n-tsgs) %>%
  dplyr::select(state, tsgs, not_tsgs) %>%
  as.data.frame %>%
  column_to_rownames("state") %>%
  fisher.test()



# tumour supressor gene status of the DEGs
degs2 <- degs %>%
  left_join(tsgs %>% mutate(tsg = 1) %>% dplyr::select(GeneSymbol, tsg), by = c("geneName" = "GeneSymbol")) %>%
  mutate(tsg = replace_na(tsg, 0))
write.table(degs2 %>% filter(tsg == 1), "./files/degs_TvsN_stomach_thyroid_tsgs.txt", sep="\t", quote=F, row.names=F)



# barplot of number of TSGs up-regulated in tumour/normal tissue by DEG type
degs2_bp <- degs2 %>%
  filter(state != "common", tsg == 1) %>%
  dplyr::select(tissue, state, logFC_males, logFC_females) %>%
  gather(-c(tissue, state), key = "gender", value = "logFC") %>%
  na.exclude() %>%
  dplyr::select(-gender) %>%
  mutate(signal = if_else(logFC > 0, "up_tumour", "up_normal")) %>%
  ggplot(mapping = aes(x = state, y = ..count.., fill = signal)) +
    geom_bar(position="dodge") +
    facet_wrap(~ tissue) +
    scale_fill_manual(values=c("#a1d76a", "#ca0020"), name="Expression", labels = c("Up-regulated normal", "Up-regulated tumour")) +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=13),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=18),
      legend.position = "bottom") +
    scale_x_discrete(labels = c("Female\nspecific", "Male\nspecifc")) +
    labs(y = "Number of genes", x = "DEGs type", title = "DEGs with TSG activity") +
    guides(fill=guide_legend(nrow=2))
ggsave(filename="thyroid_stomach_degs_tsgs_barpl.png", plot = degs2_bp, path = "./plots/thyroid_stomach_tumourVSnormal", width=5, height=5)
unlink("thyroid_stomach_degs_tsgs_barpl.png")




# fraction of DEGs that are onco-genes
degs_ocg <- degs %>%
  inner_join(ocgs %>% filter(!(OncogeneName %in% tsgs$GeneSymbol)) %>% dplyr::select(OncogeneName), by = c("geneName" = "OncogeneName")) %>%
  group_by(tissue, state) %>%
  count() %>%
  ungroup()

fraction_degs_ocg <- degs %>%
  group_by(tissue, state) %>%
  count() %>%
  ungroup() %>%
  mutate(ocgs = degs_ocg$n) %>%
  mutate(fraction = ocgs/n)


# fisher-test
# null hypothesis:
# the number of OCGs by gene group (common, female-specific, male-specific) is similar

thyroid_contingency2 <- fraction_degs_ocg %>%
  filter(tissue == "Thyroid", state != "common") %>%
  mutate(not_ocgs = n-ocgs) %>%
  dplyr::select(state, ocgs, not_ocgs) %>%
  as.data.frame %>%
  column_to_rownames("state") %>%
  fisher.test()

stomach_contingency2 <- fraction_degs_ocg %>%
  filter(tissue == "Stomach", state != "common") %>%
  mutate(not_ocgs = n-ocgs) %>%
  dplyr::select(state, ocgs, not_ocgs) %>%
  as.data.frame %>%
  column_to_rownames("state") %>%
  fisher.test()


# tumour supressor gene status of the DEGs
degs3 <- degs %>%
  left_join(ocgs %>% filter(!(OncogeneName %in% tsgs$GeneSymbol)) %>% mutate(ocg = 1) %>% dplyr::select(OncogeneName, ocg), by = c("geneName" = "OncogeneName")) %>%
  mutate(ocg = replace_na(ocg, 0))



# barplot of number of TSGs up-regulated in tumour/normal tissue by DEG type
degs3_bp <- degs3 %>%
  filter(state != "common") %>%
  filter(ocg == 1) %>%
  dplyr::select(tissue, state, logFC_males, logFC_females) %>%
  gather(-c(tissue, state), key = "gender", value = "logFC") %>%
  na.exclude() %>%
  dplyr::select(-gender) %>%
  mutate(signal = if_else(logFC > 0, "up_tumour", "up_normal")) %>%
  ggplot(mapping = aes(x = state, y = ..count.., fill = signal)) +
    geom_bar(position="dodge") +
    facet_wrap(~ tissue) +
    scale_fill_manual(values=c("#66c2a5", "#8da0cb"), name="OCG state", labels = c("Up-regulated normal", "Up-regulated tumour")) +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=13),
      axis.text.y = element_text(colour="black", size=12),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15),
      legend.position = "bottom") +
    scale_x_discrete(labels = c("Female\nspecific", "Male\nspecifc")) +
    labs(y = "Number of OCGs", x = "DEGs group") +
    guides(fill=guide_legend(nrow=2))
ggsave(filename="thyroid_stomach_degs_ocgs_barpl.png", plot = degs3_bp, path = "./plots/thyroid_stomach_tumourVSnormal", width=4, height=5)
unlink("thyroid_stomach_degs_ocgs_barpl.png")





# lncRNAs description
lncrnas <- read_tsv("./data/lncRNA_cancer/all_lnc_cancer_associations.txt") %>%
  dplyr::rename(gene_name = `LncRNA name`, cancer_type = `cancer type`,
    function_description = `function description`, pubmed_id = `pubmed id`,
    ICD_0_3 = `ICD-0-3`, ICD_0_3_1 = `ICD-0-3_1`)

degs_lncrnas <- degs %>%
  filter(geneType == "lincRNA") %>%
  dplyr::select(tissue, state, geneName) %>%
  inner_join(lncrnas, by = c("geneName" = "gene_name")) %>%
  dplyr::select(tissue, state, geneName, cancer_type) %>%
  distinct() %>%
  group_by(tissue, state, geneName) %>%
  count() %>%
  ungroup()


# sex steroid hormone receptors
horm_recpt <- read.gmt("./data/gene_lists/c5.mf.v6.2.symbols.gmt") %>%
  filter(ont == "GO_STEROID_HORMONE_RECEPTOR_ACTIVITY") %>%
  pull(gene)


degs_recp <- degs %>%
  group_by(tissue, state) %>%
  summarise(recept = length(intersect(geneName, horm_recpt))) %>%
  ungroup()


save(list=ls(), file="r_workspaces/thyroid_stomach_tcga_characterization.RData")
