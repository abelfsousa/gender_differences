# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# enrichment analysis


# -- Thyroid


library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(viridis)



# tumour networks

males_tumour_modules <- read_tsv("./files/thyroid_males_tumour_modules.txt")

males_tumour_modules_info <- males_tumour_modules %>%
  group_by(moduleL, state2) %>%
  summarise(genes=n())

females_tumour_modules <- read_tsv("./files/thyroid_females_tumour_modules.txt")

females_tumour_modules_info <- females_tumour_modules %>%
  group_by(moduleL, state2) %>%
  summarise(genes=n())




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

all_terms <- bind_rows(kegg, go_bp, onco, immuno) %>% mutate(ont = str_replace_all(ont, "_", " "))
all_terms2 <- bind_rows(kegg2, go_bp2, onco2, immuno2) %>% mutate(ont = str_replace_all(ont, "_", " "))


universe_enr <- read.table("./files/thyroid_males_tumour_fpkm_wgcna.txt", sep="\t", h=T) %>% colnames


# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))

universe_enr <- tcga.geneIDs.annot[tcga.geneIDs.annot$geneID %in% universe_enr, "geneName"]






# TUMOUR male-specific modules
black_tumour_male <- males_tumour_modules %>% filter(moduleL == "black") %>% pull(geneName)
black_tumour_male <- enrich_test(black_tumour_male, universe_enr, all_terms, all_terms2, 1, 1)
black_tumour_male@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

blue_tumour_male <- males_tumour_modules %>% filter(moduleL == "blue") %>% pull(geneName)
blue_tumour_male <- enrich_test(blue_tumour_male, universe_enr, all_terms, all_terms2, 1, 1)
blue_tumour_male@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

green_tumour_male <- males_tumour_modules %>% filter(moduleL == "green") %>% pull(geneName)
green_tumour_male <- enrich_test(green_tumour_male, universe_enr, all_terms, all_terms2, 1, 1)
green_tumour_male@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)


# TUMOUR female-specific modules
# most important: cyan greenyellow lightgreen lightcyan pink salmon tan
cyan_tumour_fem <- females_tumour_modules %>% filter(moduleL == "cyan") %>% pull(geneName)
cyan_tumour_fem <- enrich_test(cyan_tumour_fem, universe_enr, all_terms, all_terms2, 1, 1)
cyan_tumour_fem@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

greenyellow_tumour_fem <- females_tumour_modules %>% filter(moduleL == "greenyellow") %>% pull(geneName)
greenyellow_tumour_fem <- enrich_test(greenyellow_tumour_fem, universe_enr, all_terms, all_terms2, 1, 1)
greenyellow_tumour_fem@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

lightgreen_tumour_fem <- females_tumour_modules %>% filter(moduleL == "lightgreen") %>% pull(geneName)
lightgreen_tumour_fem <- enrich_test(lightgreen_tumour_fem, universe_enr, all_terms, all_terms2, 1, 1)
lightgreen_tumour_fem@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

lightcyan_tumour_fem <- females_tumour_modules %>% filter(moduleL == "lightcyan") %>% pull(geneName)
lightcyan_tumour_fem <- enrich_test(lightcyan_tumour_fem, universe_enr, all_terms, all_terms2, 1, 1)
lightcyan_tumour_fem@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

pink_tumour_fem <- females_tumour_modules %>% filter(moduleL == "pink") %>% pull(geneName)
pink_tumour_fem <- enrich_test(pink_tumour_fem, universe_enr, all_terms, all_terms2, 1, 1)
pink_tumour_fem@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

salmon_tumour_fem <- females_tumour_modules %>% filter(moduleL == "salmon") %>% pull(geneName)
salmon_tumour_fem <- enrich_test(salmon_tumour_fem, universe_enr, all_terms, all_terms2, 1, 1)
salmon_tumour_fem@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

tan_tumour_fem <- females_tumour_modules %>% filter(moduleL == "tan") %>% pull(geneName)
tan_tumour_fem <- enrich_test(tan_tumour_fem, universe_enr, all_terms, all_terms2, 1, 1)
tan_tumour_fem@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

green_tumour_fem <- females_tumour_modules %>% filter(moduleL == "green") %>% pull(geneName)
green_tumour_fem <- enrich_test(green_tumour_fem, universe_enr, all_terms, all_terms2, 1, 1)
green_tumour_fem@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)

grey60_tumour_fem <- females_tumour_modules %>% filter(moduleL == "grey60") %>% pull(geneName)
grey60_tumour_fem <- enrich_test(grey60_tumour_fem, universe_enr, all_terms, all_terms2, 1, 1)
grey60_tumour_fem@result %>% filter(p.adjust < 0.05) %>% dplyr::select(Description, ID, Count, p.adjust) %>% filter(Description == "GO_BP") %>% head(10)



# TUMOUR common modules in males and females
