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
  group_by(moduleL, state) %>%
  summarise(genes=n())

females_tumour_modules <- read_tsv("./files/thyroid_females_tumour_modules.txt")

females_tumour_modules_info <- females_tumour_modules %>%
  group_by(moduleL, state) %>%
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

  enr <- enr@result %>% dplyr::select(Description, ID, Count, p.adjust)
  return(list(enr))

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





# TUMOUR female-specific modules

females_tumour_modules_enr <- females_tumour_modules %>%
  filter(state == "female-specific") %>%
  dplyr::select(moduleL, geneName) %>%
  group_by(moduleL) %>%
  mutate(geneName = list(geneName)) %>%
  unique() %>%
  rowwise() %>%
  mutate(enr = enrich_test(geneName, universe_enr, all_terms, all_terms2, 1, 1)) %>%
  dplyr::select(-geneName) %>%
  unnest() %>%
  filter(p.adjust < 0.05) %>%
  group_by(moduleL, Description) %>% top_n(-10, p.adjust)
