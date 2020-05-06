# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Enrichment analysis of female-specific module


# -- Thyroid


library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(viridis)


# tumour networks
females_tumour_modules <- read_tsv("./files/thyroid_females_tumour_modules.txt") %>%
  mutate(tissue = "tumour", sex = "female", state2 = if_else(state == "female-specific", "female-specific", "conserved"))


# hypergeometric test
all_modules_enr <- read_tsv("./files/thyroid_modules_hypergeo_enr.txt")
all_modules_gsea <- read_tsv("./files/thyroid_modules_gsea_enr.txt")




# plotting enrichment results
greenyellow_mod_hyp <- all_modules_enr %>%
  filter(tissue == "tumour", state == "lowly-conserved", moduleL == "greenyellow", p.adjust < 0.05) %>%
  filter(Description %in% c("GO_BP", "KEGG") ) %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(Description) %>%
  top_n(5, log10_p) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count)) %>%
  filter(str_detect(ID, "DISEASE", negate = T)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_grid(Description ~ .,
    scales = "free",
    space = "free_y",
    labeller=labeller(
      state = c("female-specific" = "Female-specific"),
      Description = c("GO_BP" = "GO BP", "KEGG" = "KEGG"))) +
    theme(
      axis.title.x=element_text(colour="black", size=14),
      axis.title.y=element_blank(),
      axis.text.y=element_text(colour="black", size=12),
      axis.text.x=element_text(colour="black", size=12),
      strip.text = element_text(colour="black", size=14),
      strip.background = element_blank(),
      legend.text = element_text(colour="black", size=12),
      legend.title = element_text(colour="black", size=14)) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adjusted P-value\n(-log10)") +
  scale_y_continuous(name = "Count")
ggsave(filename="thyroid_female_spc_greenyellow_hyp.png", plot=greenyellow_mod_hyp, path="./plots/wgcna_networks_traits/", width = 9, height = 3)
ggsave(filename="thyroid_female_spc_greenyellow_hyp.pdf", plot=greenyellow_mod_hyp, path="./plots/wgcna_networks_traits/", width = 9, height = 3)
unlink("thyroid_female_spc_greenyellow_hyp.png")
unlink("thyroid_female_spc_greenyellow_hyp.pdf")



greenyellow_mod_gsea <- all_modules_gsea %>%
  filter(tissue == "tumour", state == "lowly-conserved", moduleL == "greenyellow", p.adjust < 0.05) %>%
  filter(Description %in% c("GO_BP", "KEGG") ) %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(Description) %>%
  top_n(5, enrichmentScore) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, enrichmentScore)) %>%
  filter(str_detect(ID, "DISEASE", negate = T)) %>%
  ggplot(mapping = aes(x=ID, y = enrichmentScore, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_grid(Description ~ .,
    scales = "free",
    space = "free_y",
    labeller=labeller(
      state2 = c("female-specific" = "Female-specific"),
      Description = c("GO_BP" = "GO BP", "KEGG" = "KEGG"))) +
  theme(
    axis.title.x=element_text(colour="black", size=15),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_text(colour="black", size=13),
    strip.text = element_text(colour="black", size=14),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=13),
    legend.title = element_text(colour="black", size=15)) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
  scale_y_continuous(name = "Enrichment score")
ggsave(filename="thyroid_female_spc_greenyellow_gsea.png", plot=greenyellow_mod_gsea, path="./plots/wgcna_networks_traits/", width = 9, height = 3)
ggsave(filename="thyroid_female_spc_greenyellow_gsea.pdf", plot=greenyellow_mod_gsea, path="./plots/wgcna_networks_traits/", width = 9, height = 3)
unlink("thyroid_female_spc_greenyellow_gsea.png")
unlink("thyroid_female_spc_greenyellow_gsea.pdf")



# Gene Set Enrichment Analysis on the greenyellow module

# load gene lists
kegg <- read.gmt("./data/gene_lists/c2.cp.kegg.v6.2.symbols.gmt") %>% mutate(ont = str_replace(ont, "KEGG_", ""))
kegg2 <- data.frame(ont = kegg$ont, name = "KEGG") %>% dplyr::distinct()

all_terms <- kegg %>% mutate(ont = str_replace_all(ont, "_", " "))
all_terms2 <- kegg2 %>% mutate(ont = str_replace_all(ont, "_", " "))



# gsea test (kme as sorting variable)
greenyellow_kme_sorted <- females_tumour_modules %>%
  filter(tissue == "tumour", moduleL == "greenyellow") %>%
  dplyr::arrange(desc(kme))

glist <- greenyellow_kme_sorted$kme
names(glist) <- greenyellow_kme_sorted$geneName

enr <- GSEA(
  geneList = glist,
  TERM2GENE = all_terms,
  TERM2NAME = all_terms2,
  exponent = 1,
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = FALSE,
  by = "fgsea")



greenyellow_gsea_ox_phospho <- gseaplot(enr, geneSetID = "OXIDATIVE PHOSPHORYLATION",  by = "runningScore", title = "OXIDATIVE PHOSPHORYLATION", color = "black", color.line = "darkblue", color.vline = "grey60") +
  theme_classic() +
  theme(
    plot.title=element_blank(),
    axis.title=element_text(colour="black", size=20),
    axis.text=element_text(colour="black", size=18)) +
  labs(x = "Ranked list of genes", y = "ES")
ggsave(filename="greenyellow_gsea_ox_phospho.png", plot=greenyellow_gsea_ox_phospho, path="./plots/wgcna_networks_traits/", width = 5, height = 3)
ggsave(filename="greenyellow_gsea_ox_phospho.pdf", plot=greenyellow_gsea_ox_phospho, path="./plots/wgcna_networks_traits/", width = 5, height = 3)
unlink("greenyellow_gsea_ox_phospho.png")
unlink("greenyellow_gsea_ox_phospho.pdf")
