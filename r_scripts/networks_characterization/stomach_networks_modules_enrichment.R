# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# enrichment analysis


# -- Stomach


library(tidyverse)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)

source("./r_scripts/utils.R")

# tumour networks
males_tumour_modules <- read_tsv("./files/stomach_males_tumour_modules.txt") %>%
  mutate(tissue = "tumour", sex = "male", state2 = if_else(state == "male-specific", "male-specific", "conserved"))

females_tumour_modules <- read_tsv("./files/stomach_females_tumour_modules.txt") %>%
  mutate(tissue = "tumour", sex = "female", state2 = if_else(state == "female-specific", "female-specific", "conserved"))


# normal networks
males_normal_modules <- read_tsv("./files/stomach_males_normal_modules.txt") %>%
  mutate(tissue = "normal", sex = "male", state2 = if_else(state == "male-specific", "male-specific", "conserved"))

females_normal_modules <- read_tsv("./files/stomach_females_normal_modules.txt") %>%
  mutate(tissue = "normal", sex = "female", state2 = if_else(state == "female-specific", "female-specific", "conserved"))



stomach_nets <- rbind(males_tumour_modules, females_tumour_modules, males_normal_modules, females_normal_modules)


# barplot of number of modules by network
stomach_nets_bp <- stomach_nets %>%
  group_by(tissue, sex, moduleL, state, state2) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(state = fct_relevel(state, c("highly-conserved", "moderately-conserved", "lowly-conserved", "male-specific"))) %>%
  ggplot(mapping = aes(x=sex, y = ..count.., fill = state)) +
  geom_bar(position = "dodge") +
  theme_classic() +
  facet_wrap( ~ tissue,
    labeller=labeller(tissue = c("tumour" = "Tumour", "normal" = "Normal"))) +
  theme(
    axis.title=element_text(colour="black", size=14),
    axis.text=element_text(colour="black", size=12),
    strip.text = element_text(colour="black", size=14),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=11),
    legend.title = element_text(colour="black", size=13),
    legend.position = "bottom") +
  labs(x = "Gender", y = "Count", fill = "Preservation") +
  scale_x_discrete(labels = c("Female", "Male")) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1", "#deebf7", "orange"), labels = c("highly-preserved", "moderately-preserved", "lowly-preserved", "male-specific")) +
  scale_y_continuous(limits = c(NA,15)) +
  guides(fill=guide_legend(nrow=2))
ggsave(filename="stomach_nets_number_modules.png", plot=stomach_nets_bp, path="./plots/wgcna_networks/", width = 5, height = 3)
ggsave(filename="stomach_nets_number_modules.pdf", plot=stomach_nets_bp, path="./plots/wgcna_networks/", width = 5, height = 3)
unlink("stomach_nets_number_modules.png")
unlink("stomach_nets_number_modules.pdf")



# hypergeometric test of GO terms and KEGG pathways

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


universe_enr <- read.table("./files/stomach_males_tumour_fpkm_wgcna.txt", sep="\t", h=T) %>% colnames


# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))

universe_enr <- tcga.geneIDs.annot[tcga.geneIDs.annot$geneID %in% universe_enr, "geneName"]





# enrichment analysis

# hypergeometric test
all_modules_enr <- stomach_nets %>%
  dplyr::select(tissue, sex, moduleL, state, state2, geneName) %>%
  group_by(tissue, sex, moduleL, state, state2) %>%
  summarise(geneName = list(geneName)) %>%
  ungroup() %>%
  mutate(enr = map(.x = geneName, .f = enrich_test, universe = universe_enr, terms1 = all_terms, terms2 = all_terms2, p_adj=1, q_value=1)) %>%
  dplyr::select(-geneName) %>%
  unnest(cols = c(enr))
write.table(all_modules_enr, "./files/stomach_modules_hypergeo_enr.txt", sep="\t", quote=F, row.names=F)



# gsea test (kme as sorting variable)
all_modules_gsea <- stomach_nets %>%
  dplyr::select(tissue, sex, moduleL, state, state2, kme, geneName) %>%
  group_by(tissue, sex, moduleL, state, state2) %>%
  dplyr::arrange(desc(kme), .by_group = TRUE) %>%
  group_by(tissue, sex, moduleL, state, state2) %>%
  nest() %>%
  ungroup() %>%
  mutate(enr = map(.x = data, .f = gsea_test, terms1 = all_terms, terms2 = all_terms2, p_adj=1)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(enr))
write.table(all_modules_gsea, "./files/stomach_modules_gsea_enr.txt", sep="\t", quote=F, row.names=F)



all_modules_enr_summary <- all_modules_enr %>%
  filter(p.adjust < 0.05) %>%
  group_by(tissue, sex, state, state2, Description, ID) %>%
  dplyr::summarise(mod = paste(moduleL, collapse=" ")) %>%
  ungroup()




# plotting enrichment results

# tumour_gender_sp_mod_hyp <- all_modules_enr %>%
#   filter(p.adjust < 0.05 & tissue == "tumour" & (state2 == "female-specific" | state2 == "male-specific" ) ) %>%
#   filter(Description %in% c("GO_BP", "CM", "KEGG", "POS") ) %>%
#   mutate(log10_p = -log10(p.adjust)) %>%
#   group_by(state2, Description) %>%
#   top_n(5, log10_p) %>%
#   ungroup() %>%
#   mutate_if(is.character, as.factor) %>%
#   mutate(ID = fct_reorder(ID, Count)) %>%
#   ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
#   geom_bar(stat="identity") +
#   theme_classic() +
#   facet_grid(Description ~ state2,
#     scales = "free",
#     space = "free_y",
#     labeller=labeller(
#       state2 = c("male-specific" = "Male-specific", "conserved" = "Conserved"),
#       Description = c("GO_BP" = "GO BP", "KEGG" = "KEGG", "ONCO" = "ONCOGENIC\ngene sets", "POS" = "Pos", "CM" = "Cancer\nmodules"))) +
#   theme(
#     axis.title.x=element_text(colour="black", size=16),
#     axis.title.y=element_blank(),
#     axis.text.y=element_text(colour="black", size=13),
#     axis.text.x=element_text(colour="black", size=14),
#     strip.text = element_text(colour="black", size=15),
#     strip.background = element_blank(),
#     legend.text = element_text(colour="black", size=14),
#     legend.title = element_text(colour="black", size=16)) +
#   coord_flip() +
#   scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
#   scale_y_continuous(name = "Number of genes") +
#   labs(color="Sex")
# ggsave(filename="stomach_tumour_male_female_spc_hyp.png", plot=tumour_gender_sp_mod_hyp, path="./plots/wgcna_modules_enrichment/", width = 10, height = 5)
# ggsave(filename="stomach_tumour_male_female_spc_hyp.pdf", plot=tumour_gender_sp_mod_hyp, path="./plots/wgcna_modules_enrichment/", width = 10, height = 5)
# unlink("stomach_tumour_male_female_spc_hyp.png")
# unlink("stomach_tumour_male_female_spc_hyp.pdf")
#
#
# normal_gender_sp_mod_hyp <- all_modules_enr %>%
#   filter(p.adjust < 0.05, tissue == "normal") %>%
#   mutate(log10_p = -log10(p.adjust)) %>%
#   group_by(state2, Description) %>%
#   top_n(5, log10_p) %>%
#   ungroup() %>%
#   mutate_if(is.character, as.factor) %>%
#   mutate(ID = fct_reorder(ID, Count)) %>%
#   ggplot(mapping = aes(x=ID, y = Count, fill = log10_p, color = sex)) +
#   geom_bar(stat="identity") +
#   theme_classic() +
#   facet_grid(Description ~ state2,
#     scales = "free",
#     space = "free_y",
#     labeller=labeller(
#       state2 = c("conserved" = "Conserved"),
#       Description = c("GO_BP" = "GO BP", "KEGG" = "KEGG\npathways", "ONCO" = "ONCOGENIC\ngene sets", "POS" = "Positional\ngene sets", "CM" = "Cancer\nmodules"))) +
#   theme(
#     axis.title.x=element_text(colour="black", size=16),
#     axis.title.y=element_blank(),
#     axis.text.y=element_text(colour="black", size=13),
#     axis.text.x=element_text(colour="black", size=14),
#     strip.text = element_text(colour="black", size=15),
#     strip.background = element_blank(),
#     legend.text = element_text(colour="black", size=14),
#     legend.title = element_text(colour="black", size=16)) +
#   coord_flip() +
#   scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
#   scale_y_continuous(name = "Number of genes") +
#   labs(color="Sex")
# ggsave(filename="stomach_normal_male_female_spc_hyp.png", plot=normal_gender_sp_mod_hyp, path="./plots/wgcna_modules_enrichment/", width = 14, height = 10)
# ggsave(filename="stomach_normal_male_female_spc_hyp.pdf", plot=normal_gender_sp_mod_hyp, path="./plots/wgcna_modules_enrichment/", width = 14, height = 10)
# unlink("stomach_normal_male_female_spc_hyp.png")
# unlink("stomach_normal_male_female_spc_hyp.pdf")
#
#
# tumour_gender_sp_mod_gsea <- all_modules_gsea %>%
#   filter(p.adjust < 0.05 & tissue == "tumour" & (state2 == "female-specific" | state2 == "male-specific" ) ) %>%
#   mutate(log10_p = -log10(p.adjust)) %>%
#   group_by(state2, Description) %>%
#   top_n(5, enrichmentScore) %>%
#   ungroup() %>%
#   mutate_if(is.character, as.factor) %>%
#   mutate(ID = fct_reorder(ID, enrichmentScore)) %>%
#   ggplot(mapping = aes(x=ID, y = enrichmentScore, fill = log10_p)) +
#   geom_bar(stat="identity") +
#   theme_classic() +
#   facet_grid(Description ~ state2,
#     scales = "free",
#     space = "free_y",
#     labeller=labeller(
#       state2 = c("male-specific" = "Male-specific", "conserved" = "Conserved"),
#       Description = c("IMMUNO" = "Immuno", "GO_BP" = "GO BP", "KEGG" = "KEGG\npathways", "ONCO" = "ONCOGENIC\ngene sets", "POS" = "Positional\ngene sets", "CM" = "Cancer\nmodules"))) +
#   theme(
#     axis.title.x=element_text(colour="black", size=15),
#     axis.title.y=element_blank(),
#     axis.text.y=element_text(colour="black", size=12),
#     axis.text.x=element_text(colour="black", size=13),
#     strip.text = element_text(colour="black", size=14),
#     strip.background = element_blank(),
#     legend.text = element_text(colour="black", size=13),
#     legend.title = element_text(colour="black", size=15)) +
#   coord_flip() +
#   scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
#   scale_y_continuous(name = "Enrichment score") +
#   labs(color="Sex")
# ggsave(filename="stomach_tumour_male_female_spc_gsea.png", plot=tumour_gender_sp_mod_gsea, path="./plots/wgcna_modules_enrichment/", width = 10, height = 4)
# ggsave(filename="stomach_tumour_male_female_spc_gsea.pdf", plot=tumour_gender_sp_mod_gsea, path="./plots/wgcna_modules_enrichment/", width = 10, height = 4)
# unlink("stomach_tumour_male_female_spc_gsea.png")
# unlink("stomach_tumour_male_female_spc_gsea.pdf")
#
#
# normal_gender_sp_mod_gsea <- all_modules_gsea %>%
#   filter(p.adjust < 0.05, tissue == "normal") %>%
#   mutate(log10_p = -log10(p.adjust)) %>%
#   group_by(state2, Description) %>%
#   top_n(5, enrichmentScore) %>%
#   ungroup() %>%
#   mutate_if(is.character, as.factor) %>%
#   mutate(ID = fct_reorder(ID, enrichmentScore)) %>%
#   ggplot(mapping = aes(x=ID, y = enrichmentScore, fill = log10_p, color = sex)) +
#   geom_bar(stat="identity") +
#   theme_classic() +
#   facet_grid(Description ~ state2,
#     scales = "free",
#     space = "free_y",
#     labeller=labeller(
#       state2 = c("conserved" = "Conserved"),
#       Description = c("GO_BP" = "GO BP", "KEGG" = "KEGG\npathways", "ONCO" = "ONCOGENIC\ngene sets", "POS" = "Positional\ngene sets", "CM" = "Cancer\nmodules"))) +
#   theme(
#     axis.title.x=element_text(colour="black", size=15),
#     axis.title.y=element_blank(),
#     axis.text.y=element_text(colour="black", size=12),
#     axis.text.x=element_text(colour="black", size=13),
#     strip.text = element_text(colour="black", size=14),
#     strip.background = element_blank(),
#     legend.text = element_text(colour="black", size=13),
#     legend.title = element_text(colour="black", size=15)) +
#   coord_flip() +
#   scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
#   scale_y_continuous(name = "Enrichment score") +
#   labs(color="Sex")
# ggsave(filename="stomach_normal_male_female_spc_gsea.png", plot=normal_gender_sp_mod_gsea, path="./plots/wgcna_modules_enrichment/", width = 11, height = 7)
# ggsave(filename="stomach_normal_male_female_spc_gsea.pdf", plot=normal_gender_sp_mod_gsea, path="./plots/wgcna_modules_enrichment/", width = 11, height = 7)
# unlink("stomach_normal_male_female_spc_gsea.png")
# unlink("stomach_normal_male_female_spc_gsea.pdf")





save(list=ls(), file="r_workspaces/stomach_networks_modules_enrichment.RData")
