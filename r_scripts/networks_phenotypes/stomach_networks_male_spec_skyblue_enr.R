# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Enrichment analysis of male-specific module


# -- Stomach


library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(viridis)


# tumour networks
males_tumour_modules <- read_tsv("./files/stomach_males_tumour_modules.txt") %>%
  mutate(tissue = "tumour", sex = "male", state2 = if_else(state == "male-specific", "male-specific", "conserved"))


# hypergeometric test
all_modules_enr <- read_tsv("./files/stomach_modules_hypergeo_enr.txt")
all_modules_gsea <- read_tsv("./files/stomach_modules_gsea_enr.txt")




# plotting enrichment results
skyblue_mod_hyp <- all_modules_enr %>%
  filter(tissue == "tumour", state == "male-specific", moduleL == "skyblue", p.adjust < 0.05) %>%
  filter(Description != "POS") %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  group_by(Description) %>%
  top_n(5, log10_p) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  theme_classic() +
  facet_wrap( ~ Description,
    scales = "free",
    #space = "free_y",
    labeller=labeller(Description = c("GO_BP" = "GO BP"))) +
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
ggsave(filename="stomach_male_spc_skyblue_hyp.png", plot=skyblue_mod_hyp, path="./plots/wgcna_networks_traits/", width = 6, height = 3)
ggsave(filename="stomach_male_spc_skyblue_hyp.pdf", plot=skyblue_mod_hyp, path="./plots/wgcna_networks_traits/", width = 6, height = 3)
unlink("stomach_male_spc_skyblue_hyp.png")
unlink("stomach_male_spc_skyblue_hyp.pdf")
