# Understanding Gender Differential Susceptibility in Cancer


# Enrichment comparison between gender-specific DEGs in thyroid and stomach cancers


library(tidyverse)
library(RColorBrewer)
library(viridis)



# load enrichment data
stomach_enr <- read_tsv("./files/stomach_male_female_degs_enr.txt") %>%
  filter(p.adjust < 0.05, state != "common") %>%
  dplyr::select(-c(GeneRatio, BgRatio)) %>%
  mutate(log10_p = -log10(p.adjust), tissue = "Stomach")

thyroid_enr <- read_tsv("./files/thyroid_male_female_degs_enr.txt") %>%
  filter(p.adjust < 0.05, state != "common") %>%
  dplyr::select(-c(GeneRatio, BgRatio)) %>%
  mutate(log10_p = -log10(p.adjust), tissue = "Thyroid")


common_protected <- intersect(thyroid_enr %>% filter(state == "male_specific") %>% pull(ID), stomach_enr %>% filter(state == "female_specific") %>% pull(ID))
common_susceptible <- intersect(thyroid_enr %>% filter(state == "female_specific") %>% pull(ID), stomach_enr %>% filter(state == "male_specific") %>% pull(ID))

comp_protected <- bind_rows(
  thyroid_enr %>% filter(state == "male_specific" & ID %in% common_protected),
  stomach_enr %>% filter(state == "female_specific" & ID %in% common_protected)) %>%
  arrange(ID) %>%
  mutate(suscpt = "less_affected")


comp_susceptible <- bind_rows(
  thyroid_enr %>% filter(state == "female_specific" & ID %in% common_susceptible),
  stomach_enr %>% filter(state == "male_specific" & ID %in% common_susceptible)) %>%
  arrange(ID) %>%
  mutate(suscpt = "more_affected")



comparison_n <- bind_rows(comp_protected, comp_susceptible) %>%
  filter(Description != "IMMUNO") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count)) %>%
  group_by(ID) %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  mutate(geneID = list(Reduce(intersect, geneID)) ) %>%
  ungroup() %>%
  dplyr::select(Description, ID, geneID, suscpt) %>%
  unique() %>%
  mutate(n = sapply(geneID, length))
  #group_by(ID) %>% mutate(n = length(geneID[[1]]))
  #rowwise() %>% mutate(n = length(geneID))


comparison <- bind_rows(comp_protected, comp_susceptible) %>%
  mutate(Description = str_replace_all(Description, "_", " ")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_reorder(ID, Count)) %>%
  group_by(ID) %>%
  mutate(geneID2 = str_split(geneID, "/")) %>%
  mutate(geneID2 = list(Reduce(intersect, geneID2)) ) %>%
  ungroup() %>%
  mutate(n = sapply(geneID2, length)) %>%
  ggplot(mapping = aes(x=ID, y = Count, fill = log10_p)) +
  geom_bar(stat="identity") +
  #geom_text(aes(label=Description, y=7), color="white") +
  geom_text(aes(label=n), vjust=0.5, hjust=1.5, color="white") +
  theme_classic() +
  facet_grid(suscpt ~ tissue,
    scales = "free",
    space = "free_y",
    labeller=labeller(
      suscpt = c("less_affected" = "Less", "more_affected" = "More"))) +
  theme(
    axis.title.x=element_text(colour="black", size=15),
    axis.title.y=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_text(colour="black", size=13),
    plot.title = element_blank(),
    strip.text.y = element_text(colour="black", size=12),
    strip.text.x = element_text(colour="black", size=18),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=13),
    legend.title = element_text(colour="black", size=15)) +
  coord_flip() +
  scale_fill_viridis(option="D", name="Adj p-val\n(-log10)") +
  scale_y_continuous(name = "Number of genes")
ggsave(filename="stomach_thyroid_common_degs_susceptibility.png", plot=comparison, path="./plots/thyroid_stomach_tumourVSnormal/", width = 10, height = 2)
ggsave(filename="stomach_thyroid_common_degs_susceptibility.pdf", plot=comparison, path="./plots/thyroid_stomach_tumourVSnormal/", width = 10, height = 2)
unlink("stomach_thyroid_common_degs_susceptibility.png")
unlink("stomach_thyroid_common_degs_susceptibility.pdf")



# load degs table
stomach_degs <- read_tsv("./files/stomach_males_females_signf_degs.txt") %>%
  filter(state != "common") %>%
  arrange(state) %>%
  mutate(logFC = c(na.exclude(logFC_females), na.exclude(logFC_males))) %>%
  dplyr::select(geneName, state, logFC) %>%
  mutate(tissue = "Stomach")


thyroid_degs <- read_tsv("./files/thyroid_males_females_signf_degs.txt") %>%
  filter(state != "common") %>%
  arrange(state) %>%
  mutate(logFC = c(na.exclude(logFC_females), na.exclude(logFC_males))) %>%
  dplyr::select(geneName, state, logFC) %>%
  mutate(tissue = "Thyroid")


stomach_thyroid_degs <- bind_rows(stomach_degs, thyroid_degs)


comp_distribution <- comparison$data %>%
  dplyr::select(-geneID2, -n) %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest() %>%
  dplyr::select(ID, tissue, suscpt, state, geneID) %>%
  inner_join(stomach_thyroid_degs, by=c("tissue", "state", "geneID" = "geneName")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_infreq(ID)) %>%
  ggplot(mapping = aes(x=tissue, y = logFC, fill = tissue)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  theme_classic() +
  facet_grid( ~ ID,
    scales = "free",
    labeller=labeller(
      ID = c("CELLULAR RESPONSE TO CYTOKINE STIMULUS" = "CELL RESP\nCYTOKINES",
      "CELL CELL ADHESION" = "CELL CELL\nADHESION",
      "SINGLE ORGANISM CELL ADHESION" = "SINGLE ORG\nADHESION"))) +
  theme(
    axis.title=element_text(colour="black", size=18),
    axis.text.y=element_text(colour="black", size=15),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.title = element_blank(),
    strip.text = element_text(colour="black", size=10),
    strip.background = element_blank(),
    legend.text = element_text(colour="black", size=16),
    legend.title = element_text(colour="black", size=18)) +
  scale_y_continuous(name = "log2FC") +
  scale_x_discrete(name = "Enriched term") +
  scale_fill_manual(name = "Cancer", values = c("#d95f02", "#7570b3"))
ggsave(filename="stomach_thyroid_comp_distribution.png", plot=comp_distribution, path="./plots/thyroid_stomach_tumourVSnormal/", width = 5, height = 5)
ggsave(filename="stomach_thyroid_comp_distribution.pdf", plot=comp_distribution, path="./plots/thyroid_stomach_tumourVSnormal/", width = 5, height = 5)
unlink("stomach_thyroid_comp_distribution.png")
unlink("stomach_thyroid_comp_distribution.pdf")





save(list=ls(), file="r_workspaces/comparison_enr_thca_stad.RData")
