# Understanding Gender Differential Susceptibility in Cancer


# Gene expression correlation
# Between hormone receptor genes and degs
# Characterization


library(tidyverse)
library(RColorBrewer)
library(viridis)



# load degs~receptor correlation
degs_recpt_cor <- read_tsv("./files/degs_tumourVSnormal_cor_receptors.txt") %>%
  dplyr::select(tissue, state, gene, receptor, r = estimate)



# load degs tables
stomach_degs <- read_tsv("./files/stomach_males_females_signf_degs.txt") %>%
  dplyr::select(gene = geneName, everything()) %>%
  mutate(tissue = "Stomach")


thyroid_degs <- read_tsv("./files/thyroid_males_females_signf_degs.txt") %>%
  dplyr::select(gene = geneName, everything()) %>%
  mutate(tissue = "Thyroid")


degs <- bind_rows(stomach_degs, thyroid_degs) %>%
  dplyr::select(tissue, everything())


degs_count <- degs %>%
  group_by(tissue, state) %>%
  summarise(counts = n()) %>%
  ungroup()



# number of genes with at least one receptor correlated > 0.6
degs_recpt_cor_n <- degs_recpt_cor %>%
  filter(abs(r) > 0.6) %>%
  #filter( !(tissue == "Stomach" & receptor %in% (degs %>% filter(tissue == "Stomach") %>% pull(gene))) ) %>%
  #filter( !(tissue == "Thyroid" & receptor %in% (degs %>% filter(tissue == "Thyroid") %>% pull(gene))) ) %>%
  group_by(tissue, state, gene) %>%
  summarise(counts = n()) %>%
  group_by(tissue, state) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(t = degs_count$counts) %>%
  mutate(p = n/t)



# plot of number of genes with at least one receptor correlated > 0.6
plot1 <- degs_recpt_cor_n %>%
  ggplot(mapping = aes(x = state, y = p, fill = state)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ tissue) +
    scale_fill_manual(values=c("#fdbb84", "#D7301F", "#3182bd"), labels=c("Common", "Female-specific", "Male-specific"), name = "DEGs group") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_text(colour="black", size=12),
      axis.text.y = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    scale_x_discrete(labels = c("Common", "Female-specific", "Male-specific")) +
    labs(y = "%", x = "DEGs group")
ggsave(filename="degs_recpt_cor_0_6.png", plot = plot1, path = "./plots/mechanisms", width=10, height=4)
unlink("degs_recpt_cor_0_6.png")




plot2 <- degs_recpt_cor %>%
  mutate_if(is.character, as.factor) %>%
  mutate(gene = fct_reorder(gene, as.numeric(state))) %>%
  ggplot(mapping = aes(x = gene, y = r, color = state)) +
    geom_point(size = 0.5) +
    facet_wrap(~ tissue) +
    scale_color_manual(values=c("#fdbb84", "#D7301F", "#3182bd"), labels=c("Common", "Female-specific", "Male-specific"), name = "DEGs group") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    labs(x = "DEGs", y = "Pearson's r")
ggsave(filename="degs_recpt_cor_scat.png", plot = plot2, path = "./plots/mechanisms", width=10, height=4)
unlink("degs_recpt_cor_scat.png")




plot3 <- degs_recpt_cor %>%
  mutate_if(is.character, as.factor) %>%
  mutate(gene = fct_reorder(gene, as.numeric(state))) %>%
  ggplot(mapping = aes(x = receptor, y = r, fill = state)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~ tissue) +
    scale_color_manual(values=c("#fdbb84", "#D7301F", "#3182bd"), labels=c("Common", "Female-specific", "Male-specific"), name = "DEGs group") +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=15),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(colour="black", size=13),
      legend.text=element_text(colour="black", size=13),
      legend.title=element_text(colour="black", size=15),
      strip.background = element_blank(),
      strip.text = element_text(colour="black", size=15)) +
    labs(x = "Hormone receptors", y = "Pearson's r")
ggsave(filename="degs_recpt_cor_boxpl.png", plot = plot3, path = "./plots/mechanisms", width=15, height=4)
unlink("degs_recpt_cor_boxpl.png")
