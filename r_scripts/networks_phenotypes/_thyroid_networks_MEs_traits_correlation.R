# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# correlation of module eigengenes with sample traits


# -- Thyroid


library(tidyverse)
library(RColorBrewer)
library(viridis)




males_tumour_network_MEs <- read_tsv("./files/thyroid_males_tumour_network_MEs_traits.txt")
females_tumour_network_MEs <- read_tsv("./files/thyroid_females_tumour_network_MEs_traits.txt")


males_MEs_tumour_stage <- males_tumour_network_MEs %>%
  group_by(network, moduleL, moduleN, state) %>%
  do(broom::tidy(cor.test(.$ME, .$tumor_stage, method = "pearson"))) %>%
  ungroup()

males_MEs_hist_type <- males_tumour_network_MEs %>%
  group_by(network, moduleL, moduleN, state) %>%
  do(broom::tidy(cor.test(.$ME, .$histological_type, method = "pearson"))) %>%
  ungroup()


females_MEs_tumour_stage <- females_tumour_network_MEs %>%
  group_by(network, moduleL, moduleN, state) %>%
  do(broom::tidy(cor.test(.$ME, .$tumor_stage, method = "pearson"))) %>%
  ungroup()

females_MEs_hist_type <- females_tumour_network_MEs %>%
  group_by(network, moduleL, moduleN, state) %>%
  do(broom::tidy(cor.test(.$ME, .$histological_type, method = "pearson"))) %>%
  ungroup()
