# Understanding Gender Differential Susceptibility in Cancer


# Comparison of MvsF DEGs and TvsN DEGs


library(tidyverse)
library(ggpubr)
library(cowplot)



# Thyroid

thyroid_TvsN_degs <- read_tsv("./files/thyroid_males_females_signf_degs.txt") %>% mutate(tissue = "Thyroid")
thyroid_MvsF_degs <- read_tsv("./files/thyroid_tumour_normal_signf_degs.txt") %>% mutate(tissue = "Thyroid")

common_thyroid <- inner_join(
  thyroid_MvsF_degs %>% dplyr::select(geneName, stateMvsF = state),
  thyroid_TvsN_degs %>% dplyr::select(geneName, stateTvsN = state),
  by = "geneName") %>%
  group_by(stateMvsF, stateTvsN) %>%
  count()



# Stomach

stomach_TvsN_degs <- read_tsv("./files/stomach_males_females_signf_degs.txt") %>% mutate(tissue = "Stomach")
stomach_MvsF_degs <- read_tsv("./files/stomach_tumour_normal_signf_degs.txt") %>% mutate(tissue = "Stomach")

common_stomach <- inner_join(
  stomach_MvsF_degs %>% dplyr::select(geneName, stateMvsF = state),
  stomach_TvsN_degs %>% dplyr::select(geneName, stateTvsN = state),
  by = "geneName") %>%
  group_by(stateMvsF, stateTvsN) %>%
  count()



# expression data

# TCGA gene annotation
tcga_annotation <- read_tsv("./data/annotation/geneAnnot.gencode.v22.txt") %>%
  mutate(geneID = str_replace(geneID, "\\.[0-9]+", "")) %>%
  dplyr::rename(ensemblGeneID = geneID)


# thyroid

# load fpkm and metadata
thyroid_tcga_fpkm <- read_tsv("./files/thyroid_tcga_gtex_fpkm.txt") %>%
  gather(key = "sample", value="fpkm", -gene) %>%
  filter(!str_detect(sample, "GTEX")) %>%
  mutate(sample = str_replace_all(sample, "-", "."))

thyroid_tcga_meta <- read_tsv("./files/tcga_thca_meta.txt") %>%
  dplyr::rename(sample = submitter_id) %>%
  mutate(sample = str_replace_all(sample, "-", "."))

thyroid_tcga_fpkm <- thyroid_tcga_fpkm %>%
  inner_join(thyroid_tcga_meta[, c("sample", "sample_type_id", "gender")], by = "sample") %>%
  inner_join(tcga_annotation, by = c("gene" = "ensemblGeneID")) %>%
  filter(sample_type_id != "06")


# stomach

# load fpkm and metadata
stomach_tcga_fpkm <- read_tsv("./files/stomach_tcga_gtex_fpkm.txt") %>%
  gather(key = "sample", value="fpkm", -gene) %>%
  filter(!str_detect(sample, "GTEX")) %>%
  mutate(sample = str_replace_all(sample, "-", "."))

stomach_tcga_meta <- read_tsv("./files/tcga_stad_meta.txt") %>%
  dplyr::rename(sample = submitter_id) %>%
  mutate(sample = str_replace_all(sample, "-", "."))

stomach_tcga_fpkm <- stomach_tcga_fpkm %>%
  inner_join(stomach_tcga_meta[, c("sample", "sample_type_id", "gender")], by = "sample") %>%
  inner_join(tcga_annotation, by = c("gene" = "ensemblGeneID"))




# compare expression of gender-specific DEGs between tumour and normal tissues in females and males



thca_female_specific_TvsN <- thyroid_TvsN_degs %>%
  filter(state == "female_specific") %>%
  dplyr::select(geneName, state) %>%
  inner_join(thyroid_tcga_fpkm %>% dplyr::select(geneName, sample, fpkm, gender, sample_type_id), by = c("geneName")) %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = sample_type_id)) +
  geom_boxplot() +
  facet_wrap( ~ gender,
  	labeller = labeller(gender = c("female" = "Female", "male" = "Male"))) +
  stat_compare_means() +
  theme_classic() +
  theme(
  	axis.title = element_text(color = "black", size=15),
  	axis.text = element_text(color="black", size=13),
  	strip.background = element_blank(),
  	strip.text = element_text(color="black", size=15)) +
  scale_fill_manual(values=c("#d95f02", "#7fc97f"), guide=F) +
  scale_x_discrete(labels=c("Tumour", "Normal")) +
  scale_y_continuous(trans = "log2") +
  labs(x = "Tissue", y = "FPKM (log2 scale)")




thca_male_specific_TvsN <- thyroid_TvsN_degs %>%
  filter(state == "male_specific") %>%
  dplyr::select(geneName, state) %>%
  inner_join(thyroid_tcga_fpkm %>% dplyr::select(geneName, sample, fpkm, gender, sample_type_id), by = c("geneName")) %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = sample_type_id)) +
  geom_boxplot() +
  facet_wrap( ~ gender,
    labeller = labeller(gender = c("female" = "Female", "male" = "Male"))) +
  stat_compare_means() +
  theme_classic() +
  theme(
    axis.title = element_text(color = "black", size=15),
    axis.text = element_text(color="black", size=13),
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=15)) +
  scale_fill_manual(values=c("#d95f02", "#7fc97f"), guide=F) +
  scale_x_discrete(labels=c("Tumour", "Normal")) +
  scale_y_continuous(trans = "log2") +
  labs(x = "Tissue", y = "FPKM (log2 scale)")




#stomach

stad_female_specific_TvsN <- stomach_TvsN_degs %>%
  filter(state == "female_specific") %>%
  dplyr::select(geneName, state) %>%
  inner_join(stomach_tcga_fpkm %>% dplyr::select(geneName, sample, fpkm, gender, sample_type_id), by = c("geneName")) %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = sample_type_id)) +
  geom_boxplot() +
  facet_wrap( ~ gender,
    labeller = labeller(gender = c("female" = "Female", "male" = "Male"))) +
  stat_compare_means() +
  theme_classic() +
  theme(
    axis.title = element_text(color = "black", size=15),
    axis.text = element_text(color="black", size=13),
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=15)) +
  scale_fill_manual(values=c("#d95f02", "#7fc97f"), guide=F) +
  scale_x_discrete(labels=c("Tumour", "Normal")) +
  scale_y_continuous(trans = "log2") +
  labs(x = "Tissue", y = "FPKM (log2 scale)")



stad_male_specific_TvsN <- stomach_TvsN_degs %>%
  filter(state == "male_specific") %>%
  dplyr::select(geneName, state) %>%
  inner_join(stomach_tcga_fpkm %>% dplyr::select(geneName, sample, fpkm, gender, sample_type_id), by = c("geneName")) %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = sample_type_id)) +
  geom_boxplot() +
  facet_wrap( ~ gender,
    labeller = labeller(gender = c("female" = "Female", "male" = "Male"))) +
  stat_compare_means() +
  theme_classic() +
  theme(
    axis.title = element_text(color = "black", size=15),
    axis.text = element_text(color="black", size=13),
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=15)) +
  scale_fill_manual(values=c("#d95f02", "#7fc97f"), guide=F) +
  scale_x_discrete(labels=c("Tumour", "Normal")) +
  scale_y_continuous(trans = "log2") +
  labs(x = "Tissue", y = "FPKM (log2 scale)")




gender_specific_degs_both_genders <- plot_grid(
  thca_female_specific_TvsN,
  thca_male_specific_TvsN,
  stad_female_specific_TvsN,
  stad_male_specific_TvsN,
  labels = c("THCA\nFemale-specific", "THCA\nMale-specific", "STAD\nFemale-specific", "STAD\nMale-specific"),
  nrow = 2,
  align = "v",
  label_size = 10,
  vjust = 1)
ggsave(filename="gender_specific_degs_both_genders.png", plot=gender_specific_degs_both_genders, height=10, width=10, path = "./plots/diff_expression_tumourVSnormal/")
unlink("gender_specific_degs_both_genders.png")
