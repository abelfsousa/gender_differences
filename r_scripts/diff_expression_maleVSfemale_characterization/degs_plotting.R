# Understanding Gender Differential Susceptibility in Cancer


# DEGs plotting


library(tidyverse)
library(RColorBrewer)
library(viridis)


tcga_annotation <- read_tsv("./data/annotation/geneAnnot.gencode.v22.txt") %>%
  mutate(geneID = str_replace(geneID, "\\.[0-9]+", "")) %>%
  dplyr::rename(ensemblGeneID = geneID)


# -- Stomach


# load fpkm and metadata
stomach_tcga_gtex_fpkm <- read.delim("./files/stomach_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_meta <- read.delim("./files/tcga_stad_meta.txt", row.names = c(1), stringsAsFactors=F)

stomach_tcga_fpkm <- cbind(gene = rownames(stomach_tcga_gtex_fpkm), stomach_tcga_gtex_fpkm)
stomach_tcga_fpkm <- stomach_tcga_fpkm %>%
  as.tibble() %>%
  gather(key = "sample", value="fpkm", -gene) %>%
  filter(!str_detect(sample, "GTEX"))

stomach_tcga_meta <- cbind(sample = rownames(stomach_tcga_meta), stomach_tcga_meta)
stomach_tcga_meta <- stomach_tcga_meta %>%
  as.tibble() %>%
  mutate(sample = str_replace_all(sample, "-", "."))



stomach_tcga_fpkm <- stomach_tcga_fpkm %>%
  inner_join(stomach_tcga_meta[, c("sample", "sample_type_id", "gender")], by = "sample")



# load degs and enrichment table TCGA M vs F GTEx M vs F
stomach_degs <- read_tsv("./files/stomach_tumour_normal_signf_degs.txt")
stomach_enr <- read_tsv("./files/stomach_tumour_normal_degs_enr.txt")


# load degs and enrichment table TCGA M vs F TCGA normal M vs F
stomach_degs2 <- read_tsv("./files/stomach_tumour_normal_signf_degs2.txt")
stomach_enr2 <- read_tsv("./files/stomach_tumour_normal_degs_enr2.txt")



# boxplot of tumour specific DEGs between genders (enriched in POS sets)
stomach_tumour_spec_degs <- stomach_enr2 %>%
  filter(state == "tumour_specific", p.adjust < 0.05, Description == "POS") %>%
  dplyr::select(state, Description, ID, geneID) %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest() %>%
  inner_join(tcga_annotation[, c("geneName", "ensemblGeneID")], by=c("geneID" = "geneName")) %>%
  inner_join(stomach_tcga_fpkm, by=c("ensemblGeneID" = "gene")) %>%
  mutate(sample_type_id = ifelse(sample_type_id == 1, "T", "N"))


stomach_tumour_spec_degs_boxp1 <- stomach_tumour_spec_degs %>%
  filter(ID == "chrxp11") %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = gender)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_classic() +
    facet_grid(~ geneID, scales = "free_y") +
    theme(
      axis.title = element_text(colour="black", size=14),
      axis.text = element_text(colour="black", size=10),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=14),
      strip.background = element_blank(),
      strip.text.x = element_text(colour="black", size=8)) +
    scale_y_continuous(trans = "log2") +
    labs(x = "Tissue", y = "FPKM (log2)")
ggsave(filename="stomach_tumour_spec_degs1.png", plot=stomach_tumour_spec_degs_boxp1, path = "./plots/diff_expression_maleVSfemale/", width=6, height=4)
unlink("stomach_tumour_spec_degs1.png")



stomach_tumour_spec_degs_boxp2 <- stomach_tumour_spec_degs %>%
  filter(ID == "chrxp22") %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = gender)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_classic() +
    facet_grid(~ geneID, scales = "free_y") +
    theme(
      axis.title = element_text(colour="black", size=14),
      axis.text = element_text(colour="black", size=10),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=14),
      strip.background = element_blank(),
      strip.text.x = element_text(colour="black", size=8)) +
    scale_y_continuous(trans = "log2") +
    labs(x = "Tissue", y = "FPKM (log2)")
ggsave(filename="stomach_tumour_spec_degs2.png", plot=stomach_tumour_spec_degs_boxp2, path = "./plots/diff_expression_maleVSfemale/", width=8, height=4)
unlink("stomach_tumour_spec_degs2.png")




# boxplot of normal specific DEGs between genders
stomach_normal_spec_degs <- stomach_degs2 %>%
  filter(state == "normal_specific") %>%
  dplyr::select(geneName) %>%
  inner_join(tcga_annotation[, c("geneName", "ensemblGeneID")], by="geneName") %>%
  inner_join(stomach_tcga_fpkm, by=c("ensemblGeneID" = "gene")) %>%
  mutate(sample_type_id = ifelse(sample_type_id == 1, "T", "N"))


stomach_normal_spec_degs_boxp <- stomach_normal_spec_degs %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = gender)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_classic() +
    facet_grid(~ geneName, scales = "free_y") +
    theme(
      axis.title = element_text(colour="black", size=14),
      axis.text = element_text(colour="black", size=10),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=14),
      strip.background = element_blank(),
      strip.text.x = element_text(colour="black", size=10)) +
    scale_y_continuous(trans = "log2") +
    scale_fill_discrete(name = "Gender", labels = c("Female", "Male")) +
    scale_x_discrete(labels = c("Normal", "Tumour")) +
    labs(x = "Tissue", y = "FPKM (log2)")
ggsave(filename="stomach_normal_spec_degs.png", plot=stomach_normal_spec_degs_boxp, path = "./plots/diff_expression_maleVSfemale_tcga_normal/", width=5, height=3)
unlink("stomach_normal_spec_degs.png")








# -- Thyroid


# load fpkm and metadata
thyroid_tcga_gtex_fpkm <- read.delim("./files/thyroid_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_meta <- read.delim("./files/tcga_thca_meta.txt", row.names = c(1), stringsAsFactors=F)

thyroid_tcga_fpkm <- cbind(gene = rownames(thyroid_tcga_gtex_fpkm), thyroid_tcga_gtex_fpkm)
thyroid_tcga_fpkm <- thyroid_tcga_fpkm %>%
  as.tibble() %>%
  gather(key = "sample", value="fpkm", -gene) %>%
  filter(!str_detect(sample, "GTEX"))

thyroid_tcga_meta <- cbind(sample = rownames(thyroid_tcga_meta), thyroid_tcga_meta)
thyroid_tcga_meta <- thyroid_tcga_meta %>%
  as.tibble() %>%
  mutate(sample = str_replace_all(sample, "-", "."))



thyroid_tcga_fpkm <- thyroid_tcga_fpkm %>%
  inner_join(thyroid_tcga_meta[, c("sample", "sample_type_id", "gender")], by = "sample") %>%
  filter(sample_type_id != 6)



# load degs and enrichment table TCGA M vs F GTEx M vs F
thyroid_degs <- read_tsv("./files/thyroid_tumour_normal_signf_degs.txt")
thyroid_enr <- read_tsv("./files/thyroid_tumour_normal_degs_enr.txt")


# load degs and enrichment table TCGA M vs F TCGA normal M vs F
thyroid_degs2 <- read_tsv("./files/thyroid_tumour_normal_signf_degs2.txt")
thyroid_enr2 <- read_tsv("./files/thyroid_tumour_normal_degs_enr2.txt")




# boxplot of tumour specific DEGs between genders (enriched in POS sets)
thyroid_tumour_spec_degs <- thyroid_enr2 %>%
  filter(state == "tumour_specific", p.adjust < 0.05, Description == "POS") %>%
  dplyr::select(state, Description, ID, geneID) %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest() %>%
  inner_join(tcga_annotation[, c("geneName", "ensemblGeneID")], by=c("geneID" = "geneName")) %>%
  inner_join(thyroid_tcga_fpkm, by=c("ensemblGeneID" = "gene")) %>%
  mutate(sample_type_id = ifelse(sample_type_id == 1, "T", "N"))


thyroid_tumour_spec_degs_boxp <- thyroid_tumour_spec_degs %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = gender)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_classic() +
    facet_grid(~ geneID, scales = "free_y") +
    theme(
      axis.title = element_text(colour="black", size=14),
      axis.text = element_text(colour="black", size=10),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=14),
      strip.background = element_blank(),
      strip.text.x = element_text(colour="black", size=5)) +
    scale_y_continuous(trans = "log2") +
    labs(x = "Tissue", y = "FPKM (log2)")
ggsave(filename="thyroid_tumour_spec_degs.png", plot=thyroid_tumour_spec_degs_boxp, path = "./plots/diff_expression_maleVSfemale_tcga_normal/", width=6, height=4)
unlink("thyroid_tumour_spec_degs.png")



# boxplot of normal specific DEGs between genders
thyroid_normal_spec_degs <- thyroid_degs2 %>%
  filter(state == "normal_specific") %>%
  dplyr::select(geneName) %>%
  inner_join(tcga_annotation[, c("geneName", "ensemblGeneID")], by="geneName") %>%
  inner_join(thyroid_tcga_fpkm, by=c("ensemblGeneID" = "gene")) %>%
  mutate(sample_type_id = ifelse(sample_type_id == 1, "T", "N"))


thyroid_normal_spec_degs_boxp <- thyroid_normal_spec_degs %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = gender)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_classic() +
    facet_grid(~ geneName, scales = "free_y") +
    theme(
      axis.title = element_text(colour="black", size=14),
      axis.text = element_text(colour="black", size=10),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=14),
      strip.background = element_blank(),
      strip.text.x = element_text(colour="black", size=10)) +
    scale_y_continuous(trans = "log2") +
    scale_fill_discrete(name = "Gender", labels = c("Female", "Male")) +
    scale_x_discrete(labels = c("Normal", "Tumour")) +
    labs(x = "Tissue", y = "FPKM (log2)")
ggsave(filename="thyroid_normal_spec_degs.png", plot=thyroid_normal_spec_degs_boxp, path = "./plots/diff_expression_maleVSfemale_tcga_normal/", width=5, height=3)
unlink("thyroid_normal_spec_degs.png")
