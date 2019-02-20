# Understanding Gender Differential Susceptibility in Cancer


# DEGs plotting


library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)


# load TCGA annotation
tcga_annotation <- read_tsv("./data/annotation/geneAnnot.gencode.v22.txt") %>%
  mutate(geneID = str_replace(geneID, "\\.[0-9]+", "")) %>%
  dplyr::rename(ensemblGeneID = geneID)



# load fpkm and metadata

# Thyroid
thyroid_tcga_gtex_fpkm <- read.delim("./files/thyroid_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_gtex_meta <- read.delim("./files/thyroid_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)

thyroid_tcga_gtex_fpkm <- cbind(gene = rownames(thyroid_tcga_gtex_fpkm), thyroid_tcga_gtex_fpkm)
thyroid_tcga_gtex_fpkm <- thyroid_tcga_gtex_fpkm %>%
  as.tibble() %>%
  gather(key = "sample", value="fpkm", -gene)

thyroid_tcga_gtex_meta <- cbind(sample = rownames(thyroid_tcga_gtex_meta), thyroid_tcga_gtex_meta)
thyroid_tcga_gtex_meta <- thyroid_tcga_gtex_meta %>%
  as.tibble() %>%
  mutate(sample = str_replace_all(sample, "-", "."))

thyroid_tcga_gtex_fpkm <- thyroid_tcga_gtex_fpkm %>%
  inner_join(tcga_annotation[, c("ensemblGeneID", "geneName")], by = c("gene" = "ensemblGeneID")) %>%
  inner_join(thyroid_tcga_gtex_meta[, c("sample", "sample_type", "gender")], by = "sample") %>%
  filter(!(sample_type == "TCGA_metastatic" | sample_type == "GTEx")) %>%
  mutate(tissue = "Thyroid") %>%
  dplyr::select(tissue, gene = geneName, sample_type, sample, gender, fpkm)



# Stomach
stomach_tcga_gtex_fpkm <- read.delim("./files/stomach_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_gtex_meta <- read.delim("./files/stomach_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)

stomach_tcga_gtex_fpkm <- cbind(gene = rownames(stomach_tcga_gtex_fpkm), stomach_tcga_gtex_fpkm)
stomach_tcga_gtex_fpkm <- stomach_tcga_gtex_fpkm %>%
  as.tibble() %>%
  gather(key = "sample", value="fpkm", -gene)

stomach_tcga_gtex_meta <- cbind(sample = rownames(stomach_tcga_gtex_meta), stomach_tcga_gtex_meta)
stomach_tcga_gtex_meta <- stomach_tcga_gtex_meta %>%
  as.tibble() %>%
  mutate(sample = str_replace_all(sample, "-", "."))

stomach_tcga_gtex_fpkm <- stomach_tcga_gtex_fpkm %>%
  inner_join(tcga_annotation[, c("ensemblGeneID", "geneName")], by = c("gene" = "ensemblGeneID")) %>%
  inner_join(stomach_tcga_gtex_meta[, c("sample", "sample_type", "gender")], by = "sample") %>%
  filter(!(sample_type == "GTEx")) %>%
  mutate(tissue = "Stomach") %>%
  dplyr::select(tissue, gene = geneName, sample_type, sample, gender, fpkm)



# merge fpkm from thyroid and stomach
fpkm <- bind_rows(stomach_tcga_gtex_fpkm, thyroid_tcga_gtex_fpkm) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(fpkm = log2(fpkm+1))



# PLOTTING OF HORMONE RECEPTORS DEGs



# load sex steroid hormone receptor genes
horm_recpt <- read.gmt("./data/gene_lists/c5.mf.v6.2.symbols.gmt") %>%
  filter(ont == "GO_STEROID_HORMONE_RECEPTOR_ACTIVITY") %>%
  pull(gene)


# load degs tables
stomach_degs <- read_tsv("./files/stomach_males_females_signf_degs.txt") %>%
  dplyr::select(gene = geneName, state) %>%
  mutate(tissue = "Stomach")


thyroid_degs <- read_tsv("./files/thyroid_males_females_signf_degs.txt") %>%
  dplyr::select(gene = geneName, state) %>%
  mutate(tissue = "Thyroid")


degs <- bind_rows(stomach_degs, thyroid_degs) %>%
  mutate_if(is.character, as.factor) %>%
  filter(gene %in% horm_recpt)


horm_recpt_plot1 <- fpkm %>%
  inner_join(degs, by=c("tissue", "gene")) %>%
  mutate_if(is.character, as.factor) %>%
  filter(!(state == "female_specific" & gender == "male")) %>%
  filter(!(state == "male_specific" & gender == "female")) %>%
  filter(tissue == "Thyroid") %>%
  ggplot(mapping = aes(x = gene, y = fpkm, fill = state, color = sample_type)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values=c("#fdbb84", "#D7301F", "#3182bd"), labels=c("Common", "Female-specific", "Male-specific"), name = "DEGs group") +
    scale_color_manual(values=c("#66c2a5", "#8da0cb"), name="Tissue", labels = c("Normal", "Tumour")) +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=14),
      axis.text = element_text(colour="black", size=10),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=14),
      strip.background = element_blank(),
      strip.text.x = element_text(colour="black", size=8)) +
    labs(x = "Hormone receptor DEG", y = "FPKM (log2)", title = "Thyroid")
ggsave(filename="degs_TvsN_hormone_recpt_thca.png", plot=horm_recpt_plot1, path = "./plots/diff_expression_tumourVSnormal/", width=8, height=4)
unlink("degs_TvsN_hormone_recpt_thca.png")





horm_recpt_plot2 <- fpkm %>%
  inner_join(degs, by=c("tissue", "gene")) %>%
  mutate_if(is.character, as.factor) %>%
  filter(!(state == "female_specific" & gender == "male")) %>%
  filter(!(state == "male_specific" & gender == "female")) %>%
  filter(tissue == "Stomach") %>%
  ggplot(mapping = aes(x = gene, y = fpkm, fill = state, color = sample_type)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values=c("#fdbb84", "#D7301F", "#3182bd"), labels=c("Common", "Female-specific", "Male-specific"), name = "DEGs group") +
    scale_color_manual(values=c("#66c2a5", "#8da0cb"), name="Tissue", labels = c("Normal", "Tumour")) +
    theme_classic() +
    theme(
      axis.title = element_text(colour="black", size=14),
      axis.text.x = element_text(colour="black", size=10, angle = 90),
      axis.text.y = element_text(colour="black", size=10),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=14),
      strip.background = element_blank(),
      strip.text.x = element_text(colour="black", size=8)) +
    labs(x = "Hormone receptor DEG", y = "FPKM (log2)", title = "Stomach")
ggsave(filename="degs_TvsN_hormone_recpt_stad.png", plot=horm_recpt_plot2, path = "./plots/diff_expression_tumourVSnormal/", width=8, height=4)
unlink("degs_TvsN_hormone_recpt_stad.png")
