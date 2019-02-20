# Understanding Gender Differential Susceptibility in Cancer


# PCA analysis



library(tidyverse)


# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))



# -- Thyroid


# load datasets
thyroid_tcga_gtex_counts <- read.delim("./files/thyroid_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_gtex_meta <- read.delim("./files/thyroid_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_gtex_meta <- cbind(sample=rownames(thyroid_tcga_gtex_meta), thyroid_tcga_gtex_meta)
thyroid_tcga_gtex_meta$sample <- gsub("-", ".", thyroid_tcga_gtex_meta$sample)


# reorder datasets
thyroid_tcga_gtex_counts <- thyroid_tcga_gtex_counts[, order(colnames(thyroid_tcga_gtex_counts))]
thyroid_tcga_gtex_meta <- thyroid_tcga_gtex_meta[order(thyroid_tcga_gtex_meta$sample), ]


# remove metastatic samples
thyroid_tcga_gtex_counts <- thyroid_tcga_gtex_counts[, !(colnames(thyroid_tcga_gtex_counts) %in% thyroid_tcga_gtex_meta[(thyroid_tcga_gtex_meta$sample_type == "TCGA_metastatic"), "sample"])]
thyroid_tcga_gtex_meta <- thyroid_tcga_gtex_meta[!(thyroid_tcga_gtex_meta$sample_type == "TCGA_metastatic"), ]


# add experimental group
sample_type_code <- c("GTEx" = "GTEx", "TCGA_normal" = "TCGA", "TCGA_tumour" = "TCGA")
thyroid_tcga_gtex_meta <- cbind(thyroid_tcga_gtex_meta, data.frame(data = sample_type_code[thyroid_tcga_gtex_meta$sample_type]))

# add tissue type
tissue_type_code <- c("GTEx" = "normal", "TCGA_normal" = "normal", "TCGA_tumour" = "tumour")
thyroid_tcga_gtex_meta <- cbind(thyroid_tcga_gtex_meta, data.frame(tissue_type = tissue_type_code[thyroid_tcga_gtex_meta$sample_type]))


thyroid_tcga_gtex_counts %>% colnames %>% all.equal(thyroid_tcga_gtex_meta$sample)
#TRUE




# -- Stomach


# load datasets
stomach_tcga_gtex_counts <- read.delim("./files/stomach_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_gtex_meta <- read.delim("./files/stomach_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_gtex_meta <- cbind(sample=rownames(stomach_tcga_gtex_meta), stomach_tcga_gtex_meta)
stomach_tcga_gtex_meta$sample <- gsub("-", ".", stomach_tcga_gtex_meta$sample)


# reorder datasets
stomach_tcga_gtex_counts <- stomach_tcga_gtex_counts[, order(colnames(stomach_tcga_gtex_counts))]
stomach_tcga_gtex_meta <- stomach_tcga_gtex_meta[order(stomach_tcga_gtex_meta$sample), ]


# add experimental group
sample_type_code <- c("GTEx" = "GTEx", "TCGA_normal" = "TCGA", "TCGA_tumour" = "TCGA")
stomach_tcga_gtex_meta <- cbind(stomach_tcga_gtex_meta, data.frame(data = sample_type_code[stomach_tcga_gtex_meta$sample_type]))

# add tissue type
tissue_type_code <- c("GTEx" = "normal", "TCGA_normal" = "normal", "TCGA_tumour" = "tumour")
stomach_tcga_gtex_meta <- cbind(stomach_tcga_gtex_meta, data.frame(tissue_type = tissue_type_code[stomach_tcga_gtex_meta$sample_type]))


stomach_tcga_gtex_counts %>% colnames %>% all.equal(stomach_tcga_gtex_meta$sample)
#TRUE
















# PCA thyroid and stomach
pca_thyroid <- prcomp(t(thyroid_tcga_gtex_counts), center = TRUE, scale. = TRUE)$x %>% as.data.frame()
pca_thyroid <- cbind(sample = rownames(pca_thyroid), pca_thyroid)
pca_thyroid <- pca_thyroid %>%
  as.tibble() %>%
  gather(key = "PC", value = "value", -sample) %>%
  inner_join(thyroid_tcga_gtex_meta[, c("sample", "gender", "sample_type")], by = "sample") %>%
  mutate(tissue = "Thyroid")


pca_stomach <- prcomp(t(stomach_tcga_gtex_counts), center = TRUE, scale. = TRUE)$x %>% as.data.frame()
pca_stomach <- cbind(sample = rownames(pca_stomach), pca_stomach)
pca_stomach <- pca_stomach %>%
  as.tibble() %>%
  gather(key = "PC", value = "value", -sample) %>%
  inner_join(stomach_tcga_gtex_meta[, c("sample", "gender", "sample_type")], by = "sample") %>%
  mutate(tissue = "Stomach")




# thyroid barplot
pca_thyroid_bp <- pca_thyroid %>%
  filter(PC == "PC1" | PC == "PC2") %>%
  spread(key = "PC", value = "value") %>%
  ggplot(mapping=aes(x=PC1, y=PC2, color=gender, shape=sample_type) ) +
  geom_point() +
  scale_color_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  scale_shape_manual(values=c(16, 3, 15), name="Sample", labels = c("GTEx", "TCGA normal", "TCGA tumour")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour="black", size=15),
    axis.text = element_text(colour="black", size=13),
    legend.text=element_text(colour="black", size=13),
    legend.title=element_text(colour="black", size=15),
    plot.title = element_text(colour="black", size=15, hjust = 0.5)) +
  coord_fixed() +
  labs(x = "PC1", y = "PC2", title = "Thyroid")
ggsave(filename="pca_thyroid_samples.pdf", plot=pca_thyroid_bp, path = "./plots/pre_processing/", width = 5, height = 5)
unlink("pca_thyroid_samples.pdf")



# stomach barplot
pca_stomach_bp <- pca_stomach %>%
  filter(PC == "PC1" | PC == "PC2") %>%
  spread(key = "PC", value = "value") %>%
  ggplot(mapping=aes(x=PC1, y=PC2, color=gender, shape=sample_type) ) +
  geom_point() +
  scale_color_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  scale_shape_manual(values=c(16, 3, 15), name="Sample", labels = c("GTEx", "TCGA normal", "TCGA tumour")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour="black", size=15),
    axis.text = element_text(colour="black", size=13),
    legend.text=element_text(colour="black", size=13),
    legend.title=element_text(colour="black", size=15),
    plot.title = element_text(colour="black", size=15, hjust = 0.5)) +
  coord_fixed() +
  labs(x = "PC1", y = "PC2", title = "Stomach")
ggsave(filename="pca_stomach_samples.pdf", plot=pca_stomach_bp, path = "./plots/pre_processing/", width = 5, height = 5)
unlink("pca_stomach_samples.pdf")




# thyroid and stomach barplot
pca_bp <- rbind(pca_thyroid, pca_stomach) %>%
  filter(PC == "PC1" | PC == "PC2") %>%
  spread(key = "PC", value = "value") %>%
  ggplot(mapping=aes(x=PC1, y=PC2, color=gender, shape=sample_type) ) +
  geom_point() +
  facet_wrap(~ tissue, scales = "free") +
  scale_color_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  scale_shape_manual(values=c(16, 3, 15), name="Sample", labels = c("GTEx", "TCGA normal", "TCGA tumour")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour="black", size=15),
    axis.text = element_text(colour="black", size=13),
    legend.text=element_text(colour="black", size=13),
    legend.title=element_text(colour="black", size=15),
    strip.background = element_blank(),
    strip.text = element_text(colour="black", size=15)) +
  labs(x = "PC1", y = "PC2")
ggsave(filename="pca_thyroid_stomach_samples.pdf", plot=pca_bp, path = "./plots/pre_processing/", width = 8, height = 4)
unlink("pca_thyroid_stomach_samples.pdf")






# PCA stomach gtex
stomach_gtex_meta <- read.delim("./files/gtex_stomach_meta.txt", row.names = c(1), stringsAsFactors=F)
stomach_gtex_meta <- cbind(sample=rownames(stomach_gtex_meta), stomach_gtex_meta)

pca_stomach_gtex <- prcomp(t(stomach_tcga_gtex_counts[, colnames(stomach_tcga_gtex_counts) %in% stomach_tcga_gtex_meta[stomach_tcga_gtex_meta$data == "GTEx", "sample"]]), center = TRUE, scale. = TRUE)$x %>% as.data.frame()
pca_stomach_gtex <- cbind(sample = rownames(pca_stomach_gtex), pca_stomach_gtex)
pca_stomach_gtex <- pca_stomach_gtex %>%
  as.tibble() %>%
  gather(key = "PC", value = "value", -sample) %>%
  inner_join(stomach_gtex_meta[, c("sample", "GENDER", "MHCANCERNM", "SMCENTER", "SMTSTPTREF", "SMNABTCHT")], by = "sample") %>%
  mutate(GENDER = if_else(GENDER == 1, "male", "female")) %>%
  mutate(MHCANCERNM = as.character(MHCANCERNM))



pca_stomach_gtex_bp1 <- pca_stomach_gtex %>%
  filter(PC == "PC1" | PC == "PC2" | PC == "PC3") %>%
  spread(key = "PC", value = "value") %>%
  ggplot(mapping=aes(x=PC1, y=PC2, color=GENDER, shape = SMNABTCHT) ) +
  geom_point() +
  scale_color_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  scale_shape_discrete(name="Batch", labels = c("1", "2")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour="black", size=15),
    axis.text = element_text(colour="black", size=13),
    legend.text=element_text(colour="black", size=13),
    legend.title=element_text(colour="black", size=15),
    plot.title = element_text(colour="black", size=15, hjust = 0.5)) +
  coord_fixed() +
  labs(x = "PC1", y = "PC2", title = "Stomach GTEx")
ggsave(filename="pca_stomach_gtex1.pdf", plot=pca_stomach_gtex_bp1, path = "./plots/pre_processing/", width = 4, height = 4)
unlink("pca_stomach_gtex1.pdf")


pca_stomach_gtex_bp2 <- pca_stomach_gtex %>%
  filter(PC == "PC1" | PC == "PC2" | PC == "PC3") %>%
  spread(key = "PC", value = "value") %>%
  ggplot(mapping=aes(x=PC2, y=PC3, color=GENDER, shape = SMNABTCHT) ) +
  geom_point() +
  scale_color_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  scale_shape_discrete(name="Batch", labels = c("1", "2")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour="black", size=15),
    axis.text = element_text(colour="black", size=13),
    legend.text=element_text(colour="black", size=13),
    legend.title=element_text(colour="black", size=15),
    plot.title = element_text(colour="black", size=15, hjust = 0.5)) +
  coord_fixed() +
  labs(x = "PC2", y = "PC3", title = "Stomach GTEx")
ggsave(filename="pca_stomach_gtex2.pdf", plot=pca_stomach_gtex_bp2, path = "./plots/pre_processing/", width = 4, height = 4)
unlink("pca_stomach_gtex2.pdf")




# PCA stomach and thyroid TCGA samples
pca_thyroid2 <- prcomp(t(thyroid_tcga_gtex_counts[, !str_detect(colnames(thyroid_tcga_gtex_counts), "GTEX")]), center = TRUE, scale. = TRUE)$x %>% as.data.frame()
pca_thyroid2 <- cbind(sample = rownames(pca_thyroid2), pca_thyroid2)
pca_thyroid2 <- pca_thyroid2 %>%
  as.tibble() %>%
  gather(key = "PC", value = "value", -sample) %>%
  inner_join(thyroid_tcga_gtex_meta[, c("sample", "gender", "sample_type")], by = "sample") %>%
  mutate(tissue = "Thyroid")


pca_stomach2 <- prcomp(t(stomach_tcga_gtex_counts[, !str_detect(colnames(stomach_tcga_gtex_counts), "GTEX")]), center = TRUE, scale. = TRUE)$x %>% as.data.frame()
pca_stomach2 <- cbind(sample = rownames(pca_stomach2), pca_stomach2)
pca_stomach2 <- pca_stomach2 %>%
  as.tibble() %>%
  gather(key = "PC", value = "value", -sample) %>%
  inner_join(stomach_tcga_gtex_meta[, c("sample", "gender", "sample_type")], by = "sample") %>%
  mutate(tissue = "Stomach")


# thyroid and stomach barplot
pca2_bp <- rbind(pca_thyroid2, pca_stomach2) %>%
  filter(PC == "PC1" | PC == "PC2") %>%
  spread(key = "PC", value = "value") %>%
  ggplot(mapping=aes(x=PC1, y=PC2, color=gender, shape=sample_type) ) +
  geom_point() +
  facet_wrap(~ tissue, scales = "free") +
  scale_color_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  scale_shape_manual(values=c(3, 15), name="Sample", labels = c("TCGA normal", "TCGA tumour")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour="black", size=15),
    axis.text = element_text(colour="black", size=13),
    legend.text=element_text(colour="black", size=13),
    legend.title=element_text(colour="black", size=15),
    strip.background = element_blank(),
    strip.text = element_text(colour="black", size=15)) +
  labs(x = "PC1", y = "PC2")
ggsave(filename="pca_thyroid_stomach_tcga_samples.pdf", plot=pca2_bp, path = "./plots/pre_processing/", width = 8, height = 4)
unlink("pca_thyroid_stomach_tcga_samples.pdf")





save(list=ls(), file="r_workspaces/sample_pca.RData")
