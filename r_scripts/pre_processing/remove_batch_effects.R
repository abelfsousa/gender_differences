# Understanding Gender Differential Susceptibility in Cancer


# Remove batch effects from gene expression



library(tidyverse)
library(ComplexHeatmap)

source(file = "./r_scripts/utils.R")

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy



# -- Stomach


# load datasets
stomach_tcga_gtex_counts <- read.delim("./files/stomach_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_gtex_meta <- read.delim("./files/stomach_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_gtex_meta <- cbind(sample=rownames(stomach_tcga_gtex_meta), stomach_tcga_gtex_meta)
stomach_tcga_gtex_meta$sample <- gsub("-", ".", stomach_tcga_gtex_meta$sample)

stomach_gtex_meta <- read.delim("./files/gtex_stomach_meta.txt", row.names = c(1), stringsAsFactors=F)
stomach_gtex_meta <- cbind(sample=rownames(stomach_gtex_meta), stomach_gtex_meta)


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


stomach_tcga_gtex_counts <- log2(stomach_tcga_gtex_counts+1)





# PCA stomach gtex
pca_stomach_gtex <- prcomp(t(stomach_tcga_gtex_counts[, colnames(stomach_tcga_gtex_counts) %in% stomach_tcga_gtex_meta[stomach_tcga_gtex_meta$data == "GTEx", "sample"]]), center = TRUE, scale. = TRUE)$x %>% as.data.frame()
pca_stomach_gtex <- cbind(sample = rownames(pca_stomach_gtex), pca_stomach_gtex)
pca_stomach_gtex <- pca_stomach_gtex %>%
  as.tibble() %>%
  gather(key = "PC", value = "value", -sample) %>%
  inner_join(stomach_gtex_meta[, c("sample", "SMRIN", "AGE", "ETHNCTY", "MHCANCERNM", "SMCENTER", "SMTSTPTREF", "SMNABTCHT", "SMTSISCH", "GENDER")], by = "sample") %>%
  mutate(GENDER = if_else(GENDER == 1, "male", "female")) %>%
  mutate(MHCANCERNM = as.character(MHCANCERNM))




# correlation of gtex stomach principal components to covariates
# select 10 first principal components
stomach_gtex_covars <- pca_stomach_gtex %>%
  filter(PC %in% paste0("PC", c(1:10))) %>%
  pivot_wider(names_from="PC", values_from="value")

stomach_gtex_pcs <- stomach_gtex_covars %>%
  dplyr::select(sample, starts_with("PC"))

stomach_gtex_covars <- stomach_gtex_covars %>%
    dplyr::select(-starts_with("PC")) %>%
    mutate(ETHNCTY = as.character(ETHNCTY)) %>%
    mutate(SMTSTPTREF = if_else(SMTSTPTREF == "", "unknown", SMTSTPTREF)) %>%
    mutate(SMTSISCH = if_else(is.na(SMTSISCH), round(mean(SMTSISCH, na.rm = T)), as.numeric(SMTSISCH))) %>%
    mutate(SMNABTCHT = factor(SMNABTCHT, labels = c("1", "2"))) %>%
    mutate(SMTSTPTREF = factor(SMTSTPTREF, labels = c("1", "2", "3"))) %>%
    mutate_if(is.character, as.factor)

#stomach_gtex_covars <- model.matrix( as.formula(paste0("~", paste(c(colnames(stomach_gtex_covars[, -c(1)]), "0"), collapse="+"))), stomach_gtex_covars )
stomach_gtex_covars <- model.matrix(
  object = as.formula(~ SMRIN + AGE + ETHNCTY + MHCANCERNM + SMCENTER + SMTSTPTREF + SMNABTCHT + SMTSISCH + GENDER + 0),
  data = stomach_gtex_covars,
  contrasts.arg=list(
    ETHNCTY=contrasts(stomach_gtex_covars$ETHNCTY, contrasts=F),
    MHCANCERNM=contrasts(stomach_gtex_covars$MHCANCERNM, contrasts=F),
    SMCENTER=contrasts(stomach_gtex_covars$SMCENTER, contrasts=F),
    SMTSTPTREF=contrasts(stomach_gtex_covars$SMTSTPTREF, contrasts=F),
    SMNABTCHT=contrasts(stomach_gtex_covars$SMNABTCHT, contrasts=F),
    GENDER=contrasts(stomach_gtex_covars$GENDER, contrasts=F)))

stomach_gtex_cors <- cor(stomach_gtex_pcs[, -c(1)], stomach_gtex_covars)

stomach_gtex_cors_heatmap <- Heatmap(
  matrix = stomach_gtex_cors,
  name = "Pearson's r",
  border = F,
  cluster_rows = T,
  cluster_columns = T,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 2))

pdf(file = "./plots/pre_processing/stomach_gtex_cors_pcs_covars_heatmap.pdf", width = 8, height = 5)
stomach_gtex_cors_heatmap
dev.off()

png(file = "./plots/pre_processing/stomach_gtex_cors_pcs_covars_heatmap.png", width = 8, height = 5, units = "in", res = 300)
stomach_gtex_cors_heatmap
dev.off()



# regress-out covars from gene expression

stomach_meta <- stomach_gtex_meta %>%
  as_tibble() %>%
  #select(c("sample", "SMRIN", "AGE", "ETHNCTY", "MHCANCERNM", "SMCENTER", "SMTSTPTREF", "SMNABTCHT", "SMTSISCH", "GENDER")) %>%
  select(c("sample", "SMRIN", "AGE", "ETHNCTY", "MHCANCERNM", "SMCENTER", "SMTSTPTREF", "SMNABTCHT", "SMTSISCH")) %>%
  mutate(SMTSTPTREF = if_else(SMTSTPTREF == "", "unknown", SMTSTPTREF)) %>%
  mutate(SMTSISCH = if_else(is.na(SMTSISCH), round(mean(SMTSISCH, na.rm = T)), as.numeric(SMTSISCH))) %>%
  mutate(ETHNCTY = as.character(ETHNCTY)) %>%
  mutate(MHCANCERNM = as.character(MHCANCERNM))

stomach_expr <- stomach_tcga_gtex_counts[, colnames(stomach_tcga_gtex_counts) %in% stomach_tcga_gtex_meta[stomach_tcga_gtex_meta$data == "GTEx", "sample"]]

stomach_expr <- stomach_expr %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "log2fpkm") %>%
  group_by(gene) %>%
  nest() %>%
  mutate(data2 = map(.x = data, .f = remove_batch2, covs = stomach_meta)) %>%
  select(-data) %>%
  unnest()

stomach_expr_mat <- stomach_expr %>%
  pivot_wider(names_from="sample", values_from="residual") %>%
  column_to_rownames(var = "gene")


# PCA stomach gtex
pca_stomach_gtex <- prcomp(t(stomach_expr_mat), center = TRUE, scale. = TRUE)$x %>%
  as.data.frame()
pca_stomach_gtex <- cbind(sample = rownames(pca_stomach_gtex), pca_stomach_gtex)
pca_stomach_gtex <- pca_stomach_gtex %>%
  as.tibble() %>%
  gather(key = "PC", value = "value", -sample)


stomach_gtex_pcs <- pca_stomach_gtex %>%
  filter(PC %in% paste0("PC", c(1:10))) %>%
  pivot_wider(names_from="PC", values_from="value")


stomach_gtex_cors <- cor(stomach_gtex_pcs[, -c(1)], stomach_gtex_covars)


stomach_gtex_cors_heatmap <- Heatmap(
  matrix = stomach_gtex_cors,
  name = "Pearson's r",
  border = F,
  cluster_rows = T,
  cluster_columns = T,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 2))

pdf(file = "./plots/pre_processing/stomach_gtex_cors_pcs_covars_heatmapRegOut.pdf", width = 8, height = 5)
stomach_gtex_cors_heatmap
dev.off()

png(file = "./plots/pre_processing/stomach_gtex_cors_pcs_covars_heatmapRegOut.png", width = 8, height = 5, units = "in", res = 300)
stomach_gtex_cors_heatmap
dev.off()


stomach_expr_mat <- stomach_expr %>%
  pivot_wider(names_from="sample", values_from="residual")


write_tsv(stomach_expr_mat, "./files/stomach_tcga_gtex_fpkm_GTEX_regressed_out_COVARS.txt")














# -- Thyroid


# load datasets
thyroid_tcga_gtex_counts <- read.delim("./files/thyroid_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_gtex_meta <- read.delim("./files/thyroid_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_gtex_meta <- cbind(sample=rownames(thyroid_tcga_gtex_meta), thyroid_tcga_gtex_meta)
thyroid_tcga_gtex_meta$sample <- gsub("-", ".", thyroid_tcga_gtex_meta$sample)

thyroid_gtex_meta <- read.delim("./files/gtex_thyroid_meta.txt", row.names = c(1), stringsAsFactors=F)
thyroid_gtex_meta <- cbind(sample=rownames(thyroid_gtex_meta), thyroid_gtex_meta)


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


thyroid_tcga_gtex_counts <- log2(thyroid_tcga_gtex_counts+1)




# PCA thyroid gtex
pca_thyroid_gtex <- prcomp(t(thyroid_tcga_gtex_counts[, colnames(thyroid_tcga_gtex_counts) %in% thyroid_tcga_gtex_meta[thyroid_tcga_gtex_meta$data == "GTEx", "sample"]]), center = TRUE, scale. = TRUE)$x %>% as.data.frame()
pca_thyroid_gtex <- cbind(sample = rownames(pca_thyroid_gtex), pca_thyroid_gtex)
pca_thyroid_gtex <- pca_thyroid_gtex %>%
  as.tibble() %>%
  gather(key = "PC", value = "value", -sample) %>%
  inner_join(thyroid_gtex_meta[, c("sample", "SMRIN", "AGE", "ETHNCTY", "MHCANCERNM", "SMCENTER", "SMTSTPTREF", "SMNABTCHT", "SMTSISCH", "GENDER")], by = "sample") %>%
  mutate(GENDER = if_else(GENDER == 1, "male", "female")) %>%
  mutate(MHCANCERNM = as.character(MHCANCERNM))




# correlation of gtex thyroid principal components to covariates
# select 10 first principal components
thyroid_gtex_covars <- pca_thyroid_gtex %>%
  filter(PC %in% paste0("PC", c(1:10))) %>%
  pivot_wider(names_from="PC", values_from="value")

thyroid_gtex_pcs <- thyroid_gtex_covars %>%
  dplyr::select(sample, starts_with("PC"))

thyroid_gtex_covars <- thyroid_gtex_covars %>%
    dplyr::select(-starts_with("PC")) %>%
    mutate(ETHNCTY = as.character(ETHNCTY)) %>%
    mutate(SMTSTPTREF = if_else(SMTSTPTREF == "", "unknown", SMTSTPTREF)) %>%
    mutate(SMTSISCH = if_else(is.na(SMTSISCH), round(mean(SMTSISCH, na.rm = T)), as.numeric(SMTSISCH))) %>%
    mutate(MHCANCERNM = if_else(is.na(MHCANCERNM), "2", MHCANCERNM)) %>%
    mutate(SMNABTCHT = factor(SMNABTCHT, labels = c("1", "2", "3"))) %>%
    mutate(SMTSTPTREF = factor(SMTSTPTREF, labels = c("1", "2"))) %>%
    mutate_if(is.character, as.factor)

#thyroid_gtex_covars <- model.matrix( as.formula(paste0("~", paste(c(colnames(thyroid_gtex_covars[, -c(1)]), "0"), collapse="+"))), thyroid_gtex_covars )
thyroid_gtex_covars <- model.matrix(
  object = as.formula(~ SMRIN + AGE + ETHNCTY + MHCANCERNM + SMCENTER + SMTSTPTREF + SMNABTCHT + SMTSISCH + GENDER + 0),
  data = thyroid_gtex_covars,
  contrasts.arg=list(
    ETHNCTY=contrasts(thyroid_gtex_covars$ETHNCTY, contrasts=F),
    MHCANCERNM=contrasts(thyroid_gtex_covars$MHCANCERNM, contrasts=F),
    SMCENTER=contrasts(thyroid_gtex_covars$SMCENTER, contrasts=F),
    SMTSTPTREF=contrasts(thyroid_gtex_covars$SMTSTPTREF, contrasts=F),
    SMNABTCHT=contrasts(thyroid_gtex_covars$SMNABTCHT, contrasts=F),
    GENDER=contrasts(thyroid_gtex_covars$GENDER, contrasts=F)))

thyroid_gtex_cors <- cor(thyroid_gtex_pcs[, -c(1)], thyroid_gtex_covars)

thyroid_gtex_cors_heatmap <- Heatmap(
  matrix = thyroid_gtex_cors,
  name = "Pearson's r",
  border = F,
  cluster_rows = T,
  cluster_columns = T,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 2))

pdf(file = "./plots/pre_processing/thyroid_gtex_cors_pcs_covars_heatmap.pdf", width = 8, height = 5)
thyroid_gtex_cors_heatmap
dev.off()

png(file = "./plots/pre_processing/thyroid_gtex_cors_pcs_covars_heatmap.png", width = 8, height = 5, units = "in", res = 300)
thyroid_gtex_cors_heatmap
dev.off()



# regress-out covars from gene expression

thyroid_meta <- thyroid_gtex_meta %>%
  as_tibble() %>%
  #select(c("sample", "SMRIN", "AGE", "ETHNCTY", "MHCANCERNM", "SMCENTER", "SMTSTPTREF", "SMNABTCHT", "SMTSISCH", "GENDER")) %>%
  select(c("sample", "SMRIN", "AGE", "ETHNCTY", "MHCANCERNM", "SMCENTER", "SMTSTPTREF", "SMNABTCHT", "SMTSISCH")) %>%
  mutate(MHCANCERNM = if_else(is.na(MHCANCERNM), "2", as.character(MHCANCERNM))) %>%
  mutate(SMTSISCH = if_else(is.na(SMTSISCH), round(mean(SMTSISCH, na.rm = T)), as.numeric(SMTSISCH))) %>%
  mutate(ETHNCTY = as.character(ETHNCTY))

thyroid_expr <- thyroid_tcga_gtex_counts[, colnames(thyroid_tcga_gtex_counts) %in% thyroid_tcga_gtex_meta[thyroid_tcga_gtex_meta$data == "GTEx", "sample"]]

thyroid_expr <- thyroid_expr %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "log2fpkm") %>%
  group_by(gene) %>%
  nest() %>%
  mutate(data2 = map(.x = data, .f = remove_batch2, covs = thyroid_meta)) %>%
  select(-data) %>%
  unnest()

thyroid_expr_mat <- thyroid_expr %>%
  pivot_wider(names_from="sample", values_from="residual") %>%
  column_to_rownames(var = "gene")


# PCA thyroid gtex
pca_thyroid_gtex <- prcomp(t(thyroid_expr_mat), center = TRUE, scale. = TRUE)$x %>%
  as.data.frame()
pca_thyroid_gtex <- cbind(sample = rownames(pca_thyroid_gtex), pca_thyroid_gtex)
pca_thyroid_gtex <- pca_thyroid_gtex %>%
  as.tibble() %>%
  gather(key = "PC", value = "value", -sample)


thyroid_gtex_pcs <- pca_thyroid_gtex %>%
  filter(PC %in% paste0("PC", c(1:10))) %>%
  pivot_wider(names_from="PC", values_from="value")


thyroid_gtex_cors <- cor(thyroid_gtex_pcs[, -c(1)], thyroid_gtex_covars)


thyroid_gtex_cors_heatmap <- Heatmap(
  matrix = thyroid_gtex_cors,
  name = "Pearson's r",
  border = F,
  cluster_rows = T,
  cluster_columns = T,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 2))

pdf(file = "./plots/pre_processing/thyroid_gtex_cors_pcs_covars_heatmapRegOut.pdf", width = 8, height = 5)
thyroid_gtex_cors_heatmap
dev.off()

png(file = "./plots/pre_processing/thyroid_gtex_cors_pcs_covars_heatmapRegOut.png", width = 8, height = 5, units = "in", res = 300)
thyroid_gtex_cors_heatmap
dev.off()


thyroid_expr_mat <- thyroid_expr %>%
  pivot_wider(names_from="sample", values_from="residual")


write_tsv(thyroid_expr_mat, "./files/thyroid_tcga_gtex_fpkm_GTEX_regressed_out_COVARS.txt")
