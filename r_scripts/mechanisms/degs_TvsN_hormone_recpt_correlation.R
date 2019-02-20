# Understanding Gender Differential Susceptibility in Cancer


# Gene expression correlation
# Between hormone receptor genes and degs


library(tidyverse)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)




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
  mutate(receptor = list(horm_recpt)) %>%
  unnest() %>%
  dplyr::select(tissue, everything()) %>%
  mutate_if(is.character, as.factor)



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
  filter(sample_type == "TCGA_tumour") %>%
  mutate(tissue = "Thyroid") %>%
  dplyr::select(tissue, gene = geneName, sample, gender, fpkm)



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
  filter(sample_type == "TCGA_tumour") %>%
  mutate(tissue = "Stomach") %>%
  dplyr::select(tissue, gene = geneName, sample, gender, fpkm)



# merge fpkm from thyroid and stomach
fpkm <- bind_rows(stomach_tcga_gtex_fpkm, thyroid_tcga_gtex_fpkm) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(fpkm = log2(fpkm+1))



# merge fpkm data with degs/receptor table
degs_recpt_cor <- degs %>%
  inner_join(fpkm, by = c("tissue", "gene")) %>%
  dplyr::rename(fpkm_deg = fpkm) %>%
  inner_join(fpkm, by = c("tissue", "receptor" = "gene", "sample", "gender")) %>%
  dplyr::rename(fpkm_recpt = fpkm) %>%
  filter(!(state == "female_specific" & gender == "male")) %>%
  filter(!(state == "male_specific" & gender == "female")) %>%
  filter(!(gene == receptor))



# calculate Pearson r
degs_recpt_cor <- degs_recpt_cor %>%
  group_by(tissue, state, gene, receptor) %>%
  do(broom::tidy(cor.test(.$fpkm_deg, .$fpkm_recpt, method = "pearson"))) %>%
  ungroup()
write.table(degs_recpt_cor, "./files/degs_tumourVSnormal_cor_receptors.txt", sep="\t", quote=F, row.names=F)
