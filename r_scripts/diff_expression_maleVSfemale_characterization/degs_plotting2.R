# Understanding Gender Differential Susceptibility in Cancer


# DEGs plotting


library(tidyverse)
library(RColorBrewer)
library(viridis)
library(data.table)


tcga_annotation <- read_tsv("./data/annotation/geneAnnot.gencode.v22.txt") %>%
  mutate(geneID = str_replace(geneID, "\\.[0-9]+", "")) %>%
  dplyr::rename(ensemblGeneID = geneID)





# -- Thyroid


# load fpkm and metadata
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
  inner_join(thyroid_tcga_gtex_meta[, c("sample", "sample_type", "gender")], by = "sample") %>%
  filter(!(sample_type == "TCGA_metastatic" | sample_type == "TCGA_normal"))



# load degs
mf_compCl <- read_tsv("./files/thyroid_male_female_compCluster_GO.txt") %>%
  filter(Cluster == "up normal female") %>%
  filter(Description == "regulation of cell activation" | Description == "lymphocyte activation" | Description == "leukocyte cell-cell adhesion") %>%
  dplyr::select(Cluster, Description, geneID) %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest()

tn_degs <- read_tsv("./files/thyroid_tumour_normal_degs_enr.txt") %>%
  filter(state == "normal_specific") %>%
  filter(ID == "REGULATION OF CELL ACTIVATION" | ID == "LYMPHOCYTE ACTIVATION" | ID == "LEUKOCYTE CELL CELL ADHESION") %>%
  dplyr::select(state, ID, geneID) %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest()

# [1] "LEUKOCYTE MIGRATION"
# [2] "LYMPHOCYTE ACTIVATION"
# [3] "REGULATION OF CELL ACTIVATION"
# [4] "POSITIVE REGULATION OF IMMUNE RESPONSE"
# [5] "LEUKOCYTE CHEMOTAXIS"
# [6] "ACTIVATION OF IMMUNE RESPONSE"
# [7] "ADAPTIVE IMMUNE RESPONSE"
# [8] "CELL CHEMOTAXIS"
# [9] "POSITIVE REGULATION OF CELL ACTIVATION"
# [10] "CELLULAR EXTRAVASATION"                




thyroid_degs <- inner_join(mf_compCl, tn_degs, by=c("geneID")) %>%
  dplyr::select(geneID) %>%
  unique() %>%
  inner_join(tcga_annotation[, c("geneName", "ensemblGeneID")], by=c("geneID" = "geneName")) %>%
  inner_join(thyroid_tcga_gtex_fpkm, by=c("ensemblGeneID" = "gene"))


thyroid_degs_boxp <- thyroid_degs %>%
  ggplot(mapping = aes(x = sample_type, y = fpkm, fill = gender)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_classic() +
    facet_grid(~ geneID, scales = "free_y") +
    theme(
      axis.title = element_text(colour="black", size=14),
      axis.text = element_text(colour="black", size=8),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=14),
      strip.background = element_blank(),
      strip.text.x = element_text(colour="black", size=10)) +
    scale_y_continuous(trans = "log2") +
    labs(x = "Tissue", y = "FPKM (log2)")
ggsave(filename="thyroid_degs.png", plot=thyroid_degs_boxp, path = "./plots/diff_expression_maleVSfemale/", width=11, height=4)
unlink("thyroid_degs.png")
