# Understanding Gender Differential Susceptibility in Cancer


# Comparison of the methylation status of differentially expressed genes between sample groups


library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(cowplot)


# load datasets


# differentially expressed genes
thca_TvsN <- read_tsv("./files/thyroid_males_females_signf_degs.txt")
thca_MvsF <- read_tsv("./files/thyroid_tumour_normal_signf_degs.txt")

stad_TvsN <- read_tsv("./files/stomach_males_females_signf_degs.txt")
stad_MvsF <- read_tsv("./files/stomach_tumour_normal_signf_degs.txt")



# methylation data

# thyroid
thca_meth <- read_tsv("./data/tcga/firebrowser/gdac.broadinstitute.org_THCA.Methylation_Preprocess.Level_3.2016012800.0.0/THCA.meth.by_min_expr_corr.data.txt")
thca_meth <- thca_meth[-c(1), ]

thca_meth <- thca_meth %>%
	dplyr::select(-c(X2, X3, X4)) %>%
	dplyr::rename(gene = `Hybridization REF`) %>%
	gather(key = "sample", value = "beta_value", -gene) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(beta_value = as.numeric(beta_value))

thca_meth2 <- read_tsv("./data/tcga/cbioportal/tcga_data_provisional/thca_tcga/data_methylation_hm450.txt")

thca_meth2 <- thca_meth2 %>%
	dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(gene = Hugo_Symbol) %>%
  gather(key = "sample", value = "beta_value", -gene) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(beta_value = as.numeric(beta_value))


# correlate both datasets
thca_meth_cor <- inner_join(thca_meth, thca_meth2, by = c("gene", "sample")) %>%
	na.exclude() %>%
	do(broom::tidy(cor.test(.$beta_value.x, .$beta_value.y, method = "pearson")))
#0.835


# stomach
stad_meth <- read_tsv("./data/tcga/firebrowser/gdac.broadinstitute.org_STAD.Methylation_Preprocess.Level_3.2016012800.0.0/STAD.meth.by_min_expr_corr.data.txt")
stad_meth <- stad_meth[-c(1), ]

stad_meth <- stad_meth %>%
  dplyr::select(-c(X2, X3, X4)) %>%
  dplyr::rename(gene = `Hybridization REF`) %>%
  gather(key = "sample", value = "beta_value", -gene) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(beta_value = as.numeric(beta_value))


stad_meth2_450 <- read_tsv("./data/tcga/cbioportal/tcga_data_provisional/stad_tcga/data_methylation_hm450.txt")

stad_meth2_450 <- stad_meth2_450 %>%
	dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(gene = Hugo_Symbol) %>%
  gather(key = "sample", value = "beta_value", -gene) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(beta_value = as.numeric(beta_value))


stad_meth2_450_normals <- read_tsv("./data/tcga/cbioportal/tcga_data_provisional/stad_tcga/data_methylation_hm450_normals.txt")

stad_meth2_450_normals <- stad_meth2_450_normals %>%
	dplyr::select(-Entrez_Gene_Id) %>%
	dplyr::rename(gene = Hugo_Symbol) %>%
	gather(key = "sample", value = "beta_value", -gene) %>%
	mutate(sample = str_replace_all(sample, "-", ".")) %>%
	mutate(beta_value = as.numeric(beta_value))


stad_meth2_27 <- read_tsv("./data/tcga/cbioportal/tcga_data_provisional/stad_tcga/data_methylation_hm27.txt")

stad_meth2_27 <- stad_meth2_27 %>%
	dplyr::select(-Entrez_Gene_Id) %>%
	dplyr::rename(gene = Hugo_Symbol) %>%
	gather(key = "sample", value = "beta_value", -gene) %>%
	mutate(sample = str_replace_all(sample, "-", ".")) %>%
	mutate(beta_value = as.numeric(beta_value))


stad_meth2_27_normals <- read_tsv("./data/tcga/cbioportal/tcga_data_provisional/stad_tcga/data_methylation_hm27_normals.txt")

stad_meth2_27_normals <- stad_meth2_27_normals %>%
	dplyr::select(-Entrez_Gene_Id) %>%
	dplyr::rename(gene = Hugo_Symbol) %>%
	gather(key = "sample", value = "beta_value", -gene) %>%
	mutate(sample = str_replace_all(sample, "-", ".")) %>%
	mutate(beta_value = as.numeric(beta_value))


stad_meth2 <- bind_rows(stad_meth2_450, stad_meth2_450_normals, stad_meth2_27, stad_meth2_27_normals)


# correlate both datasets
stad_meth_cor <- inner_join(stad_meth, stad_meth2, by = c("gene", "sample")) %>%
	na.exclude() %>%
	do(broom::tidy(cor.test(.$beta_value.x, .$beta_value.y, method = "pearson")))
#1



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







# thyroid DEGs

# tumour vs normal

thca_TvsN_meth <- thca_TvsN %>%
  filter(state != "common") %>%
  dplyr::select(genes, geneName, state, logFC_males, logFC_females) %>%
  gather(key = "foldchange", value = "logFC", -c(genes, geneName, state)) %>%
  filter(!is.na(logFC)) %>%
  dplyr::select(-foldchange) %>%
  inner_join(thyroid_tcga_fpkm %>% dplyr::select(gene, geneName, sample, fpkm, gender, sample_type_id), by = c("genes" = "gene", "geneName")) %>%
  filter(!(state == "male_specific" & gender == "female")) %>%
  filter(!(state == "female_specific" & gender == "male")) %>%
  inner_join(thca_meth, by = c("geneName" = "gene", "sample"))
#thca_TvsN_meth %>% group_by(state, sample_type_id) %>% summarise(median_beta = median(beta_value, na.rm=T)) %>% ungroup()
#thca_TvsN_meth %>% group_by(state, sample_type_id) %>% summarise(count = length(unique(sample)))


thca_TvsN_meth_boxp <- thca_TvsN_meth %>%
  ggplot(mapping = aes(x = sample_type_id, y = beta_value, fill = sample_type_id)) +
  geom_boxplot() +
  facet_wrap( ~ state,
    labeller = labeller(state = c("female_specific" = "Female specific", "male_specific" = "Male specific"))) +
  stat_compare_means() +
  theme_classic() +
  theme(
    axis.title = element_text(color = "black", size=15),
    axis.text = element_text(color="black", size=13),
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=15)) +
  scale_fill_manual(values=c("#d95f02", "#7fc97f"), guide=F) +
  scale_x_discrete(labels=c("Tumour", "Normal")) +
  labs(x = "Tissue", y = "Beta value")

thca_TvsN_fpkm_boxp <- thca_TvsN_meth %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = sample_type_id)) +
  geom_boxplot() +
  facet_wrap( ~ state,
    labeller = labeller(state = c("female_specific" = "Female specific", "male_specific" = "Male specific"))) +
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

thca_TvsN_plot <- plot_grid(thca_TvsN_fpkm_boxp, thca_TvsN_meth_boxp, labels = c("", ""), nrow = 2, align = "v")
ggsave(filename="thca_TvsN_rna_meth_plot.png", plot=thca_TvsN_plot, height=12, width=6, path = "./plots/methylation_analysis/")
unlink("thca_TvsN_rna_meth_plot.png")



# DEGs that are TSGs
tsgs <- read_tsv("./files/degs_TvsN_stomach_thyroid_tsgs.txt") %>%
	filter(tissue == "Thyroid" & state != "common")

thca_TvsN_meth_boxp_tsgs <- thca_TvsN_meth %>%
	filter(geneName %in% tsgs$geneName) %>%
	ggplot(mapping = aes(x = sample_type_id, y = beta_value, fill = sample_type_id)) +
	geom_boxplot() +
	facet_wrap( ~ state,
		labeller = labeller(state = c("female_specific" = "Female-specific", "male_specific" = "Male-specific"))) +
	stat_compare_means() +
	theme_classic() +
	theme(
		axis.title = element_text(color = "black", size=18),
	  axis.text = element_text(color="black", size=14),
	  strip.background = element_blank(),
	  strip.text = element_text(color="black", size=18),
		plot.title = element_text(color="black", size=18)) +
	scale_fill_manual(values=c("#ca0020", "#a1d76a"), guide=F) +
	scale_x_discrete(labels=c("Tumour", "Normal")) +
	labs(x = "", y = "Beta value", title = "Methylation")

thca_TvsN_fpkm_boxp_tsgs <- thca_TvsN_meth %>%
	filter(geneName %in% tsgs$geneName) %>%
	ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = sample_type_id)) +
	geom_boxplot() +
	facet_wrap( ~ state,
		labeller = labeller(state = c("female_specific" = "Female-specific", "male_specific" = "Male-specific"))) +
	stat_compare_means() +
	theme_classic() +
	theme(
		axis.title = element_text(color = "black", size=18),
		axis.text = element_text(color="black", size=14),
		strip.background = element_blank(),
		strip.text = element_text(color="black", size=18),
		plot.title = element_text(color="black", size=18)) +
	scale_fill_manual(values=c("#ca0020", "#a1d76a"), guide=F) +
	scale_x_discrete(labels=c("Tumour", "Normal")) +
	scale_y_continuous(trans = "log2") +
	labs(x = "Tissue", y = "FPKM (log2 scale)", title = "Expression")

thca_TvsN_plot_tsgs <- plot_grid(thca_TvsN_meth_boxp_tsgs, thca_TvsN_fpkm_boxp_tsgs, labels = c("", ""), nrow = 2, align = "v")
ggsave(filename="thca_TvsN_rna_meth_plot_tsgs.png", plot=thca_TvsN_plot_tsgs, height=12, width=6, path = "./plots/methylation_analysis/")
unlink("thca_TvsN_rna_meth_plot_tsgs.png")


# compare methylation and expression of female-specific DEGs (TSGs) between tumour and normal tissues in females and males
thca_female_specific_TvsN_tsgs <- tsgs %>%
	filter(state == "female_specific") %>%
	dplyr::select(tissue, geneName, state) %>%
	inner_join(thyroid_tcga_fpkm %>% dplyr::select(geneName, sample, fpkm, gender, sample_type_id), by = c("geneName")) %>%
	inner_join(thca_meth, by = c("geneName" = "gene", "sample"))


foo1 <- thca_female_specific_TvsN_tsgs %>%
	ggplot(mapping = aes(x = sample_type_id, y = beta_value, fill = sample_type_id)) +
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
	labs(x = "Tissue", y = "Beta value")

foo2 <- thca_female_specific_TvsN_tsgs %>%
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

foo3 <- plot_grid(foo2, foo1, labels = c("", ""), nrow = 2, align = "v")
ggsave(filename="thca_female_specific_TvsN_tsgs_male_female.png", plot=foo3, height=12, width=6, path = "./plots/methylation_analysis/")
unlink("thca_female_specific_TvsN_tsgs_male_female.png")






# male vs female in tumours

thca_MvsF_meth <- thca_MvsF %>%
  filter(state == "tumour_specific") %>%
  dplyr::select(genes, geneName, state, logFC = log2FC_tumour) %>%
  inner_join(thyroid_tcga_fpkm %>% dplyr::select(gene, geneName, sample, fpkm, gender, sample_type_id), by = c("genes" = "gene", "geneName")) %>%
  filter(!(sample_type_id == "11")) %>%
  inner_join(thca_meth2, by = c("geneName" = "gene", "sample"))
#thca_MvsF_meth %>% group_by(state, gender) %>% summarise(median_beta = median(beta_value, na.rm=T)) %>% ungroup()
#thca_MvsF_meth %>% group_by(state, gender) %>% summarise(count = length(unique(sample))) %>% ungroup()


thca_MvsF_meth_boxp <- thca_MvsF_meth %>%
  ggplot(mapping = aes(x = gender, y = beta_value, fill = gender)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_classic() +
  theme(
    axis.title = element_text(color = "black", size=15),
    axis.text = element_text(color="black", size=13)) +
  scale_fill_manual(values=c("pink", "blue"), guide=F) +
  scale_x_discrete(labels=c("Female", "Male")) +
  labs(x = "Gender", y = "Beta value")

thca_MvsF_fpkm_boxp <- thca_MvsF_meth %>%
  ggplot(mapping = aes(x = gender, y = fpkm, fill = gender)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_classic() +
  theme(
    axis.title = element_text(color = "black", size=15),
    axis.text = element_text(color="black", size=13)) +
  scale_fill_manual(values=c("pink", "blue"), guide=F) +
  scale_x_discrete(labels=c("Female", "Male")) +
  scale_y_continuous(trans = "log2") +
  labs(x = "Gender", y = "FPKM (log2 scale)")

thca_MvsF_plot <- plot_grid(thca_MvsF_fpkm_boxp, thca_MvsF_meth_boxp, labels = c("", ""), nrow = 2, align = "v")
ggsave(filename="thca_MvsF_tumours_rna_meth_plot.png", plot=thca_MvsF_plot, height=12, width=6, path = "./plots/methylation_analysis/")
unlink("thca_MvsF_tumours_rna_meth_plot.png")



# stomach DEGs

# tumour vs normal

stad_TvsN_meth <- stad_TvsN %>%
  filter(state != "common") %>%
  dplyr::select(genes, geneName, state, logFC_males, logFC_females) %>%
  gather(key = "foldchange", value = "logFC", -c(genes, geneName, state)) %>%
  filter(!is.na(logFC)) %>%
  dplyr::select(-foldchange) %>%
  inner_join(stomach_tcga_fpkm %>% dplyr::select(gene, geneName, sample, fpkm, gender, sample_type_id), by = c("genes" = "gene", "geneName")) %>%
  filter(!(state == "male_specific" & gender == "female")) %>%
  filter(!(state == "female_specific" & gender == "male"))
  #inner_join(stad_meth2, by = c("geneName" = "gene", "sample"))
#stad_TvsN_meth %>% group_by(state, sample_type_id) %>% summarise(median_beta = median(beta_value, na.rm=T)) %>% ungroup()
#stad_TvsN_meth %>% group_by(state, sample_type_id) %>% summarise(count = length(unique(sample)))

stad_TvsN_fpkm_boxp <- stad_TvsN_meth %>%
  ggplot(mapping = aes(x = sample_type_id, y = fpkm, fill = sample_type_id)) +
  geom_boxplot() +
  facet_wrap( ~ state,
    labeller = labeller(state = c("female_specific" = "Female specific", "male_specific" = "Male specific"))) +
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
ggsave(filename="stad_TvsN_fpkm_plot.png", plot=stad_TvsN_fpkm_boxp, height=6, width=6, path = "./plots/methylation_analysis/")
unlink("stad_TvsN_fpkm_plot.png")



# male vs female in tumours

stad_MvsF_meth <- stad_MvsF %>%
  filter(state == "tumour_specific") %>%
  dplyr::select(genes, geneName, state, logFC = log2FC_tumour) %>%
  inner_join(stomach_tcga_fpkm %>% dplyr::select(gene, geneName, sample, fpkm, gender, sample_type_id), by = c("genes" = "gene", "geneName")) %>%
  filter(!(state == "tumour_specific" & sample_type_id == "11")) %>%
  inner_join(stad_meth2, by = c("geneName" = "gene", "sample"))
#stad_MvsF_meth %>% group_by(state, gender) %>% summarise(median_beta = median(beta_value, na.rm=T)) %>% ungroup()
#stad_MvsF_meth %>% group_by(state, gender) %>% summarise(count = length(unique(sample))) %>% ungroup()


stad_MvsF_meth_boxp <- stad_MvsF_meth %>%
  ggplot(mapping = aes(x = gender, y = beta_value, fill = gender)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_classic() +
  theme(
    axis.title = element_text(color = "black", size=15),
    axis.text = element_text(color="black", size=13)) +
  scale_fill_manual(values=c("pink", "blue"), guide=F) +
  scale_x_discrete(labels=c("Female", "Male")) +
  labs(x = "Gender", y = "Beta value")

stad_MvsF_fpkm_boxp <- stad_MvsF_meth %>%
  ggplot(mapping = aes(x = gender, y = fpkm, fill = gender)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_classic() +
  theme(
    axis.title = element_text(color = "black", size=15),
    axis.text = element_text(color="black", size=13)) +
  scale_fill_manual(values=c("pink", "blue"), guide=F) +
  scale_x_discrete(labels=c("Female", "Male")) +
  scale_y_continuous(trans = "log2") +
  labs(x = "Gender", y = "FPKM (log2 scale)")

stad_MvsF_plot <- plot_grid(stad_MvsF_fpkm_boxp, stad_MvsF_meth_boxp, labels = c("", ""), nrow = 2, align = "v")
ggsave(filename="stad_MvsF_tumours_rna_meth_plot.png", plot=stad_MvsF_plot, height=12, width=6, path = "./plots/methylation_analysis/")
unlink("stad_MvsF_tumours_rna_meth_plot.png")
