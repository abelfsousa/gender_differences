# Understanding Gender Differential Susceptibility in Cancer


# Comparison of the methylation status of differentially expressed genes between sample groups


library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(cowplot)



# load datasets

thyroid_tcga_meta <- read_tsv("./files/tcga_thca_meta.txt") %>%
  dplyr::rename(sample = submitter_id) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  filter(sample_type_id != "06")

stomach_tcga_meta <- read_tsv("./files/tcga_stad_meta.txt") %>%
  dplyr::rename(sample = submitter_id) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  filter(sample_type_id != "06")


# differentially expressed genes
thca_TvsN <- read_tsv("./files/thyroid_males_females_signf_degs.txt")
thca_MvsF <- read_tsv("./files/thyroid_tumour_normal_signf_degs.txt")

stad_TvsN <- read_tsv("./files/stomach_males_females_signf_degs.txt")
stad_MvsF <- read_tsv("./files/stomach_tumour_normal_signf_degs.txt")


# methylation data

# thyroid
thca_meth <- read_tsv("./data/tcga/firebrowser/gdac.broadinstitute.org_THCA.Methylation_Preprocess.Level_3.2016012800.0.0/THCA_methylation450_mean_beta_promoters.txt.gz") %>%
  rename_all(.funs = ~ str_sub(.x, 1, 15)) %>%
  pivot_longer(-Gene_Symbol, names_to = "sample", values_to = "beta") %>%
  rename(gene = Gene_Symbol) %>%
  filter(sample %in% thyroid_tcga_meta$sample) %>%
  filter(str_detect(gene, ";", negate = T))


# stomach
stad_meth <- read_tsv("./data/tcga/firebrowser/gdac.broadinstitute.org_STAD.Methylation_Preprocess.Level_3.2016012800.0.0/STAD_methylation450_mean_beta_promoters.txt.gz") %>%
  rename_all(.funs = ~ str_sub(.x, 1, 15)) %>%
  pivot_longer(-Gene_Symbol, names_to = "sample", values_to = "beta") %>%
  rename(gene = Gene_Symbol) %>%
  filter(sample %in% stomach_tcga_meta$sample) %>%
  filter(str_detect(gene, ";", negate = T))



# compare beta distribution between DEGs / SBGs

# Thyroid T vs N DEGs
comparison_boxplot <- thca_meth %>%
  left_join(thca_TvsN[, c("geneName", "state")], by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background")) %>%
  ggplot(mapping = aes(x = state, y = beta, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3), c(1,4))) +
  scale_x_discrete(labels = c("background", "common", "female-specific", "male-specific")) +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black")) +
  labs(x = "DEGs", y = "Beta (mean values across promoter probes)", title = "Thyroid")

ggsave(filename="thca_TvsN_beta_distribution.png", plot=comparison_boxplot, height=6, width=4, path = "./plots/methylation_analysis_promoters_distributions/")
unlink("thca_TvsN_beta_distribution.png")

comparison_density <- thca_meth %>%
  left_join(thca_TvsN[, c("geneName", "state")], by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background")) %>%
  ggplot(mapping = aes(x = beta, color = state)) +
  geom_density() +
  scale_color_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black"))


# Stomach T vs N DEGs
comparison_boxplot <- stad_meth %>%
  left_join(stad_TvsN[, c("geneName", "state")], by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background")) %>%
  ggplot(mapping = aes(x = state, y = beta, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3), c(1,4)), label.y = c(0.4, 0.45, 0.5), tip.length = 0.01) +
  scale_x_discrete(labels = c("background", "common", "female-specific", "male-specific")) +
  scale_y_continuous(limits = c(NA, 0.5)) +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black")) +
  labs(x = "DEGs", y = "Beta (mean values across promoter probes)", title = "Stomach")

ggsave(filename="stad_TvsN_beta_distribution.png", plot=comparison_boxplot, height=6, width=4, path = "./plots/methylation_analysis_promoters_distributions/")
unlink("stad_TvsN_beta_distribution.png")

comparison_density <- stad_meth %>%
  left_join(stad_TvsN[, c("geneName", "state")], by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background")) %>%
  ggplot(mapping = aes(x = beta, color = state)) +
  geom_density() +
  scale_color_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black"))



# Thyroid/Stomach T vs N DEGs
thyroid <- thca_meth %>%
  left_join(thca_TvsN[, c("geneName", "state")], by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background"), tissue = "Thyroid")

stomach <- stad_meth %>%
  left_join(stad_TvsN[, c("geneName", "state")], by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background"), tissue = "Stomach")


comparison_boxplot <- thyroid %>%
  bind_rows(stomach) %>%
  ggplot(mapping = aes(x = state, y = beta, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  facet_wrap(~ tissue) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3), c(1,4))) +
  scale_x_discrete(labels = c("background", "common", "female\nspecific", "male\nspecific")) +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, color = "black")) +
  labs(x = "DEGs", y = "Beta (mean values across promoter probes)")

ggsave(filename="thca_stad_TvsN_beta_distribution.png", plot=comparison_boxplot, height=6, width=6, path = "./plots/methylation_analysis_promoters_distributions/")
unlink("thca_stad_TvsN_beta_distribution.png")



# Thyroid/Stomach M vs F SBGs
thyroid <- thca_meth %>%
  mutate(sample_type = str_sub(sample, 14, 15)) %>%
  filter(sample_type != "11") %>%
  left_join(filter(thca_MvsF, str_detect(state, "normal", negate=T)) %>% select(geneName, state), by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background"), tissue = "Thyroid")

stomach <- stad_meth %>%
  mutate(sample_type = str_sub(sample, 14, 15)) %>%
  filter(sample_type != "11") %>%
  left_join(filter(stad_MvsF, str_detect(state, "normal", negate=T)) %>% select(geneName, state), by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background"), tissue = "Stomach")


comparison_boxplot <- thyroid %>%
  bind_rows(stomach) %>%
  ggplot(mapping = aes(x = state, y = beta, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  facet_wrap(~ tissue) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1,3)), label.y = 0.5, tip.length = 0.01) +
  scale_x_discrete(labels = c("background", "common", "tumour\nspecific")) +
  scale_y_continuous(limits = c(NA, 0.5)) +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, color = "black")) +
  labs(x = "SBGs", y = "Beta (mean values across promoter probes)")

ggsave(filename="thca_stad_MvsF_beta_distribution.png", plot=comparison_boxplot, height=6, width=6, path = "./plots/methylation_analysis_promoters_distributions/")
unlink("thca_stad_MvsF_beta_distribution.png")


comparison_density <- thyroid %>%
  bind_rows(stomach) %>%
  ggplot(mapping = aes(x = beta, color = state)) +
  geom_density() +
  facet_wrap(~ tissue) +
  scale_color_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black"))



# DEGs that are TSGs

thca_TvsN_meth <- thca_TvsN %>%
  filter(state != "common") %>%
  dplyr::select(genes, geneName, state, logFC_males, logFC_females) %>%
  gather(key = "foldchange", value = "logFC", -c(genes, geneName, state)) %>%
  filter(!is.na(logFC)) %>%
  dplyr::select(-foldchange) %>%
  inner_join(thca_meth, by = c("geneName" = "gene")) %>%
  inner_join(thyroid_tcga_meta[, c("sample", "gender", "sample_type_id")], by = c("sample")) %>%
  filter((state == "male_specific" & gender == "male") | (state == "female_specific" & gender == "female"))

tsgs <- read_tsv("./files/degs_TvsN_stomach_thyroid_tsgs.txt") %>%
	filter(tissue == "Thyroid" & state != "common") %>%
	gather(key = "logFC", value = "value", -c(genes, geneName, chrom, state, tissue, geneType, tsg)) %>%
	na.exclude() %>%
	dplyr::select(-logFC) %>%
	dplyr::rename(logFC = value) %>%
	filter(logFC < 0)

thca_TvsN_meth_boxp_tsgs <- thca_TvsN_meth %>%
	filter(geneName %in% tsgs$geneName) %>%
	ggplot(mapping = aes(x = sample_type_id, y = beta, fill = sample_type_id)) +
	geom_boxplot() +
	facet_wrap( ~ state,
		labeller = labeller(state = c("female_specific" = "Female-specific", "male_specific" = "Male-specific"))) +
	stat_compare_means(size = 5) +
	theme_classic() +
	theme(
		axis.title = element_text(color = "black", size=18),
	  axis.text = element_text(color="black", size=14),
	  strip.background = element_blank(),
	  strip.text = element_text(color="black", size=18),
		plot.title = element_text(color="black", size=18),
		legend.title = element_text(color="black", size=18),
		legend.text = element_text(color="black", size=14)) +
	scale_fill_manual(values=c("#ca0020", "#a1d76a"), labels = c("Tumour", "Normal"), name = "Tissue") +
	scale_x_discrete(labels=c("Tumour", "Normal")) +
	labs(x = "Tissue", y = "Beta value", title = "")



# compare variance distribution between DEGs / SBGs

# tumour vs normal comparison

thca_meth_var <- thca_meth %>%
  group_by(gene) %>%
  summarise(beta_var = var(beta, na.rm=T)) %>%
  ungroup() %>%
  filter(!is.na(beta_var))


stad_meth_var <- stad_meth %>%
  group_by(gene) %>%
  summarise(beta_var = var(beta, na.rm=T)) %>%
  ungroup() %>%
  filter(!is.na(beta_var))



comparison <- thca_meth_var %>%
  left_join(thca_TvsN[, c("geneName", "state")], by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background")) %>%
  ggplot(mapping = aes(x = state, y = beta_var, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3), c(1,4)), label.y = c(0.005, 0.006, 0.007), tip.length = 0.001) +
  scale_x_discrete(labels = c("background", "common", "female-specific", "male-specific")) +
  scale_y_continuous(limits = c(NA, 0.0075)) +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black")) +
  labs(x = "DEGs", y = "Variance of beta values", title = "Thyroid")

ggsave(filename="thca_TvsN_beta_var_distribution.png", plot=comparison, height=6, width=4, path = "./plots/methylation_analysis_promoters_distributions/")
unlink("thca_TvsN_beta_var_distribution.png")


comparison <- thca_meth_var %>%
  left_join(thca_TvsN[, c("geneName", "state")], by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background")) %>%
  ggplot(mapping = aes(x = beta_var, color = state)) +
  geom_density(alpha = 0.8) +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black")) +
  labs(x = "Variance", y = "Density")
#comparison


comparison <- stad_meth_var %>%
  left_join(stad_TvsN[, c("geneName", "state")], by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background")) %>%
  ggplot(mapping = aes(x = state, y = beta_var, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3), c(1,4)), label.y = c(0.03, 0.035, 0.04), tip.length = 0.01) +
  scale_x_discrete(labels = c("background", "common", "female-specific", "male-specific")) +
  scale_y_continuous(limits = c(NA, 0.04)) +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black")) +
  labs(x = "DEGs", y = "Variance of beta values", title = "Stomach")

ggsave(filename="stad_TvsN_beta_var_distribution.png", plot=comparison, height=6, width=4, path = "./plots/methylation_analysis_promoters_distributions/")
unlink("stad_TvsN_beta_var_distribution.png")




# male vs female comparison
thca_meth_var <- thca_meth %>%
  mutate(sample_type = str_sub(sample, 14, 15)) %>%
  filter(sample_type != "11") %>%
  group_by(gene) %>%
  summarise(beta_var = var(beta, na.rm=T)) %>%
  ungroup() %>%
  filter(!is.na(beta_var))


stad_meth_var <- stad_meth %>%
  mutate(sample_type = str_sub(sample, 14, 15)) %>%
  filter(sample_type != "11") %>%
  group_by(gene) %>%
  summarise(beta_var = var(beta, na.rm=T)) %>%
  ungroup() %>%
  filter(!is.na(beta_var))


comparison <- thca_meth_var %>%
  left_join(filter(thca_MvsF, str_detect(state, "normal", negate=T)) %>% select(geneName, state), by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background")) %>%
  ggplot(mapping = aes(x = state, y = beta_var, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3)), label.y = c(0.005, 0.006), tip.length = 0.001) +
  scale_x_discrete(labels = c("background", "common", "tumour-specific")) +
  scale_y_continuous(limits = c(NA, 0.0075)) +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black")) +
  labs(x = "DEGs", y = "Variance of beta values", title = "Thyroid")

ggsave(filename="thca_MvsF_beta_var_distribution.png", plot=comparison, height=6, width=4, path = "./plots/methylation_analysis_promoters_distributions/")
unlink("thca_MvsF_beta_var_distribution.png")


comparison <- stad_meth_var %>%
  left_join(filter(stad_MvsF, str_detect(state, "normal", negate=T)) %>% select(geneName, state), by = c("gene" = "geneName")) %>%
  mutate(state = replace_na(state, "background")) %>%
  ggplot(mapping = aes(x = state, y = beta_var, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2), c(1,3)), label.y = c(0.03, 0.035), tip.length = 0.01) +
  scale_x_discrete(labels = c("background", "common", "tumour-specific")) +
  scale_y_continuous(limits = c(NA, 0.04)) +
  scale_fill_brewer(type = "div", palette = "BrBG") +
  theme(
    axis.text = element_text(size = 8, color = "black")) +
  labs(x = "DEGs", y = "Variance of beta values", title = "Stomach")

ggsave(filename="stad_MvsF_beta_var_distribution.png", plot=comparison, height=6, width=4, path = "./plots/methylation_analysis_promoters_distributions/")
unlink("stad_MvsF_beta_var_distribution.png")




save(list=ls(), file="r_workspaces/methylation_analysis_promoters_distributions.RData")
