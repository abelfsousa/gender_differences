# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data

# Samples correlation



library(WGCNA)
library(tidyverse)
library(gplots)

source("utils.R")




# -- Stomach

# load R workspace
load("./r_workspaces/tcga_gtex_stomach_wgcna_networks.RData")




# build correlation matrixes
males_tumour_corMat <- cor(t(males_tumour))
males_normal_corMat <- cor(t(males_normal))
females_tumour_corMat <- cor(t(females_tumour))
females_normal_corMat <- cor(t(females_normal))



males_tumour_corMat2 <- males_tumour_corMat %>% get_upper_tri %>% as.numeric %>% na.exclude %>% as.numeric

males_normal_corMat2 <- males_normal_corMat %>% get_upper_tri %>% as.numeric %>% na.exclude %>% as.numeric

females_tumour_corMat2 <- females_tumour_corMat %>% get_upper_tri %>% as.numeric %>% na.exclude %>% as.numeric

females_normal_corMat2 <- females_normal_corMat %>% get_upper_tri %>% as.numeric %>% na.exclude %>% as.numeric


stomach_correlations <- tibble(
	data = c("males_tumor", "males_normal", "females_tumor", "females_normal"),
	pearson = list(males_tumour_corMat2, males_normal_corMat2, females_tumour_corMat2, females_normal_corMat2)) %>%
	unnest()



stomach_correlations_plot <- ggplot(data=stomach_correlations, mapping=aes(x=pearson, colour=data, fill=data)) +
    geom_density(alpha=0.5) +
    theme(axis.text = element_text(colour="black", size=13),
        axis.title = element_text(colour="black", size=15),
        plot.title = element_text(colour="black", size=15)) +
    scale_fill_discrete(guide=F) +
    labs(x="Pearson correlation coefficient", y="Density", colour="Data", title="Stomach")
ggsave("stomach_wgcna_samples_correlations_density_plot.png", plot=stomach_correlations_plot, path="./plots/")
unlink("stomach_wgcna_samples_correlations_density_plot.png")




