# Understanding Gender Differential Susceptibility in Cancer


# Descriptive plots about number and type of samples, between others.


library(tidyverse)
library(RColorBrewer)



# -- Thyroid

sample_type_code <- c("GTEx" = "GTEx", "TCGA_normal" = "TCGA", "TCGA_tumour" = "TCGA")
tissue_type_code <- c("GTEx" = "normal", "TCGA_normal" = "normal", "TCGA_tumour" = "tumour")


# load metadata
thyroid_tcga_gtex_meta <- read.delim("./files/thyroid_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_gtex_meta <- thyroid_tcga_gtex_meta %>%
	mutate(sample = gsub("-", ".", rownames(thyroid_tcga_gtex_meta))) %>%
	as.tibble() %>%
	# remove metastatic samples
	filter(!(sample_type == "TCGA_metastatic")) %>%
	# add experimental group
	mutate(data = sample_type_code[sample_type]) %>%
	# add tissue type
	mutate(tissue_type = tissue_type_code[sample_type]) %>%
	mutate(tissue = "Thyroid") %>%
	dplyr::select(sample, everything())


stomach_tcga_gtex_meta <- read.delim("./files/stomach_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_gtex_meta <- stomach_tcga_gtex_meta %>%
	mutate(sample = gsub("-", ".", rownames(stomach_tcga_gtex_meta))) %>%
	as.tibble() %>%
	# remove metastatic samples
	filter(!(sample_type == "TCGA_metastatic")) %>%
	# add experimental group
	mutate(data = sample_type_code[sample_type]) %>%
	# add tissue type
	mutate(tissue_type = tissue_type_code[sample_type]) %>%
	mutate(tissue = "Stomach") %>%
	dplyr::select(sample, everything())


thyroid_stomach_metadata <- bind_rows(thyroid_tcga_gtex_meta, stomach_tcga_gtex_meta)


ggplot(data=thyroid_stomach_metadata, mapping=aes(x=tissue, y=..count.., fill=tissue_type)) +
	geom_bar(position="dodge") +
	scale_fill_manual(values=c("#99d594", "#fc8d59"), name = "Tissue type") +
	theme_classic() +
	theme(axis.title.y = element_text(colour="black", size=15),
		axis.title.x = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        legend.text=element_text(colour="black", size=13),
        legend.title=element_text(colour="black", size=15)) +
    labs(y = "Number of samples", x = "Tissue", title="")
ggsave(filename="thyroid_stomach_samples_tissue_type.png", path = "./plots/summary_plots")
unlink("thyroid_stomach_samples_tissue_type.png")


ggplot(data=thyroid_stomach_metadata, mapping=aes(x=tissue, y=..count.., fill=gender)) +
	geom_bar(position="dodge") +
	scale_fill_manual(values=c("pink", "blue"), name="Gender") +
	theme_classic() +
	theme(axis.title.y = element_text(colour="black", size=15),
		axis.title.x = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        legend.text=element_text(colour="black", size=13),
        legend.title=element_text(colour="black", size=15)) +
  labs(y = "Number of samples", x = "Tissue")
ggsave(filename="thyroid_stomach_samples_gender.png", path = "./plots/summary_plots", width=6, height=6)
unlink("thyroid_stomach_samples_gender.png")


ggplot(data=thyroid_stomach_metadata, mapping=aes(x=tissue, y=..count.., fill=sample_type)) +
	geom_bar(position="dodge") +
	scale_fill_manual(values=c("#99d594", "#7fbf7b", "#fc8d59"), name = "Sample type") +
	theme_classic() +
	theme(axis.title.y = element_text(colour="black", size=15),
		axis.title.x = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        legend.text=element_text(colour="black", size=13),
        legend.title=element_text(colour="black", size=15)) +
    labs(y = "Number of samples", x = "Tissue", title="")
ggsave(filename="thyroid_stomach_samples_data_type.png", path = "./plots/summary_plots")
unlink("thyroid_stomach_samples_data_type.png")


ggplot(data=thyroid_stomach_metadata, mapping=aes(x=tissue, y=..count.., fill=sample_type, color=gender)) +
	geom_bar(position="dodge", size = 1) +
	scale_fill_manual(values=c("#99d594", "#7fbf7b", "#fc8d59"), labels=c("GTEx", "TCGA normal", "TCGA tumour"), name="Sample type") +
	scale_color_manual(values=c("pink", "blue"), name="Gender") +
	theme_classic() +
	theme(axis.title.y = element_text(colour="black", size=15),
		axis.title.x = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13),
        legend.text=element_text(colour="black", size=13),
        legend.title=element_text(colour="black", size=15)) +
	scale_y_continuous(limits = c(NA, 400)) +
    labs(y = "Number of samples", x = "Tissue")
ggsave(filename="thyroid_stomach_samples_all.png", path = "./plots/summary_plots", width=6, height=6)
unlink("thyroid_stomach_samples_all.png")


save(list=ls(), file="r_workspaces/summary_plots.RData")
