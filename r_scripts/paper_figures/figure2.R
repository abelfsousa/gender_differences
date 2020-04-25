# figure 2 for paper

library(tidyverse)
library(cowplot)


#load plots
volcano_stomach <- read_rds("./r_objects/plots/figure2/diff_expr_stomach_all_tcga_maleVSfemale_volcano_plot.rds")
volcano_thyroid <- read_rds("./r_objects/plots/figure2/diff_expr_thyroid_all_tcga_maleVSfemale_volcano_plot.rds")

venn_stomach <- read_rds("./r_objects/plots/figure2/stomach_tumour_normal_signf_degs_maleVsfemale_venn_ggplot.rds")
venn_thyroid <- read_rds("./r_objects/plots/figure2/thyroid_tumour_normal_signf_degs_maleVsfemale_venn_ggplot.rds")

chr_stomach <- read_rds("./r_objects/plots/figure2/stomach_degs_chr_barplot.rds")
chr_thyroid <- read_rds("./r_objects/plots/figure2/thyroid_degs_chr_barplot.rds")

enr_stomach <- read_rds("./r_objects/plots/figure2/stomach_tumour_normal_enrichment_barplot_normal_specific.rds")
enr_thyroid <- read_rds("./r_objects/plots/figure2/thyroid_tumour_normal_enrichment_barplot_normal_specific.rds")

#enr_stomach <- enr_stomach + theme(axis.text.y = element_text(size = 5))
#enr_thyroid <- enr_thyroid + theme(axis.text.y = element_text(size = 5))


first_row = plot_grid(volcano_stomach, volcano_thyroid, labels = c("A", "B"))
second_row = plot_grid(venn_stomach, chr_stomach, venn_thyroid, chr_stomach, labels = c("C", "D", "E", "F"), nrow=1)
third_row = plot_grid(enr_stomach, enr_thyroid, labels = c("G", "H"), nrow=1)


gg_all = plot_grid(first_row, second_row, third_row, labels=NULL, nrow=3)
ggsave(filename="figure2.png", plot = gg_all, path = "./plots/paper_figures/", width=20, height=30, units = "cm")
unlink("figure2.png")
