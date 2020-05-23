library(tidyverse)


annotation <- read_tsv("./data/annotation/geneAnnot.gencode.v22.txt") %>%
  mutate(gene = str_split(geneID, "\\.", simplify=T)[,1]) %>%
  dplyr::select(-geneID) %>%
  dplyr::select(gene, everything())



# TCGA thyroid

tcga_thyroid <- read_tsv(file = "./files/tcga_thca_fpkm.txt")

tcga_thyroid <- tcga_thyroid %>%
  pivot_longer(-gene, names_to = "sample", values_to = "fpkm") %>%
  mutate(fpkm = log2(fpkm+1)) %>%
  group_by(gene) %>%
  summarise(median_fpkm = median(fpkm), sd_fpkm = sd(fpkm)) %>%
  ungroup() %>%
  inner_join(annotation[, c("gene", "chrom")], by = "gene") %>%
  filter(chrom != "chrM") %>%
  mutate(chrom_type = map_chr(.x = chrom, .f = ~ if(.x %in% c("chrX", "chrY")){.x}else{"autosome"}))

per <- tcga_thyroid %>%
  summarise(n = n(), nn = sum(median_fpkm > 1), p = (nn/n)*100)

tcga_thyroid_plot <- tcga_thyroid %>%
  mutate(gene = fct_reorder(gene, median_fpkm, .fun = function(x) max(x))) %>%
  ggplot() +
  geom_point(mapping = aes(x = gene, y = median_fpkm), size = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.3) +
  annotate(geom = "text", x = 5000, y = 9, label = str_c("median expression > 1 = ", round(per$p,0), "%")) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("grey", "red", "blue")) +
  labs(y = "Median expression (log2 FPKM + 1)", color = "chr", title = "TCGA - Thyroid")

ggsave(filename="median_expression_tcga_thyroid.png", plot = tcga_thyroid_plot, path = "./plots/summary_plots", width=5, height=2)
ggsave(filename="median_expression_tcga_thyroid.pdf", plot = tcga_thyroid_plot, path = "./plots/summary_plots", width=5, height=2)
unlink("median_expression_tcga_thyroid.png")
unlink("median_expression_tcga_thyroid.pdf")



# TCGA stomach

tcga_stomach <- read_tsv(file = "./files/tcga_stad_fpkm.txt")

tcga_stomach <- tcga_stomach %>%
  pivot_longer(-gene, names_to = "sample", values_to = "fpkm") %>%
  mutate(fpkm = log2(fpkm+1)) %>%
  group_by(gene) %>%
  summarise(median_fpkm = median(fpkm), sd_fpkm = sd(fpkm)) %>%
  ungroup() %>%
  inner_join(annotation[, c("gene", "chrom")], by = "gene") %>%
  filter(chrom != "chrM") %>%
  mutate(chrom_type = map_chr(.x = chrom, .f = ~ if(.x %in% c("chrX", "chrY")){.x}else{"autosome"}))

per <- tcga_stomach %>%
  summarise(n = n(), nn = sum(median_fpkm > 1), p = (nn/n)*100)

tcga_stomach_plot <- tcga_stomach %>%
  mutate(gene = fct_reorder(gene, median_fpkm, .fun = function(x) max(x))) %>%
  ggplot() +
  geom_point(mapping = aes(x = gene, y = median_fpkm), size = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.3) +
  annotate(geom = "text", x = 5000, y = 9, label = str_c("median expression > 1 = ", round(per$p,0), "%")) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("grey", "red", "blue")) +
  labs(y = "Median expression (log2 FPKM + 1)", color = "chr", title = "TCGA - Stomach")

ggsave(filename="median_expression_tcga_stomach.png", plot = tcga_stomach_plot, path = "./plots/summary_plots", width=5, height=2)
ggsave(filename="median_expression_tcga_stomach.pdf", plot = tcga_stomach_plot, path = "./plots/summary_plots", width=5, height=2)
unlink("median_expression_tcga_stomach.png")
unlink("median_expression_tcga_stomach.pdf")




# GTEx thyroid

gtex_thyroid <- read_tsv(file = "./files/gtex_thyroid_fpkm.txt")

gtex_thyroid <- gtex_thyroid %>%
  pivot_longer(-gene, names_to = "sample", values_to = "fpkm") %>%
  mutate(fpkm = log2(fpkm+1)) %>%
  group_by(gene) %>%
  summarise(median_fpkm = median(fpkm), sd_fpkm = sd(fpkm)) %>%
  ungroup() %>%
  inner_join(annotation[, c("gene", "chrom")], by = "gene") %>%
  filter(chrom != "chrM") %>%
  mutate(chrom_type = map_chr(.x = chrom, .f = ~ if(.x %in% c("chrX", "chrY")){.x}else{"autosome"}))

per <- gtex_thyroid %>%
  summarise(n = n(), nn = sum(median_fpkm > 1), p = (nn/n)*100)

gtex_thyroid_plot <- gtex_thyroid %>%
  mutate(gene = fct_reorder(gene, median_fpkm, .fun = function(x) max(x))) %>%
  ggplot() +
  geom_point(mapping = aes(x = gene, y = median_fpkm), size = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.3) +
  annotate(geom = "text", x = 5000, y = 9, label = str_c("median expression > 1 = ", round(per$p,0), "%")) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("grey", "red", "blue")) +
  labs(y = "Median expression (log2 FPKM + 1)", color = "chr", title = "GTEx - Thyroid")

ggsave(filename="median_expression_gtex_thyroid.png", plot = gtex_thyroid_plot, path = "./plots/summary_plots", width=5, height=2)
ggsave(filename="median_expression_gtex_thyroid.pdf", plot = gtex_thyroid_plot, path = "./plots/summary_plots", width=5, height=2)
unlink("median_expression_gtex_thyroid.png")
unlink("median_expression_gtex_thyroid.pdf")




# GTEx stomach

gtex_stomach <- read_tsv(file = "./files/gtex_stomach_fpkm.txt")

gtex_stomach <- gtex_stomach %>%
  pivot_longer(-gene, names_to = "sample", values_to = "fpkm") %>%
  mutate(fpkm = log2(fpkm+1)) %>%
  group_by(gene) %>%
  summarise(median_fpkm = median(fpkm), sd_fpkm = sd(fpkm)) %>%
  ungroup() %>%
  inner_join(annotation[, c("gene", "chrom")], by = "gene") %>%
  filter(chrom != "chrM") %>%
  mutate(chrom_type = map_chr(.x = chrom, .f = ~ if(.x %in% c("chrX", "chrY")){.x}else{"autosome"}))

per <- gtex_stomach %>%
  summarise(n = n(), nn = sum(median_fpkm > 1), p = (nn/n)*100)

gtex_stomach_plot <- gtex_stomach %>%
  mutate(gene = fct_reorder(gene, median_fpkm, .fun = function(x) max(x))) %>%
  ggplot() +
  geom_point(mapping = aes(x = gene, y = median_fpkm), size = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.3) +
  annotate(geom = "text", x = 5000, y = 9, label = str_c("median expression > 1 = ", round(per$p,0), "%")) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 10)) +
  scale_color_manual(values = c("grey", "red", "blue")) +
  labs(y = "Median expression (log2 FPKM + 1)", color = "chr", title = "GTEx - Stomach")

ggsave(filename="median_expression_gtex_stomach.png", plot = gtex_stomach_plot, path = "./plots/summary_plots", width=5, height=2)
ggsave(filename="median_expression_gtex_stomach.pdf", plot = gtex_stomach_plot, path = "./plots/summary_plots", width=5, height=2)
unlink("median_expression_gtex_stomach.png")
unlink("median_expression_gtex_stomach.pdf")
