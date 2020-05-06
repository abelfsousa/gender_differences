# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Visualization female-specific module
# greenyellow


# -- Thyroid


library(tidyverse)
library(WGCNA)



# load R workspace
load("./r_workspaces/tcga_gtex_thyroid_wgcna_networks.RData")


# load female tumour modules
females_tumour_modules <- read_tsv("./files/thyroid_females_tumour_modules.txt")


# select greenyellow module
greenyellow <- females_tumour_modules %>%
  filter(moduleL == "greenyellow")


# calculate mean expression for females and males
females_tumour_mean <- females_tumour %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as.tibble() %>%
  gather(key = "sample", value = "fpkm", -gene) %>%
  filter(gene %in% greenyellow$genes_ens) %>%
  group_by(gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  inner_join(greenyellow[, c("genes_ens", "geneName")], by=c("gene" = "genes_ens")) %>%
  dplyr::select(geneName, fpkm)
write.table(females_tumour_mean, "./files/thyroid_greenyellow_females_mean_expression.txt", col.names=F, sep = "\t", row.names=F, quote=F)

males_tumour_mean <- males_tumour %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as.tibble() %>%
  gather(key = "sample", value = "fpkm", -gene) %>%
  filter(gene %in% greenyellow$genes_ens) %>%
  group_by(gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  inner_join(greenyellow[, c("genes_ens", "geneName")], by=c("gene" = "genes_ens")) %>%
  dplyr::select(geneName, fpkm)
write.table(males_tumour_mean, "./files/thyroid_greenyellow_males_mean_expression.txt", col.names=F, sep = "\t", row.names=F, quote=F)



# calculate gene-gene correlations
females_tumour_cor <- females_tumour[, greenyellow$genes_ens] %>% cor() %>% abs()
males_tumour_cor <- males_tumour[, greenyellow$genes_ens] %>% cor() %>% abs()



# export the network into edge and node list files Cytoscape can read
cyt_females = exportNetworkToCytoscape(females_tumour_cor,
  edgeFile = "files/thyroid_greenyellow_module_females_edges.txt",
  nodeFile = "files/thyroid_greenyellow_module_females_nodes.txt",
  weighted = TRUE,
  threshold = 0,
  nodeNames = greenyellow$genes_ens,
  altNodeNames = greenyellow$geneName)

cyt_males = exportNetworkToCytoscape(males_tumour_cor,
  edgeFile = "files/thyroid_greenyellow_module_males_edges.txt",
  nodeFile = "files/thyroid_greenyellow_module_males_nodes.txt",
  weighted = TRUE,
  threshold = 0,
  nodeNames = greenyellow$genes_ens,
  altNodeNames = greenyellow$geneName)


cyt_females2 <- cyt_females$edgeData %>%
  mutate(pair = paste(fromAltName, toAltName, sep="_")) %>%
  filter(weight > 0.7)
write.table(cyt_females2, "./files/thyroid_greenyellow_module_females_edges.txt", col.names=T, sep = "\t", row.names=F, quote=F)

females_genes <- unique(c(as.character(cyt_females2$fromAltName), as.character(cyt_females2$toAltName)))

cyt_males2 <- cyt_males$edgeData %>%
  mutate(pair = paste(fromAltName, toAltName, sep="_")) %>%
  filter(pair %in% cyt_females2$pair)
  #filter(fromAltName %in% females_genes | toAltName %in% females_genes)
write.table(cyt_males2, "./files/thyroid_greenyellow_module_males_edges.txt", col.names=T, sep = "\t", row.names=F, quote=F)



# compare distribution of correlations between genders

wilcox.test(cyt_males$edgeData$weight, cyt_females$edgeData$weight)

comparing_distribution <- bind_rows(
  cyt_females$edgeData %>% dplyr::select(fromAltName, toAltName, weight) %>% mutate(gender = "females"),
  cyt_males$edgeData %>% dplyr::select(fromAltName, toAltName, weight) %>% mutate(gender = "males"),
)

comparing_distribution_stats <- comparing_distribution %>%
  group_by(gender) %>%
  summarise(median = median(weight))

comparing_distribution_boxpl <- ggplot(data = comparing_distribution, mapping = aes(x = weight, fill = gender, color = gender)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = comparing_distribution_stats, mapping=aes(xintercept=median, group=gender, color=gender), linetype="dashed", size=0.8, show.legend=F) +
  annotate("text", x = 0.85, y = 2.6, label = "Wilcoxon,\nP-value < 2.2e-16", size = 6.5) +
  theme_classic() +
  scale_fill_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  scale_color_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  theme(
    axis.title = element_text(colour="black", size=20),
    axis.text = element_text(colour="black", size=18),
    legend.text = element_text(colour="black", size=18),
    legend.title = element_text(colour="black", size=20),
    legend.position = "bottom") +
  labs(x = "Pearson's r (absolute pairwise gene correlation)", y = "Density") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 3)) +
  guides(color = guide_legend(nrow=1))
ggsave(filename="thyroid_greenyellow_module_cor_dist.png", plot = comparing_distribution_boxpl, path = "./plots/wgcna_networks_traits", width=7, height=4)
ggsave(filename="thyroid_greenyellow_module_cor_dist.pdf", plot = comparing_distribution_boxpl, path = "./plots/wgcna_networks_traits", width=7, height=4)
unlink("thyroid_greenyellow_module_cor_dist.png")
unlink("thyroid_greenyellow_module_cor_dist.pdf")
