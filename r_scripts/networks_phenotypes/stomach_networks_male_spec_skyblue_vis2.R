# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Visualization male-specific module
# skyblue


# -- Stomach


library(tidyverse)
library(WGCNA)



# load R workspace
load("./r_workspaces/tcga_gtex_stomach_wgcna_networks.RData")


# load male tumour modules
males_tumour_modules <- read_tsv("./files/stomach_males_tumour_modules.txt")


# select skyblue module
skyblue <- males_tumour_modules %>%
  filter(moduleL == "skyblue")


# calculate mean expression for females and males
females_tumour_mean <- females_tumour %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as.tibble() %>%
  gather(key = "sample", value = "fpkm", -gene) %>%
  filter(gene %in% skyblue$genes_ens) %>%
  group_by(gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  inner_join(skyblue[, c("genes_ens", "geneName")], by=c("gene" = "genes_ens")) %>%
  dplyr::select(geneName, fpkm)
write.table(females_tumour_mean, "./files/stomach_skyblue_females_mean_expression.txt", col.names=F, sep = "\t", row.names=F, quote=F)

males_tumour_mean <- males_tumour %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as.tibble() %>%
  gather(key = "sample", value = "fpkm", -gene) %>%
  filter(gene %in% skyblue$genes_ens) %>%
  group_by(gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  inner_join(skyblue[, c("genes_ens", "geneName")], by=c("gene" = "genes_ens")) %>%
  dplyr::select(geneName, fpkm)
write.table(males_tumour_mean, "./files/stomach_skyblue_males_mean_expression.txt", col.names=F, sep = "\t", row.names=F, quote=F)



# calculate gene-gene correlations
females_tumour_cor <- females_tumour[, skyblue$genes_ens] %>% cor() %>% abs()
males_tumour_cor <- males_tumour[, skyblue$genes_ens] %>% cor() %>% abs()



# export the network into edge and node list files Cytoscape can read
cyt_females = exportNetworkToCytoscape(females_tumour_cor,
  edgeFile = "files/stomach_skyblue_module_females_edges.txt",
  nodeFile = "files/stomach_skyblue_module_females_nodes.txt",
  weighted = TRUE,
  threshold = 0,
  nodeNames = skyblue$genes_ens,
  altNodeNames = skyblue$geneName)

cyt_males = exportNetworkToCytoscape(males_tumour_cor,
  edgeFile = "files/stomach_skyblue_module_males_edges.txt",
  nodeFile = "files/stomach_skyblue_module_males_nodes.txt",
  weighted = TRUE,
  threshold = 0,
  nodeNames = skyblue$genes_ens,
  altNodeNames = skyblue$geneName)


cyt_males2 <- cyt_males$edgeData %>%
  mutate(pair = paste(fromAltName, toAltName, sep="_")) %>%
  filter(weight > 0.7)
write.table(cyt_males2, "./files/stomach_skyblue_module_males_edges.txt", col.names=T, sep = "\t", row.names=F, quote=F)

males_genes <- unique(c(as.character(cyt_males2$fromAltName), as.character(cyt_males2$toAltName)))

cyt_females2 <- cyt_females$edgeData %>%
  mutate(pair = paste(fromAltName, toAltName, sep="_")) %>%
  filter(pair %in% cyt_males2$pair)
  #filter(fromAltName %in% males_genes | toAltName %in% males_genes)
write.table(cyt_females2, "./files/stomach_skyblue_module_females_edges.txt", col.names=T, sep = "\t", row.names=F, quote=F)



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
  geom_density(alpha = 0.5) +
  geom_vline(data = comparing_distribution_stats, mapping=aes(xintercept=median, group=gender, color=gender), linetype="dashed", size=0.8) +
  annotate("text", x = 0.15, y = 2, label = "Wilcoxon,\np < 2.2e-16") +
  theme_classic() +
  scale_fill_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  scale_color_manual(values=c("#fbb4b9", "#74a9cf"), name="Gender", labels = c("Female", "Male")) +
  theme(
    axis.title.y = element_text(colour="black", size=15),
    axis.title.x = element_text(colour="black", size=15),
    axis.text.x = element_text(colour="black", size=12),
    axis.text.y = element_text(colour="black", size=12),
    legend.text = element_text(colour="black", size=12),
    legend.title = element_text(colour="black", size=15)) +
  labs(x = "Pearson's r", y = "Density")
ggsave(filename="stomach_skyblue_module_cor_dist.png", plot = comparing_distribution_boxpl, path = "./plots/wgcna_networks_traits", width=6, height=4)
unlink("stomach_skyblue_module_cor_dist.png")
