# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data




library(WGCNA)
library(tidyverse)
library(gplots)

options(bitmapType = "cairo")

source("utils.R")




# -- Thyroid

# load R workspace
load("./r_workspaces/tcga_gtex_thyroid_wgcna_networks.RData")



# build correlation matrixes
males_tumour_corMat <- cor(males_tumour)
males_normal_corMat <- cor(males_normal)
females_tumour_corMat <- cor(females_tumour)
females_normal_corMat <- cor(females_normal)


males_tumour_corMat2 <- males_tumour_corMat %>% get_upper_tri %>% as.numeric %>% na.exclude %>% as.numeric

males_normal_corMat2 <- males_normal_corMat %>% get_upper_tri %>% as.numeric %>% na.exclude %>% as.numeric

females_tumour_corMat2 <- females_tumour_corMat %>% get_upper_tri %>% as.numeric %>% na.exclude %>% as.numeric

females_normal_corMat2 <- females_normal_corMat %>% get_upper_tri %>% as.numeric %>% na.exclude %>% as.numeric



thyroid_correlations <- data.frame(
    males_tumor = males_tumour_corMat2,
    males_normal = males_normal_corMat2,
    females_tumor = females_tumour_corMat2,
    females_normal = females_normal_corMat2) %>%
    as.tibble() %>%
    gather(key="data", value="pearson")


thyroid_correlations_plot <- ggplot(data=thyroid_correlations, mapping=aes(x=data, y=pearson, fill=data)) +
    geom_boxplot(outlier.size = 0.5) +
    theme(axis.text.y = element_text(colour="black", size=13),
        axis.title.y = element_text(colour="black", size=15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
    labs(x="", y="Pearson correlation coefficient", fill="tissue / gender")
ggsave("thyroid_correlations_plot.png", plot=thyroid_correlations_plot, path="./plots/")
unlink("thyroid_correlations_plot.png")


thyroid_correlations_densityplot <- ggplot(data=thyroid_correlations, mapping=aes(x=pearson, fill=data, colour=data)) +
    geom_density(alpha = 0.2) +
    geom_vline(xintercept=0, color="black", linetype=2, size = 0.3) +
    theme(axis.text = element_text(colour="black", size=13), axis.title = element_text(colour="black", size=15)) +
    labs(x="Pearson correlation coefficient", y="Density", fill="tissue / gender")
ggsave("thyroid_correlations_densityplot.png", plot=thyroid_correlations_densityplot, path="./plots/")
unlink("thyroid_correlations_densityplot.png")


# distance matrix and hierarchical clustering
males_tumour_corMat_dist <- as.dist(1-males_tumour_corMat)
males_tumour_corMat_hc <- hclust(males_tumour_corMat_dist, method="average")

males_normal_corMat_dist <- as.dist(1-males_normal_corMat)
males_normal_corMat_hc <- hclust(males_normal_corMat_dist, method="average")

females_tumour_corMat_dist <- as.dist(1-females_tumour_corMat)
females_tumour_corMat_hc <- hclust(females_tumour_corMat_dist, method="average")

females_normal_corMat_dist <- as.dist(1-females_normal_corMat)
females_normal_corMat_hc <- hclust(females_normal_corMat_dist, method="average")




png(file="./plots/thyroid_males_tumour_wgcna_modules_heatmap.png")
heatmap.2(x=males_tumour_corMat,
    Rowv=as.dendrogram(males_tumour_corMat_hc),
    Colv="Rowv",
    dendrogram="none",
    symm = T,
    distfun = NULL,
    hclustfun = NULL,
    scale="none",
    trace="none",
    key=TRUE,
    keysize=1.2,
    symkey=TRUE,
    key.title = NA,
    key.xlab = NA,
    key.ylab = NA,
    col=colorpanel(50, "blue", "white", "red"),
    labRow=NA,
    labCol=NA,
    ColSideColors = labels2colors(males_tumour_network$colors),
    RowSideColors = labels2colors(males_tumour_network$colors))
dev.off()


png(file="./plots/thyroid_males_normal_wgcna_modules_heatmap.png")
heatmap.2(x=males_normal_corMat,
    Rowv=as.dendrogram(males_normal_corMat_hc),
    Colv="Rowv",
    dendrogram="none",
    symm = T,
    distfun = NULL,
    hclustfun = NULL,
    scale="none",
    trace="none",
    key=TRUE,
    keysize=1.2,
    symkey=TRUE,
    key.title = NA,
    key.xlab = NA,
    key.ylab = NA,
    col=colorpanel(50, "blue", "white", "red"),
    labRow=NA,
    labCol=NA,
    ColSideColors = labels2colors(males_normal_network$colors),
    RowSideColors = labels2colors(males_normal_network$colors))
dev.off()


png(file="./plots/thyroid_females_tumour_wgcna_modules_heatmap.png")
heatmap.2(x=females_tumour_corMat,
    Rowv=as.dendrogram(females_tumour_corMat_hc),
    Colv="Rowv",
    dendrogram="none",
    symm = T,
    distfun = NULL,
    hclustfun = NULL,
    scale="none",
    trace="none",
    key=TRUE,
    keysize=1.2,
    symkey=TRUE,
    key.title = NA,
    key.xlab = NA,
    key.ylab = NA,
    col=colorpanel(50, "blue", "white", "red"),
    labRow=NA,
    labCol=NA,
    ColSideColors = labels2colors(females_tumour_network$colors),
    RowSideColors = labels2colors(females_tumour_network$colors))
dev.off()


png(file="./plots/thyroid_females_normal_wgcna_modules_heatmap.png")
heatmap.2(x=females_normal_corMat,
    Rowv=as.dendrogram(females_normal_corMat_hc),
    Colv="Rowv",
    dendrogram="none",
    symm = T,
    distfun = NULL,
    hclustfun = NULL,
    scale="none",
    trace="none",
    key=TRUE,
    keysize=1.2,
    symkey=TRUE,
    key.title = NA,
    key.xlab = NA,
    key.ylab = NA,
    col=colorpanel(50, "blue", "white", "red"),
    labRow=NA,
    labCol=NA,
    ColSideColors = labels2colors(females_normal_network$colors),
    RowSideColors = labels2colors(females_normal_network$colors))
dev.off()




