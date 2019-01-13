# Understanding Gender Differential Susceptibility in Cancer


# Correlation networks
# TCGA + GTEx data


library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(limma)

library(gplots)
library(dynamicTreeCut)
library(moduleColor)

source("utils.R")





# -- Thyroid

# load datasets
thyroid_tcga_gtex_counts <- read.delim("./files/thyroid_tcga_gtex_counts.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_gtex_meta <- read.delim("./files/thyroid_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)
rownames(thyroid_tcga_gtex_meta) <- gsub("-", ".", rownames(thyroid_tcga_gtex_meta))


# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))


# reorder datasets
thyroid_tcga_gtex_counts <- thyroid_tcga_gtex_counts[, order(colnames(thyroid_tcga_gtex_counts))]
thyroid_tcga_gtex_meta <- thyroid_tcga_gtex_meta[order(rownames(thyroid_tcga_gtex_meta)), ]


# remove metastatic samples
thyroid_tcga_gtex_counts <- thyroid_tcga_gtex_counts[, !colnames(thyroid_tcga_gtex_counts) %in% rownames(thyroid_tcga_gtex_meta[(thyroid_tcga_gtex_meta$sample_type == "TCGA_metastatic"), ])]
thyroid_tcga_gtex_meta <- thyroid_tcga_gtex_meta[!(thyroid_tcga_gtex_meta$sample_type == "TCGA_metastatic"), ]


#remove normal TCGA samples
thyroid_tcga_gtex_counts <- thyroid_tcga_gtex_counts[, !colnames(thyroid_tcga_gtex_counts) %in% rownames(thyroid_tcga_gtex_meta[(thyroid_tcga_gtex_meta$sample_type == "TCGA_normal"), ])]
thyroid_tcga_gtex_meta <- thyroid_tcga_gtex_meta[!(thyroid_tcga_gtex_meta$sample_type == "TCGA_normal"), ]



thyroid_tcga_gtex_counts %>% colnames %>% all.equal(rownames(thyroid_tcga_gtex_meta))
#TRUE





# -- Normalize counts


# define design matrix
design <- model.matrix(~ gender + ethnicity + age + sample_type, data = thyroid_tcga_gtex_meta)



# DGEList object using edgeR
thyroid_tcga_gtex_dgelist <- DGEList(counts=thyroid_tcga_gtex_counts)

# scale normalization using TMM method
thyroid_tcga_gtex_dgelist <- calcNormFactors(thyroid_tcga_gtex_dgelist)

# transform the counts using voom
thyroid_tcga_gtex_dgelist <- voom(thyroid_tcga_gtex_dgelist, design = design)$E



# pca betweewn GTEx and TCGA samples
pca1 <- prcomp(t(thyroid_tcga_gtex_dgelist), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca1 <- pca1 %>%
    mutate(data = as.character(thyroid_tcga_gtex_meta$sample_type)) %>%
    mutate(gender = as.character(thyroid_tcga_gtex_meta$gender))


#scatterplot
ggplot( data=pca1, mapping=aes(x=PC1, y=PC2, color=data, shape=gender) ) +
    geom_point() +
    scale_shape_manual(values=c(19, 15)) +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15))
ggsave(filename="pca1_thyroid_tcga_gtex_samples_voom.pdf", path = "./plots/")
unlink("pca1_thyroid_tcga_gtex_samples_voom.pdf")




# regress-out batch effect between GTEx and TCGA samples
# use sample_type (TCGA_tumour or GTEx) as covariate 
thyroid_tcga_gtex_dgelist_cor <- sapply( rownames(thyroid_tcga_gtex_dgelist), 
	function(x) remove_batch(x, thyroid_tcga_gtex_dgelist, thyroid_tcga_gtex_meta$sample_type) )

thyroid_tcga_gtex_dgelist_cor <- thyroid_tcga_gtex_dgelist_cor %>% t() %>% as.data.frame()



# pca betweewn GTEx and TCGA samples
pca2 <- prcomp(t(thyroid_tcga_gtex_dgelist_cor), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca2 <- pca2 %>%
    mutate(data = as.character(thyroid_tcga_gtex_meta$sample_type)) %>%
    mutate(gender = as.character(thyroid_tcga_gtex_meta$gender))


#scatterplot
ggplot( data=pca2, mapping=aes(x=PC1, y=PC2, color=data, shape=gender) ) +
    geom_point() +
    scale_shape_manual(values=c(19, 15)) +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15))
ggsave(filename="pca2_thyroid_tcga_gtex_samples_voom.pdf", path = "./plots/")
unlink("pca2_thyroid_tcga_gtex_samples_voom.pdf")






save(list=ls(), file="r_workspaces/tcga_gtex_thyroid_geneCorrelation.RData")





# -- Males tumour

males_tumour <- thyroid_tcga_gtex_dgelist_cor[, colnames(thyroid_tcga_gtex_dgelist_cor) %in% rownames(thyroid_tcga_gtex_meta[(thyroid_tcga_gtex_meta$gender == "male" & thyroid_tcga_gtex_meta$sample_type == "TCGA_tumour"), ])]


# build correlation matrix
males_tumour_corMat <- cor(t(males_tumour))


# select genes with pearson r higher than 0.1 in 40% of the gene correlations
#males_tumour_genes <- rowSums(abs(males_tumour_corMat) > 0.1) > ncol(males_tumour_corMat)*0.7
#males_tumour_corMat <- males_tumour_corMat[males_tumour_genes, males_tumour_genes]


# distance matrix and hierarchical clustering
males_tumour_corMat_dist <- as.dist(1-males_tumour_corMat)
males_tumour_corMat_hc <- hclust(males_tumour_corMat_dist, method="average")


# detect clusters
males_tumour_clusters1 <- cutreeDynamic(males_tumour_corMat_hc, method = "hybrid", distM = as.matrix(males_tumour_corMat_dist), minClusterSize = 50)
males_tumour_clusters2 <- cutreeDynamic(males_tumour_corMat_hc, method = "tree", distM = NULL, minClusterSize = 50)





png(file="./plots/thyroid_males_tumour_heatmap.png")
heatmap.2(x=males_tumour_corMat,
    Rowv=as.dendrogram(males_tumour_corMat_hc),
    Colv="Rowv",
    dendrogram="both",
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
    ColSideColors = labels2colors(males_tumour_clusters2),
    RowSideColors = labels2colors(males_tumour_clusters2))
dev.off()






# -- Males normal

males_normal <- thyroid_tcga_gtex_dgelist_cor[, colnames(thyroid_tcga_gtex_dgelist_cor) %in% rownames(thyroid_tcga_gtex_meta[(thyroid_tcga_gtex_meta$gender == "male" & thyroid_tcga_gtex_meta$sample_type == "GTEx"), ])]


# build correlation matrix
males_normal_corMat <- cor(t(males_normal))


# select genes with pearson r higher than 0.1 in 40% of the gene correlations
#males_normal_genes <- rowSums(abs(males_normal_corMat) > 0.1) > ncol(males_normal_corMat)*0.7
#males_normal_corMat <- males_normal_corMat[males_normal_genes, males_normal_genes]


# distance matrix and hierarchical clustering
males_normal_corMat_dist <- as.dist(1-males_normal_corMat)
males_normal_corMat_hc <- hclust(males_normal_corMat_dist, method="average")


# detect clusters
males_normal_clusters1 <- cutreeDynamic(males_normal_corMat_hc, method = "hybrid", distM = as.matrix(males_normal_corMat_dist), minClusterSize = 50)
males_normal_clusters2 <- cutreeDynamic(males_normal_corMat_hc, method = "tree", distM = NULL, minClusterSize = 50)



png(file="./plots/thyroid_males_normal_heatmap.png")
heatmap.2(x=males_normal_corMat,
    Rowv=as.dendrogram(males_normal_corMat_hc),
    Colv="Rowv",
    dendrogram="both",
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
    ColSideColors = labels2colors(males_normal_clusters2),
    RowSideColors = labels2colors(males_normal_clusters2))
dev.off()






# -- Females tumour

females_tumour <- thyroid_tcga_gtex_dgelist_cor[, colnames(thyroid_tcga_gtex_dgelist_cor) %in% rownames(thyroid_tcga_gtex_meta[(thyroid_tcga_gtex_meta$gender == "female" & thyroid_tcga_gtex_meta$sample_type == "TCGA_tumour"), ])]


# build correlation matrix
females_tumour_corMat <- cor(t(females_tumour))


# select genes with pearson r higher than 0.1 in 40% of the gene correlations
#males_tumour_genes <- rowSums(abs(males_tumour_corMat) > 0.1) > ncol(males_tumour_corMat)*0.7
#males_tumour_corMat <- males_tumour_corMat[males_tumour_genes, males_tumour_genes]


# distance matrix and hierarchical clustering
females_tumour_corMat_dist <- as.dist(1-females_tumour_corMat)
females_tumour_corMat_hc <- hclust(females_tumour_corMat_dist, method="average")


png(file="./plots/thyroid_females_tumour_heatmap.png")
heatmap.2(x=females_tumour_corMat,
    Rowv=as.dendrogram(females_tumour_corMat_hc),
    Colv="Rowv",
    dendrogram="both",
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
    labCol=NA)
dev.off()






# -- Females normal

females_normal <- thyroid_tcga_gtex_dgelist_cor[, colnames(thyroid_tcga_gtex_dgelist_cor) %in% rownames(thyroid_tcga_gtex_meta[(thyroid_tcga_gtex_meta$gender == "female" & thyroid_tcga_gtex_meta$sample_type == "GTEx"), ])]


# build correlation matrix
females_normal_corMat <- cor(t(females_normal))


# select genes with pearson r higher than 0.1 in 40% of the gene correlations
#males_normal_genes <- rowSums(abs(males_normal_corMat) > 0.1) > ncol(males_normal_corMat)*0.7
#males_normal_corMat <- males_normal_corMat[males_normal_genes, males_normal_genes]


# distance matrix and hierarchical clustering
females_normal_corMat_dist <- as.dist(1-females_normal_corMat)
females_normal_corMat_hc <- hclust(females_normal_corMat_dist, method="average")


png(file="./plots/thyroid_females_normal_heatmap.png")
heatmap.2(x=females_normal_corMat,
    Rowv=as.dendrogram(females_normal_corMat_hc),
    Colv="Rowv",
    dendrogram="both",
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
    labCol=NA)
dev.off()






# select upper triangle
#males_tumour_corMat <- get_upper_tri(males_tumour_corMat)



# gather dataset
#males_tumour_corMat <- males_tumour_corMat %>%
#    mutate(gene = rownames(males_tumour_corMat)) %>%
#    as.tibble() %>%
#    gather(key="gene_pair", value="pearson_r", -gene)



# heatmap
#ggplot(males_tumour_corMat, aes(x = gene, y = gene_pair, fill = pearson_r)) +
#  geom_tile() +
#  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson r") +
#  theme_classic() +
#  coord_fixed() +
#  theme(axis.text = element_blank(),
#        axis.ticks = element_blank(),
#        axis.title = element_text(colour="black", size=12)) +
#  labs(x = "Genes", y="Genes") +
#ggsave(filename="thyroid_males_tumour_heatmap.png", path = "./plots/")
#unlink("thyroid_males_tumour_heatmap.png")





