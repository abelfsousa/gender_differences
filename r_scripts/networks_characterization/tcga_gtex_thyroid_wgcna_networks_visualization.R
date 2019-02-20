# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data

# TOM plot




library(WGCNA)
library(tidyverse)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(3)



# -- Thyroid

# load R workspace
load("./r_workspaces/tcga_gtex_thyroid_wgcna_networks.RData")




# Males


## tumour
colors1 <- labels2colors(males_tumour_network$colors)

restGenes1 <- (colors1 != "grey")
diss_mat1 <- 1 - TOMsimilarityFromExpr( males_tumour[, restGenes1], power = sft_males_tumour$powerEstimate )
hier_clust1 <- hclust( as.dist(diss_mat1), method="average" )
diag(diss_mat1) = NA



png(file="./plots/wgcna_networks/thyroid_males_tumour_wgcna_modules_tomplot.png")
TOMplot( diss_mat1^10, hier_clust1, as.character(colors1[restGenes1]), main = "" )
dev.off()



## normal
colors2 <- labels2colors(males_normal_network$colors)

restGenes2 <- (colors2 != "grey")
diss_mat2 <- 1 - TOMsimilarityFromExpr( males_normal[, restGenes2], power = sft_males_normal$powerEstimate )
hier_clust2 <- hclust( as.dist(diss_mat2), method="average" )
diag(diss_mat2) = NA



png(file="./plots/wgcna_networks/thyroid_males_normal_wgcna_modules_tomplot.png")
TOMplot( diss_mat2^10, hier_clust2, as.character(colors2[restGenes2]), main = "" )
dev.off()








# Females


## tumour
colors3 <- labels2colors(females_tumour_network$colors)

restGenes3 <- (colors3 != "grey")
diss_mat3 <- 1 - TOMsimilarityFromExpr( females_tumour[, restGenes3], power = sft_females_tumour$powerEstimate )
hier_clust3 <- hclust( as.dist(diss_mat3), method="average" )
diag(diss_mat3) = NA



png(file="./plots/wgcna_networks/thyroid_females_tumour_wgcna_modules_tomplot.png")
TOMplot( diss_mat3^10, hier_clust3, as.character(colors3[restGenes3]), main = "" )
dev.off()



## normal
colors4 <- labels2colors(females_normal_network$colors)

restGenes4 <- (colors4 != "grey")
diss_mat4 <- 1 - TOMsimilarityFromExpr( females_normal[, restGenes4], power = sft_females_normal$powerEstimate )
hier_clust4 <- hclust( as.dist(diss_mat4), method="average" )
diag(diss_mat4) = NA



png(file="./plots/wgcna_networks/thyroid_females_normal_wgcna_modules_tomplot.png")
TOMplot( diss_mat4^10, hier_clust4, as.character(colors4[restGenes4]), main = "" )
dev.off()
