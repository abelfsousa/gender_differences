# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data




library(WGCNA)
library(tidyverse)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(3)




# -- Stomach

# load datasets
stomach_tcga_gtex_fpkm <- read.delim("./files/stomach_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_gtex_meta <- read.delim("./files/stomach_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)
rownames(stomach_tcga_gtex_meta) <- gsub("-", ".", rownames(stomach_tcga_gtex_meta))


# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))


# reorder datasets
stomach_tcga_gtex_fpkm <- stomach_tcga_gtex_fpkm[, order(colnames(stomach_tcga_gtex_fpkm))]
stomach_tcga_gtex_meta <- stomach_tcga_gtex_meta[order(rownames(stomach_tcga_gtex_meta)), ]


#remove normal TCGA samples
stomach_tcga_gtex_fpkm <- stomach_tcga_gtex_fpkm[, !colnames(stomach_tcga_gtex_fpkm) %in% rownames(stomach_tcga_gtex_meta[(stomach_tcga_gtex_meta$sample_type == "TCGA_normal"), ])]
stomach_tcga_gtex_meta <- stomach_tcga_gtex_meta[!(stomach_tcga_gtex_meta$sample_type == "TCGA_normal"), ]



stomach_tcga_gtex_fpkm %>% colnames %>% all.equal(rownames(stomach_tcga_gtex_meta))
#TRUE





# split data
males_tumour <- stomach_tcga_gtex_fpkm[, colnames(stomach_tcga_gtex_fpkm) %in% rownames(stomach_tcga_gtex_meta[(stomach_tcga_gtex_meta$gender == "male" & stomach_tcga_gtex_meta$sample_type == "TCGA_tumour"), ])]
males_normal <- stomach_tcga_gtex_fpkm[, colnames(stomach_tcga_gtex_fpkm) %in% rownames(stomach_tcga_gtex_meta[(stomach_tcga_gtex_meta$gender == "male" & stomach_tcga_gtex_meta$sample_type == "GTEx"), ])]


females_tumour <- stomach_tcga_gtex_fpkm[, colnames(stomach_tcga_gtex_fpkm) %in% rownames(stomach_tcga_gtex_meta[(stomach_tcga_gtex_meta$gender == "female" & stomach_tcga_gtex_meta$sample_type == "TCGA_tumour"), ])]
females_normal <- stomach_tcga_gtex_fpkm[, colnames(stomach_tcga_gtex_fpkm) %in% rownames(stomach_tcga_gtex_meta[(stomach_tcga_gtex_meta$gender == "female" & stomach_tcga_gtex_meta$sample_type == "GTEx"), ])]


#log-transform and transpose
males_tumour <- t(log2(males_tumour + 1))
males_normal <- t(log2(males_normal + 1))

females_tumour <- t(log2(females_tumour + 1))
females_normal <- t(log2(females_normal + 1))




#removing samples and genes with too many missing values
gsg_males_tumour <- goodSamplesGenes(males_tumour, verbose = 0)
gsg_males_tumour$allOK
#TRUE

gsg_males_normal <- goodSamplesGenes(males_normal, verbose = 0)
gsg_males_normal$allOK
#TRUE

gsg_females_tumour <- goodSamplesGenes(females_tumour, verbose = 0)
gsg_females_tumour$allOK
#TRUE

gsg_females_normal <- goodSamplesGenes(females_normal, verbose = 0)
gsg_females_normal$allOK
#TRUE



#cluster samples to remove outliers
pdf(file = "./plots/wgcna_networks/stomach_hclust_tcga_gtex_samples_1.pdf")
plot( hclust(dist(males_tumour), method = "average"), main = "Sample clustering to detect outliers", sub="Males TCGA tumours", labels=F, xlab="" )
plot( hclust(dist(males_normal), method = "average"), main = "Sample clustering to detect outliers", sub="Males GTEx", labels=F, xlab="" )
plot( hclust(dist(females_tumour), method = "average"), main = "Sample clustering to detect outliers", sub="Females TCGA tumours", labels=F, xlab="" )
plot( hclust(dist(females_normal), method = "average"), main = "Sample clustering to detect outliers", sub="Females GTEx", labels=F, xlab="" )
dev.off()


#detect samples above certain distance
cutTree_males_tumour <- cutreeStatic(hclust(dist(males_tumour), method = "average"), cutHeight = 123, minSize = 10)
cutTree_males_tumour <- (cutTree_males_tumour == 1)
cutTree_females_tumour <- cutreeStatic(hclust(dist(females_tumour), method = "average"), cutHeight = 117, minSize = 10)
cutTree_females_tumour <- (cutTree_females_tumour == 1)


#remove samples
males_tumour <- males_tumour[cutTree_males_tumour, ]
females_tumour <- females_tumour[cutTree_females_tumour, ]


pdf(file = "./plots/wgcna_networks/stomach_hclust_tcga_gtex_samples_2.pdf")
plot( hclust(dist(males_tumour), method = "average"), main = "Sample clustering to detect outliers", sub="Males TCGA tumours", labels=F, xlab="" )
plot( hclust(dist(males_normal), method = "average"), main = "Sample clustering to detect outliers", sub="Males GTEx", labels=F, xlab="" )
plot( hclust(dist(females_tumour), method = "average"), main = "Sample clustering to detect outliers", sub="Females TCGA tumours", labels=F, xlab="" )
plot( hclust(dist(females_normal), method = "average"), main = "Sample clustering to detect outliers", sub="Females GTEx", labels=F, xlab="" )
dev.off()






# choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))

sft_males_tumour <- pickSoftThreshold(males_tumour, powerVector = powers, verbose = 5)
sft_males_normal <- pickSoftThreshold(males_normal, powerVector = powers, verbose = 5)
sft_females_tumour <- pickSoftThreshold(females_tumour, powerVector = powers, verbose = 5)
sft_females_normal <- pickSoftThreshold(females_normal, powerVector = powers, verbose = 5)




# networks construction
males_tumour_network <- blockwiseModules(males_tumour, maxBlockSize = 20000, power = sft_males_tumour$powerEstimate, TOMType = "unsigned", minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "./r_workspaces/stomach_males_tumour", verbose = 3)
males_normal_network <- blockwiseModules(males_normal, maxBlockSize = 20000, power = sft_males_normal$powerEstimate, TOMType = "unsigned", minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "./r_workspaces/stomach_males_normal", verbose = 3)
females_tumour_network <- blockwiseModules(females_tumour, maxBlockSize = 20000, power = sft_females_tumour$powerEstimate, TOMType = "unsigned", minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "./r_workspaces/stomach_females_tumour", verbose = 3)
females_normal_network <- blockwiseModules(females_normal, maxBlockSize = 20000, power = sft_females_normal$powerEstimate, TOMType = "unsigned", minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "./r_workspaces/stomach_females_normal", verbose = 3)





# plot network dendrogram and the module colors underneath

net <- males_tumour_network
mergedColors <- labels2colors(males_tumour_network$colors)
pdf("./plots/wgcna_networks/stomach_males_tumour_wgcna_modules_dendrogram.pdf", w=12, h=10)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], groupLabels = "Modules", main= "Stomach Tumour Tissue (TCGA) - Men", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


net <- males_normal_network
mergedColors <- labels2colors(males_normal_network$colors)
pdf("./plots/wgcna_networks/stomach_males_normal_wgcna_modules_dendrogram.pdf", w=12, h=10)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], groupLabels = "Modules", main= "Stomach Normal Tissue (GTEx) - Men", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


net <- females_tumour_network
mergedColors <- labels2colors(females_tumour_network$colors)
pdf("./plots/wgcna_networks/stomach_females_tumour_wgcna_modules_dendrogram.pdf", w=12, h=10)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], groupLabels = "Modules", main= "Stomach Tumour Tissue (TCGA) - Women", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


net <- females_normal_network
mergedColors <- labels2colors(females_normal_network$colors)
pdf("./plots/wgcna_networks/stomach_females_normal_wgcna_modules_dendrogram.pdf", w=12, h=10)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], groupLabels = "Modules", main= "Stomach Normal Tissue (GTEx) - Women", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()



write.table(males_tumour, "./files/stomach_males_tumour_fpkm_wgcna.txt", quote = F, sep="\t")
write.table(males_normal, "./files/stomach_males_normal_fpkm_wgcna.txt", quote = F, sep="\t")
write.table(females_tumour, "./files/stomach_females_tumour_fpkm_wgcna.txt", quote = F, sep="\t")
write.table(females_normal, "./files/stomach_females_normal_fpkm_wgcna.txt", quote = F, sep="\t")


save(list=ls(), file="./r_workspaces/tcga_gtex_stomach_wgcna_networks.RData")
