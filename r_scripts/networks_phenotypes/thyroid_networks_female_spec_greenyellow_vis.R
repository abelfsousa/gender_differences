# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Visualization female-specific module
# greenyellow


# -- Thyroid


library(tidyverse)
library(WGCNA)



# load TOM object
load("./r_workspaces/thyroid_females_tumour-block.1.RData")


# load R workspace
load("./r_workspaces/tcga_gtex_thyroid_wgcna_networks.RData")


TOM_females_tumour <- as.matrix(TOM)
dimnames(TOM_females_tumour) <- list(colnames(females_tumour), colnames(females_tumour))
rm(TOM)


# load female tumour modules with kmE
females_tumour_modules <- read_tsv("./files/thyroid_females_tumour_modules.txt")


# select greenyellow module
greenyellow <- females_tumour_modules %>%
  filter(moduleL == "greenyellow") %>%
  mutate(status = if_else(kme > 0.8, "hub", "not_hub"))



# select the corresponding Topological Overlap
greenyellow_tom = TOM_females_tumour[greenyellow$genes_ens, greenyellow$genes_ens]



# export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(greenyellow_tom,
  edgeFile = paste("./files/CytoscapeInput-edges-", paste("greenyellow", collapse="-"), ".txt", sep=""),
  nodeFile = paste("./files/CytoscapeInput-nodes-", paste("greenyellow", collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0,
  nodeNames = greenyellow$geneName,
  altNodeNames = greenyellow$genes_ens,
  nodeAttr = greenyellow$status)
