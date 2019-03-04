# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data


# Enrichment analysis of female-specific module


# -- Stomach


library(tidyverse)
library(WGCNA)



# load TOM object
load("./r_workspaces/stomach_males_tumour-block.1.RData")


# load R workspace
load("./r_workspaces/tcga_gtex_stomach_wgcna_networks.RData")


TOM_males_tumour <- as.matrix(TOM)
dimnames(TOM_males_tumour) <- list(colnames(males_tumour), colnames(males_tumour))
rm(TOM)


# load male tumour modules with kmE
males_tumour_modules <- read_tsv("./files/stomach_males_tumour_modules.txt")


# select skyblue module
skyblue <- males_tumour_modules %>%
  filter(moduleL == "skyblue") %>%
  mutate(status = if_else(kme > 0.8, "hub", "not_hub"))



# select the corresponding Topological Overlap
skyblue_tom = TOM_males_tumour[skyblue$genes_ens, skyblue$genes_ens]



# export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(skyblue_tom,
  edgeFile = paste("./files/CytoscapeInput-edges-", paste("skyblue", collapse="-"), ".txt", sep=""),
  nodeFile = paste("./files/CytoscapeInput-nodes-", paste("skyblue", collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0,
  nodeNames = skyblue$geneName,
  altNodeNames = skyblue$genes_ens,
  nodeAttr = skyblue$status)

hub_genes <- skyblue %>% filter(status == "hub") %>% pull(geneName)

cyt[[1]] <- cyt[[1]] %>%
  mutate(edgeAtt = if_else((fromNode %in% hub_genes) & (toNode %in% hub_genes), "hub", "not_hub"))
write.table(cyt[[1]], "./files/CytoscapeInput-edges-skyblue2.txt", quote = F, row.names=F, sep="\t")
