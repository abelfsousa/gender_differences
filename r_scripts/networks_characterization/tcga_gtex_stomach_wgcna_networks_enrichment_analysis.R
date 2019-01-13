# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data




library(WGCNA)
library(tidyverse)
library(gplots)
library(clusterProfiler)
library(org.Hs.eg.db)




# -- Stomach

# load R workspace
load("./r_workspaces/tcga_gtex_stomach_wgcna_networks.RData")



all_genes <- colnames(males_tumour)



enrichment_2modules_go <- function(all_genes, network1, module1, network2, module2){
    
    module_genes1 <- all_genes[which(labels2colors(network1$colors) == module1)]
    module_genes2 <- all_genes[which(labels2colors(network2$colors) == module2)]

    genes_modules1_2 <- intersect(module_genes1, module_genes2)


    go_enrichment <- enrichGO(
        gene = genes_modules1_2,
        universe = all_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500,
        pool=TRUE)


    return(go_enrichment)

}



enrichment_1module_go <- function(all_genes, network, module){
    
    module_genes <- all_genes[which(labels2colors(network$colors) == module)]


    go_enrichment <- enrichGO(
        gene = module_genes,
        universe = all_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500,
        pool=TRUE)

    return(go_enrichment)

}




enrichment_2modules_kegg <- function(all_genes, network1, module1, network2, module2){
    
    module_genes1 <- all_genes[which(labels2colors(network1$colors) == module1)]
    module_genes2 <- all_genes[which(labels2colors(network2$colors) == module2)]

    genes_modules1_2 <- intersect(module_genes1, module_genes2)


    kegg_enrichment <- enrichKEGG(
        gene = bitr(genes_modules1_2, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
        organism = "hsa",
        keyType = "ncbi-geneid",
        universe = bitr(all_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500,
        use_internal_data = FALSE)

    return(kegg_enrichment)

}




enrichment_1module_kegg <- function(all_genes, network, module){
    
    module_genes <- all_genes[which(labels2colors(network$colors) == module)]


    kegg_enrichment <- enrichKEGG(
        gene = bitr(module_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
        organism = "hsa",
        keyType = "ncbi-geneid",
        universe = bitr(all_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500,
        use_internal_data = FALSE)

    return(kegg_enrichment)

}





###############################
#Males tumour vs females tumour
###############################



# turquoise module is conserved in both networks
turquoise_module_females_males_tumourGO <- enrichment_2modules_go(all_genes, males_tumour_network, "turquoise", females_tumour_network, "turquoise")
turquoise_module_females_males_tumourkeggPath <- enrichment_2modules_kegg(all_genes, males_tumour_network, "turquoise", females_tumour_network, "turquoise")

turquoise_module_females_males_tumourGO_dot <- dotplot(turquoise_module_females_males_tumourGO, x = "GeneRatio", color = "p.adjust", showCategory = 10, split = NULL, font.size = 20, title = "Turquoise module\nGO terms: BP")
ggsave(filename="stomach_tcga_turquoise_module_females_males_tumourGO_dot.png", plot=turquoise_module_females_males_tumourGO_dot, path = "./plots/wgcna_modules_enrichment/", width=8, height=10)
unlink("stomach_tcga_turquoise_module_females_males_tumourGO_dot.png")


turquoise_module_females_males_tumourkeggPath_dot <- dotplot(turquoise_module_females_males_tumourkeggPath, x = "GeneRatio", color = "p.adjust", showCategory = 10, font.size = 20, title = "Turquoise module\nKEGG pathways")
ggsave(filename="stomach_tcga_turquoise_module_females_males_tumourkeggPath_dot.png", plot=turquoise_module_females_males_tumourkeggPath_dot, path = "./plots/wgcna_modules_enrichment/", width=12, height=8)
unlink("stomach_tcga_turquoise_module_females_males_tumourkeggPath_dot.png")



# green module is conserved in both networks
green_module_females_males_tumourGO <- enrichment_2modules_go(all_genes, males_tumour_network, "green", females_tumour_network, "green")
green_module_females_males_tumourkeggPath <- enrichment_2modules_kegg(all_genes, males_tumour_network, "green", females_tumour_network, "green")

green_module_females_males_tumourGO_dot <- dotplot(green_module_females_males_tumourGO, x = "GeneRatio", color = "p.adjust", showCategory = 10, split = NULL, font.size = 20, title = "Green module\nGO terms: BP")
ggsave(filename="stomach_tcga_green_module_females_males_tumourGO_dot.png", plot=green_module_females_males_tumourGO_dot, path = "./plots/wgcna_modules_enrichment/", width=12, height=8)
unlink("stomach_tcga_green_module_females_males_tumourGO_dot.png")


green_module_females_males_tumourkeggPath_dot <- dotplot(green_module_females_males_tumourkeggPath, x = "GeneRatio", color = "p.adjust", showCategory = 10, font.size = 20, title = "Green module\nKEGG pathways")
ggsave(filename="stomach_tcga_green_module_females_males_tumourkeggPath_dot.png", plot=green_module_females_males_tumourkeggPath_dot, path = "./plots/wgcna_modules_enrichment/", width=10, height=8)
unlink("stomach_tcga_green_module_females_males_tumourkeggPath_dot.png")



# yellow(females)/brown(males) module is conserved in both networks
yellow_brown_module_females_males_tumourGO <- enrichment_2modules_go(all_genes, males_tumour_network, "brown", females_tumour_network, "yellow")
yellow_brown_module_females_males_tumourkeggPath <- enrichment_2modules_kegg(all_genes, males_tumour_network, "brown", females_tumour_network, "yellow")


yellow_brown_module_females_males_tumourGO_dot <- dotplot(yellow_brown_module_females_males_tumourGO, x = "GeneRatio", color = "p.adjust", showCategory = 10, split = NULL, font.size = 20, title = "Yellow/Brown module\nGO terms: BP")
ggsave(filename="stomach_tcga_yellow_brown_module_females_males_tumourGO_dot.png", plot=yellow_brown_module_females_males_tumourGO_dot, path = "./plots/wgcna_modules_enrichment/", width=10, height=8)
unlink("stomach_tcga_yellow_brown_module_females_males_tumourGO_dot.png")


yellow_brown_module_females_males_tumourkeggPath_dot <- dotplot(yellow_brown_module_females_males_tumourkeggPath, x = "GeneRatio", color = "p.adjust", showCategory = 10, font.size = 20, title = "Yellow/Brown module\nKEGG pathways")
ggsave(filename="stomach_tcga_yellow_brown_module_females_males_tumourkeggPath_dot.png", plot=yellow_brown_module_females_males_tumourkeggPath_dot, path = "./plots/wgcna_modules_enrichment/", width=9, height=6)
unlink("stomach_tcga_yellow_brown_module_females_males_tumourkeggPath_dot.png")







save(list=ls(), file="./r_workspaces/tcga_gtex_stomach_wgcna_networks_enrichment_analysis.RData")