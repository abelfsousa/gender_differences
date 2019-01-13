# Understanding Gender Differential Susceptibility in Cancer


# Tumour vs Normal differential expression by gender
# TCGA + GTEx data

# Differential expression using limma R package


library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(limma)

source("utils.R")





# -- Stomach

# load datasets
stomach_tcga_gtex_counts <- read.delim("./files/stomach_tcga_gtex_counts.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_gtex_meta <- read.delim("./files/stomach_tcga_gtex_meta.txt", row.names = c(1), stringsAsFactors=F)
rownames(stomach_tcga_gtex_meta) <- gsub("-", ".", rownames(stomach_tcga_gtex_meta))

# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))



# reorder datasets
stomach_tcga_gtex_counts <- stomach_tcga_gtex_counts[, order(colnames(stomach_tcga_gtex_counts))]
stomach_tcga_gtex_meta <- stomach_tcga_gtex_meta[order(rownames(stomach_tcga_gtex_meta)), ]


# add experimental group
sample_type_code <- c("GTEx" = "GTEx", "TCGA_normal" = "TCGA", "TCGA_tumour" = "TCGA")
stomach_tcga_gtex_meta <- cbind(stomach_tcga_gtex_meta, data.frame(data = sample_type_code[stomach_tcga_gtex_meta$sample_type]))


# add tissue type
tissue_type_code <- c("GTEx" = "normal", "TCGA_normal" = "normal", "TCGA_tumour" = "tumour")
stomach_tcga_gtex_meta <- cbind(stomach_tcga_gtex_meta, data.frame(tissue_type = tissue_type_code[stomach_tcga_gtex_meta$sample_type]))


stomach_tcga_gtex_counts %>% colnames %>% all.equal(rownames(stomach_tcga_gtex_meta))
#TRUE





# -- Males

stomach_males_meta <- stomach_tcga_gtex_meta[(stomach_tcga_gtex_meta$gender == "male"), ]
stomach_males_counts <- stomach_tcga_gtex_counts[ ,colnames(stomach_tcga_gtex_counts) %in% rownames(stomach_males_meta)]

stomach_males_counts %>% colnames %>% all.equal(rownames(stomach_males_meta))
#TRUE


# define design matrix
design_males <- model.matrix(~ ethnicity + age + data + tissue_type, data = stomach_males_meta)


# DGEList object using edgeR
stomach_males_dgelist <- DGEList(counts=stomach_males_counts)

# scale normalization using TMM method
stomach_males_dgelist <- calcNormFactors(stomach_males_dgelist)

# transform the counts using voom
stomach_males_dgelist <- voom(stomach_males_dgelist, design = design_males)




# pca betweewn GTEx and TCGA samples
pca1 <- prcomp(t(stomach_males_dgelist$E), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca1 <- pca1 %>%
    mutate(sample_type = as.character(stomach_males_meta$sample_type))


#scatterplot
pca1_males_plot <- ggplot( data=pca1, mapping=aes(x=PC1, y=PC2, color=as.factor(sample_type)) ) +
    geom_point() +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15)) +
    theme(legend.text=element_text(colour="black", size=13), legend.title=element_text(colour="black", size=15)) +
    scale_y_continuous(limits = c(-200, 200)) +
    scale_x_continuous(limits = c(-200, 200)) +
    labs(title="Males") +
    coord_fixed()
ggsave(filename="pca1_stomach_males_tcga_gtex_voom.pdf", plot=pca1_males_plot, path = "./plots/")
unlink("pca1_stomach_males_tcga_gtex_voom.pdf")




# regress-out batch effect between GTEx and TCGA samples
# use experimental group (TCGA or GTEx) as covariate 
stomach_males_counts_cor <- removeBatchEffect(stomach_males_dgelist, stomach_males_meta$data)

stomach_males_counts_cor <- stomach_males_counts_cor %>% as.data.frame()



# pca betweewn GTEx and TCGA samples
pca2 <- prcomp(t(stomach_males_counts_cor), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca2 <- pca2 %>%
    mutate(sample_type = as.character(stomach_males_meta$sample_type))


#scatterplot
pca2_males_plot <- ggplot( data=pca2, mapping=aes(x=PC1, y=PC2, color=as.factor(sample_type)) ) +
    geom_point() +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15)) +
    theme(legend.text=element_text(colour="black", size=13), legend.title=element_text(colour="black", size=15)) +
    scale_y_continuous(limits = c(-200, 200)) +
    scale_x_continuous(limits = c(-200, 200)) +
    labs(title="Males") +
    coord_fixed()
ggsave(filename="pca2_stomach_males_tcga_gtex_voom.pdf", plot=pca2_males_plot, path = "./plots/")
unlink("pca2_stomach_males_tcga_gtex_voom.pdf")



# fit the linear model for each gene
fit_males <- lmFit(stomach_males_dgelist, design_males)

#empirical bayes statistics for differential expression
fit_males <- eBayes(fit_males)


# extract summary and gene table with statistics
diff_expr_males_summary <- summary(decideTests(fit_males))
diff_expr_males_table <- topTable(fit_males, coef = "tissue_typetumour", n = Inf, sort.by = "p", p = 1)

diff_expr_males_table <- diff_expr_males_table %>%
    mutate(genes = rownames(diff_expr_males_table)) %>%
    # add gene annotations
    left_join(tcga.geneIDs.annot, by=c("genes" = "geneID")) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))
write.table(diff_expr_males_table, file="./files/diff_expr_stomach_males_tcga_gtex_tumourVSnormal.txt", sep = "\t", quote=F, row.names=F)
write.table(diff_expr_males_table[(diff_expr_males_table$adj.P.Val < 0.05), ], file="./files/signf_diff_expr_stomach_males_tcga_gtex_tumourVSnormal.txt", sep = "\t", quote=F, row.names=F)



diff_expr_males_table_wilcoxon <- run.wilcoxon.diffExpr(
    rownames(stomach_males_meta[stomach_males_meta$tissue_type=="tumour", ]),
    rownames(stomach_males_meta[stomach_males_meta$tissue_type=="normal", ]),
    stomach_males_dgelist$E
)
write.table(diff_expr_males_table_wilcoxon, file="./files/diff_expr_stomach_males_tcga_gtex_tumourVSnormal_wilcoxon.txt", sep = "\t", quote=F, row.names=F)






# -- Females


stomach_females_meta <- stomach_tcga_gtex_meta[(stomach_tcga_gtex_meta$gender == "female"), ]
stomach_females_counts <- stomach_tcga_gtex_counts[ ,colnames(stomach_tcga_gtex_counts) %in% rownames(stomach_females_meta)]

stomach_females_counts %>% colnames %>% all.equal(rownames(stomach_females_meta))
#TRUE


# define design matrix
design_females <- model.matrix(~ ethnicity + age + data + tissue_type, data = stomach_females_meta)


# DGEList object using edgeR
stomach_females_dgelist <- DGEList(counts=stomach_females_counts)

# scale normalization using TMM method
stomach_females_dgelist <- calcNormFactors(stomach_females_dgelist)

# transform the counts using voom
stomach_females_dgelist <- voom(stomach_females_dgelist, design = design_females)




# pca betweewn GTEx and TCGA samples
pca1_females <- prcomp(t(stomach_females_dgelist$E), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca1_females <- pca1_females %>%
    mutate(sample_type = as.character(stomach_females_meta$sample_type))


#scatterplot
pca1_females_plot <- ggplot( data=pca1_females, mapping=aes(x=PC1, y=PC2, color=as.factor(sample_type)) ) +
    geom_point() +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15)) +
    theme(legend.text=element_text(colour="black", size=13), legend.title=element_text(colour="black", size=15)) +
    scale_y_continuous(limits = c(-200, 200)) +
    scale_x_continuous(limits = c(-200, 200)) +
    labs(title="Females") +
    coord_fixed()
ggsave(filename="pca1_stomach_females_tcga_gtex_voom.pdf", plot=pca1_females_plot, path = "./plots/")
unlink("pca1_stomach_females_tcga_gtex_voom.pdf")




# regress-out batch effect between GTEx and TCGA samples
# use experimental group (TCGA or GTEx) as covariate
stomach_females_counts_cor <- removeBatchEffect(stomach_females_dgelist, stomach_females_meta$data)

stomach_females_counts_cor <- stomach_females_counts_cor %>% as.data.frame()



# pca betweewn GTEx and TCGA samples
pca2_females <- prcomp(t(stomach_females_counts_cor), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca2_females <- pca2_females %>%
    mutate(sample_type = as.character(stomach_females_meta$sample_type))


#scatterplot
pca2_females_plot <- ggplot( data=pca2_females, mapping=aes(x=PC1, y=PC2, color=as.factor(sample_type)) ) +
    geom_point() +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15)) +
    theme(legend.text=element_text(colour="black", size=13), legend.title=element_text(colour="black", size=15)) +
    scale_y_continuous(limits = c(-200, 200)) +
    scale_x_continuous(limits = c(-200, 200)) +
    labs(title="Females") +
    coord_fixed()
ggsave(filename="pca2_stomach_females_tcga_gtex_voom.pdf", plot=pca2_females_plot, path = "./plots/")
unlink("pca2_stomach_females_tcga_gtex_voom.pdf")



# fit the linear model for each gene
fit_females <- lmFit(stomach_females_dgelist, design_females)

#empirical bayes statistics for differential expression
fit_females <- eBayes(fit_females)


# extract summary and gene table with statistics
diff_expr_females_summary <- summary(decideTests(fit_females))
diff_expr_females_table <- topTable(fit_females, coef = "tissue_typetumour", n = Inf, sort.by = "p", p = 1)

diff_expr_females_table <- diff_expr_females_table %>%
    mutate(genes = rownames(diff_expr_females_table)) %>%
    # add gene annotations
    left_join(tcga.geneIDs.annot, by=c("genes" = "geneID")) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))
write.table(diff_expr_females_table, file="./files/diff_expr_stomach_females_tcga_gtex_tumourVSnormal.txt", sep = "\t", quote=F, row.names=F)
write.table(diff_expr_females_table[(diff_expr_females_table$adj.P.Val < 0.05), ], file="./files/signf_diff_expr_stomach_females_tcga_gtex_tumourVSnormal.txt", sep = "\t", quote=F, row.names=F)



diff_expr_females_table_wilcoxon <- run.wilcoxon.diffExpr(
    rownames(stomach_females_meta[stomach_females_meta$tissue_type=="tumour", ]),
    rownames(stomach_females_meta[stomach_females_meta$tissue_type=="normal", ]),
    stomach_females_dgelist$E
)
write.table(diff_expr_females_table_wilcoxon, file="./files/diff_expr_stomach_females_tcga_gtex_tumourVSnormal_wilcoxon.txt", sep = "\t", quote=F, row.names=F)






save(list=ls(), file="r_workspaces/tcga_gtex_stomach_diffExpr_tumourVSnormal.RData")


