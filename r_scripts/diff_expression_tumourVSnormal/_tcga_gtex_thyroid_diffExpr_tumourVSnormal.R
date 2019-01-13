# Understanding Gender Differential Susceptibility in Cancer


# Tumour vs Normal differential expression by gender
# TCGA + GTEx data

# Differential expression using limma R package


library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(limma)

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


# add experimental group
sample_type_code <- c("GTEx" = "GTEx", "TCGA_normal" = "TCGA", "TCGA_tumour" = "TCGA")
thyroid_tcga_gtex_meta <- cbind(thyroid_tcga_gtex_meta, data.frame(data = sample_type_code[thyroid_tcga_gtex_meta$sample_type]))

# add tissue type
tissue_type_code <- c("GTEx" = "normal", "TCGA_normal" = "normal", "TCGA_tumour" = "tumour")
thyroid_tcga_gtex_meta <- cbind(thyroid_tcga_gtex_meta, data.frame(tissue_type = tissue_type_code[thyroid_tcga_gtex_meta$sample_type]))


thyroid_tcga_gtex_counts %>% colnames %>% all.equal(rownames(thyroid_tcga_gtex_meta))
#TRUE





# -- Males

thyroid_males_meta <- thyroid_tcga_gtex_meta[(thyroid_tcga_gtex_meta$gender == "male"), ]
thyroid_males_counts <- thyroid_tcga_gtex_counts[ ,colnames(thyroid_tcga_gtex_counts) %in% rownames(thyroid_males_meta)]

thyroid_males_counts %>% colnames %>% all.equal(rownames(thyroid_males_meta))
#TRUE


# define design matrix
design_males <- model.matrix(~ ethnicity + age + data + tissue_type, data = thyroid_males_meta)


# DGEList object using edgeR
thyroid_males_dgelist <- DGEList(counts=thyroid_males_counts)

# scale normalization using TMM method
thyroid_males_dgelist <- calcNormFactors(thyroid_males_dgelist)

# transform the counts using voom
thyroid_males_dgelist <- voom(thyroid_males_dgelist, design = design_males)




# pca betweewn GTEx and TCGA samples
pca1 <- prcomp(t(thyroid_males_dgelist$E), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca1 <- pca1 %>%
    mutate(sample_type = as.character(thyroid_males_meta$sample_type))


#scatterplot
pca1_males_plot <- ggplot( data=pca1, mapping=aes(x=PC1, y=PC2, color=as.factor(sample_type)) ) +
    geom_point() +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15)) +
    theme(legend.text=element_text(colour="black", size=13), legend.title=element_text(colour="black", size=15)) +
    scale_y_continuous(limits = c(-210, 210)) +
    scale_x_continuous(limits = c(-210, 210)) +
    labs(title="Males") +
    coord_fixed()
ggsave(filename="pca1_thyroid_males_tcga_gtex_voom.pdf", plot=pca1_males_plot, path = "./plots/")
unlink("pca1_thyroid_males_tcga_gtex_voom.pdf")




# regress-out batch effect between GTEx and TCGA samples
# use experimental group (TCGA or GTEx) as covariate 
thyroid_males_counts_cor <- removeBatchEffect(thyroid_males_dgelist, thyroid_males_meta$data)

thyroid_males_counts_cor <- thyroid_males_counts_cor %>% as.data.frame()



# pca betweewn GTEx and TCGA samples
pca2 <- prcomp(t(thyroid_males_counts_cor), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca2 <- pca2 %>%
    mutate(sample_type = as.character(thyroid_males_meta$sample_type))


#scatterplot
pca2_males_plot <- ggplot( data=pca2, mapping=aes(x=PC1, y=PC2, color=as.factor(sample_type)) ) +
    geom_point() +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15)) +
    theme(legend.text=element_text(colour="black", size=13), legend.title=element_text(colour="black", size=15)) +
    scale_y_continuous(limits = c(-210, 210)) +
    scale_x_continuous(limits = c(-210, 210)) +
    labs(title="Males") +
    coord_fixed()
ggsave(filename="pca2_thyroid_males_tcga_gtex_voom.pdf", plot=pca2_males_plot, path = "./plots/")
unlink("pca2_thyroid_males_tcga_gtex_voom.pdf")



# fit the linear model for each gene
fit_males <- lmFit(thyroid_males_dgelist, design_males)

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
write.table(diff_expr_males_table, file="./files/diff_expr_thyroid_males_tcga_gtex_tumourVSnormal.txt", sep = "\t", quote=F, row.names=F)
write.table(diff_expr_males_table[(diff_expr_males_table$adj.P.Val < 0.05), ], file="./files/signf_diff_expr_thyroid_males_tcga_gtex_tumourVSnormal.txt", sep = "\t", quote=F, row.names=F)



diff_expr_males_table_wilcoxon <- run.wilcoxon.diffExpr(
    rownames(thyroid_males_meta[thyroid_males_meta$tissue_type=="tumour", ]),
    rownames(thyroid_males_meta[thyroid_males_meta$tissue_type=="normal", ]),
    thyroid_males_dgelist$E
)
write.table(diff_expr_males_table_wilcoxon, file="./files/diff_expr_thyroid_males_tcga_gtex_tumourVSnormal_wilcoxon.txt", sep = "\t", quote=F, row.names=F)





# -- Females


thyroid_females_meta <- thyroid_tcga_gtex_meta[(thyroid_tcga_gtex_meta$gender == "female"), ]
thyroid_females_counts <- thyroid_tcga_gtex_counts[ ,colnames(thyroid_tcga_gtex_counts) %in% rownames(thyroid_females_meta)]

thyroid_females_counts %>% colnames %>% all.equal(rownames(thyroid_females_meta))
#TRUE


# define design matrix
design_females <- model.matrix(~ ethnicity + age + data + tissue_type, data = thyroid_females_meta)


# DGEList object using edgeR
thyroid_females_dgelist <- DGEList(counts=thyroid_females_counts)

# scale normalization using TMM method
thyroid_females_dgelist <- calcNormFactors(thyroid_females_dgelist)

# transform the counts using voom
thyroid_females_dgelist <- voom(thyroid_females_dgelist, design = design_females)




# pca betweewn GTEx and TCGA samples
pca1_females <- prcomp(t(thyroid_females_dgelist$E), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca1_females <- pca1_females %>%
    mutate(sample_type = as.character(thyroid_females_meta$sample_type))


#scatterplot
pca1_females_plot <- ggplot( data=pca1_females, mapping=aes(x=PC1, y=PC2, color=as.factor(sample_type)) ) +
    geom_point() +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15)) +
    theme(legend.text=element_text(colour="black", size=13), legend.title=element_text(colour="black", size=15)) +
    scale_y_continuous(limits = c(-210, 210)) +
    scale_x_continuous(limits = c(-210, 210)) +
    labs(title="Females") +
    coord_fixed()
ggsave(filename="pca1_thyroid_females_tcga_gtex_voom.pdf", plot=pca1_females_plot, path = "./plots/")
unlink("pca1_thyroid_females_tcga_gtex_voom.pdf")




# regress-out batch effect between GTEx and TCGA samples
# use experimental group (TCGA or GTEx) as covariate
thyroid_females_counts_cor <- removeBatchEffect(thyroid_females_dgelist, thyroid_females_meta$data)

thyroid_females_counts_cor <- thyroid_females_counts_cor %>% as.data.frame()



# pca betweewn GTEx and TCGA samples
pca2_females <- prcomp(t(thyroid_females_counts_cor), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca2_females <- pca2_females %>%
    mutate(sample_type = as.character(thyroid_females_meta$sample_type))


#scatterplot
pca2_females_plot <- ggplot( data=pca2_females, mapping=aes(x=PC1, y=PC2, color=as.factor(sample_type)) ) +
    geom_point() +
    theme(plot.title=element_text(size=15)) +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=15), axis.text.x=element_text(colour="black", size=15)) +
    theme(legend.text=element_text(colour="black", size=13), legend.title=element_text(colour="black", size=15)) +
    scale_y_continuous(limits = c(-300, 300)) +
    scale_x_continuous(limits = c(-300, 300)) +
    labs(title="Females") +
    coord_fixed()
ggsave(filename="pca2_thyroid_females_tcga_gtex_voom.pdf", plot=pca2_females_plot, path = "./plots/")
unlink("pca2_thyroid_females_tcga_gtex_voom.pdf")



# fit the linear model for each gene
fit_females <- lmFit(thyroid_females_dgelist, design_females)

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
write.table(diff_expr_females_table, file="./files/diff_expr_thyroid_females_tcga_gtex_tumourVSnormal.txt", sep = "\t", quote=F, row.names=F)
write.table(diff_expr_females_table[(diff_expr_females_table$adj.P.Val < 0.05), ], file="./files/signf_diff_expr_thyroid_females_tcga_gtex_tumourVSnormal.txt", sep = "\t", quote=F, row.names=F)


diff_expr_females_table_wilcoxon <- run.wilcoxon.diffExpr(
    rownames(thyroid_females_meta[thyroid_females_meta$tissue_type=="tumour", ]),
    rownames(thyroid_females_meta[thyroid_females_meta$tissue_type=="normal", ]),
    thyroid_females_dgelist$E
)
write.table(diff_expr_females_table_wilcoxon, file="./files/diff_expr_thyroid_females_tcga_gtex_tumourVSnormal_wilcoxon.txt", sep = "\t", quote=F, row.names=F)






save(list=ls(), file="r_workspaces/tcga_gtex_thyroid_diffExpr_tumourVSnormal.RData")


