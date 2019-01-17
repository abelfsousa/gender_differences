# Understanding Gender Differential Susceptibility in Cancer


# Tumour vs Normal differential expression by gender
# TCGA + GTEx data

# Differential expression using limma R package


library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(limma)



edgeR_diff_expression <- function(design, factor, expression.data, gene.annot){

	#differential expression analysis by edgeR

	#adjust for covariates when defining design matrix:
	#model.matrix(~ cov1 + cov2 + cov..., data)

	library(edgeR)

	#create a DGEList object
	y <- DGEList( counts=expression.data, group=factor, genes=rownames(expression.data) )

	#normalization TMM
	y <- calcNormFactors(y)
	#estimating NB dispertions
	y <- estimateDisp(y, design)
	#fitting the negative binomial GLM
	fit <- glmFit(y, design)
	#likelihood ratio test
	lrt <- glmLRT(fit)

	degs <- cbind(lrt$genes, lrt$table, p.adjust(lrt$table$PValue, method="BH"))
	colnames(degs)[6] <- "FDR"
	degs <- merge(degs, gene.annot, by.x = "genes", by.y = "geneID")
	degs <- degs[,c("genes","geneName","geneType","chrom","logFC","logCPM","LR","PValue","FDR")]
	degs <- degs[order(degs$FDR,decreasing=F),]

	return(degs)
}






# -- Thyroid


# load datasets
thyroid_tcga_counts <- read.delim("./files/thyroid_tcga_gtex_counts.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_meta <- read.delim("./files/tcga_thca_meta.txt", row.names = c(1), stringsAsFactors=F)

rownames(thyroid_tcga_meta) <- gsub("-", ".", rownames(thyroid_tcga_meta))

sample_type_code <- c("1" = "tumour", "11" = "normal", "6" = "metastatic")

thyroid_tcga_meta$sample_type_id <- as.character(thyroid_tcga_meta$sample_type_id)
thyroid_tcga_meta$sample_type_id <- sample_type_code[thyroid_tcga_meta$sample_type_id]


# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))


# remove gtex samples from counts matrix
thyroid_tcga_counts <- thyroid_tcga_counts[, substring(colnames(thyroid_tcga_counts), 1, 4) != "GTEX"]


# reorder datasets
thyroid_tcga_counts <- thyroid_tcga_counts[, order(colnames(thyroid_tcga_counts))]
thyroid_tcga_meta <- thyroid_tcga_meta[order(rownames(thyroid_tcga_meta)), ]


# remove metastatic samples
thyroid_tcga_counts <- thyroid_tcga_counts[, !colnames(thyroid_tcga_counts) %in% rownames(thyroid_tcga_meta[(thyroid_tcga_meta$sample_type_id == "metastatic"), ])]
thyroid_tcga_meta <- thyroid_tcga_meta[!(thyroid_tcga_meta$sample_type_id == "metastatic"), ]



thyroid_tcga_counts %>% colnames %>% all.equal(rownames(thyroid_tcga_meta))
#TRUE





# -- Males

thyroid_males_meta <- thyroid_tcga_meta[(thyroid_tcga_meta$gender == "male"), ]
thyroid_males_counts <- thyroid_tcga_counts[ ,colnames(thyroid_tcga_counts) %in% rownames(thyroid_males_meta)]

thyroid_males_counts %>% colnames %>% all.equal(rownames(thyroid_males_meta))
#TRUE


# define design matrix
# adjust for covariates
design_males <- model.matrix(~ race + ethnicity + age_at_diagnosis + tss + portion + plate + sample_type_id, data = thyroid_males_meta)


# DGEList object using edgeR
thyroid_males_dgelist <- DGEList(counts=thyroid_males_counts)

# scale normalization using TMM method
thyroid_males_dgelist <- calcNormFactors(thyroid_males_dgelist)

# transform the counts using voom
thyroid_males_dgelist <- voom(thyroid_males_dgelist, design = design_males)

# fit the linear model for each gene
fit_males <- lmFit(thyroid_males_dgelist, design_males)

#empirical bayes statistics for differential expression
fit_males <- eBayes(fit_males)


# extract summary and gene table with statistics
diff_expr_males_summary <- summary(decideTests(fit_males))
diff_expr_males_table <- topTable(fit_males, coef = "sample_type_idtumour", n = Inf, sort.by = "p", p = 1)

diff_expr_males_table <- diff_expr_males_table %>%
    mutate(genes = rownames(diff_expr_males_table)) %>%
    # add gene annotations
    left_join(tcga.geneIDs.annot, by=c("genes" = "geneID")) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))
write.table(diff_expr_males_table, file="./files/diff_expr_thyroid_males_tcga_tumourVSnormal_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% nrow
#1471


diff_expr_males_table2 <- edgeR_diff_expression( design_males, thyroid_males_meta$sample_type_id, thyroid_males_counts, tcga.geneIDs.annot)
write.table(diff_expr_males_table2, file="./files/diff_expr_thyroid_males_tcga_tumourVSnormal_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow
#1504


length(intersect(diff_expr_males_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% pull(genes), diff_expr_males_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% pull(genes)))
#1379


# -- Females


thyroid_females_meta <- thyroid_tcga_meta[(thyroid_tcga_meta$gender == "female"), ]
thyroid_females_counts <- thyroid_tcga_counts[ ,colnames(thyroid_tcga_counts) %in% rownames(thyroid_females_meta)]

thyroid_females_counts %>% colnames %>% all.equal(rownames(thyroid_females_meta))
#TRUE


# define design matrix
# adjust for covariates
design_females <- model.matrix(~ race + ethnicity + age_at_diagnosis + tss + portion + plate + sample_type_id, data = thyroid_females_meta)


# DGEList object using edgeR
thyroid_females_dgelist <- DGEList(counts=thyroid_females_counts)

# scale normalization using TMM method
thyroid_females_dgelist <- calcNormFactors(thyroid_females_dgelist)

# transform the counts using voom
thyroid_females_dgelist <- voom(thyroid_females_dgelist, design = design_females)

# fit the linear model for each gene
fit_females <- lmFit(thyroid_females_dgelist, design_females)

#empirical bayes statistics for differential expression
fit_females <- eBayes(fit_females)


# extract summary and gene table with statistics
diff_expr_females_summary <- summary(decideTests(fit_females))
diff_expr_females_table <- topTable(fit_females, coef = "sample_type_idtumour", n = Inf, sort.by = "p", p = 1)

diff_expr_females_table <- diff_expr_females_table %>%
    mutate(genes = rownames(diff_expr_females_table)) %>%
    # add gene annotations
    left_join(tcga.geneIDs.annot, by=c("genes" = "geneID")) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))
write.table(diff_expr_females_table, file="./files/diff_expr_thyroid_females_tcga_tumourVSnormal_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% nrow
#1202



diff_expr_females_table2 <- edgeR_diff_expression( design_females, thyroid_females_meta$sample_type_id, thyroid_females_counts, tcga.geneIDs.annot)
write.table(diff_expr_females_table2, file="./files/diff_expr_thyroid_females_tcga_tumourVSnormal_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow
#1148


length(intersect(diff_expr_females_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% pull(genes), diff_expr_females_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% pull(genes)))
#1047


save(list=ls(), file="r_workspaces/tcga_thyroid_diffExpr_tumourVSnormal.RData")
