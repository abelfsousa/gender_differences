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






# -- Stomach


# load datasets
stomach_tcga_counts <- read.delim("./files/stomach_tcga_gtex_counts.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_meta <- read.delim("./files/tcga_stad_meta.txt", row.names = c(1), stringsAsFactors=F)

rownames(stomach_tcga_meta) <- gsub("-", ".", rownames(stomach_tcga_meta))

sample_type_code <- c("1" = "tumour", "11" = "normal", "6" = "metastatic")

stomach_tcga_meta$sample_type_id <- as.character(stomach_tcga_meta$sample_type_id)
stomach_tcga_meta$sample_type_id <- sample_type_code[stomach_tcga_meta$sample_type_id]


# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))


# remove gtex samples from counts matrix
stomach_tcga_counts <- stomach_tcga_counts[, substring(colnames(stomach_tcga_counts), 1, 4) != "GTEX"]


# reorder datasets
stomach_tcga_counts <- stomach_tcga_counts[, order(colnames(stomach_tcga_counts))]
stomach_tcga_meta <- stomach_tcga_meta[order(rownames(stomach_tcga_meta)), ]



stomach_tcga_counts %>% colnames %>% all.equal(rownames(stomach_tcga_meta))
#TRUE





# -- Males

stomach_males_meta <- stomach_tcga_meta[(stomach_tcga_meta$gender == "male"), ]
stomach_males_counts <- stomach_tcga_counts[ ,colnames(stomach_tcga_counts) %in% rownames(stomach_males_meta)]

stomach_males_counts %>% colnames %>% all.equal(rownames(stomach_males_meta))
#TRUE


# define design matrix
# adjust for covariates
design_males <- model.matrix(~ race + ethnicity + age_at_diagnosis + portion + plate + sample_type_id, data = stomach_males_meta)


# DGEList object using edgeR
stomach_males_dgelist <- DGEList(counts=stomach_males_counts)

# scale normalization using TMM method
stomach_males_dgelist <- calcNormFactors(stomach_males_dgelist)

# transform the counts using voom
stomach_males_dgelist <- voom(stomach_males_dgelist, design = design_males)

# fit the linear model for each gene
fit_males <- lmFit(stomach_males_dgelist, design_males)

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
write.table(diff_expr_males_table, file="./files/diff_expr_stomach_males_tcga_tumourVSnormal_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% nrow
#2367


diff_expr_males_table2 <- edgeR_diff_expression( design_males, stomach_males_meta$sample_type_id, stomach_males_counts, tcga.geneIDs.annot)
write.table(diff_expr_males_table2, file="./files/diff_expr_stomach_males_tcga_tumourVSnormal_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow
#2294


length(intersect(diff_expr_males_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% pull(genes), diff_expr_males_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% pull(genes)))
#2019


# -- Females


stomach_females_meta <- stomach_tcga_meta[(stomach_tcga_meta$gender == "female"), ]
stomach_females_counts <- stomach_tcga_counts[ ,colnames(stomach_tcga_counts) %in% rownames(stomach_females_meta)]

stomach_females_counts %>% colnames %>% all.equal(rownames(stomach_females_meta))
#TRUE


# define design matrix
# adjust for covariates
design_females <- model.matrix(~ race + ethnicity + age_at_diagnosis + portion + plate + sample_type_id, data = stomach_females_meta)


# DGEList object using edgeR
stomach_females_dgelist <- DGEList(counts=stomach_females_counts)

# scale normalization using TMM method
stomach_females_dgelist <- calcNormFactors(stomach_females_dgelist)

# transform the counts using voom
stomach_females_dgelist <- voom(stomach_females_dgelist, design = design_females)

# fit the linear model for each gene
fit_females <- lmFit(stomach_females_dgelist, design_females)

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
write.table(diff_expr_females_table, file="./files/diff_expr_stomach_females_tcga_tumourVSnormal_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% nrow
#1841



diff_expr_females_table2 <- edgeR_diff_expression( design_females, stomach_females_meta$sample_type_id, stomach_females_counts, tcga.geneIDs.annot)
write.table(diff_expr_females_table2, file="./files/diff_expr_stomach_females_tcga_tumourVSnormal_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow
#1856


length(intersect(diff_expr_females_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% pull(genes), diff_expr_females_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% pull(genes)))
#1586


save(list=ls(), file="r_workspaces/tcga_stomach_diffExpr_tumourVSnormal.RData")
