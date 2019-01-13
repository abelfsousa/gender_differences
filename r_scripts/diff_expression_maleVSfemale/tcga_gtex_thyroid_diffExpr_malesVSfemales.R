# Understanding Gender Differential Susceptibility in Cancer


# Male vs Female differential expression by tissue (normal and tumour)
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
thyroid_tcga_gtex_counts <- read.delim("./files/thyroid_tcga_gtex_counts.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_meta <- read.delim("./files/tcga_thca_meta.txt", row.names = c(1), stringsAsFactors=F)
thyroid_gtex_meta <- read.delim("./files/gtex_thyroid_meta.txt", row.names = c(1), stringsAsFactors=F)

rownames(thyroid_tcga_meta) <- gsub("-", ".", rownames(thyroid_tcga_meta))

sample_type_code <- c("1" = "tumour", "11" = "normal", "6" = "metastatic")

thyroid_tcga_meta$sample_type_id <- as.character(thyroid_tcga_meta$sample_type_id)
thyroid_tcga_meta$sample_type_id <- sample_type_code[thyroid_tcga_meta$sample_type_id]


# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))

gtex.geneIDs.annot <- read.table("./data/gtex.thyroid.stomach/geneAnnot.gencode.v19.txt", h=T, sep="\t", stringsAsFactors=F)
gtex.geneIDs.annot$geneID <- sapply(strsplit(gtex.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))




# -- Tumour

thyroid_tcga_tumour_meta <- thyroid_tcga_meta[thyroid_tcga_meta$sample_type_id == "tumour", ]
thyroid_tcga_tumour_counts <- thyroid_tcga_gtex_counts[, colnames(thyroid_tcga_gtex_counts) %in% rownames(thyroid_tcga_tumour_meta)]


# reorder datasets
thyroid_tcga_tumour_counts <- thyroid_tcga_tumour_counts[, order(colnames(thyroid_tcga_tumour_counts))]
thyroid_tcga_tumour_meta <- thyroid_tcga_tumour_meta[order(rownames(thyroid_tcga_tumour_meta)), ]


thyroid_tcga_tumour_counts %>% colnames %>% all.equal(rownames(thyroid_tcga_tumour_meta))
#TRUE


# define design matrix
# adjust for covariates
design_tcga_tumour <- model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + histological_type + tss + portion + plate + gender, data = thyroid_tcga_tumour_meta)


# DGEList object using edgeR
tcga_tumour_dgelist <- DGEList(counts=thyroid_tcga_tumour_counts)

# scale normalization using TMM method
tcga_tumour_dgelist <- calcNormFactors(tcga_tumour_dgelist)

# transform the counts using voom
tcga_tumour_dgelist <- voom(tcga_tumour_dgelist, design = design_tcga_tumour)

# fit the linear model for each gene
fit_tcga_tumour <- lmFit(tcga_tumour_dgelist, design_tcga_tumour)

#empirical bayes statistics for differential expression
fit_tcga_tumour <- eBayes(fit_tcga_tumour)


# extract summary and gene table with statistics
diff_expr_tcga_tumour_summary <- summary(decideTests(fit_tcga_tumour))
diff_expr_tcga_tumour_table <- topTable(fit_tcga_tumour, coef = "gendermale", n = Inf, sort.by = "p", p = 1)

diff_expr_tcga_tumour_table <- diff_expr_tcga_tumour_table %>%
    mutate(genes = rownames(diff_expr_tcga_tumour_table)) %>%
    # add gene annotations
    left_join(tcga.geneIDs.annot, by=c("genes" = "geneID")) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))
write.table(diff_expr_tcga_tumour_table, file="./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_tumour_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% nrow
#11


diff_expr_tcga_tumour_table2 <- edgeR_diff_expression( design_tcga_tumour, thyroid_tcga_tumour_meta$gender, thyroid_tcga_tumour_counts, tcga.geneIDs.annot)
write.table(diff_expr_tcga_tumour_table2, file="./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_tumour_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow
#12


length(intersect(diff_expr_tcga_tumour_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% pull(genes), diff_expr_tcga_tumour_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% pull(genes)))
#10




# -- Normal - TCGA

thyroid_tcga_normal_meta <- thyroid_tcga_meta[thyroid_tcga_meta$sample_type_id == "normal", ]
thyroid_tcga_normal_counts <- thyroid_tcga_gtex_counts[, colnames(thyroid_tcga_gtex_counts) %in% rownames(thyroid_tcga_normal_meta)]


# reorder datasets
thyroid_tcga_normal_counts <- thyroid_tcga_normal_counts[, order(colnames(thyroid_tcga_normal_counts))]
thyroid_tcga_normal_meta <- thyroid_tcga_normal_meta[order(rownames(thyroid_tcga_normal_meta)), ]


thyroid_tcga_normal_counts %>% colnames %>% all.equal(rownames(thyroid_tcga_normal_meta))
#TRUE


# define design matrix
# adjust for covariates
design_tcga_normal <- model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + portion + gender, data = thyroid_tcga_normal_meta)


# DGEList object using edgeR
tcga_normal_dgelist <- DGEList(counts=thyroid_tcga_normal_counts)

# scale normalization using TMM method
tcga_normal_dgelist <- calcNormFactors(tcga_normal_dgelist)

# transform the counts using voom
tcga_normal_dgelist <- voom(tcga_normal_dgelist, design = design_tcga_normal)

# fit the linear model for each gene
fit_tcga_normal <- lmFit(tcga_normal_dgelist, design_tcga_normal)

#empirical bayes statistics for differential expression
fit_tcga_normal <- eBayes(fit_tcga_normal)


# extract summary and gene table with statistics
diff_expr_tcga_normal_summary <- summary(decideTests(fit_tcga_normal))
diff_expr_tcga_normal_table <- topTable(fit_tcga_normal, coef = "gendermale", n = Inf, sort.by = "p", p = 1)

diff_expr_tcga_normal_table <- diff_expr_tcga_normal_table %>%
    mutate(genes = rownames(diff_expr_tcga_normal_table)) %>%
    # add gene annotations
    left_join(tcga.geneIDs.annot, by=c("genes" = "geneID")) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))
write.table(diff_expr_tcga_normal_table, file="./files/diff_expr_thyroid_normal_tcga_maleVSfemale_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_normal_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% nrow
#9


diff_expr_tcga_normal_table2 <- edgeR_diff_expression( design_tcga_normal, thyroid_tcga_normal_meta$gender, thyroid_tcga_normal_counts, tcga.geneIDs.annot)
write.table(diff_expr_tcga_normal_table2, file="./files/diff_expr_thyroid_normal_tcga_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_normal_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow
#10


length(intersect(diff_expr_tcga_normal_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% pull(genes), diff_expr_tcga_normal_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% pull(genes)))
#9




# -- Normal - GTEx

thyroid_gtex_meta <- as.data.frame(na.exclude(thyroid_gtex_meta))
thyroid_gtex_counts <- thyroid_tcga_gtex_counts[, colnames(thyroid_tcga_gtex_counts) %in% rownames(thyroid_gtex_meta)]


# reorder datasets
thyroid_gtex_counts <- thyroid_gtex_counts[, order(colnames(thyroid_gtex_counts))]
thyroid_gtex_meta <- thyroid_gtex_meta[order(rownames(thyroid_gtex_meta)), ]


thyroid_gtex_counts %>% colnames %>% all.equal(rownames(thyroid_gtex_meta))
#TRUE


gender_code <- c("1" = "male", "2" = "female")
thyroid_gtex_meta$GENDER <- as.character(thyroid_gtex_meta$GENDER)
thyroid_gtex_meta$GENDER <- gender_code[thyroid_gtex_meta$GENDER]
thyroid_gtex_meta$ETHNCTY <- as.character(thyroid_gtex_meta$ETHNCTY)
thyroid_gtex_meta$MHCANCERNM <- as.character(thyroid_gtex_meta$MHCANCERNM)


# define design matrix
# adjust for covariates
design_gtex <- model.matrix(~ SMRIN + AGE + ETHNCTY + MHCANCERNM + SMCENTER + SMTSTPTREF + SMNABTCHT + GENDER, data = thyroid_gtex_meta)


# DGEList object using edgeR
gtex_dgelist <- DGEList(counts=thyroid_gtex_counts)

# scale normalization using TMM method
gtex_dgelist <- calcNormFactors(gtex_dgelist)

# transform the counts using voom
gtex_dgelist <- voom(gtex_dgelist, design = design_gtex)

# fit the linear model for each gene
fit_gtex <- lmFit(gtex_dgelist, design_gtex)

#empirical bayes statistics for differential expression
fit_gtex <- eBayes(fit_gtex)


# extract summary and gene table with statistics
diff_expr_gtex_summary <- summary(decideTests(fit_gtex))
diff_expr_gtex_table <- topTable(fit_gtex, coef = "GENDERmale", n = Inf, sort.by = "p", p = 1)

diff_expr_gtex_table <- diff_expr_gtex_table %>%
    mutate(genes = rownames(diff_expr_gtex_table)) %>%
    # add gene annotations
    left_join(tcga.geneIDs.annot, by=c("genes" = "geneID")) %>%
    mutate(fdr = if_else(adj.P.Val > 0.05, "FDR > 0.05", if_else(adj.P.Val <= 0.01, "FDR <= 0.01", "FDR <= 0.05")))
write.table(diff_expr_gtex_table, file="./files/diff_expr_thyroid_normal_gtex_maleVSfemale_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_gtex_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% nrow
#11


diff_expr_gtex_table2 <- edgeR_diff_expression( design_gtex, thyroid_gtex_meta$gender, thyroid_gtex_counts, tcga.geneIDs.annot)
write.table(diff_expr_gtex_table2, file="./files/diff_expr_thyroid_normal_gtex_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_gtex_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow
#13


length(intersect(diff_expr_gtex_table2 %>% filter(FDR < 0.05, abs(logFC)>1) %>% pull(genes), diff_expr_gtex_table %>% filter(adj.P.Val < 0.05, abs(logFC)>1) %>% pull(genes)))
#10




save(list=ls(), file="r_workspaces/tcga_gtex_thyroid_diffExpr_malesVSfemales.RData")

