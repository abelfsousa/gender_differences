# Understanding Gender Differential Susceptibility in Cancer


# Male vs Female differential expression by tissue (normal and tumour)
# TCGA + GTEx data


library(tidyverse)


source(file = "./r_scripts/utils.R")



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

diff_expr_tcga_tumour_table <- edgeR_diff_expression( design_tcga_tumour, thyroid_tcga_tumour_meta$gender, thyroid_tcga_tumour_counts, tcga.geneIDs.annot)
write.table(diff_expr_tcga_tumour_table, file="./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_tumour_table %>% filter(FDR < 0.05) %>% nrow()
#99


#diff_expr_tcga_tumour_table <- sam_diff_expression(as.matrix(thyroid_tcga_tumour_counts), thyroid_tcga_tumour_meta$gender, tcga.geneIDs.annot)



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
design_tcga_normal <- model.matrix(~ race + ethnicity + age_at_diagnosis + portion + gender, data = thyroid_tcga_normal_meta)

diff_expr_tcga_normal_table <- edgeR_diff_expression( design_tcga_normal, thyroid_tcga_normal_meta$gender, thyroid_tcga_normal_counts, tcga.geneIDs.annot)
write.table(diff_expr_tcga_normal_table, file="./files/diff_expr_thyroid_normal_tcga_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_normal_table %>% filter(FDR < 0.05) %>% nrow()
#15


#diff_expr_tcga_normal_table <- sam_diff_expression(thyroid_tcga_normal_counts, thyroid_tcga_normal_meta$gender, tcga.geneIDs.annot)



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

diff_expr_gtex_table <- edgeR_diff_expression( design_gtex, thyroid_gtex_meta$gender, thyroid_gtex_counts, tcga.geneIDs.annot)
write.table(diff_expr_gtex_table, file="./files/diff_expr_thyroid_normal_gtex_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_gtex_table %>% filter(FDR < 0.05) %>% nrow()
#754


#diff_expr_gtex_table <- sam_diff_expression(thyroid_gtex_counts, thyroid_gtex_meta$gender, tcga.geneIDs.annot)


save(list=ls(), file="r_workspaces/tcga_gtex_thyroid_diffExpr_malesVSfemales.RData")
