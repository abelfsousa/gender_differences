# Understanding Gender Differential Susceptibility in Cancer


# Male vs Female differential expression by tissue (normal and tumour)
# TCGA + GTEx data


library(tidyverse)


source(file = "./r_scripts/utils.R")



# -- Stomach


# load datasets
stomach_tcga_gtex_counts <- read.delim("./files/stomach_tcga_gtex_counts.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_meta <- read.delim("./files/tcga_stad_meta.txt", row.names = c(1), stringsAsFactors=F)
stomach_gtex_meta <- read.delim("./files/gtex_stomach_meta.txt", row.names = c(1), stringsAsFactors=F)

rownames(stomach_tcga_meta) <- gsub("-", ".", rownames(stomach_tcga_meta))

sample_type_code <- c("1" = "tumour", "11" = "normal", "6" = "metastatic")

stomach_tcga_meta$sample_type_id <- as.character(stomach_tcga_meta$sample_type_id)
stomach_tcga_meta$sample_type_id <- sample_type_code[stomach_tcga_meta$sample_type_id]


# load annotation
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)
tcga.geneIDs.annot$geneID <- sapply(strsplit(tcga.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))

gtex.geneIDs.annot <- read.table("./data/gtex.thyroid.stomach/geneAnnot.gencode.v19.txt", h=T, sep="\t", stringsAsFactors=F)
gtex.geneIDs.annot$geneID <- sapply(strsplit(gtex.geneIDs.annot$geneID, split=".", fixed=T), function(x)(x[1]))




# -- Tumour

stomach_tcga_tumour_meta <- stomach_tcga_meta[stomach_tcga_meta$sample_type_id == "tumour", ]
stomach_tcga_tumour_counts <- stomach_tcga_gtex_counts[, colnames(stomach_tcga_gtex_counts) %in% rownames(stomach_tcga_tumour_meta)]


# reorder datasets
stomach_tcga_tumour_counts <- stomach_tcga_tumour_counts[, order(colnames(stomach_tcga_tumour_counts))]
stomach_tcga_tumour_meta <- stomach_tcga_tumour_meta[order(rownames(stomach_tcga_tumour_meta)), ]


stomach_tcga_tumour_counts %>% colnames %>% all.equal(rownames(stomach_tcga_tumour_meta))
#TRUE


# define design matrix
# adjust for covariates
design_tcga_tumour <- model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + Lauren_Class + portion + plate + gender, data = stomach_tcga_tumour_meta)

diff_expr_tcga_tumour_table <- edgeR_diff_expression(design_tcga_tumour, stomach_tcga_tumour_meta$gender, stomach_tcga_tumour_counts, tcga.geneIDs.annot)
write.table(diff_expr_tcga_tumour_table, file="./files/diff_expr_stomach_tumour_tcga_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_tumour_table %>% filter(FDR < 0.05) %>% nrow()
#44




# -- Normal - TCGA

stomach_tcga_normal_meta <- stomach_tcga_meta[stomach_tcga_meta$sample_type_id == "normal", ]
stomach_tcga_normal_counts <- stomach_tcga_gtex_counts[, colnames(stomach_tcga_gtex_counts) %in% rownames(stomach_tcga_normal_meta)]


# reorder datasets
stomach_tcga_normal_counts <- stomach_tcga_normal_counts[, order(colnames(stomach_tcga_normal_counts))]
stomach_tcga_normal_meta <- stomach_tcga_normal_meta[order(rownames(stomach_tcga_normal_meta)), ]


stomach_tcga_normal_counts %>% colnames %>% all.equal(rownames(stomach_tcga_normal_meta))
#TRUE


# define design matrix
# adjust for covariates
design_tcga_normal <- model.matrix(~ race + age_at_diagnosis + portion + gender, data = stomach_tcga_normal_meta)

diff_expr_tcga_normal_table <- edgeR_diff_expression(design_tcga_normal, stomach_tcga_normal_meta$gender, stomach_tcga_normal_counts, tcga.geneIDs.annot)
write.table(diff_expr_tcga_normal_table, file="./files/diff_expr_stomach_normal_tcga_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_normal_table %>% filter(FDR < 0.05) %>% nrow()
#2




# -- Normal - GTEx

stomach_gtex_meta <- as.data.frame(na.exclude(stomach_gtex_meta))
stomach_gtex_counts <- stomach_tcga_gtex_counts[, colnames(stomach_tcga_gtex_counts) %in% rownames(stomach_gtex_meta)]


# reorder datasets
stomach_gtex_counts <- stomach_gtex_counts[, order(colnames(stomach_gtex_counts))]
stomach_gtex_meta <- stomach_gtex_meta[order(rownames(stomach_gtex_meta)), ]


stomach_gtex_counts %>% colnames %>% all.equal(rownames(stomach_gtex_meta))
#TRUE


gender_code <- c("1" = "male", "2" = "female")
stomach_gtex_meta$GENDER <- as.character(stomach_gtex_meta$GENDER)
stomach_gtex_meta$GENDER <- gender_code[stomach_gtex_meta$GENDER]
stomach_gtex_meta$ETHNCTY <- as.character(stomach_gtex_meta$ETHNCTY)
stomach_gtex_meta$MHCANCERNM <- as.character(stomach_gtex_meta$MHCANCERNM)


# define design matrix
# adjust for covariates
design_gtex <- model.matrix(~ SMRIN + AGE + ETHNCTY + MHCANCERNM + SMCENTER + SMTSTPTREF + SMNABTCHT + GENDER, data = stomach_gtex_meta)

diff_expr_gtex_table <- edgeR_diff_expression(design_gtex, stomach_gtex_meta$gender, stomach_gtex_counts, tcga.geneIDs.annot)
write.table(diff_expr_gtex_table, file="./files/diff_expr_stomach_normal_gtex_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_gtex_table %>% filter(FDR < 0.05) %>% nrow()
#61



save(list=ls(), file="r_workspaces/tcga_gtex_stomach_diffExpr_malesVSfemales.RData")
