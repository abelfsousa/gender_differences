# Understanding Gender Differential Susceptibility in Cancer


# Male vs Female differential expression by tissue (normal and tumour)
# TCGA + GTEx data


library(tidyverse)


source(file = "./r_scripts/utils.R")



# -- Thyroid


# load datasets
thyroid_tcga_gtex_counts <- read.delim("./files/thyroid_tcga_gtex_counts.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_gtex_fpkm <- read.delim("./files/thyroid_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
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
thyroid_tcga_tumour_fpkm <- thyroid_tcga_gtex_fpkm[, colnames(thyroid_tcga_gtex_fpkm) %in% rownames(thyroid_tcga_tumour_meta)]

# reorder datasets
thyroid_tcga_tumour_counts <- thyroid_tcga_tumour_counts[, order(colnames(thyroid_tcga_tumour_counts))]
thyroid_tcga_tumour_fpkm <- thyroid_tcga_tumour_fpkm[, order(colnames(thyroid_tcga_tumour_fpkm))]
thyroid_tcga_tumour_meta <- thyroid_tcga_tumour_meta[order(rownames(thyroid_tcga_tumour_meta)), ]


thyroid_tcga_tumour_counts %>% colnames %>% all.equal(rownames(thyroid_tcga_tumour_meta))
thyroid_tcga_tumour_fpkm %>% colnames %>% all.equal(rownames(thyroid_tcga_tumour_meta))
#TRUE


# define design matrix
# adjust for covariates
design_tcga_tumour <- model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + histological_type + tss + portion + plate + gender, data = thyroid_tcga_tumour_meta)

diff_expr_tcga_tumour_table <- edgeR_diff_expression( design_tcga_tumour, thyroid_tcga_tumour_meta$gender, thyroid_tcga_tumour_counts, tcga.geneIDs.annot)
write.table(diff_expr_tcga_tumour_table, file="./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_tumour_table %>% filter(FDR < 0.05) %>% nrow()
#128

diff_expr_tcga_tumour_table <- noiseq_diff_expression(thyroid_tcga_tumour_counts, thyroid_tcga_tumour_meta[, "gender", drop=F], "gender", tcga.geneIDs.annot)
write.table(diff_expr_tcga_tumour_table, file="./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_noiseqbio.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_tumour_table %>% filter(FDR < 0.05) %>% nrow()
#25

diff_expr_tcga_tumour_table <- wilcoxon_diff_expression(rownames(thyroid_tcga_tumour_meta[thyroid_tcga_tumour_meta$gender == "female", ]), rownames(thyroid_tcga_tumour_meta[thyroid_tcga_tumour_meta$gender == "male", ]), log2(thyroid_tcga_tumour_fpkm+1), tcga.geneIDs.annot)
write.table(diff_expr_tcga_tumour_table, file="./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_wilcoxon.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_tumour_table %>% filter(FDR < 0.05) %>% nrow()
#60

diff_expr_tcga_tumour_table <- limma_diff_expression(design_tcga_tumour, thyroid_tcga_tumour_meta$gender, thyroid_tcga_tumour_counts, tcga.geneIDs.annot)
write.table(diff_expr_tcga_tumour_table, file="./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_tcga_tumour_table %>% filter(FDR < 0.05) %>% nrow()
#81




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
#26




# -- Normal - GTEx

thyroid_gtex_meta <- as.data.frame(na.exclude(thyroid_gtex_meta))
thyroid_gtex_counts <- thyroid_tcga_gtex_counts[, colnames(thyroid_tcga_gtex_counts) %in% rownames(thyroid_gtex_meta)]
thyroid_gtex_fpkm <- thyroid_tcga_gtex_fpkm[, colnames(thyroid_tcga_gtex_fpkm) %in% rownames(thyroid_gtex_meta)]


# reorder datasets
thyroid_gtex_counts <- thyroid_gtex_counts[, order(colnames(thyroid_gtex_counts))]
thyroid_gtex_fpkm <- thyroid_gtex_fpkm[, order(colnames(thyroid_gtex_fpkm))]
thyroid_gtex_meta <- thyroid_gtex_meta[order(rownames(thyroid_gtex_meta)), ]


thyroid_gtex_counts %>% colnames %>% all.equal(rownames(thyroid_gtex_meta))
thyroid_gtex_fpkm %>% colnames %>% all.equal(rownames(thyroid_gtex_meta))
#TRUE


gender_code <- c("1" = "male", "2" = "female")
thyroid_gtex_meta$GENDER <- as.character(thyroid_gtex_meta$GENDER)
thyroid_gtex_meta$GENDER <- gender_code[thyroid_gtex_meta$GENDER]
thyroid_gtex_meta$ETHNCTY <- as.character(thyroid_gtex_meta$ETHNCTY)
thyroid_gtex_meta$MHCANCERNM <- as.character(thyroid_gtex_meta$MHCANCERNM)


# define design matrix
# adjust for covariates
design_gtex <- model.matrix(~ SMRIN + AGE + ETHNCTY + MHCANCERNM + SMCENTER + SMTSTPTREF + SMNABTCHT + SMTSISCH + GENDER, data = thyroid_gtex_meta)

diff_expr_gtex_table <- edgeR_diff_expression(design_gtex, thyroid_gtex_meta$gender, thyroid_gtex_counts, tcga.geneIDs.annot)
write.table(diff_expr_gtex_table, file="./files/diff_expr_thyroid_normal_gtex_maleVSfemale_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_gtex_table %>% filter(FDR < 0.05) %>% nrow()
#691

diff_expr_gtex_table <- noiseq_diff_expression(thyroid_gtex_counts, thyroid_gtex_meta[, "GENDER", drop=F], "GENDER", tcga.geneIDs.annot)
write.table(diff_expr_gtex_table, file="./files/diff_expr_thyroid_normal_gtex_maleVSfemale_noiseqbio.txt", sep = "\t", quote=F, row.names=F)

diff_expr_gtex_table %>% filter(FDR < 0.05) %>% nrow()
#57

diff_expr_gtex_table <- wilcoxon_diff_expression(rownames(thyroid_gtex_meta[thyroid_gtex_meta$GENDER == "female", ]), rownames(thyroid_gtex_meta[thyroid_gtex_meta$GENDER == "male", ]), log2(thyroid_gtex_fpkm+1), tcga.geneIDs.annot)
write.table(diff_expr_gtex_table, file="./files/diff_expr_thyroid_normal_gtex_maleVSfemale_wilcoxon.txt", sep = "\t", quote=F, row.names=F)

diff_expr_gtex_table %>% filter(FDR < 0.05) %>% nrow()
#869

diff_expr_gtex_table <- limma_diff_expression(design_gtex, thyroid_gtex_meta$gender, thyroid_gtex_counts, tcga.geneIDs.annot)
write.table(diff_expr_gtex_table, file="./files/diff_expr_thyroid_normal_gtex_maleVSfemale_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_gtex_table %>% filter(FDR < 0.05) %>% nrow()
#525


# write datasets
datasets <- thyroid_gtex_fpkm %>%
	rownames_to_column(var = "geneID") %>%
	as_tibble() %>%
	inner_join(tcga.geneIDs.annot[, c("geneID", "geneName")]) %>%
	select(geneID, geneName, everything()) %>%
	select(-geneID)

females <- datasets %>%
	select_at(.vars = c("geneName", rownames(thyroid_gtex_meta[thyroid_gtex_meta$GENDER == "female", ])))

males <- datasets %>%
	select_at(.vars = c("geneName", rownames(thyroid_gtex_meta[thyroid_gtex_meta$GENDER == "male", ])))

write_tsv(females, "./files/thyroid_gtex_females_fpkm.txt")
write_tsv(males, "./files/thyroid_gtex_males_fpkm.txt")


save(list=ls(), file="r_workspaces/tcga_gtex_thyroid_diffExpr_malesVSfemales.RData")
