# Understanding Gender Differential Susceptibility in Cancer


# Tumour vs Normal differential expression by gender
# TCGA + GTEx data


library(tidyverse)


source(file = "./r_scripts/utils.R")



# -- Thyroid


# load datasets
thyroid_tcga_counts <- read.delim("./files/thyroid_tcga_gtex_counts.txt", row.names = c(1), stringsAsFactors=F)
thyroid_tcga_fpkm <- read.delim("./files/thyroid_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
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
thyroid_tcga_fpkm <- thyroid_tcga_fpkm[, substring(colnames(thyroid_tcga_fpkm), 1, 4) != "GTEX"]


# reorder datasets
thyroid_tcga_counts <- thyroid_tcga_counts[, order(colnames(thyroid_tcga_counts))]
thyroid_tcga_fpkm <- thyroid_tcga_fpkm[, order(colnames(thyroid_tcga_fpkm))]
thyroid_tcga_meta <- thyroid_tcga_meta[order(rownames(thyroid_tcga_meta)), ]


# remove metastatic samples
thyroid_tcga_counts <- thyroid_tcga_counts[, !colnames(thyroid_tcga_counts) %in% rownames(thyroid_tcga_meta[(thyroid_tcga_meta$sample_type_id == "metastatic"), ])]
thyroid_tcga_fpkm <- thyroid_tcga_fpkm[, !colnames(thyroid_tcga_fpkm) %in% rownames(thyroid_tcga_meta[(thyroid_tcga_meta$sample_type_id == "metastatic"), ])]
thyroid_tcga_meta <- thyroid_tcga_meta[!(thyroid_tcga_meta$sample_type_id == "metastatic"), ]



thyroid_tcga_counts %>% colnames %>% all.equal(rownames(thyroid_tcga_meta))
thyroid_tcga_fpkm %>% colnames %>% all.equal(rownames(thyroid_tcga_meta))
#TRUE





# -- Males

thyroid_males_meta <- thyroid_tcga_meta[(thyroid_tcga_meta$gender == "male"), ]
thyroid_males_counts <- thyroid_tcga_counts[ ,colnames(thyroid_tcga_counts) %in% rownames(thyroid_males_meta)]
thyroid_males_fpkm <- thyroid_tcga_fpkm[ ,colnames(thyroid_tcga_fpkm) %in% rownames(thyroid_males_meta)]

thyroid_males_counts %>% colnames %>% all.equal(rownames(thyroid_males_meta))
thyroid_males_fpkm %>% colnames %>% all.equal(rownames(thyroid_males_meta))
#TRUE


# define design matrix
# adjust for covariates
design_males <- model.matrix(~ race + ethnicity + age_at_diagnosis + tss + portion + plate + sample_type_id, data = thyroid_males_meta)

diff_expr_males_table <- edgeR_diff_expression(design_males, thyroid_males_meta$sample_type_id, thyroid_males_counts, tcga.geneIDs.annot)
write.table(diff_expr_males_table, file="./files/diff_expr_thyroid_males_tcga_tumourVSnormal_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow()
#1504

diff_expr_males_table <- noiseq_diff_expression(thyroid_males_counts, thyroid_males_meta[, "sample_type_id", drop=F], "sample_type_id", tcga.geneIDs.annot)
write.table(diff_expr_males_table, file="./files/diff_expr_thyroid_males_tcga_tumourVSnormal_noiseqbio.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(FDR < 0.05, abs(log2FC)>1) %>% nrow()
#1402

diff_expr_males_table <- wilcoxon_diff_expression(rownames(thyroid_males_meta[thyroid_males_meta$sample_type_id == "tumour", ]), rownames(thyroid_males_meta[thyroid_males_meta$sample_type_id == "normal", ]), log2(thyroid_males_fpkm+1), tcga.geneIDs.annot)
write.table(diff_expr_males_table, file="./files/diff_expr_thyroid_males_tcga_tumourVSnormal_wilcoxon.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(FDR < 0.05, abs(log2FC)>1) %>% nrow()
#568

diff_expr_males_table <- limma_diff_expression(design_males, thyroid_males_meta$sample_type_id, thyroid_males_counts, tcga.geneIDs.annot)
write.table(diff_expr_males_table, file="./files/diff_expr_thyroid_males_tcga_tumourVSnormal_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow()
#1471




# -- Females


thyroid_females_meta <- thyroid_tcga_meta[(thyroid_tcga_meta$gender == "female"), ]
thyroid_females_counts <- thyroid_tcga_counts[ ,colnames(thyroid_tcga_counts) %in% rownames(thyroid_females_meta)]
thyroid_females_fpkm <- thyroid_tcga_fpkm[ ,colnames(thyroid_tcga_fpkm) %in% rownames(thyroid_females_meta)]

thyroid_females_counts %>% colnames %>% all.equal(rownames(thyroid_females_meta))
thyroid_females_fpkm %>% colnames %>% all.equal(rownames(thyroid_females_meta))
#TRUE


# define design matrix
# adjust for covariates
design_females <- model.matrix(~ race + ethnicity + age_at_diagnosis + tss + portion + plate + sample_type_id, data = thyroid_females_meta)

diff_expr_females_table <- edgeR_diff_expression(design_females, thyroid_females_meta$sample_type_id, thyroid_females_counts, tcga.geneIDs.annot)
write.table(diff_expr_females_table, file="./files/diff_expr_thyroid_females_tcga_tumourVSnormal_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow()
#1148

diff_expr_females_table <- noiseq_diff_expression(thyroid_females_counts, thyroid_females_meta[, "sample_type_id", drop=F], "sample_type_id", tcga.geneIDs.annot)
write.table(diff_expr_females_table, file="./files/diff_expr_thyroid_females_tcga_tumourVSnormal_noiseqbio.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(FDR < 0.05, abs(log2FC)>1) %>% nrow()
#1136

diff_expr_females_table <- wilcoxon_diff_expression(rownames(thyroid_females_meta[thyroid_females_meta$sample_type_id == "tumour", ]), rownames(thyroid_females_meta[thyroid_females_meta$sample_type_id == "normal", ]), log2(thyroid_females_fpkm+1), tcga.geneIDs.annot)
write.table(diff_expr_females_table, file="./files/diff_expr_thyroid_females_tcga_tumourVSnormal_wilcoxon.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(FDR < 0.05, abs(log2FC)>1) %>% nrow()
#520

diff_expr_females_table <- limma_diff_expression(design_females, thyroid_females_meta$sample_type_id, thyroid_females_counts, tcga.geneIDs.annot)
write.table(diff_expr_females_table, file="./files/diff_expr_thyroid_females_tcga_tumourVSnormal_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow()
#1202




save(list=ls(), file="./r_workspaces/tcga_thyroid_diffExpr_tumourVSnormal.RData")
