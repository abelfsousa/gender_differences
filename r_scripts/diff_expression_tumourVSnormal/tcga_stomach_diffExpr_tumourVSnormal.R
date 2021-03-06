# Understanding Gender Differential Susceptibility in Cancer


# Tumour vs Normal differential expression by gender
# TCGA + GTEx data


library(tidyverse)


source(file = "./r_scripts/utils.R")



# -- Stomach


# load datasets
stomach_tcga_counts <- read.delim("./files/stomach_tcga_gtex_counts.txt", row.names = c(1), stringsAsFactors=F)
stomach_tcga_fpkm <- read.delim("./files/stomach_tcga_gtex_fpkm.txt", row.names = c(1), stringsAsFactors=F)
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
stomach_tcga_fpkm <- stomach_tcga_fpkm[, substring(colnames(stomach_tcga_fpkm), 1, 4) != "GTEX"]


# reorder datasets
stomach_tcga_counts <- stomach_tcga_counts[, order(colnames(stomach_tcga_counts))]
stomach_tcga_fpkm <- stomach_tcga_fpkm[, order(colnames(stomach_tcga_fpkm))]
stomach_tcga_meta <- stomach_tcga_meta[order(rownames(stomach_tcga_meta)), ]



stomach_tcga_counts %>% colnames %>% all.equal(rownames(stomach_tcga_meta))
stomach_tcga_fpkm %>% colnames %>% all.equal(rownames(stomach_tcga_meta))
#TRUE





# -- Males

stomach_males_meta <- stomach_tcga_meta[(stomach_tcga_meta$gender == "male"), ]
stomach_males_counts <- stomach_tcga_counts[ ,colnames(stomach_tcga_counts) %in% rownames(stomach_males_meta)]
stomach_males_fpkm <- stomach_tcga_fpkm[ ,colnames(stomach_tcga_fpkm) %in% rownames(stomach_males_meta)]

stomach_males_counts %>% colnames %>% all.equal(rownames(stomach_males_meta))
stomach_males_fpkm %>% colnames %>% all.equal(rownames(stomach_males_meta))
#TRUE


# define design matrix
# adjust for covariates
design_males <- model.matrix(~ race + ethnicity + age_at_diagnosis + portion + plate + sample_type_id, data = stomach_males_meta)

diff_expr_males_table <- edgeR_diff_expression(design_males, stomach_males_meta$sample_type_id, stomach_males_counts, tcga.geneIDs.annot)
write.table(diff_expr_males_table, file="./files/diff_expr_stomach_males_tcga_tumourVSnormal_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow()
#2294

diff_expr_males_table <- noiseq_diff_expression(stomach_males_counts, stomach_males_meta[, "sample_type_id", drop=F], "sample_type_id", tcga.geneIDs.annot)
write.table(diff_expr_males_table, file="./files/diff_expr_stomach_males_tcga_tumourVSnormal_noiseqbio.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(FDR < 0.05, abs(log2FC)>1) %>% nrow()
#1655

diff_expr_males_table <- wilcoxon_diff_expression(rownames(stomach_males_meta[stomach_males_meta$sample_type_id == "tumour", ]), rownames(stomach_males_meta[stomach_males_meta$sample_type_id == "normal", ]), log2(stomach_males_fpkm+1), tcga.geneIDs.annot)
write.table(diff_expr_males_table, file="./files/diff_expr_stomach_males_tcga_tumourVSnormal_wilcoxon.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(FDR < 0.05, abs(log2FC)>1) %>% nrow()
#495

diff_expr_males_table <- limma_diff_expression(design_males, stomach_males_meta$sample_type_id, stomach_males_counts, tcga.geneIDs.annot)
write.table(diff_expr_males_table, file="./files/diff_expr_stomach_males_tcga_tumourVSnormal_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_males_table %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow()
#2367




# -- Females


stomach_females_meta <- stomach_tcga_meta[(stomach_tcga_meta$gender == "female"), ]
stomach_females_counts <- stomach_tcga_counts[ ,colnames(stomach_tcga_counts) %in% rownames(stomach_females_meta)]
stomach_females_fpkm <- stomach_tcga_fpkm[ ,colnames(stomach_tcga_fpkm) %in% rownames(stomach_females_meta)]

stomach_females_counts %>% colnames %>% all.equal(rownames(stomach_females_meta))
stomach_females_fpkm %>% colnames %>% all.equal(rownames(stomach_females_meta))
#TRUE


# define design matrix
# adjust for covariates
design_females <- model.matrix(~ race + ethnicity + age_at_diagnosis + portion + plate + sample_type_id, data = stomach_females_meta)

diff_expr_females_table <- edgeR_diff_expression(design_females, stomach_females_meta$sample_type_id, stomach_females_counts, tcga.geneIDs.annot)
write.table(diff_expr_females_table, file="./files/diff_expr_stomach_females_tcga_tumourVSnormal_edgeR.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow()
#1856

diff_expr_females_table <- noiseq_diff_expression(stomach_females_counts, stomach_females_meta[, "sample_type_id", drop=F], "sample_type_id", tcga.geneIDs.annot)
write.table(diff_expr_females_table, file="./files/diff_expr_stomach_females_tcga_tumourVSnormal_noiseqbio.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(FDR < 0.05, abs(log2FC)>1) %>% nrow()
#836

diff_expr_females_table <- wilcoxon_diff_expression(rownames(stomach_females_meta[stomach_females_meta$sample_type_id == "tumour", ]), rownames(stomach_females_meta[stomach_females_meta$sample_type_id == "normal", ]), log2(stomach_females_fpkm+1), tcga.geneIDs.annot)
write.table(diff_expr_females_table, file="./files/diff_expr_stomach_females_tcga_tumourVSnormal_wilcoxon.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(FDR < 0.05, abs(log2FC)>1) %>% nrow()
#566

diff_expr_females_table <- limma_diff_expression(design_females, stomach_females_meta$sample_type_id, stomach_females_counts, tcga.geneIDs.annot)
write.table(diff_expr_females_table, file="./files/diff_expr_stomach_females_tcga_tumourVSnormal_limma.txt", sep = "\t", quote=F, row.names=F)

diff_expr_females_table %>% filter(FDR < 0.05, abs(logFC)>1) %>% nrow()
#1841





save(list=ls(), file="./r_workspaces/tcga_stomach_diffExpr_tumourVSnormal.RData")
