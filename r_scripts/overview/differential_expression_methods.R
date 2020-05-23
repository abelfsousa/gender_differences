library(tidyverse)


#Stomach TCGA Male vs Female
edger <- read_tsv("./files/diff_expr_stomach_tumour_tcga_maleVSfemale_edgeR.txt") %>%
  filter(FDR < 0.05)
limma <- read_tsv("./files/diff_expr_stomach_tumour_tcga_maleVSfemale_limma.txt") %>%
  filter(FDR < 0.05)

(length(intersect(edger$genes, limma$genes))/length(union(edger$genes, limma$genes)))*100
(length(intersect(edger$genes, limma$genes))/nrow(edger))*100
(length(intersect(edger$genes, limma$genes))/nrow(limma))*100 #95%


#Stomach GTEx Male vs Female
edger <- read_tsv("./files/diff_expr_stomach_normal_gtex_maleVSfemale_edgeR.txt") %>%
  filter(FDR < 0.05)
limma <- read_tsv("./files/diff_expr_stomach_normal_gtex_maleVSfemale_limma.txt") %>%
  filter(FDR < 0.05)

(length(intersect(edger$genes, limma$genes))/length(union(edger$genes, limma$genes)))*100
(length(intersect(edger$genes, limma$genes))/nrow(edger))*100
(length(intersect(edger$genes, limma$genes))/nrow(limma))*100 #97%


#Thyroid TCGA Male vs Female
edger <- read_tsv("./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_edgeR.txt") %>%
  filter(FDR < 0.05)
limma <- read_tsv("./files/diff_expr_thyroid_tumour_tcga_maleVSfemale_limma.txt") %>%
  filter(FDR < 0.05)

(length(intersect(edger$genes, limma$genes))/length(union(edger$genes, limma$genes)))*100
(length(intersect(edger$genes, limma$genes))/nrow(edger))*100
(length(intersect(edger$genes, limma$genes))/nrow(limma))*100 #86%


#Thyroid GTEx Male vs Female
edger <- read_tsv("./files/diff_expr_thyroid_normal_gtex_maleVSfemale_edgeR.txt") %>%
  filter(FDR < 0.05)
limma <- read_tsv("./files/diff_expr_thyroid_normal_gtex_maleVSfemale_limma.txt") %>%
  filter(FDR < 0.05)

(length(intersect(edger$genes, limma$genes))/length(union(edger$genes, limma$genes)))*100
(length(intersect(edger$genes, limma$genes))/nrow(edger))*100
(length(intersect(edger$genes, limma$genes))/nrow(limma))*100 #92%


#Stomach males Tumour vs Normal
edger <- read_tsv("./files/diff_expr_stomach_males_tcga_tumourVSnormal_edgeR.txt") %>%
  filter(FDR < 0.05 & abs(logFC)>1)
limma <- read_tsv("./files/diff_expr_stomach_males_tcga_tumourVSnormal_limma.txt") %>%
  filter(FDR < 0.05 & abs(logFC)>1)

(length(intersect(edger$genes, limma$genes))/length(union(edger$genes, limma$genes)))*100
(length(intersect(edger$genes, limma$genes))/nrow(edger))*100
(length(intersect(edger$genes, limma$genes))/nrow(limma))*100 #85%


#Stomach females Tumour vs Normal
edger <- read_tsv("./files/diff_expr_stomach_females_tcga_tumourVSnormal_edgeR.txt") %>%
  filter(FDR < 0.05 & abs(logFC)>1)
limma <- read_tsv("./files/diff_expr_stomach_females_tcga_tumourVSnormal_limma.txt") %>%
  filter(FDR < 0.05 & abs(logFC)>1)

(length(intersect(edger$genes, limma$genes))/length(union(edger$genes, limma$genes)))*100
(length(intersect(edger$genes, limma$genes))/nrow(edger))*100
(length(intersect(edger$genes, limma$genes))/nrow(limma))*100 #86%


#Thyroid males Tumour vs Normal
edger <- read_tsv("./files/diff_expr_thyroid_males_tcga_tumourVSnormal_edgeR.txt") %>%
  filter(FDR < 0.05 & abs(logFC)>1)
limma <- read_tsv("./files/diff_expr_thyroid_males_tcga_tumourVSnormal_limma.txt") %>%
  filter(FDR < 0.05 & abs(logFC)>1)

(length(intersect(edger$genes, limma$genes))/length(union(edger$genes, limma$genes)))*100
(length(intersect(edger$genes, limma$genes))/nrow(edger))*100
(length(intersect(edger$genes, limma$genes))/nrow(limma))*100 #94%


#Thyroid females Tumour vs Normal
edger <- read_tsv("./files/diff_expr_thyroid_females_tcga_tumourVSnormal_edgeR.txt") %>%
  filter(FDR < 0.05 & abs(logFC)>1)
limma <- read_tsv("./files/diff_expr_thyroid_females_tcga_tumourVSnormal_limma.txt") %>%
  filter(FDR < 0.05 & abs(logFC)>1)

(length(intersect(edger$genes, limma$genes))/length(union(edger$genes, limma$genes)))*100
(length(intersect(edger$genes, limma$genes))/nrow(edger))*100
(length(intersect(edger$genes, limma$genes))/nrow(limma))*100 #87%
