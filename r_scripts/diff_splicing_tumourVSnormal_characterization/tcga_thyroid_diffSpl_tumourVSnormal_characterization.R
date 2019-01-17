# Understanding Gender Differential Susceptibility in Cancer


# Characterization of differentially spliced exons between tumour and normal samples in Gender




library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)
library(viridis)



# -- Thyroid

# load datasets
diff_spl_males_exons <- read_tsv("./files/diff_splicing/THCA/THCA.male.TvsN.limma.diff.splicing.exon-level.txt")
diff_spl_males_genes <- read_tsv("./files/diff_splicing/THCA/THCA.male.TvsN.limma.diff.splicing.gene-level-Simes.txt")

diff_spl_females_exons <- read_tsv("./files/diff_splicing/THCA/THCA.female.TvsN.limma.diff.splicing.exon-level.txt")
diff_spl_females_genes <- read_tsv("./files/diff_splicing/THCA/THCA.female.TvsN.limma.diff.splicing.gene-level-Simes.txt")


# significant DSGs between males and females
males_signf_dsgs <- diff_spl_males_genes %>%
    filter(FDR <= 0.05)

females_signf_dsgs <- diff_spl_females_genes %>%
    filter(FDR <= 0.05)

males_signf_dsgs <- males_signf_dsgs %>%
    mutate(state = if_else(geneName %in% females_signf_dsgs$geneName, "common", "male_specific"))

females_signf_dsgs <- females_signf_dsgs %>%
    mutate(state = if_else(geneName %in% males_signf_dsgs$geneName, "common", "female_specific"))


#diff_spl_males_exons %>% filter(geneName %in% males_signf_dsgs$geneName) %>% filter(FDR < 0.05)
#diff_spl_females_exons %>% filter(geneName %in% females_signf_dsgs$geneName) %>% filter(FDR < 0.05)


# all significant DEGs between males and females
males_females_signf_dsgs <- full_join(
    males_signf_dsgs, females_signf_dsgs,
    by=c("gene", "geneName", "chrom", "NExons", "P.Value", "FDR", "state"))
#    dplyr::rename(P.Value.males = P.Value.x, P.Value.females = P.Value.y, FDR.males = FDR.x, FDR.females = FDR.y)
write.table(males_females_signf_dsgs, "./files/thyroid_males_females_signf_dsgs.txt", sep="\t", quote=F, row.names=F)



# Venn Diagram
males_females_signf_dsgs_venn <- draw.pairwise.venn(
    4,
    3,
    0,
    scaled = T,
    category = c("", ""),
    lty = rep("blank", 2),
    fill = c("light blue", "pink"),
    alpha = rep(0.5, 2),
    cex = rep(2, 3),
    cat.cex = rep(1.5, 2))
ggsave(filename="thyroid_males_females_signf_dsgs_venn.png", plot=males_females_signf_dsgs_venn, path = "./plots/diff_splicing_tumourVSnormal/", width=3, height=3)
ggsave(filename="thyroid_males_females_signf_dsgs_venn.pdf", plot=males_females_signf_dsgs_venn, path = "./plots/diff_splicing_tumourVSnormal/", width=3, height=3)
unlink("thyroid_males_females_signf_dsgs_venn.png")
unlink("thyroid_males_females_signf_dsgs_venn.pdf")





# load stomach DEGs between tissue type
thyroid_males_females_signf_degs <- read_tsv("./files/thyroid_males_females_signf_degs.txt")

#males_females_signf_dsgs %>% filter(geneName %in% thyroid_males_females_signf_degs$geneName)
#thyroid_males_females_signf_degs %>% filter(geneName %in% males_females_signf_dsgs$geneName)
