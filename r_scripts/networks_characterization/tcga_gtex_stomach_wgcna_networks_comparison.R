# Understanding Gender Differential Susceptibility in Cancer


# GENE-CORRELATION NETWORK ANALYSES USING WGCNA PACKAGE
# TCGA + GTEx data




# -- Stomach

# load R workspace
load("./r_workspaces/tcga_gtex_stomach_wgcna_networks.RData")

source("./r_scripts/utils.R")



library(WGCNA)
library(tidyverse)
library(gplots)





#Males tumour vs females tumour
males_females_tumour_overlapTable <- overlap_table_wgcna(
    males_tumour_network,
    females_tumour_network,
    "Stomach Tumours (TCGA) - Men",
    "Stomach Tumours (TCGA) - Women",
    "./plots/wgcna_networks/tcga_stomach_tumours_wgcna_networks_comparison_malesVSfemales.pdf", 20, 20)

males_females_tumour_comparison <- networks_comparison(
    males_tumour,
    females_tumour,
    males_tumour_network,
    females_tumour_network,
    males_females_tumour_overlapTable)


pdf(file="./plots/wgcna_networks/tcga_stomach_tumours_wgcna_networks_comparison_malesVSfemales2.pdf", w=12, h=8)
heatmap.2( x=as.matrix(males_females_tumour_comparison$corTable2),
    Rowv=T,
    Colv=T,
    dendrogram="both",
    scale = "none",
    trace = "none",
    key = FALSE,
    margins=c(7,7),
    col = c("grey", "blue", "green"),
    main = "Stomach Tumours (TCGA): Men modules vs Women modules",
    xlab = "Women modules",
    ylab = "Men modules")
legend("topleft", fill = c("grey", "blue", "green"), legend = c("Pair of modules not conserved", "Pair of modules conserved with\nloss of topology", "Pair of modules conserved with\npreservation of topology"), y.intersp = 2 )
dev.off()


males_females_tumour_table1 <- -log10(males_females_tumour_overlapTable$pTable[rownames(males_females_tumour_overlapTable$pTable) != "grey", colnames(males_females_tumour_overlapTable$pTable) != "grey"])
males_females_tumour_table1[is.infinite(males_females_tumour_table1)] <- 400
males_females_tumour_table1 <- males_females_tumour_table1 %>%
    #rowSums() %>%
    apply(., 1, max) %>%
    as.data.frame() %>%
    setNames(., "p_value")
males_females_tumour_table1 <- bind_cols(module = rownames(males_females_tumour_table1), males_females_tumour_table1, network = rep("males", 34)) %>%
    arrange(p_value) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(module = fct_reorder(module, order(p_value)))


males_females_tumour_table2 <- -log10(males_females_tumour_overlapTable$pTable[rownames(males_females_tumour_overlapTable$pTable) != "grey", colnames(males_females_tumour_overlapTable$pTable) != "grey"])
males_females_tumour_table2[is.infinite(males_females_tumour_table2)] <- 400
males_females_tumour_table2 <- males_females_tumour_table2 %>%
    #colSums() %>%
    apply(., 2, max) %>%
    as.data.frame() %>%
    setNames(., "p_value")
males_females_tumour_table2 <- bind_cols(module = rownames(males_females_tumour_table2), males_females_tumour_table2, network = rep("females", 22)) %>%
    arrange(p_value) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(module = fct_reorder(module, order(p_value)))


males_females_tumour_table <- bind_rows(males_females_tumour_table1, males_females_tumour_table2) %>%
    group_by(network) %>%
    arrange(p_value) %>%
    ungroup()


males_females_tumour_table_plot1 <- ggplot(data=males_females_tumour_table1, mapping=aes(x=module, y=p_value, fill=module)) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13, angle=45, vjust = 1, hjust=1),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    scale_fill_manual(values = as.character(males_females_tumour_table1$module)) +
    labs(x = "Module", y = "Maximum p-value (-log10)", title="Males")
ggsave(filename="tcga_stomach_tumours_wgcna_networks_comparison_malesVSfemales_bp1.png", plot=males_females_tumour_table_plot1, path = "./plots/wgcna_networks/")
unlink("tcga_stomach_tumours_wgcna_networks_comparison_malesVSfemales_bp1.png")



males_females_tumour_table_plot2 <- ggplot(data=males_females_tumour_table2, mapping=aes(x=module, y=p_value, fill=module)) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13, angle=45, vjust = 1, hjust=1),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    scale_fill_manual(values = as.character(males_females_tumour_table2$module)) +
    labs(x = "Module", y = "Maximum p-value (-log10)", title="Females")
ggsave(filename="tcga_stomach_tumours_wgcna_networks_comparison_malesVSfemales_bp2.png", plot=males_females_tumour_table_plot2, path = "./plots/wgcna_networks/")
unlink("tcga_stomach_tumours_wgcna_networks_comparison_malesVSfemales_bp2.png")





#Males normal vs females normal
males_females_normal_overlapTable <- overlap_table_wgcna(
    males_normal_network,
    females_normal_network,
    "Stomach Normal (GTEx) - Men",
    "Stomach Normal (GTEx) - Women",
    "./plots/wgcna_networks/gtex_stomach_wgcna_networks_comparison_malesVSfemales.pdf", 12, 10)

males_females_normal_comparison <- networks_comparison(
    males_normal,
    females_normal,
    males_normal_network,
    females_normal_network,
    males_females_normal_overlapTable)


pdf(file="./plots/wgcna_networks/gtex_stomach_wgcna_networks_comparison_malesVSfemales2.pdf", w=12, h=8)
heatmap.2( x=as.matrix(males_females_normal_comparison$corTable2),
    Rowv=T,
    Colv=T,
    dendrogram="both",
    scale = "none",
    trace = "none",
    key = FALSE,
    margins=c(7,7),
    col = c("grey", "blue", "green"),
    main = "Stomach Normal (GTEx): Males modules vs Females modules",
    xlab = "Women modules",
    ylab = "Men modules")
legend("topleft", fill = c("grey", "blue", "green"), legend = c("Pair of modules not conserved", "Pair of modules conserved with\nloss of topology", "Pair of modules conserved with\npreservation of topology"), y.intersp = 2 )
dev.off()



males_females_normal_table1 <- -log10(males_females_normal_overlapTable$pTable[rownames(males_females_normal_overlapTable$pTable) != "grey", colnames(males_females_normal_overlapTable$pTable) != "grey"])
males_females_normal_table1[is.infinite(males_females_normal_table1)] <- 400
males_females_normal_table1 <- males_females_normal_table1 %>%
    #rowSums() %>%
    apply(., 1, max) %>%
    as.data.frame() %>%
    setNames(., "p_value")
males_females_normal_table1 <- bind_cols(module = rownames(males_females_normal_table1), males_females_normal_table1, network = rep("males", 8)) %>%
    arrange(p_value) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(module = fct_reorder(module, order(p_value)))


males_females_normal_table2 <- -log10(males_females_normal_overlapTable$pTable[rownames(males_females_normal_overlapTable$pTable) != "grey", colnames(males_females_normal_overlapTable$pTable) != "grey"])
males_females_normal_table2[is.infinite(males_females_normal_table2)] <- 400
males_females_normal_table2 <- males_females_normal_table2 %>%
    #colSums() %>%
    apply(., 2, max) %>%
    as.data.frame() %>%
    setNames(., "p_value")
males_females_normal_table2 <- bind_cols(module = rownames(males_females_normal_table2), males_females_normal_table2, network = rep("females", 12)) %>%
    arrange(p_value) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(module = fct_reorder(module, order(p_value)))


males_females_normal_table <- bind_rows(males_females_normal_table1, males_females_normal_table2) %>%
    group_by(network) %>%
    arrange(p_value) %>%
    ungroup()


males_females_normal_table_plot1 <- ggplot(data=males_females_normal_table1, mapping=aes(x=module, y=p_value, fill=module)) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13, angle=45, vjust = 1, hjust=1),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    scale_fill_manual(values = as.character(males_females_normal_table1$module)) +
    labs(x = "Module", y = "Maximum p-value (-log10)", title="Males")
ggsave(filename="gtex_stomach_wgcna_networks_comparison_malesVSfemales_bp1.png", plot=males_females_normal_table_plot1, path = "./plots/wgcna_networks/")
unlink("gtex_stomach_wgcna_networks_comparison_malesVSfemales_bp1.png")



males_females_normal_table_plot2 <- ggplot(data=males_females_normal_table2, mapping=aes(x=module, y=p_value, fill=module)) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(axis.title = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=13),
        axis.text.x = element_text(colour="black", size=13, angle=45, vjust = 1, hjust=1),
        plot.title = element_text(colour="black", size=15),
        legend.position = "none") +
    scale_fill_manual(values = as.character(males_females_normal_table2$module)) +
    labs(x = "Module", y = "Maximum p-value (-log10)", title="Females")
ggsave(filename="gtex_stomach_wgcna_networks_comparison_malesVSfemales_bp2.png", plot=males_females_normal_table_plot2, path = "./plots/wgcna_networks/")
unlink("gtex_stomach_wgcna_networks_comparison_malesVSfemales_bp2.png")