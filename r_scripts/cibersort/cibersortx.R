library(tidyverse)

# thyroid immune cell fraction between males and females

males <- read_tsv("./files/thyroid_gtex_males_cibersort_cell_fractions.txt")
females <- read_tsv("./files/thyroid_gtex_females_cibersort_cell_fractions.txt")

males <- males %>%
	pivot_longer(-Mixture, names_to="cell", values_to="fraction") %>%
	rename(sample=Mixture) %>%
	mutate(gender = "males")

females <- females %>%
	pivot_longer(-Mixture, names_to="cell", values_to="fraction") %>%
	rename(sample=Mixture) %>%
	mutate(gender = "females")

dat <- bind_rows(males, females) %>%
	filter(!cell %in% c("Correlation", "P-value", "RMSE"))

plot <- dat %>%
	ggplot(mapping=aes(x = cell, y = fraction, fill = gender)) +
	geom_boxplot() +
	coord_flip() +
	ggpubr::stat_compare_means(method = "wilcox.test", label.y=0.75) +
	scale_y_continuous(limits = c(NA, 1))


# thyroid immune cell gene expression between males and females

readF <- function(file, folder){

	open_f <- str_c(folder, file, sep = "/")

	cell <- str_extract(file, "[A-Za-z0-9]+_Window40.txt")
	cell <- str_replace(cell, "_Window40.txt", "")

	dat <- read_tsv(open_f) %>%
		pivot_longer(-GeneSymbol, names_to="sample", values_to="expression") %>%
		mutate(cell = cell)

	dat
}

wilcx <- function(df){

	g1 <- df %>%
		filter(gender == "females") %>%
		pull(expression)

	g2 <- df %>%
		filter(gender == "males") %>%
		pull(expression)

	if(all(g1 == 1) | all(is.na(g1)) | all(g2 == 1) | all(is.na(g2))){
		res <- tibble(mvar = NA, fvar = NA, p.val = NA, log2FC = NA)
	} else {
		p.val <- wilcox.test(g2, g1)$p.value
		log2FC <- log2((median(g1)+1)/(median(g2)+1))
		res <- tibble(mvar = var(g2), fvar = var(g1), p.val = p.val, log2FC = log2FC)
	}
	res
}


female_files <- c("CIBERSORTxHiRes_Job9_Bcells_Window40.txt", "CIBERSORTxHiRes_Job9_Dendriticcells_Window40.txt", "CIBERSORTxHiRes_Job9_Eosinophils_Window40.txt", "CIBERSORTxHiRes_Job9_Mastcells_Window40.txt", "CIBERSORTxHiRes_Job9_Monocytes_Window40.txt", "CIBERSORTxHiRes_Job9_Neutrophils_Window40.txt", "CIBERSORTxHiRes_Job9_NKcells_Window40.txt", "CIBERSORTxHiRes_Job9_Plasmacells_Window40.txt", "CIBERSORTxHiRes_Job9_TcellsCD4_Window40.txt", "CIBERSORTxHiRes_Job9_TcellsCD8_Window40.txt")
male_files <- c("CIBERSORTxHiRes_Job10_Bcells_Window40.txt", "CIBERSORTxHiRes_Job10_Dendriticcells_Window40.txt", "CIBERSORTxHiRes_Job10_Eosinophils_Window40.txt", "CIBERSORTxHiRes_Job10_Mastcells_Window40.txt", "CIBERSORTxHiRes_Job10_Monocytes_Window40.txt", "CIBERSORTxHiRes_Job10_Neutrophils_Window40.txt", "CIBERSORTxHiRes_Job10_NKcells_Window40.txt", "CIBERSORTxHiRes_Job10_Plasmacells_Window40.txt", "CIBERSORTxHiRes_Job10_TcellsCD4_Window40.txt", "CIBERSORTxHiRes_Job10_TcellsCD8_Window40.txt")


females <- map_dfr(.x = female_files, .f = readF, folder = "./files/thyroid_gtex_females_cibersort_immune_expression") %>%
	mutate(gender = "females")

males <- map_dfr(.x = male_files, .f = readF, folder = "./files/thyroid_gtex_males_cibersort_immune_expression") %>%
	mutate(gender = "males")


dat <- bind_rows(females, males) %>%
	group_by(GeneSymbol, cell) %>%
	nest() %>%
	ungroup() %>%
	mutate(p_values = map(.x = data, .f = wilcx)) %>%
	unnest(p_values) %>%
	filter(!is.na(p.val)) %>%
	mutate(p_adjust = p.adjust(p.val, method = "BH"))

plot <- dat %>%
	ggplot(mapping=aes(x = log2FC, y = -log10(p_adjust))) +
		geom_point() +
		geom_hline(yintercept=-log10(0.05), linetype = "dashed") +
		geom_vline(xintercept=0, linetype="dashed") +
		facet_wrap(~ cell)




# stomach immune cell gene expression between males and females

male_files <- c("CIBERSORTxHiRes_Job11_Bcells_Window40.txt", "CIBERSORTxHiRes_Job11_Dendriticcells_Window40.txt", "CIBERSORTxHiRes_Job11_Eosinophils_Window40.txt", "CIBERSORTxHiRes_Job11_Mastcells_Window40.txt", "CIBERSORTxHiRes_Job11_Monocytes_Window40.txt", "CIBERSORTxHiRes_Job11_Neutrophils_Window40.txt", "CIBERSORTxHiRes_Job11_NKcells_Window40.txt", "CIBERSORTxHiRes_Job11_Plasmacells_Window40.txt", "CIBERSORTxHiRes_Job11_TcellsCD4_Window40.txt", "CIBERSORTxHiRes_Job11_TcellsCD8_Window40.txt")
female_files <- c("CIBERSORTxHiRes_Job12_Bcells_Window40.txt", "CIBERSORTxHiRes_Job12_Dendriticcells_Window40.txt", "CIBERSORTxHiRes_Job12_Eosinophils_Window40.txt", "CIBERSORTxHiRes_Job12_Mastcells_Window40.txt", "CIBERSORTxHiRes_Job12_Monocytes_Window40.txt", "CIBERSORTxHiRes_Job12_Neutrophils_Window40.txt", "CIBERSORTxHiRes_Job12_NKcells_Window40.txt", "CIBERSORTxHiRes_Job12_Plasmacells_Window40.txt", "CIBERSORTxHiRes_Job12_TcellsCD4_Window40.txt", "CIBERSORTxHiRes_Job12_TcellsCD8_Window40.txt")


females <- map_dfr(.x = female_files, .f = readF, folder = "./files/stomach_gtex_females_cibersort_immune_expression") %>%
	mutate(gender = "females")

males <- map_dfr(.x = male_files, .f = readF, folder = "./files/stomach_gtex_males_cibersort_immune_expression") %>%
	mutate(gender = "males")


dat <- bind_rows(females, males) %>%
	group_by(GeneSymbol, cell) %>%
	nest() %>%
	ungroup() %>%
	mutate(p_values = map(.x = data, .f = wilcx)) %>%
	unnest(p_values) %>%
	filter(!is.na(p.val)) %>%
	mutate(p_adjust = p.adjust(p.val, method = "BH"))

plot <- dat %>%
	ggplot(mapping=aes(x = log2FC, y = -log10(p_adjust))) +
		geom_point() +
		geom_hline(yintercept=-log10(0.05), linetype = "dashed") +
		geom_vline(xintercept=0, linetype="dashed") +
		facet_wrap(~ cell)
