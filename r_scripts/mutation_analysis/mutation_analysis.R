# Understanding Gender Differential Susceptibility in Cancer


# Comparison of mutation status of differentially expressed genes between sample groups


library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(data.table)


# load datasets


# differentially expressed genes
thca_TvsN <- read_tsv("./files/thyroid_males_females_signf_degs.txt")
thca_MvsF <- read_tsv("./files/thyroid_tumour_normal_signf_degs.txt")

stad_TvsN <- read_tsv("./files/stomach_males_females_signf_degs.txt")
stad_MvsF <- read_tsv("./files/stomach_tumour_normal_signf_degs.txt")



# tumour vs normal DEGs tsgs

# thyroid

# female-specific
thca_female_tsgs <- read_tsv("./files/degs_TvsN_stomach_thyroid_tsgs.txt") %>%
	filter(tissue == "Thyroid" & state == "female_specific")

# male-specific
thca_male_tsgs <- read_tsv("./files/degs_TvsN_stomach_thyroid_tsgs.txt") %>%
  filter(tissue == "Thyroid" & state == "male_specific")


# stomach

# female-specific
stad_female_tsgs <- read_tsv("./files/degs_TvsN_stomach_thyroid_tsgs.txt") %>%
filter(tissue == "Stomach" & state == "female_specific")

# male-specific
stad_male_tsgs <- read_tsv("./files/degs_TvsN_stomach_thyroid_tsgs.txt") %>%
filter(tissue == "Stomach" & state == "male_specific")






# expression data

# TCGA gene annotation
tcga_annotation <- read_tsv("./data/annotation/geneAnnot.gencode.v22.txt") %>%
  mutate(geneID = str_replace(geneID, "\\.[0-9]+", "")) %>%
  dplyr::rename(ensemblGeneID = geneID)


# thyroid

# load fpkm and metadata
thyroid_tcga_fpkm <- read_tsv("./files/thyroid_tcga_gtex_fpkm.txt") %>%
  gather(key = "sample", value="fpkm", -gene) %>%
  filter(!str_detect(sample, "GTEX")) %>%
  mutate(sample = str_replace_all(sample, "-", "."))

thyroid_tcga_meta <- read_tsv("./files/tcga_thca_meta.txt") %>%
  dplyr::rename(sample = submitter_id) %>%
  mutate(sample = str_replace_all(sample, "-", "."))

thyroid_tcga_fpkm <- thyroid_tcga_fpkm %>%
  inner_join(thyroid_tcga_meta[, c("sample", "sample_type_id", "gender")], by = "sample") %>%
  inner_join(tcga_annotation, by = c("gene" = "ensemblGeneID")) %>%
  filter(sample_type_id != "06")


# stomach

# load fpkm and metadata
stomach_tcga_fpkm <- read_tsv("./files/stomach_tcga_gtex_fpkm.txt") %>%
  gather(key = "sample", value="fpkm", -gene) %>%
  filter(!str_detect(sample, "GTEX")) %>%
  mutate(sample = str_replace_all(sample, "-", "."))

stomach_tcga_meta <- read_tsv("./files/tcga_stad_meta.txt") %>%
  dplyr::rename(sample = submitter_id) %>%
  mutate(sample = str_replace_all(sample, "-", "."))

stomach_tcga_fpkm <- stomach_tcga_fpkm %>%
  inner_join(stomach_tcga_meta[, c("sample", "sample_type_id", "gender")], by = "sample") %>%
  inner_join(tcga_annotation, by = c("gene" = "ensemblGeneID"))



# mutation data

# thyroid
thca_mut <- fread("./data/tcga/cbioportal/thca_tcga_pan_can_atlas_2018/data_mutations_extended.txt") %>%
  as_tibble() %>%
  dplyr::select(sample=Tumor_Sample_Barcode, gene=Hugo_Symbol, Consequence, Variant_Classification, Variant_Type, ref_allele=Reference_Allele, tumor_allele1=Tumor_Seq_Allele1, tumor_allele2=Tumor_Seq_Allele2, Start_Position, End_Position) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  inner_join(thyroid_tcga_meta %>% filter(sample_type_id == "01") %>% dplyr::select(sample, gender), by = "sample") %>%
  dplyr::select(sample, gender, everything()) %>%
  filter(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")) %>%
  group_by(sample, gene) %>%
  filter(!duplicated(Start_Position)) %>%
  ungroup() %>%
  filter(gene != "TTN")

thca_mut_n <- thca_mut %>%
  group_by(sample, gender, gene) %>%
  summarise(times_mut = n()) %>%
  ungroup() %>%
  arrange(desc(times_mut))


# stomach
stad_mut <- fread("./data/tcga/cbioportal/stad_tcga_pan_can_atlas_2018/data_mutations_extended.txt") %>%
  as_tibble() %>%
  dplyr::select(sample=Tumor_Sample_Barcode, gene=Hugo_Symbol, Consequence, Variant_Classification, Variant_Type, ref_allele=Reference_Allele, tumor_allele1=Tumor_Seq_Allele1, tumor_allele2=Tumor_Seq_Allele2, Start_Position, End_Position) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  inner_join(stomach_tcga_meta %>% filter(sample_type_id == "01") %>% dplyr::select(sample, gender), by = "sample") %>%
  dplyr::select(sample, gender, everything()) %>%
  filter(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")) %>%
  group_by(sample, gene) %>%
  filter(!duplicated(Start_Position)) %>%
  ungroup() %>%
  filter(gene != "TTN")

stad_mut_n <- stad_mut %>%
  group_by(sample, gender, gene) %>%
  summarise(times_mut = n()) %>%
  ungroup() %>%
  arrange(desc(times_mut))


# mutation frequence

# thyroid

# total
thca_mut_freq <- thca_mut_n %>%
  group_by(gene) %>%
  count() %>%
  mutate(freq = n/length(unique(thca_mut_n$sample))) %>%
  arrange(desc(freq))

# males
thca_mut_freq1 <- thca_mut_n %>%
  filter(gender == "male") %>%
  group_by(gene) %>%
  count() %>%
  ungroup() %>%
  mutate(freq = n/length(unique(thca_mut_n[thca_mut_n$gender == "male", c("sample"), drop=T]))) %>%
  mutate(gender = "male") %>%
  arrange(freq) %>%
  mutate(gene = factor(gene, unique(gene))) %>%
  mutate(male_tsgs = as.character(if_else(gene %in% thca_male_tsgs$geneName, 1, 0))) %>%
  mutate(female_tsgs = as.character(if_else(gene %in% thca_female_tsgs$geneName, 1, 0)))


# females
thca_mut_freq2 <- thca_mut_n %>%
  filter(gender == "female") %>%
  group_by(gene) %>%
  count() %>%
  ungroup() %>%
  mutate(freq = n/length(unique(thca_mut_n[thca_mut_n$gender == "female", c("sample"), drop=T]))) %>%
  mutate(gender = "female") %>%
  arrange(freq) %>%
  mutate(gene = factor(gene, unique(gene))) %>%
  mutate(male_tsgs = as.character(if_else(gene %in% thca_male_tsgs$geneName, 1, 0))) %>%
  mutate(female_tsgs = as.character(if_else(gene %in% thca_female_tsgs$geneName, 1, 0)))

thca_mut_freq_gender <- bind_rows(thca_mut_freq1, thca_mut_freq2)



# female-specific DEGs TSGs along male and female mutation freq distribution
#thca_mut_n %>% filter(gene %in% thca_female_tsgs$geneName)

thca_mut_freq_males_plot <- thca_mut_freq1 %>%
  ggplot(mapping = aes(x = gene, y = freq, color = female_tsgs)) +
  geom_segment(aes(x=gene, xend=gene, y=0, yend=freq, color = female_tsgs), size = 0.1) +
  geom_point(size = 0.5) +
  geom_text_repel(data = thca_mut_freq1 %>% filter(female_tsgs == "1"), mapping=aes(x = gene, y = freq, label = gene), color = "deeppink", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +
  scale_color_manual(values = c("#74a9cf", "deeppink"), guide = F) +

  #geom_segment(data = thca_mut_freq1 %>% tail %>% slice(5:6), aes(x=gene, xend=gene, y=0, yend=freq), size = 0.5, color = "orange") +
  #geom_point(data = thca_mut_freq1 %>% tail %>% slice(5:6), aes(x = gene, y = freq), color = "orange", size = 1) +
  geom_text(data = thca_mut_freq1 %>% tail %>% slice(5:6), mapping=aes(x = gene, y = freq, label = gene), size = 2, color = "black", vjust = -0.5, hjust = 1) +

  annotate("text",
    x = 200,
    y = 0.5,
    label = paste(paste("genes: ", nrow(thca_mut_freq1)), paste("samples: ", length(unique(thca_mut_n[thca_mut_n$gender == "male", c("sample"), drop=T]))), sep="\n"),
    size = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
		plot.title = element_text(size = 7, color="black")) +
  labs(x = "", y = "Frequency", title = "Males mutational distribution\nFemale-specific TSGs differentially expressed T-N in pink")



thca_mut_freq_females_plot <- thca_mut_freq2 %>%
  ggplot(mapping = aes(x = gene, y = freq, color = female_tsgs)) +
  geom_segment(aes(x=gene, xend=gene, y=0, yend=freq, color = female_tsgs), size = 0.1) +
  geom_point(size = 0.5) +
	geom_text_repel(data = thca_mut_freq2 %>% filter(female_tsgs == "1"), mapping=aes(x = gene, y = freq, label = gene), color = "deeppink", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +
  scale_color_manual(values = c("#fbb4b9", "deeppink"), guide = F) +

  #geom_segment(data = thca_mut_freq2 %>% tail %>% slice(5:6), aes(x=gene, xend=gene, y=0, yend=freq), size = 0.5, color = "orange") +
  #geom_point(data = thca_mut_freq2 %>% tail %>% slice(5:6), aes(x = gene, y = freq), color = "orange", size = 1) +
  geom_text(data = thca_mut_freq2 %>% tail %>% slice(5:6), mapping=aes(x = gene, y = freq, label = gene), size = 2, color = "black", vjust = -0.5, hjust = 1) +

  annotate("text",
    x = 431,
    y = 0.5,
    label = paste(paste("genes: ", nrow(thca_mut_freq2)), paste("samples: ", length(unique(thca_mut_n[thca_mut_n$gender == "female", c("sample"), drop=T]))), sep="\n"),
    size = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
	plot.title = element_text(size = 7, color="black")) +
  labs(x = "Gene", y = "Frequency", title = "Females mutational distribution\nFemale-specific TSGs differentially expressed T-N in pink")

thca_mut_freq_genders_plot <- plot_grid(thca_mut_freq_males_plot, thca_mut_freq_females_plot, labels = c("", ""), nrow = 2, align = "v")
ggsave(filename="thca_mut_freq_genders_tsgs_females.png", plot=thca_mut_freq_genders_plot, height=5, width=5, path = "./plots/mutation_analysis")
unlink("thca_mut_freq_genders_tsgs_females.png")



# male-specific DEGs TSGs along male and female mutation freq distribution
#thca_mut_n %>% filter(gene %in% thca_male_tsgs$geneName)

thca_mut_freq_males_plot <- thca_mut_freq1 %>%
  ggplot(mapping = aes(x = gene, y = freq, color = male_tsgs)) +
  geom_segment(aes(x=gene, xend=gene, y=0, yend=freq, color = male_tsgs), size = 0.1) +
  geom_point(size = 0.5) +
	geom_text_repel(data = thca_mut_freq1 %>% filter(male_tsgs == "1"), mapping=aes(x = gene, y = freq, label = gene), color = "blue", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +
  scale_color_manual(values = c("#74a9cf", "blue"), guide = F) +

  #geom_segment(data = thca_mut_freq1 %>% tail %>% slice(5:6), aes(x=gene, xend=gene, y=0, yend=freq), size = 0.5, color = "orange") +
  #geom_point(data = thca_mut_freq1 %>% tail %>% slice(5:6), aes(x = gene, y = freq), color = "orange", size = 1) +
  geom_text(data = thca_mut_freq1 %>% tail %>% slice(5:6), mapping=aes(x = gene, y = freq, label = gene), size = 2, color = "black", vjust = -0.5, hjust = 1) +

  annotate("text",
    x = 200,
    y = 0.5,
    label = paste(paste("genes: ", nrow(thca_mut_freq1)), paste("samples: ", length(unique(thca_mut_n[thca_mut_n$gender == "male", c("sample"), drop=T]))), sep="\n"),
    size = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
		plot.title = element_text(size = 7, color="black")) +
  labs(x = "", y = "Frequency", title = "Males mutational distribution\nMale-specific TSGs differentially expressed T-N in blue")



thca_mut_freq_females_plot <- thca_mut_freq2 %>%
  ggplot(mapping = aes(x = gene, y = freq, color = male_tsgs)) +
  geom_segment(aes(x=gene, xend=gene, y=0, yend=freq, color = male_tsgs), size = 0.1) +
  geom_point(size = 0.5) +
	geom_text_repel(data = thca_mut_freq2 %>% filter(male_tsgs == "1"), mapping=aes(x = gene, y = freq, label = gene), color = "blue", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +
  scale_color_manual(values = c("#fbb4b9", "blue"), guide = F) +

  #geom_segment(data = thca_mut_freq2 %>% tail %>% slice(5:6), aes(x=gene, xend=gene, y=0, yend=freq), size = 0.5, color = "orange") +
  #geom_point(data = thca_mut_freq2 %>% tail %>% slice(5:6), aes(x = gene, y = freq), color = "orange", size = 1) +
  geom_text(data = thca_mut_freq2 %>% tail %>% slice(5:6), mapping=aes(x = gene, y = freq, label = gene), size = 2, color = "black", vjust = -0.5, hjust = 1) +

  annotate("text",
    x = 431,
    y = 0.5,
    label = paste(paste("genes: ", nrow(thca_mut_freq2)), paste("samples: ", length(unique(thca_mut_n[thca_mut_n$gender == "female", c("sample"), drop=T]))), sep="\n"),
    size = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
		plot.title = element_text(size = 7, color="black")) +
  labs(x = "Gene", y = "Frequency", title = "Females mutational distribution\nMale-specific TSGs differentially expressed T-N in blue")

thca_mut_freq_genders_plot <- plot_grid(thca_mut_freq_males_plot, thca_mut_freq_females_plot, labels = c("", ""), nrow = 2, align = "v")
ggsave(filename="thca_mut_freq_genders_tsgs_males.png", plot=thca_mut_freq_genders_plot, height=5, width=5, path = "./plots/mutation_analysis")
unlink("thca_mut_freq_genders_tsgs_males.png")







# stomach

# total
stad_mut_freq <- stad_mut_n %>%
  group_by(gene) %>%
  count() %>%
  mutate(freq = n/length(unique(stad_mut_n$sample))) %>%
  arrange(desc(freq))

# males
stad_mut_freq1 <- stad_mut_n %>%
  filter(gender == "male") %>%
  group_by(gene) %>%
  count() %>%
  ungroup() %>%
  mutate(freq = n/length(unique(stad_mut_n[stad_mut_n$gender == "male", c("sample"), drop=T]))) %>%
  mutate(gender = "male") %>%
  arrange(freq) %>%
  mutate(gene = factor(gene, unique(gene))) %>%
  mutate(male_tsgs = as.character(if_else(gene %in% stad_male_tsgs$geneName, 1, 0))) %>%
  mutate(female_tsgs = as.character(if_else(gene %in% stad_female_tsgs$geneName, 1, 0)))

# females
stad_mut_freq2 <- stad_mut_n %>%
  filter(gender == "female") %>%
  group_by(gene) %>%
  count() %>%
  ungroup() %>%
  mutate(freq = n/length(unique(stad_mut_n[stad_mut_n$gender == "female", c("sample"), drop=T]))) %>%
  mutate(gender = "female") %>%
  arrange(freq) %>%
  mutate(gene = factor(gene, unique(gene))) %>%
  mutate(male_tsgs = as.character(if_else(gene %in% stad_male_tsgs$geneName, 1, 0))) %>%
  mutate(female_tsgs = as.character(if_else(gene %in% stad_female_tsgs$geneName, 1, 0)))


stad_mut_freq_gender <- bind_rows(stad_mut_freq1, stad_mut_freq2)


# female-specific DEGs TSGs along male and female mutation freq distribution
#stad_mut_n %>% filter(gene %in% stad_female_tsgs$geneName)

stad_mut_freq_males_plot <- stad_mut_freq1 %>%
  ggplot(mapping = aes(x = gene, y = freq, color = female_tsgs)) +
  geom_segment(aes(x=gene, xend=gene, y=0, yend=freq, color = female_tsgs), size = 0.1) +
  geom_point(size = 0.5) +
	geom_text_repel(data = stad_mut_freq1 %>% filter(female_tsgs == "1", gene %in% (stad_mut_freq2 %>% filter(female_tsgs == "1") %>% arrange(desc(n)) %>% head(10) %>% pull(gene))),
		mapping=aes(x = gene, y = freq, label = gene), color = "deeppink", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +
  scale_color_manual(values = c("#74a9cf", "deeppink"), guide = F) +

  #geom_segment(data = stad_mut_freq1 %>% tail, aes(x=gene, xend=gene, y=0, yend=freq), size = 0.5, color = "orange") +
  #geom_point(data = stad_mut_freq1 %>% tail, aes(x = gene, y = freq), color = "orange", size = 1) +
  geom_text_repel(data = stad_mut_freq1 %>% tail, mapping=aes(x = gene, y = freq, label = gene), color = "black", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +

  annotate("text",
    x = 2000,
    y = 0.45,
    label = paste(paste("genes: ", nrow(stad_mut_freq1)), paste("samples: ", length(unique(stad_mut_n[stad_mut_n$gender == "male", c("sample"), drop=T]))), sep="\n"),
    size = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
		plot.title = element_text(size = 7, color="black")) +
  labs(x = "", y = "Frequency", title = "Males mutational distribution\nFemale-specific TSGs differentially expressed T-N in pink")



stad_mut_freq_females_plot <- stad_mut_freq2 %>%
  ggplot(mapping = aes(x = gene, y = freq, color = female_tsgs)) +
  geom_segment(aes(x=gene, xend=gene, y=0, yend=freq, color = female_tsgs), size = 0.1) +
  geom_point(size = 0.5) +
	geom_text_repel(data = stad_mut_freq2 %>% filter(female_tsgs == "1") %>% arrange(desc(n)) %>% head(10),
		mapping=aes(x = gene, y = freq, label = gene), color = "deeppink", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +
  scale_color_manual(values = c("#fbb4b9", "deeppink"), guide = F) +

  #geom_segment(data = stad_mut_freq2 %>% tail, aes(x=gene, xend=gene, y=0, yend=freq), size = 0.5, color = "orange") +
  #geom_point(data = stad_mut_freq2 %>% tail, aes(x = gene, y = freq), color = "orange", size = 1) +
	geom_text_repel(data = stad_mut_freq2 %>% tail, mapping=aes(x = gene, y = freq, label = gene), color = "black", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +

  annotate("text",
    x = 2000,
    y = 0.45,
    label = paste(paste("genes: ", nrow(stad_mut_freq2)), paste("samples: ", length(unique(stad_mut_n[stad_mut_n$gender == "female", c("sample"), drop=T]))), sep="\n"),
    size = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
		plot.title = element_text(size = 7, color="black")) +
  labs(x = "Gene", y = "Frequency", title = "Females mutational distribution\nFemale-specific TSGs differentially expressed T-N in pink")

stad_mut_freq_genders_plot <- plot_grid(stad_mut_freq_males_plot, stad_mut_freq_females_plot, labels = c("", ""), nrow = 2, align = "v")
ggsave(filename="stad_mut_freq_genders_tsgs_females.png", plot=stad_mut_freq_genders_plot, height=5, width=5, path = "./plots/mutation_analysis")
unlink("stad_mut_freq_genders_tsgs_females.png")


# male-specific DEGs TSGs along male and female mutation freq distribution
#thca_mut_n %>% filter(gene %in% thca_male_tsgs$geneName)

stad_mut_freq_males_plot <- stad_mut_freq1 %>%
  ggplot() +
  geom_segment(mapping = aes(x=gene, xend=gene, y=0, yend=freq), size = 0.1, color = "#74a9cf") +
  geom_point(mapping = aes(x = gene, y = freq), size = 0.5, color = "#74a9cf") +

	geom_segment(data = stad_mut_freq1 %>% filter(male_tsgs == "1") %>% arrange(desc(freq)) %>% head(10), aes(x=gene, xend=gene, y=0, yend=freq), size = 0.5, color = "black") +
	geom_point(data = stad_mut_freq1 %>% filter(male_tsgs == "1") %>% arrange(desc(freq)) %>% head(10), aes(x = gene, y = freq), color = "black", size = 1) +
	geom_text_repel(data = stad_mut_freq1 %>% filter(male_tsgs == "1") %>% arrange(desc(freq)) %>% head(10),
		mapping=aes(x = gene, y = freq, label = gene), color = "black", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +

  annotate("text",
    x = 2000,
    y = 0.45,
    label = paste(paste("genes: ", nrow(stad_mut_freq1)), paste("samples: ", length(unique(stad_mut_n[stad_mut_n$gender == "male", c("sample"), drop=T]))), sep="\n"),
    size = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_text(color = "black", size = 15),
		plot.title = element_text(size = 7, color="black")) +
	scale_y_continuous(labels = scales::percent, limits = c(NA, 0.6)) +
  labs(x = "", y = "Frequency", title = "Males: mutational distribution\nMale-specific TSGs differentially expressed T-N highlighted (top 10 mutated)")



stad_mut_freq_females_plot <- stad_mut_freq2 %>%
  ggplot() +
  geom_segment(mapping = aes(x=gene, xend=gene, y=0, yend=freq), size = 0.1, color = "#fbb4b9") +
  geom_point(mapping = aes(x = gene, y = freq), size = 0.5, color = "#fbb4b9") +

	geom_segment(data = stad_mut_freq2 %>% filter(gene %in% (stad_mut_freq1 %>% filter(male_tsgs == "1") %>% arrange(desc(freq)) %>% head(10) %>% pull(gene) %>% as.character())),
		aes(x=gene, xend=gene, y=0, yend=freq), size = 0.5, color = "black") +
	geom_point(data = stad_mut_freq2 %>% filter(gene %in% (stad_mut_freq1 %>% filter(male_tsgs == "1") %>% arrange(desc(freq)) %>% head(10) %>% pull(gene) %>% as.character())),
		aes(x = gene, y = freq), color = "black", size = 1) +
	geom_text_repel(data = stad_mut_freq2 %>% filter(gene %in% (stad_mut_freq1 %>% filter(male_tsgs == "1") %>% arrange(desc(freq)) %>% head(10) %>% pull(gene) %>% as.character())),
		mapping=aes(x = gene, y = freq, label = gene), color = "black", size = 2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines")) +

  annotate("text",
    x = 2000,
    y = 0.45,
    label = paste(paste("genes: ", nrow(stad_mut_freq2)), paste("samples: ", length(unique(stad_mut_n[stad_mut_n$gender == "female", c("sample"), drop=T]))), sep="\n"),
    size = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
		axis.title = element_text(color = "black", size = 15),
		plot.title = element_text(size = 7, color="black")) +
	scale_y_continuous(labels = scales::percent, limits = c(NA, 0.6)) +
  labs(x = "Gene", y = "Frequency", title = "Females: mutational distribution\nMale-specific TSGs differentially expressed T-N highlighted (top 10 mutated)")

stad_mut_freq_genders_plot <- plot_grid(stad_mut_freq_males_plot, stad_mut_freq_females_plot, labels = c("", ""), nrow = 2, align = "v")
ggsave(filename="stad_mut_freq_genders_tsgs_males.png", plot=stad_mut_freq_genders_plot, height=5, width=5, path = "./plots/mutation_analysis")
unlink("stad_mut_freq_genders_tsgs_males.png")


# compare male and female mutational frequency
gender_frequency <- stad_mut_freq1 %>%
	filter(male_tsgs == "1") %>% arrange(desc(freq)) %>% head(10) %>% dplyr::select(gene, freq, n, gender) %>%
	bind_rows(stad_mut_freq2 %>% filter(gene %in% (stad_mut_freq1 %>% filter(male_tsgs == "1") %>% arrange(desc(freq)) %>% head(10) %>% pull(gene) %>% as.character())) %>% dplyr::select(gene, freq, n, gender)) %>%
	mutate(gene = fct_inorder(gene)) %>%
	ggplot() +
	geom_bar(mapping = aes(x = gene, y = freq, fill = gender), stat = "identity", position = "dodge") +
	theme(
		axis.text.x = element_text(color = "black", size = 15, angle = 90, hjust = 1),
		axis.title.x = element_blank(),
		axis.title.y = element_text(color = "black", size = 15),
		axis.text.y = element_text(color = "black", size = 13),
		plot.title = element_text(size = 7, color="black")) +
	labs(x = "", y = "Frequency", title = "") +
	scale_y_continuous(labels = scales::percent, limits = c(NA, 0.08)) +
	scale_fill_manual(name = "Gender", labels = c("Female", "Male"), values = c("#fbb4b9", "#74a9cf"))
ggsave(filename="stad_mut_freq_genders_tsgs_males_freq.png", plot=gender_frequency, height=2, width=5, path = "./plots/mutation_analysis")
unlink("stad_mut_freq_genders_tsgs_males_freq.png")










#stad_mut_freq_males_plot$data %>% filter(female_tsgs == "1") %>% as.data.frame
#stad_mut_freq_females_plot$data %>% filter(female_tsgs == "1") %>% as.data.frame
#intersect(stad_mut_freq_males_plot$data %>% filter(female_tsgs == "1") %>% pull(gene), stad_mut_freq_females_plot$data %>% filter(female_tsgs == "1") %>% pull(gene))

#stad_mut_freq_males_plot$data %>% filter(male_tsgs == "1") %>% as.data.frame
#stad_mut_freq_females_plot$data %>% filter(male_tsgs == "1") %>% as.data.frame
#intersect(stad_mut_freq_males_plot$data %>% filter(male_tsgs == "1") %>% pull(gene), stad_mut_freq_females_plot$data %>% filter(male_tsgs == "1") %>% pull(gene))



# male vs female DEGs in tumours


# thyroid

thca_MvsF_mut_freq_gender <- thca_MvsF %>%
  filter(state == "tumour_specific") %>%
  dplyr::select(genes, geneName, state, logFC = log2FC_tumour) %>%
  inner_join(thca_mut_freq_gender, by = c("geneName" = "gene")) %>%
  arrange(geneName) %>%
  as.data.frame()

#thca_MvsF_mut_freq_gender %>% group_by(gender) %>% summarise(freq = median(freq))
# higher freq in men: probably because they are less

#thca_MvsF_mut_freq_gender %>% group_by(gender) %>% summarise(count = sum(n))
# female: 16 mutations
# male: 8 mutations
# females have more mutations

#thca_MvsF_mut_freq_gender %>% group_by(geneName) %>% count()
#mutations are dispersed





# stomach

stad_MvsF_mut_freq_gender <- stad_MvsF %>%
  filter(state == "tumour_specific") %>%
  dplyr::select(genes, geneName, state, logFC = log2FC_tumour) %>%
  inner_join(stad_mut_freq_gender, by = c("geneName" = "gene")) %>%
  arrange(geneName) %>%
  as.data.frame()
#stad_MvsF_mut_freq_gender %>% group_by(gender) %>% summarise(freq = median(freq))
# higher freq in women: probably because they are less

#stad_MvsF_mut_freq_gender %>% group_by(gender) %>% summarise(count = sum(n))
# female: 73 mutations
# male: 74 mutations

# similar number of mutations: much higher than in thyroid though



# pesquisar por genetic signatures (BGI paper)
