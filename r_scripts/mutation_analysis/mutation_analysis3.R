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
  filter(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation")) %>%
	#filter(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")) %>%
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
  filter(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation")) %>%
	#filter(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")) %>%
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
	ungroup() %>%
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




# stomach

# total
stad_mut_freq <- stad_mut_n %>%
  group_by(gene) %>%
  count() %>%
	ungroup() %>%
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




# stomach

stad_mut_freq_male_tsgs_down_tumor <- stad_mut_freq1 %>% filter(gene %in% (stad_male_tsgs %>% filter(logFC_males < 0) %>% pull(geneName)))
stad_mut_freq_male_tsgs_up_tumor <- stad_mut_freq1 %>% filter(gene %in% (stad_male_tsgs %>% filter(logFC_males > 0) %>% pull(geneName)))

stad_mut_freq_female_tsgs_down_tumor <- stad_mut_freq2 %>% filter(gene %in% (stad_female_tsgs %>% filter(logFC_females < 0) %>% pull(geneName)))
stad_mut_freq_female_tsgs_up_tumor <- stad_mut_freq2 %>% filter(gene %in% (stad_female_tsgs %>% filter(logFC_females > 0) %>% pull(geneName)))




# thyroid

thca_mut_freq_male_tsgs_down_tumor <- thca_mut_freq1 %>% filter(gene %in% (thca_male_tsgs %>% filter(logFC_males < 0) %>% pull(geneName)))
thca_mut_freq_male_tsgs_up_tumor <- thca_mut_freq1 %>% filter(gene %in% (thca_male_tsgs %>% filter(logFC_males > 0) %>% pull(geneName)))

thca_mut_freq_female_tsgs_down_tumor <- thca_mut_freq2 %>% filter(gene %in% (thca_female_tsgs %>% filter(logFC_females < 0) %>% pull(geneName)))
thca_mut_freq_female_tsgs_up_tumor <- thca_mut_freq2 %>% filter(gene %in% (thca_female_tsgs %>% filter(logFC_females > 0) %>% pull(geneName)))
