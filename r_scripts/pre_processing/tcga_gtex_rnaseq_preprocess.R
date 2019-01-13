#Understanding Gender Differential Susceptibility in Cancer

#What are the transcriptomic dissimilarities between genders that may contribute or underlie a differential cancer incidence?


#Thyroid and Stomach cancers


#https://cancergenome.nih.gov/cancersselected/thyroid
#https://cancergenome.nih.gov/cancersselected/stomachcancer
#https://tcga-data.nci.nih.gov/docs/publications/thca_2014/
#https://tcga-data.nci.nih.gov/docs/publications/stad_2014/


#==============================================================
#
#	RNA-seq AND METADATA PRE-PROCESSING
#
#==============================================================


#RNA-seq data from TCGA GDC Data Portal:
#https://portal.gdc.cancer.gov/

#HARMONIZED data
#https://gdc.cancer.gov/about-data/data-harmonization-and-generation/genomic-data-harmonization-0#Overview

#RNA-seq and normalization pipeline
#https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/


#important info about the samples IDs:
#https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
#https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables
#https://gdc.cancer.gov/about-data/data-harmonization-and-generation/clinical-data-harmonization


#removing expression files from original directories
#for dir in `ls`; do cd $dir; cp `ls | grep -v "logs"` ../../files; cd ../; done
#for dir in `ls -d */ | egrep "^.{37}"`; do cd $dir; cp `ls | grep -v "logs"` ../files; cd ../; done
#`ls | grep -v "logs" | grep -v ".tar.gz" | grep -v "annotations.txt"`



library(tidyverse)
library(edgeR)
library(rjson)



get.metaData <- function(metadata.json, flag){

	#metadata processing
	#returns list containing two data.frames:
	#1 information by expression file
	#2 information by sample (person)

	metadata <- fromJSON(file = metadata.json)
	exprFiles.info <- c()
	samplesClinical.info <- c()
	for(i in 1:length(metadata)){
		file_name <- metadata[[i]]$file_name
		file_type <- strsplit(metadata[[i]]$submitter_id, "_", fixed=T)[[1]][2]
		submitter_id <- substring(metadata[[i]]$cases[[1]]$samples[[1]]$submitter_id, 1, 15)
		sample_type_id <- metadata[[i]]$cases[[1]]$samples[[1]]$sample_type_id
		sample_type <- metadata[[i]]$cases[[1]]$samples[[1]]$sample_type
		
		exprFiles.info <- rbind(exprFiles.info, c(file_name, file_type, submitter_id, sample_type_id, sample_type))
		
		if(metadata[[i]]$cases[[1]]$demographic$gender == "not reported"){
			gender <- NA
		}
		else{
			gender <- metadata[[i]]$cases[[1]]$demographic$gender
		}
		if(metadata[[i]]$cases[[1]]$demographic$race == "not reported"){
			race <- NA
		}
		else{
			race <- metadata[[i]]$cases[[1]]$demographic$race
		}
		if(metadata[[i]]$cases[[1]]$demographic$ethnicity == "not reported"){
			ethnicity <- NA
		}
		else{
			ethnicity <- metadata[[i]]$cases[[1]]$demographic$ethnicity
		}
		if(class(metadata[[i]]$cases[[1]]$demographic$year_of_birth) == "NULL"){
			age_actual <- NA
		}
		else{
			age_actual <- as.numeric(substring(Sys.Date(),1,4))-metadata[[i]]$cases[[1]]$demographic$year_of_birth
		}
		if(class(metadata[[i]]$cases[[1]]$diagnoses[[1]]$age_at_diagnosis) == "NULL"){
			age_at_diagnosis <- NA
		}
		else{
			age_at_diagnosis <- round(metadata[[i]]$cases[[1]]$diagnoses[[1]]$age_at_diagnosis/365)
		}
		if(metadata[[i]]$cases[[1]]$diagnoses[[1]]$tumor_stage == "not reported"){
			tumor_stage <- NA
		}
		else{
			tumor_stage <- metadata[[i]]$cases[[1]]$diagnoses[[1]]$tumor_stage
		}
		if(metadata[[i]]$cases[[1]]$diagnoses[[1]]$vital_status == "not reported"){
			vital_status <- NA
		}
		else{
			vital_status <- metadata[[i]]$cases[[1]]$diagnoses[[1]]$vital_status
		}

		tss <- substring(metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 6, 7)
		vial <- substring(metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 16, 16)
		portion <- substring(metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 18, 19)
		analyte <- substring(metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 20, 20)
		plate <- substring(metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 22, 25)
		center <- substring(metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 27, 28)
		
		samplesClinical.info <- rbind(samplesClinical.info, c(submitter_id, sample_type_id, gender, race, ethnicity, age_actual, age_at_diagnosis, tumor_stage, vital_status, tss, vial, portion, analyte, plate, center))
	}
	
	exprFiles.info <- as.data.frame(exprFiles.info, stringsAsFactors=F)
	colnames(exprFiles.info) <- c("file_name", "file_type", "submitter_id", "sample_type_id", "sample_type")
	
	samplesClinical.info <- as.data.frame(samplesClinical.info, stringsAsFactors=F)
	colnames(samplesClinical.info) <- c("submitter_id", "sample_type_id", "gender", "race", "ethnicity", "age_actual", "age_at_diagnosis", "tumor_stage", "vital_status", "tss", "vial", "portion", "analyte", "plate", "center")
	samplesClinical.info <- unique(samplesClinical.info)
	rownames(samplesClinical.info) <- samplesClinical.info$submitter_id
	
	write.table(exprFiles.info, paste(flag, "exprFiles_info.txt", sep = "."), quote=F, sep="\t", row.names=F)
	write.table(samplesClinical.info, paste(flag, "samplesClinical_info.txt", sep = "."), quote=F, sep="\t", row.names=F)
	
	meta <- list(exprFiles.info = exprFiles.info, samplesClinical.info = samplesClinical.info)
	return(meta)
}


get.expressionData <- function(exprFiles.info, files.dir, geneIDs){

	#counts and FPKM processing
	#returns list containing two data.frames:
	#1 FPKM for all genes (rows) and samples (columns)
	#2 counts for all genes (rows) and samples (columns)

	fpkm <- data.frame(row.names = geneIDs)
	counts <- data.frame(row.names = geneIDs)
	for(i in 1:nrow(exprFiles.info)){
		if(exprFiles.info[i,"file_type"] == "fpkm"){
			expr.file <- read.table(gzfile(paste(files.dir, as.character(exprFiles.info[i,"file_name"]), sep="")), header=F)
			rownames(expr.file) <- expr.file$V1
			expr.file <- expr.file[geneIDs, , drop=F]
			fpkm <- cbind(fpkm, expr.file$V2)
			colnames(fpkm)[ncol(fpkm)] <- as.character(exprFiles.info[i,"submitter_id"])
		}
		else if(exprFiles.info[i,"file_type"] == "count"){
			expr.file <- read.table(gzfile(paste(files.dir, as.character(exprFiles.info[i,"file_name"]), sep="")), header=F)
			expr.file <- expr.file[-c(60484,60485,60486,60487,60488), ,drop=F]
			counts <- cbind(counts, expr.file$V2)
			colnames(counts)[ncol(counts)] <- as.character(exprFiles.info[i,"submitter_id"])
		}
	}
	expression.data <- list(fpkm = fpkm, counts = counts)
	return(expression.data)
}


filter.genes.tcga <- function(expression.counts, geneIDs, cpm.cutO, p.cutO){

	#filter out the low-expressed genes
	#returns a vector containing the ensembl gene IDs of the genes that passed the filter


	tumour.data <- expression.counts[,substring(colnames(expression.counts), 14, 15) %in% c("01")]
	normal.data <- expression.counts[,substring(colnames(expression.counts), 14, 15) %in% c("11")]


	#select genes that have at least # reads per million mapped reads (CPM) in at least #% of tumour OR normal tissue
	keep <- ( rowSums( cpm(tumour.data) > cpm.cutO ) >= ncol(tumour.data)*p.cutO ) | ( rowSums( cpm(normal.data) > cpm.cutO ) >= ncol(normal.data)*p.cutO )
	keep.ensids <- geneIDs[keep]


	return(keep.ensids)
}

filter.genes.gtex <- function(expression.counts, geneIDs, cpm.cutO, p.cutO){

	#filter out the low-expressed genes
	#returns a vector containing the ensembl gene IDs of the genes that passed the filter


	#select genes that have at least # reads per million mapped reads (CPM) in at least #% of samples
	keep <- rowSums( cpm(expression.counts) > cpm.cutO ) >= ncol(expression.counts)*p.cutO
	keep.ensids <- geneIDs[keep]


	return(keep.ensids)
}




#========================
#	read annotation files
#========================


tcga.geneIDs <- read.table("./data/annotation/geneIDs_from_counts.txt", stringsAsFactors=F)$V1
tcga.geneIDs.annot <- read.table("./data/annotation/geneAnnot.gencode.v22.txt", sep="\t", h=T, stringsAsFactors=F)

tcga.geneIDs.annot.2 <- tcga.geneIDs.annot
tcga.geneIDs.annot.2$geneID <- sapply(strsplit(tcga.geneIDs.annot.2$geneID, split=".", fixed=T), function(x)(x[1]))



#=====================================
#	THCA - Papillary Thyroid Carcinoma
#=====================================


meta.file.thca <- "/Volumes/toshiba/Projects/Gender_differences/Gender_cancer_raw_data/TCGA.GDC/harmonized_data/rnaseq/THCA/metadata.cart.2016-10-14T14-04-10.146974.json"
files.dir.thca <- "/Volumes/toshiba/Projects/Gender_differences/Gender_cancer_raw_data/TCGA.GDC/harmonized_data/rnaseq/THCA/files/"


#read metadata
metadata.thca <- get.metaData(meta.file.thca, "./files/THCA")

thca_paper.info <- read.csv("./data/from_paper/thca_additional_info.csv", sep = ";", h=T, stringsAsFactors=F, na.strings=c("","NA"))
thca_paper.info <- thca_paper.info[, -c(11)]
thca_paper.info$sample <- substring(thca_paper.info$sample, 1, 15)
metadata.thca$samplesClinical.info <- merge(metadata.thca$samplesClinical.info, thca_paper.info, by.x = "submitter_id", by.y = "sample")
rownames(metadata.thca$samplesClinical.info) <- metadata.thca$samplesClinical.info$submitter_id


#remove NAs from clinical info
metadata.thca$samplesClinical.info[is.na(metadata.thca$samplesClinical.info$age_actual), c("age_actual")] <- as.character(round(mean(as.numeric(metadata.thca$samplesClinical.info[!is.na(metadata.thca$samplesClinical.info$age_actual), c("age_actual")]))))
metadata.thca$samplesClinical.info[is.na(metadata.thca$samplesClinical.info$age_at_diagnosis), c("age_at_diagnosis")] <- as.character(round(mean(as.numeric(metadata.thca$samplesClinical.info[!is.na(metadata.thca$samplesClinical.info$age_at_diagnosis), c("age_at_diagnosis")]))))
metadata.thca$samplesClinical.info[is.na(metadata.thca$samplesClinical.info)] <- "unknown"
metadata.thca$samplesClinical.info$histological_type[metadata.thca$samplesClinical.info$histological_type == "Other"] <- "unknown"

#read expression data (counts and FPKM)
expressiondata.thca <- get.expressionData(metadata.thca$exprFiles.info, files.dir.thca, tcga.geneIDs)




#================================
#	STAD - Stomach Adenocarcinoma
#================================


meta.file.stad <- "/Volumes/toshiba/Projects/Gender_differences/Gender_cancer_raw_data/TCGA.GDC/harmonized_data/rnaseq/STAD/metadata.cart.2016-10-17T11-05-37.830972.json"
files.dir.stad <- "/Volumes/toshiba/Projects/Gender_differences/Gender_cancer_raw_data/TCGA.GDC/harmonized_data/rnaseq/STAD/files/"


#read metadata
metadata.stad <- get.metaData(meta.file.stad, "./files/STAD")

stad_paper.info <- read.csv("./data/from_paper/stad_additional_info.csv", sep = ";", h=T, stringsAsFactors=F)
metadata.stad$samplesClinical.info <- merge(metadata.stad$samplesClinical.info, stad_paper.info, by.x = "submitter_id", by.y = "TCGA_barcode")
rownames(metadata.stad$samplesClinical.info) <- metadata.stad$samplesClinical.info$submitter_id


#remove NAs from clinical info
metadata.stad$samplesClinical.info[is.na(metadata.stad$samplesClinical.info$age_actual), c("age_actual")] <- as.character(round(mean(as.numeric(metadata.stad$samplesClinical.info[!is.na(metadata.stad$samplesClinical.info$age_actual), c("age_actual")]))))
metadata.stad$samplesClinical.info[is.na(metadata.stad$samplesClinical.info$age_at_diagnosis), c("age_at_diagnosis")] <- as.character(round(mean(as.numeric(metadata.stad$samplesClinical.info[!is.na(metadata.stad$samplesClinical.info$age_at_diagnosis), c("age_at_diagnosis")]))))
metadata.stad$samplesClinical.info[is.na(metadata.stad$samplesClinical.info)] <- "unknown"


#read expression data (counts and FPKM)
expressiondata.stad <- get.expressionData(metadata.stad$exprFiles.info, files.dir.stad, tcga.geneIDs)



#===============
#	filter genes
#===============


#select protein-coding and lincRNA genes
#https://www.gencodegenes.org/gencode_biotypes.html
#http://www.ensembl.org/info/genome/genebuild/ncrna.html
accepted.biotypes <- c("protein_coding", "lincRNA")


accepted.genes <- sort(tcga.geneIDs.annot[ ( tcga.geneIDs.annot$geneType %in% accepted.biotypes ), c("geneID") ])
#27470


#select genes that have at least 5 reads per million mapped reads (CPM) in at least 20% of tumour OR normal samples
keep.ensids.thca <- filter.genes.tcga(expressiondata.thca$counts[accepted.genes,], accepted.genes, 5, 0.2)
#12960
keep.ensids.stad <- filter.genes.tcga(expressiondata.stad$counts[accepted.genes,], accepted.genes, 5, 0.2)
#13674



#==========
#split data
#==========


#thca - thyroid
tcga_thca_counts <- expressiondata.thca$counts[keep.ensids.thca, ]
rownames(tcga_thca_counts) <- sapply(strsplit(rownames(tcga_thca_counts), split=".", fixed=T), function(x)(x[1]))
tcga_thca_fpkm <- expressiondata.thca$fpkm[keep.ensids.thca, ]
rownames(tcga_thca_fpkm) <- sapply(strsplit(rownames(tcga_thca_fpkm), split=".", fixed=T), function(x)(x[1]))
tcga_thca_meta <- metadata.thca$samplesClinical.info

tcga_thca_counts <- tcga_thca_counts %>%
	mutate(gene = rownames(tcga_thca_counts)) %>%
	dplyr::select(gene, everything())

tcga_thca_fpkm <- tcga_thca_fpkm %>%
	mutate(gene = rownames(tcga_thca_fpkm)) %>%
	dplyr::select(gene, everything())

write.table(tcga_thca_counts, "./files/tcga_thca_counts.txt", quote=F, sep="\t", row.names=F)
write.table(tcga_thca_fpkm, "./files/tcga_thca_fpkm.txt", quote=F, sep="\t", row.names=F)
write.table(tcga_thca_meta, "./files/tcga_thca_meta.txt", quote=F, sep="\t", row.names=F)



#stad - stomach
tcga_stad_counts <- expressiondata.stad$counts[keep.ensids.stad, ]
rownames(tcga_stad_counts) <- sapply(strsplit(rownames(tcga_stad_counts), split=".", fixed=T), function(x)(x[1]))
tcga_stad_fpkm <- expressiondata.stad$fpkm[keep.ensids.stad, ]
rownames(tcga_stad_fpkm) <- sapply(strsplit(rownames(tcga_stad_fpkm), split=".", fixed=T), function(x)(x[1]))
tcga_stad_meta <- metadata.stad$samplesClinical.info

tcga_stad_counts <- tcga_stad_counts %>%
	mutate(gene = rownames(tcga_stad_counts)) %>%
	dplyr::select(gene, everything())

tcga_stad_fpkm <- tcga_stad_fpkm %>%
	mutate(gene = rownames(tcga_stad_fpkm)) %>%
	dplyr::select(gene, everything())

write.table(tcga_stad_counts, "./files/tcga_stad_counts.txt", quote=F, sep="\t", row.names=F)
write.table(tcga_stad_fpkm, "./files/tcga_stad_fpkm.txt", quote=F, sep="\t", row.names=F)
write.table(tcga_stad_meta, "./files/tcga_stad_meta.txt", quote=F, sep="\t", row.names=F)








#=========
#GTEx data
#=========

#read GTEx data for thyroid and stomach tissues


get.GTEx.data <- function(readcnts.file, rpkm.file, metadata.file){


	readcnts <- read.table(readcnts.file, sep = "\t", h = T, stringsAsFactors = F)
	rownames(readcnts) <- readcnts$Name
	#print(all.equal(geneIDs, readcnts$Name))
	readcnts <- readcnts[,-c(1,2)]
	
	rpkm <- read.table(rpkm.file, sep = "\t", h = T, stringsAsFactors = F)
	rownames(rpkm) <- rpkm$Name
	#print(all.equal(geneIDs, rpkm$Name))
	rpkm <- rpkm[,-c(1,2)]

	metadata <- read.table(metadata.file, sep = "\t", h = T, stringsAsFactors = F)

	thyroid.counts <- readcnts[ , colnames(readcnts) %in% rownames(metadata[(metadata$SMTSD == "Thyroid"), ]) ]
	thyroid.rpkm <- rpkm[ , colnames(rpkm) %in% rownames(metadata[(metadata$SMTSD == "Thyroid"), ]) ]
	thyroid.meta <- metadata[(metadata$SMTSD == "Thyroid"), ]

	stomach.counts <- readcnts[ , colnames(readcnts) %in% rownames(metadata[(metadata$SMTSD == "Stomach"), ]) ]
	stomach.rpkm <- rpkm[ , colnames(rpkm) %in% rownames(metadata[(metadata$SMTSD == "Stomach"), ]) ]
	stomach.meta <- metadata[(metadata$SMTSD == "Stomach"), ]

	data <- list(thyroid = NULL, stomach = NULL)
	data$thyroid <- list(counts = thyroid.counts, rpkm = thyroid.rpkm, metadata = thyroid.meta)
	data$stomach <- list(counts = stomach.counts, rpkm = stomach.rpkm, metadata = stomach.meta)

	return(data)
}




#cat readcnt.thyroid_stomach.tab | awk '{print $2"\t"$3}' | sed "1d" > geneIDs.txt
gtex.geneIDs <- read.table("./data/gtex.thyroid.stomach/geneIDs.txt", h=F, sep="\t", stringsAsFactors=F)
gtex.geneIDs.annot <- read.table("./data/gtex.thyroid.stomach/geneAnnot.gencode.v19.txt", h=T, sep="\t", stringsAsFactors=F)
gtex.geneIDs.annot.2 <- gtex.geneIDs.annot
gtex.geneIDs.annot.2$geneID <- sapply( strsplit(gtex.geneIDs.annot.2$geneID, split = ".", fixed = T), function(x)(x[1]) )

gtex.readcnts.file <- "./data/gtex.thyroid.stomach/readcnt.thyroid_stomach.tab"
gtex.rpkm.file <- "./data/gtex.thyroid.stomach/rpkm.thyroid_stomach.tab"
gtex.metadata.file <- "./data/gtex.thyroid.stomach/metadata.tab"

gtex.data <- get.GTEx.data(gtex.readcnts.file, gtex.rpkm.file, gtex.metadata.file)


#===============
#	filter genes
#===============


#select protein-coding and lincRNA genes
#https://www.gencodegenes.org/gencode_biotypes.html
#http://www.ensembl.org/info/genome/genebuild/ncrna.html
accepted.biotypes <- c("protein_coding", "lincRNA")


accepted.genes.gtex <- sort(gtex.geneIDs.annot[ ( gtex.geneIDs.annot$geneType %in% accepted.biotypes ), c("geneID") ])
#27459


#select genes that have at least 5 reads per million mapped reads (CPM) in at least 20% of samples
keep.ensids.thyroid <- filter.genes.gtex(gtex.data$thyroid$counts[intersect(accepted.genes.gtex, rownames(gtex.data$thyroid$counts)), ], intersect(accepted.genes.gtex, rownames(gtex.data$thyroid$counts)), 5, 0.2)
#12501

keep.ensids.stomach <- filter.genes.gtex(gtex.data$stomach$counts[intersect(accepted.genes.gtex, rownames(gtex.data$stomach$counts)), ], intersect(accepted.genes.gtex, rownames(gtex.data$stomach$counts)), 5, 0.2)
#12371



#==========
#split data
#==========


#thyroid
gtex_thyroid_counts <- gtex.data$thyroid$counts[keep.ensids.thyroid, ]
rownames(gtex_thyroid_counts) <- sapply(strsplit(rownames(gtex_thyroid_counts), split=".", fixed=T), function(x)(x[1]))
gtex_thyroid_fpkm <- gtex.data$thyroid$rpkm[keep.ensids.thyroid, ]
rownames(gtex_thyroid_fpkm) <- sapply(strsplit(rownames(gtex_thyroid_fpkm), split=".", fixed=T), function(x)(x[1]))
gtex_thyroid_meta <- gtex.data$thyroid$meta

gtex_thyroid_counts <- gtex_thyroid_counts %>%
	mutate(gene = rownames(gtex_thyroid_counts)) %>%
	dplyr::select(gene, everything())

gtex_thyroid_fpkm <- gtex_thyroid_fpkm %>%
	mutate(gene = rownames(gtex_thyroid_fpkm)) %>%
	dplyr::select(gene, everything())

write.table(gtex_thyroid_counts, "./files/gtex_thyroid_counts.txt", quote=F, sep="\t", row.names=F)
write.table(gtex_thyroid_fpkm, "./files/gtex_thyroid_fpkm.txt", quote=F, sep="\t", row.names=F)
write.table(gtex_thyroid_meta, "./files/gtex_thyroid_meta.txt", quote=F, sep="\t", row.names=F)



#stomach
gtex_stomach_counts <- gtex.data$stomach$counts[keep.ensids.stomach, ]
rownames(gtex_stomach_counts) <- sapply(strsplit(rownames(gtex_stomach_counts), split=".", fixed=T), function(x)(x[1]))
gtex_stomach_fpkm <- gtex.data$stomach$rpkm[keep.ensids.stomach, ]
rownames(gtex_stomach_fpkm) <- sapply(strsplit(rownames(gtex_stomach_fpkm), split=".", fixed=T), function(x)(x[1]))
gtex_stomach_meta <- gtex.data$stomach$meta

gtex_stomach_counts <- gtex_stomach_counts %>%
	mutate(gene = rownames(gtex_stomach_counts)) %>%
	dplyr::select(gene, everything())

gtex_stomach_fpkm <- gtex_stomach_fpkm %>%
	mutate(gene = rownames(gtex_stomach_fpkm)) %>%
	dplyr::select(gene, everything())

write.table(gtex_stomach_counts, "./files/gtex_stomach_counts.txt", quote=F, sep="\t", row.names=F)
write.table(gtex_stomach_fpkm, "./files/gtex_stomach_fpkm.txt", quote=F, sep="\t", row.names=F)
write.table(gtex_stomach_meta, "./files/gtex_stomach_meta.txt", quote=F, sep="\t", row.names=F)





#========================
#merge tcga and gtex data
#========================


# thyroid
thyroid_tcga_gtex_counts <- inner_join(tcga_thca_counts, gtex_thyroid_counts, by="gene")
thyroid_tcga_gtex_fpkm <- inner_join(tcga_thca_fpkm, gtex_thyroid_fpkm, by="gene")

gtex_thyroid_meta <- gtex_thyroid_meta %>%
	mutate(sample_type = "GTEx")

thyroid_tcga_gtex_meta <- rbind(setNames(tcga_thca_meta[, c("submitter_id", "gender", "ethnicity", "age_at_diagnosis", "sample_type_id")], c("sample", "gender", "ethnicity", "age", "sample_type")),
	setNames(gtex_thyroid_meta[, c("sel.samples.sel.tissues.iddots", "GENDER", "ETHNCTY", "AGE", "sample_type")], c("sample", "gender", "ethnicity", "age", "sample_type")))

gender_code <- c(female = "female", male = "male", "1" = "male", "2" = "female")
ethnicity_code <- c(`hispanic or latino` = "hispanic or latino", `not hispanic or latino` = "not hispanic or latino", unknown = "unknown", "0" = "not hispanic or latino", "1" = "hispanic or latino", "99" = "unknown", "98" = "unknown")
sample_type_code <- c(GTEx = "GTEx", "01" = "TCGA_tumour", "11" = "TCGA_normal", "06" = "TCGA_metastatic")

thyroid_tcga_gtex_meta <- thyroid_tcga_gtex_meta %>%
	mutate(gender = gender_code[thyroid_tcga_gtex_meta$gender]) %>%
	mutate(ethnicity = ethnicity_code[thyroid_tcga_gtex_meta$ethnicity]) %>%
	mutate(sample_type = sample_type_code[thyroid_tcga_gtex_meta$sample_type])

write.table(thyroid_tcga_gtex_counts, "./files/thyroid_tcga_gtex_counts.txt", quote=F, sep="\t", row.names=F)
write.table(thyroid_tcga_gtex_fpkm, "./files/thyroid_tcga_gtex_fpkm.txt", quote=F, sep="\t", row.names=F)
write.table(thyroid_tcga_gtex_meta, "./files/thyroid_tcga_gtex_meta.txt", quote=F, sep="\t", row.names=F)



# stomach
stomach_tcga_gtex_counts <- inner_join(tcga_stad_counts, gtex_stomach_counts, by="gene")
stomach_tcga_gtex_fpkm <- inner_join(tcga_stad_fpkm, gtex_stomach_fpkm, by="gene")

gtex_stomach_meta <- gtex_stomach_meta %>%
	mutate(sample_type = "GTEx")

stomach_tcga_gtex_meta <- rbind(setNames(tcga_stad_meta[, c("submitter_id", "gender", "ethnicity", "age_at_diagnosis", "sample_type_id")], c("sample", "gender", "ethnicity", "age", "sample_type")),
	setNames(gtex_stomach_meta[, c("sel.samples.sel.tissues.iddots", "GENDER", "ETHNCTY", "AGE", "sample_type")], c("sample", "gender", "ethnicity", "age", "sample_type")))

stomach_tcga_gtex_meta <- stomach_tcga_gtex_meta %>%
	mutate(gender = gender_code[stomach_tcga_gtex_meta$gender]) %>%
	mutate(ethnicity = ethnicity_code[stomach_tcga_gtex_meta$ethnicity]) %>%
	mutate(sample_type = sample_type_code[stomach_tcga_gtex_meta$sample_type])

write.table(stomach_tcga_gtex_counts, "./files/stomach_tcga_gtex_counts.txt", quote=F, sep="\t", row.names=F)
write.table(stomach_tcga_gtex_fpkm, "./files/stomach_tcga_gtex_fpkm.txt", quote=F, sep="\t", row.names=F)
write.table(stomach_tcga_gtex_meta, "./files/stomach_tcga_gtex_meta.txt", quote=F, sep="\t", row.names=F)





save(list=ls(), file="./r_workspaces/tcga_gtex_rnaseq_preprocess.RData")

