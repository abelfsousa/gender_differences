#Understanding Gender Differential Susceptibility in Cancer

#GDC Data Portal

#**LEGACY TCGA DATA**#

#metadata processing
#returns list containing two data.frames:
#1)information by expression file
#2)information by sample (person)
retrieve.metaData <- function(metadata.json, flag){
	library("rjson")
	metadata <- fromJSON(file = metadata.json)
	exprFiles.info <- c()
	samplesClinical.info <- c()
	for(i in 1:length(metadata)){
		file_name <- metadata[[i]]$file_name
		if("unnormalized" %in% metadata[[i]]$tags){ file_type <- "gene.unnormalized" }
		else if("normalized" %in% metadata[[i]]$tags){ file_type <- "gene.normalized" }
		else if("exon" %in% metadata[[i]]$tags & !"junction" %in% metadata[[i]]$tags){ file_type <- "exon.quantification" }
		else if("exon" %in% metadata[[i]]$tags & "junction" %in% metadata[[i]]$tags){ file_type <- "exon.junction_quantification" }
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
			age <- NA
		}
		else{
			age <- as.numeric(substring(Sys.Date(),1,4))-metadata[[i]]$cases[[1]]$demographic$year_of_birth
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
		portion <- substring(metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 18, 19)
		plate <- substring(metadata[[i]]$associated_entities[[1]]$entity_submitter_id, 22, 25)
				
		samplesClinical.info <- rbind(samplesClinical.info, c(submitter_id, sample_type_id, gender, race, ethnicity, age, age_at_diagnosis, tumor_stage, vital_status, tss, portion, plate))
	}
	
	exprFiles.info <- as.data.frame(exprFiles.info, stringsAsFactors=F)
	colnames(exprFiles.info) <- c("file_name", "file_type", "submitter_id", "sample_type_id", "sample_type")
	
	samplesClinical.info <- as.data.frame(samplesClinical.info, stringsAsFactors=F)
	colnames(samplesClinical.info) <- c("submitter_id", "sample_type_id", "gender", "race", "ethnicity", "age", "age_at_diagnosis", "tumor_stage", "vital_status", "tss", "portion", "plate")
	samplesClinical.info <- unique(samplesClinical.info)
	rownames(samplesClinical.info) <- samplesClinical.info$submitter_id
	
	write.table(exprFiles.info, paste(flag, "exprFiles_info.txt", sep = "."), quote=F, sep="\t", row.names=F)
	write.table(samplesClinical.info, paste(flag, "samplesClinical_info.txt", sep = "."), quote=F, sep="\t", row.names=F)
	
	meta <- list(exprFiles.info = exprFiles.info, samplesClinical.info = samplesClinical.info)
	return(meta)
}


#expression data processing
#returns list containing three data.frames:
#1)raw counts for all genes (rows) and samples (columns)
#2)normalized counts for all genes (rows) and samples (columns)
#3)raw counts for all exons (rows) and samples (columns)
retrieve.expressionData <- function(exprFiles.info, files.dir, geneIDs, exonIDs, junctionIDs){
	gene.counts <- data.frame(row.names = geneIDs)
	gene.norm <- data.frame(row.names = geneIDs)
	exon.counts <- data.frame(row.names = exonIDs)
	exon.rpkm <- data.frame(row.names = exonIDs)
	exon.junction <- data.frame(row.names = unique(junctionIDs))
	for(i in 1:nrow(exprFiles.info)){
		print(i)
		if(exprFiles.info[i,"file_type"] == "gene.unnormalized"){
			expr.file <- read.table(paste(files.dir, as.character(exprFiles.info[i,"file_name"]), sep=""), header=T)
			gene.counts <- cbind(gene.counts, round(expr.file$raw_count))
			colnames(gene.counts)[ncol(gene.counts)] <- as.character(exprFiles.info[i,"submitter_id"])
		}
		else if(exprFiles.info[i,"file_type"] == "gene.normalized"){
			expr.file <- read.table(paste(files.dir, as.character(exprFiles.info[i,"file_name"]), sep=""), header=T)
			gene.norm <- cbind(gene.norm, expr.file$normalized_count)
			colnames(gene.norm)[ncol(gene.norm)] <- as.character(exprFiles.info[i,"submitter_id"])
		}
		else if(exprFiles.info[i,"file_type"] == "exon.quantification"){
			expr.file <- read.table(paste(files.dir, as.character(exprFiles.info[i,"file_name"]), sep=""), header=T)
			exon.counts <- cbind(exon.counts, round(expr.file$raw_counts))
			colnames(exon.counts)[ncol(exon.counts)] <- as.character(exprFiles.info[i,"submitter_id"])
			exon.rpkm <- cbind(exon.rpkm, expr.file$RPKM)
			colnames(exon.rpkm)[ncol(exon.rpkm)] <- as.character(exprFiles.info[i,"submitter_id"])
		}
		else if(exprFiles.info[i,"file_type"] == "exon.junction_quantification"){
			expr.file <- read.table(paste(files.dir, as.character(exprFiles.info[i,"file_name"]), sep=""), header=T)
			expr.file <- expr.file[match(unique(as.character(expr.file$junction)), as.character(expr.file$junction)), ]
			exon.junction <- cbind(exon.junction, round(expr.file$raw_counts))
			colnames(exon.junction)[ncol(exon.junction)] <- as.character(exprFiles.info[i,"submitter_id"])
		}
	}
	expression.data <- list(gene.counts = gene.counts, gene.norm = gene.norm, exon.counts = exon.counts, exon.rpkm = exon.rpkm, junction.counts = exon.junction)
	return(expression.data)
}


#filter out the low-expressed features
#check the Levels of metaData$exprFiles.info$sample_type_id to get the types of samples available (ex: 01/11/06)
#returns a vector containing the feature IDs that passed the filter
filter.feature <- function(expression.counts, samplesClinical.info, IDs, cpm.cutO, p.cutO){
	tumour.data <- expression.counts[,substring(colnames(expression.counts), 14, 15) %in% c("01")]
	normal.data <- expression.counts[,substring(colnames(expression.counts), 14, 15) %in% c("11")]
	
	m.tumour.data <- tumour.data[,substring(colnames(tumour.data),1,12) %in% rownames(samplesClinical.info[samplesClinical.info$gender == "male",])]
	f.tumour.data <- tumour.data[,substring(colnames(tumour.data),1,12) %in% rownames(samplesClinical.info[samplesClinical.info$gender == "female",])]

	m.normal.data <- normal.data[,substring(colnames(normal.data),1,12) %in% rownames(samplesClinical.info[samplesClinical.info$gender == "male",])]
	f.normal.data <- normal.data[,substring(colnames(normal.data),1,12) %in% rownames(samplesClinical.info[samplesClinical.info$gender == "female",])]

	#select genes that have at least # reads per million mapped reads (# counts per million)
	#in at least #% of males OR females in tumour AND normal tissue
	library(edgeR)
	keep <- rowSums(cpm(m.tumour.data) > cpm.cutO) >= ncol(m.tumour.data)*p.cutO | rowSums(cpm(f.tumour.data) > cpm.cutO) >= ncol(f.tumour.data)*p.cutO & rowSums(cpm(m.normal.data) > cpm.cutO) >= ncol(m.normal.data)*p.cutO | rowSums(cpm(f.normal.data) > cpm.cutO) >= ncol(f.normal.data)*p.cutO
	#keep <- rowSums(cpm(m.tumour.data) > cpm.cutO) >= ncol(m.tumour.data)*p.cutO | rowSums(cpm(f.tumour.data) > cpm.cutO) >= ncol(f.tumour.data)*p.cutO
	#keep <- rowSums(cpm(m.normal.data) > cpm.cutO) >= ncol(m.normal.data)*p.cutO | rowSums(cpm(f.normal.data) > cpm.cutO) >= ncol(f.normal.data)*p.cutO

	keep.ids <- IDs[keep]
	return(keep.ids)
}


#returns the expression and clinical data for tumor and normal datasets, and also for counts and normalized data
#the list is organized in the following manner:
#gene/exon data -> counts/norm data -> tumor/normal data -> expression/clinical data
split.data <- function(gene.counts, gene.norm, exon.counts, exon.rpkm, junction.counts, samplesClinical.info){
	
	gene.counts.tumour.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="01" , ]
	gene.counts.normal.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="11", ]
	gene.counts.tumour<-gene.counts[,colnames(gene.counts) %in% gene.counts.tumour.clinical$submitter_id]
	gene.counts.normal<-gene.counts[,colnames(gene.counts) %in% gene.counts.normal.clinical$submitter_id]
	gene.counts.tumour.normal.clinical <- rbind(gene.counts.tumour.clinical,gene.counts.normal.clinical)
	sample_type_code <- c("01" = "tumour","11" = "normal" )
	gene.counts.tumour.normal.clinical<-cbind(gene.counts.tumour.normal.clinical,data.frame(sample_type=sample_type_code[gene.counts.tumour.normal.clinical$sample_type_id]))
	colnames(gene.counts.tumour) <- substring(colnames(gene.counts.tumour),1,12)
	colnames(gene.counts.normal) <- substring(colnames(gene.counts.normal),1,12)

	gene.norm.tumour.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="01" , ]
	gene.norm.normal.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="11", ]
	gene.norm.tumour<-gene.norm[,colnames(gene.norm) %in% gene.norm.tumour.clinical$submitter_id]
	gene.norm.normal<-gene.norm[,colnames(gene.norm) %in% gene.norm.normal.clinical$submitter_id]
	gene.norm.tumour.normal.clinical <- rbind(gene.norm.tumour.clinical,gene.norm.normal.clinical)
	sample_type_code <- c("01" = "tumour","11" = "normal" )
	gene.norm.tumour.normal.clinical<-cbind(gene.norm.tumour.normal.clinical,data.frame(sample_type=sample_type_code[gene.norm.tumour.normal.clinical$sample_type_id]))
	colnames(gene.norm.tumour) <- substring(colnames(gene.norm.tumour),1,12)
	colnames(gene.norm.normal) <- substring(colnames(gene.norm.normal),1,12)
	rownames(gene.norm.tumour.clinical) <- substring(rownames(gene.norm.tumour.clinical),1,12)
	rownames(gene.norm.normal.clinical) <- substring(rownames(gene.norm.normal.clinical),1,12)

	exon.counts.tumour.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="01" , ]
	exon.counts.normal.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="11", ]
	exon.counts.tumour<-exon.counts[,colnames(exon.counts) %in% exon.counts.tumour.clinical$submitter_id]
	exon.counts.normal<-exon.counts[,colnames(exon.counts) %in% exon.counts.normal.clinical$submitter_id]
	exon.counts.tumour.normal.clinical <- rbind(exon.counts.tumour.clinical,exon.counts.normal.clinical)
	sample_type_code <- c("01" = "tumour","11" = "normal" )
	exon.counts.tumour.normal.clinical<-cbind(exon.counts.tumour.normal.clinical,data.frame(sample_type=sample_type_code[exon.counts.tumour.normal.clinical$sample_type_id]))
	colnames(exon.counts.tumour) <- substring(colnames(exon.counts.tumour),1,12)
	colnames(exon.counts.normal) <- substring(colnames(exon.counts.normal),1,12)
	rownames(exon.counts.tumour.clinical) <- substring(rownames(exon.counts.tumour.clinical),1,12)
	rownames(exon.counts.normal.clinical) <- substring(rownames(exon.counts.normal.clinical),1,12)

	exon.rpkm.tumour.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="01" , ]
	exon.rpkm.normal.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="11", ]
	exon.rpkm.tumour<-exon.rpkm[,colnames(exon.rpkm) %in% exon.rpkm.tumour.clinical$submitter_id]
	exon.rpkm.normal<-exon.rpkm[,colnames(exon.rpkm) %in% exon.rpkm.normal.clinical$submitter_id]
	exon.rpkm.tumour.normal.clinical <- rbind(exon.rpkm.tumour.clinical,exon.rpkm.normal.clinical)
	sample_type_code <- c("01" = "tumour","11" = "normal" )
	exon.rpkm.tumour.normal.clinical<-cbind(exon.rpkm.tumour.normal.clinical,data.frame(sample_type=sample_type_code[exon.rpkm.tumour.normal.clinical$sample_type_id]))
	colnames(exon.rpkm.tumour) <- substring(colnames(exon.rpkm.tumour),1,12)
	colnames(exon.rpkm.normal) <- substring(colnames(exon.rpkm.normal),1,12)
	rownames(exon.rpkm.tumour.clinical) <- substring(rownames(exon.rpkm.tumour.clinical),1,12)
	rownames(exon.rpkm.normal.clinical) <- substring(rownames(exon.rpkm.normal.clinical),1,12)

	
	junction.counts.tumour.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="01" , ]
	junction.counts.normal.clinical <- samplesClinical.info[samplesClinical.info$sample_type_id=="11", ]
	junction.counts.tumour<-junction.counts[,colnames(junction.counts) %in% junction.counts.tumour.clinical$submitter_id]
	junction.counts.normal<-junction.counts[,colnames(junction.counts) %in% junction.counts.normal.clinical$submitter_id]
	junction.counts.tumour.normal.clinical <- rbind(junction.counts.tumour.clinical,junction.counts.normal.clinical)
	sample_type_code <- c("01" = "tumour","11" = "normal" )
	junction.counts.tumour.normal.clinical<-cbind(junction.counts.tumour.normal.clinical,data.frame(sample_type=sample_type_code[junction.counts.tumour.normal.clinical$sample_type_id]))
	colnames(junction.counts.tumour) <- substring(colnames(junction.counts.tumour),1,12)
	colnames(junction.counts.normal) <- substring(colnames(junction.counts.normal),1,12)
	rownames(junction.counts.tumour.clinical) <- substring(rownames(junction.counts.tumour.clinical),1,12)
	rownames(junction.counts.normal.clinical) <- substring(rownames(junction.counts.normal.clinical),1,12)

	data <- list(gene = NULL, exon = NULL, junction = NULL)

	data$gene <- list(counts = NULL, norm = NULL)
	data$exon <- list(counts = NULL, rpkm = NULL)
	data$junction <- list(counts = NULL)

	data$gene$counts <- list(tumour = NULL, normal = NULL, tumournormal= NULL)
	data$gene$counts$tumour <- list(expr = gene.counts.tumour, clin = gene.counts.tumour.clinical)
	data$gene$counts$normal <- list(expr = gene.counts.normal, clin = gene.counts.normal.clinical)
	data$gene$counts$tumour.normal <- list(expr = gene.counts, clin = gene.counts.tumour.normal.clinical)

	data$gene$norm <- list(tumour = NULL, normal = NULL, tumour.normal= NULL)
	data$gene$norm$tumour <- list(expr = gene.norm.tumour, clin = gene.norm.tumour.clinical)
	data$gene$norm$normal <- list(expr = gene.norm.normal, clin = gene.norm.normal.clinical)
	data$gene$norm$tumour.normal <- list(expr = gene.norm, clin = gene.norm.tumour.normal.clinical)
	
	data$exon$counts <- list(tumour = NULL, normal = NULL, tumour.normal= NULL)
	data$exon$counts$tumour <- list(expr = exon.counts.tumour, clin = exon.counts.tumour.clinical)
	data$exon$counts$normal <- list(expr = exon.counts.normal, clin = exon.counts.normal.clinical)
	data$exon$counts$tumour.normal <- list(expr = exon.counts, clin = exon.counts.tumour.normal.clinical)

	data$exon$rpkm <- list(tumour = NULL, normal = NULL, tumour.normal= NULL)
	data$exon$rpkm$tumour <- list(expr = exon.rpkm.tumour, clin = exon.rpkm.tumour.clinical)
	data$exon$rpkm$normal <- list(expr = exon.rpkm.normal, clin = exon.rpkm.normal.clinical)
	data$exon$rpkm$tumour.normal <- list(expr = exon.rpkm, clin = exon.rpkm.tumour.normal.clinical)
	
	data$junction$counts <- list(tumour = NULL, normal = NULL, tumour.normal= NULL)
	data$junction$counts$tumour <- list(expr = junction.counts.tumour, clin = junction.counts.tumour.clinical)
	data$junction$counts$normal <- list(expr = junction.counts.normal, clin = junction.counts.normal.clinical)
	data$junction$counts$tumour.normal <- list(expr = junction.counts, clin = junction.counts.tumour.normal.clinical)

	return(data)
}


run.limma.diff.spl <- function(expression.data, design, group,exon.annot, flag, method){
	library(edgeR)
	library(limma)

	exon.annot <- exon.annot[exon.annot$exon %in% rownames(expression.data),]
	exon.annot <- cbind(exon.annot, sapply(strsplit(exon.annot$exon, split=":"), function(x){x[1]}))
	colnames(exon.annot)[4] <- "chrom"

	#create a DGEList object
	y <- DGEList(counts=expression.data, group=group, genes=exon.annot)


	#normalization TMM
	y <- calcNormFactors(y)
	v <- voom(y,design,plot=FALSE)
	#estimating NB dispertions
	#fitting the negative binomial GLM at the exon level
	#fit <- glmQLFit(y, design, robust = TRUE)
	fit <- lmFit(v, design)
	fit.de <- eBayes(fit)
	#test for differentially expressed exons between genders
	#note that we can omit the coef argument because the gender effect is the last coefficient in the model
	
	deTags <- topTable(fit.de)
	#dees = differentially expressed exons
	dees <- deTags[[1]]
	summary(decideTests(fit.de))

	#alternative splicing
	#test for differential exon usage for each gene between genders
	#testing for differential exon usage is equivalent to testing whether the exons in each gene have the same log-fold-changes as
    #the other exons in the same gene
	#note that we can omit the coef argument because the gender effect is the last coefficient in the model
	sp <- diffSplice(fit, geneid = "gene", exonid = "exon")

    #at exon-level
    #the exon-level tests test for differences between each exon and all the exons for the same gene
    #the log-fold-change of each exon is compared to the log-fold-change of the entire gene which contains that exon
    #dsel = differential splicing exon-level
	dsel <- topSplice(sp, test="t", number=Inf, FDR=Inf)
	
	#at gene-level
	#the gene-level tests test for any differences in exon usage between experimental conditions
	#converting exon-level p-values to gene-level p-values by the Simes method
	#the Simes method processes the exon-level p-values to give an overall call of differential splicing for each gene. It returns the minimum Simesusted p-values for each gene
	#the Simes p-values is likely to be more powerful when only a minority of the exons for a gene are differentially spliced (otherwise F-tests ("gene") are likely to be more powerful)
	#dsgl = differential splicing gene-level
	dsgl <- topSplice(sp, test="simes", number=Inf, FDR=Inf)


	file.name <- paste(flag, "limma", "diff.exon.expr.txt", sep = ".")
	write.table(dees, file.name, row.names=F, col.names=T, quote=F, sep="\t")

	file.name <- paste(flag, "limma","diff.splicing.exon-level.txt", sep = ".")
	write.table(dsel, file.name, row.names=F, col.names=T, quote=F, sep="\t")

	file.name <- paste(flag, "limma", "diff.splicing.gene-level-Simes.txt", sep = ".")
	write.table(dsgl, file.name, row.names=F, col.names=T, quote=F, sep="\t")

	
	data <- list(diff.ex.expr = list(dees), diff.ex.usage = list(dsel, dsgl, sp), y, fit)

	return(data)
}


#read required files and define expression files directory
geneIDs <- read.table("/home/mferreira/Abel_project/splicing_data/geneIDs.txt", sep="\t", h=F, stringsAsFactors=F)$V1
exonIDs <- read.table("/home/mferreira/Abel_project/splicing_data/exonIDs.txt", sep="\t", h=F, stringsAsFactors=F)$V1
junctionIDs <- read.table("/home/mferreira/Abel_project/splicing_data/junctionIDs.txt", sep="\t", h=F, stringsAsFactors=F)$V1
dup.junctions <- junctionIDs[duplicated(junctionIDs)]
geneIDs.annot <- read.table("/home/mferreira/Abel_project/splicing_data/geneAnnot.TCGA.hg19.June2011.txt", sep="\t", h=T, stringsAsFactors=F)
junctionIDs.annot <- read.table("/home/mferreira/Abel_project/splicing_data/junctionAnnot.TCGA.hg19.June2011.txt", sep="\t", h=T, stringsAsFactors=F)
exonIDs.annot <- read.table("/home/mferreira/Abel_project/splicing_data/exonAnnot.TCGA.hg19.June2011_compositeExon.txt", sep="\t", h=T, stringsAsFactors=F)


meta.file.STAD <- "/home/mferreira/Abel_project/splicing_data/STAD/metadata.cart.2017-01-14T12-41-03.153700.json"
files.dir.STAD <- "/home/mferreira/Abel_project/splicing_data/STAD/files/"



metaData.STAD <- retrieve.metaData(meta.file.STAD, "STAD")
stad_paper.info <- read.csv("/home/mferreira/Abel_project/New_analises_abel/data/from_paper/stad_additional_info.csv", sep = ";", h=T, stringsAsFactors=F)
metaData.STAD$samplesClinical.info <- merge(metaData.STAD$samplesClinical.info, stad_paper.info, by.x = "submitter_id", by.y = "TCGA_barcode")
rownames(metaData.STAD$samplesClinical.info) <- metaData.STAD$samplesClinical.info$submitter_id
metaData.STAD$samplesClinical.info<-metaData.STAD$samplesClinical.info[,c(1:14)]


any(is.na(metaData.STAD$samplesClinical.info))
#[1] TRUE

metaData.STAD$samplesClinical.info[is.na(metaData.STAD$samplesClinical.info[,"race"]),"race"] <- "unknown"
metaData.STAD$samplesClinical.info[is.na(metaData.STAD$samplesClinical.info[,"ethnicity"]),"ethnicity"] <- "unknown"
metaData.STAD$samplesClinical.info[is.na(metaData.STAD$samplesClinical.info[,"tumor_stage"]),"tumor_stage"] <- "unknown"
metaData.STAD$samplesClinical.info[is.na(metaData.STAD$samplesClinical.info[,"age"]),"age"] <- as.character(round(mean(as.numeric(metaData.STAD$samplesClinical.info[!is.na(metaData.STAD$samplesClinical.info[,"age"]),"age"]))))
metaData.STAD$samplesClinical.info[is.na(metaData.STAD$samplesClinical.info[,"age_at_diagnosis"]),"age_at_diagnosis"] <- "unknown"
metaData.STAD$samplesClinical.info[is.na(metaData.STAD$samplesClinical.info[,"tumor_stage"]),"tumor_stage"] <- "unknown"
metaData.STAD$samplesClinical.info[is.na(metaData.STAD$samplesClinical.info[,"vital_status"]),"vital_status"] <- "unknown"
metaData.STAD$samplesClinical.info[is.na(metaData.STAD$samplesClinical.info[,"Lauren_Class"]),"Lauren_Class"] <- "unknown"
metaData.STAD$samplesClinical.info$Lauren_Class[metaData.STAD$samplesClinical.info$Lauren_Class=="None"] <- "unknown"
metaData.STAD$samplesClinical.info[is.na(metaData.STAD$samplesClinical.info[,"Molecular_Subtype"]),"Molecular_Subtype"] <- "unknown"
metaData.STAD$samplesClinical.info$race<-as.character(metaData.STAD$samplesClinical.info$race)
metaData.STAD$samplesClinical.info$ethnicity<-as.character(metaData.STAD$samplesClinical.info$ethnicity)
metaData.STAD$samplesClinical.info$age_at_diagnosis<-as.numeric(metaData.STAD$samplesClinical.info$age_at_diagnosis)
metaData.STAD$samplesClinical.info$tumor_stage<-as.character(metaData.STAD$samplesClinical.info$tumor_stage)
metaData.STAD$samplesClinical.info$Lauren_Class<-as.character(metaData.STAD$samplesClinical.info$Lauren_Class)
metaData.STAD$samplesClinical.info$tss<-as.character(metaData.STAD$samplesClinical.info$tss)
metaData.STAD$samplesClinical.info$portion<-as.character(metaData.STAD$samplesClinical.info$portion)
metaData.STAD$samplesClinical.info$plate<-as.character(metaData.STAD$samplesClinical.info$plate)
metaData.STAD$samplesClinical.info$gender<-as.character(metaData.STAD$samplesClinical.info$gender)

any(is.na(metaData.STAD$samplesClinical.info))
#[1] FALSE



# filtered genes, exons and junxtion annotatios by genes present in diferential expression only

DEGS.STAD<-read.table("genes_dges_STAD.txt",header=F,sep="\t")
head(DEGS.STAD)
ALL.genes<-read.table("/home/mferreira/Abel_project/splicing_data/all_genes.txt",header=T,sep="\t")
head(ALL.genes)

DEGS.annot<-ALL.genes[ALL.genes$geneID %in% DEGS.STAD[,1],]
head(DEGS.annot)

exonIDs.annot$geneName<-strsplit(exonIDs.annot$gene,"\\|.*")


exonIDs.annot_new<-exonIDs.annot[exonIDs.annot$geneName %in% DEGS.annot$geneName,]
exonIDs.annot_new$geneName<-as.character(exonIDs.annot_new$geneName)
dim(exonIDs.annot)
dim(exonIDs.annot_new)
exonIDs_new<-exonIDs[exonIDs %in% exonIDs.annot_new$exon]
length(exonIDs)
length(exonIDs_new)

geneIDs.annot$geneName<-strsplit(geneIDs.annot$gene,"\\|.*")
geneIDs.annot_new<-geneIDs.annot[geneIDs.annot$geneName %in% DEGS.annot$geneName,]
geneIDs.annot_new$geneName<-as.character(geneIDs.annot_new$geneName)
geneIDs_new<-geneIDs[geneIDs %in% geneIDs.annot_new$gene]
dim(geneIDs.annot)
dim(geneIDs.annot_new)
length(geneIDs)
length(geneIDs_new)

junctionIDs.annot_new<-junctionIDs.annot[junctionIDs.annot$gene %in% geneIDs.annot_new$gene,]
junctionIDs_new<-junctionIDs[junctionIDs %in% junctionIDs.annot_new$junction]
dim(junctionIDs.annot)
dim(junctionIDs.annot_new)
length(junctionIDs)
length(junctionIDs_new)

#read expression data
expressionData.STAD <- retrieve.expressionData(metaData.STAD$exprFiles.info, files.dir.STAD, geneIDs, exonIDs, junctionIDs)
expressionData.STAD_new <- expressionData.STAD
expressionData.STAD_new$gene.counts <- expressionData.STAD_new$gene.counts[rownames(expressionData.STAD_new$gene.counts) %in% geneIDs_new,]
expressionData.STAD_new$gene.norm <- expressionData.STAD_new$gene.norm[rownames(expressionData.STAD_new$gene.norm) %in% geneIDs_new,]
expressionData.STAD_new$exon.counts <- expressionData.STAD_new$exon.counts[rownames(expressionData.STAD_new$exon.counts) %in% exonIDs_new,]
expressionData.STAD_new$exon.rpkm <- expressionData.STAD_new$gene.counts[rownames(expressionData.STAD_new$exon.rpkm) %in% exonIDs_new,]
expressionData.STAD_new$junction.counts  <- expressionData.STAD_new$junction.counts[rownames(expressionData.STAD_new$junction.counts) %in% junctionIDs_new,]
dim(expressionData.STAD$gene.counts)
dim(expressionData.STAD$gene.norm)
dim(expressionData.STAD$exon.counts)
dim(expressionData.STAD$exon.rpkm )
dim(expressionData.STAD$junction.counts)
dim(expressionData.STAD_new$gene.counts)
dim(expressionData.STAD_new$gene.norm )
dim(expressionData.STAD_new$exon.counts)
dim(expressionData.STAD_new$exon.rpkm )
dim(expressionData.STAD_new$junction.counts)

#filter out low expressed genes
#always confirm how many genes are removed before proceed with the differential expression analyses
 
#20531->15603
keep.geneids <- filter.feature(expressionData.STAD_new$gene.counts, metaData.STAD$samplesClinical.info, geneIDs_new, 1, 0.1)
gene.counts <- expressionData.STAD_new$gene.counts[keep.geneids, ]
gene.norm <- expressionData.STAD_new$gene.norm[keep.geneids, ]


#filter out low expressed exons

#always confirm how many exons are removed before proceed with the differential splicing analyses
keep.exonids <- filter.feature(expressionData.STAD_new$exon.counts, metaData.STAD$samplesClinical.info, exonIDs_new, 1, 0.1)
#keep only the exons with gene information

keep.exonids <- intersect(keep.exonids, exonIDs.annot_new$exon)
exon.counts <- expressionData.STAD_new$exon.counts[keep.exonids, ]
exon.rpkm <- expressionData.STAD_new$exon.rpkm[keep.exonids, ]
#length(unique(exonIDs.annot[exonIDs.annot$exon %in% keep.exonids, "gene"]))


#prepare data for differential expression/splicing analysis
allData.STAD <- split.data(gene.counts, gene.norm, exon.counts, exon.rpkm, expressionData.STAD_new$junction.counts, metaData.STAD$samplesClinical.info)

rm(gene.counts, gene.norm)
rm(exon.counts, exon.rpkm)




#====================================
#differential exon splicing/inclusion
#====================================
#Tumor: Male vs Female

STAD.tumour.MvsF.edgeR <- run.edgeR.diff.spl(allData.STAD$exon$counts$tumour$expr,design.tumour.MvsF ,allData.STAD$exon$counts$tumour$clin$gender, exonIDs.annot_new, "STAD.Tumor.MvsF")

#limma adjusted
design.tumour.MvsF <- model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + Lauren_Class + portion + plate + gender, allData.STAD$exon$counts$tumour$clin)
STAD.tumour.MvsF.limma <- run.limma.diff.spl(allData.STAD$exon$counts$tumour$expr,design.tumour.MvsF ,allData.STAD$exon$counts$tumour$clin$gender, exonIDs.annot_new, "STAD.Tumor.MvsF")
STAD.tumour.MvsF.limma.degs.fil <- STAD.tumour.MvsF.limma$diff.ex.usage[[2]][STAD.tumour.MvsF.limma$diff.ex.usage[[2]]$FDR<0.05, ]
dim(STAD.tumour.MvsF.limma.degs.fil )
#[1] 1 6
STAD.tumour.MvsF.limma.degs.fil
#               gene geneName chrom NExons      P.Value         FDR
#171939 PHLPP2|23035   PHLPP2 chr16     22 8.505962e-07 0.009274051

png("STAD.tumour.MvsF.limma.PHLPP2.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.tumour.MvsF.limma$diff.ex.usage[[3]], geneid="PHLPP2", genecol="geneName")
dev.off()

#Normal: Male vs Female

design.normal.MvsF <- model.matrix(~ race + age_at_diagnosis + tumor_stage + portion + gender, allData.STAD$exon$counts$normal$clin)
STAD.normal.MvsF.limma <- run.limma.diff.spl(allData.STAD$exon$counts$normal$expr, design.normal.MvsF ,allData.STAD$exon$counts$normal$clin$gender, exonIDs.annot_new, "STAD.Normal.MvsF")
STAD.normal.MvsF.limma.degs.fil <- STAD.normal.MvsF.limma$diff.ex.usage[[2]][STAD.normal.MvsF.limma$diff.ex.usage[[2]]$FDR<0.05, ]
dim(STAD.normal.MvsF.limma.degs.fil)
#[1] 11  6
STAD.normal.MvsF.limma.degs.fil
#                 gene geneName chrom NExons      P.Value          FDR
# 50908     PSMD2|5708    PSMD2  chr3     22 2.672066e-08 0.0002913354
# 4317   FAM76A|199870   FAM76A  chr1     10 4.083415e-07 0.0022260737
# 17480      DHX9|1660     DHX9  chr1     29 7.040089e-07 0.0025586030
# 70900    SIRT5|23408    SIRT5  chr6     12 1.411773e-06 0.0031741692
# 39945    SLC6A6|6533   SLC6A6  chr3     16 1.455640e-06 0.0031741692
# 197142   HNRNPL|3191   HNRNPL chr19     15 4.664501e-06 0.0084761749
# 166583   NOMO1|23420    NOMO1 chr16     31 6.263430e-06 0.0092575647
# 113081    VDAC2|7417    VDAC2 chr10     12 6.792673e-06 0.0092575647
# 46883     IQCB1|9657    IQCB1  chr3     15 1.567615e-05 0.0189907833
# 24752  SLC30A6|55676  SLC30A6  chr2     15 1.904579e-05 0.0207656210
# 67230        IK|3550       IK  chr5     20 3.672351e-05 0.0363996731

png("STAD.normal.MvsF.limma.PSMD2.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="PSMD2", genecol="geneName")
dev.off()
png("STAD.normal.MvsF.limma.FAM76A.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="FAM76A", genecol="geneName")
dev.off()
png("STAD.normal.MvsF.limma.DHX9.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="DHX9", genecol="geneName")
dev.off()
png("STAD.normal.MvsF.limma.SIRT5.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="SIRT5", genecol="geneName")
dev.off()
png("STAD.normal.MvsF.limma.SLC6A6.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="SLC6A6", genecol="geneName")
dev.off()
png("STAD.normal.MvsF.limma.HNRNPL.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="HNRNPL", genecol="geneName")
dev.off()
png("STAD.normal.MvsF.limma.NOMO1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="NOMO1", genecol="geneName")
dev.off()
png("STAD.normal.MvsF.limma.IQCB1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="IQCB1", genecol="geneName")
dev.off()
png("STAD.normal.MvsF.limma.SLC30A6.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="SLC30A6", genecol="geneName")
dev.off()
png("STAD.normal.MvsF.limma.IK.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.normal.MvsF.limma$diff.ex.usage[[3]], geneid="IK", genecol="geneName")
dev.off()
#Tumor vs Normal

#Tumor vs Normal:Female
exon_counts_tumour_normal_female_clin <- subset(allData.STAD$exon$counts$tumour.normal$clin, allData.STAD$exon$counts$tumour.normal$clin$gender=="female")
exon_counts_tumour_normal_female <- allData.STAD$exon$counts$tumour.normal$expr[,colnames(allData.STAD$exon$counts$tumour.normal$expr) %in% rownames(exon_counts_tumour_normal_female_clin)]

design.female.TvsN <-model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + portion + plate + sample_type_id, exon_counts_tumour_normal_female_clin)
STAD.female.TvsN.limma <- run.limma.diff.spl(exon_counts_tumour_normal_female, design.female.TvsN, exon_counts_tumour_normal_female_clin$sample_type_id, exonIDs.annot_new, "STAD.female.TvsN")
STAD.female.TvsN.limma.degs.fil <- STAD.female.TvsN.limma$diff.ex.usage[[2]][STAD.female.TvsN.limma$diff.ex.usage[[2]]$FDR<0.05, ]
dim(STAD.female.TvsN.limma.degs.fil)
#[1] 15  6
STAD.female.TvsN.limma.degs.fil
#                  gene geneName chrom NExons      P.Value          FDR
# 72896      ATF6B|1388    ATF6B  chr6     19 8.180488e-14 8.919186e-10
# 166718     NDE1|54820     NDE1 chr16     11 2.662850e-12 1.451653e-08
# 181717      SPOP|8405     SPOP chr17     15 2.202948e-08 8.006247e-05
# 152574       AKT1|207     AKT1 chr14     16 2.579118e-07 7.030031e-04
# 109768   NSUN6|221078    NSUN6 chr10     15 9.993643e-07 2.179214e-03
# 185316     CYTH1|9267    CYTH1 chr17     16 2.478532e-06 4.503906e-03
# 76282     RIMS1|22999    RIMS1  chr6     39 2.916247e-06 4.542263e-03
# 8387      DOCK7|85440    DOCK7  chr1     51 7.335314e-06 9.329857e-03
# 210429   MICAL3|57553   MICAL3 chr22     40 8.376792e-06 9.329857e-03
# 3927       STMN1|3925    STMN1  chr1      9 8.557146e-06 9.329857e-03
# 135880      NACA|4666     NACA chr12      9 1.179853e-05 1.167057e-02
# 52530       RNF4|6047     RNF4  chr4     10 1.284480e-05 1.167057e-02
# 213526 TMEM184B|25829 TMEM184B chr22     11 2.282640e-05 1.914433e-02
# 81916    RNF216|54476   RNF216  chr7     21 3.613614e-05 2.814231e-02
# 29503      UXS1|80146     UXS1  chr2     16 3.901862e-05 2.836133e-02

png("STAD.female.TvsN.limma.ATF6B.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="ATF6B", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.NDE1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="NDE1", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.SPOP.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="SPOP", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.AKT1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="AKT1", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.NSUN6.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="NSUN6", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.CYTH1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="CYTH1", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.RIMS1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="RIMS1", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.DOCK7.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="DOCK7", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.MICAL3.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="MICAL3", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.STMN1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="STMN1", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.NACA.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="NACA", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.RNF4.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="RNF4", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.TMEM184B.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="TMEM184B", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.RNF216.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="RNF216", genecol="geneName")
dev.off()
png("STAD.female.TvsN.limma.UXS1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.female.TvsN.limma$diff.ex.usage[[3]], geneid="UXS1", genecol="geneName")
dev.off()

#Tumor vs Normal:Male
exon_counts_tumour_normal_male_clin <- subset(allData.STAD$exon$counts$tumour.normal$clin, allData.STAD$exon$counts$tumour.normal$clin$gender=="male")
exon_counts_tumour_normal_male <- allData.STAD$exon$counts$tumour.normal$expr[,colnames(allData.STAD$exon$counts$tumour.normal$expr) %in% rownames(exon_counts_tumour_normal_male_clin)]

design.male.TvsN <-model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + portion + plate + sample_type_id, exon_counts_tumour_normal_male_clin)
STAD.male.TvsN.limma <- run.limma.diff.spl(exon_counts_tumour_normal_male, design.male.TvsN, exon_counts_tumour_normal_male_clin$sample_type_id, exonIDs.annot_new, "STAD.male.TvsN")
STAD.male.TvsN.limma.degs.fil <- STAD.male.TvsN.limma$diff.ex.usage[[2]][STAD.male.TvsN.limma$diff.ex.usage[[2]]$FDR<0.05, ]
dim(STAD.male.TvsN.limma.degs.fil)
#[1] 5 6
STAD.male.TvsN.limma.degs.fil
#                 gene geneName chrom NExons      P.Value         FDR
# 106649  FAM73B|84895   FAM73B  chr9     18 1.565568e-06 0.009765054
# 2218   DNAJC16|23341  DNAJC16  chr1     20 1.791260e-06 0.009765054
# 83756    VPS41|27072    VPS41  chr7     30 3.114260e-06 0.011318257
# 149952   DCAF4|26094    DCAF4 chr14     14 5.538409e-06 0.015096317
# 1048     ACOT7|11332    ACOT7  chr1     13 1.928978e-05 0.042063287

png("STAD.male.TvsN.limma.FAM73B.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.male.TvsN.limma$diff.ex.usage[[3]], geneid="FAM73B", genecol="geneName")
dev.off()
png("STAD.male.TvsN.limma.DNAJC16.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.male.TvsN.limma$diff.ex.usage[[3]], geneid="DNAJC16", genecol="geneName")
dev.off()
png("STAD.male.TvsN.limma.VPS41.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.male.TvsN.limma$diff.ex.usage[[3]], geneid="VPS41", genecol="geneName")
dev.off()
png("STAD.male.TvsN.limma.DCAF4.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.male.TvsN.limma$diff.ex.usage[[3]], geneid="DCAF4", genecol="geneName")
dev.off()
png("STAD.male.TvsN.limma.ACOT7.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(STAD.male.TvsN.limma$diff.ex.usage[[3]], geneid="ACOT7", genecol="geneName")
dev.off()


save(list=ls(), file="STAD.limma.spl.RData")
