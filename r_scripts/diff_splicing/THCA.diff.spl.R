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

run.limma.diff.spl <- function(expression.data, design, group,exon.annot, flag){
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


meta.file.THCA <- "/home/mferreira/Abel_project/splicing_data/THCA/metadata.cart.2017-01-14T12-35-35.540581.json"
files.dir.THCA <- "/home/mferreira/Abel_project/splicing_data/THCA/files/"



metaData.THCA <- retrieve.metaData(meta.file.THCA, "THCA")

metaData.THCA$samplesClinical.info<-metaData.THCA$samplesClinical.info[!(metaData.THCA$samplesClinical.info$sample_type_id=="06"),]
dim(metaData.THCA$samplesClinical.info)
#[1] 564  10
metaData.THCA$exprFiles.info<-metaData.THCA$exprFiles.info[!(metaData.THCA$exprFiles.info$sample_type_id=="06"),]
dim(metaData.THCA$exprFiles.info)
#[1] 2256    5


thca_paper.info <- read.csv("/home/mferreira/Abel_project/splicing_data/THCA/thca_additional_info.csv", sep = ";", h=T, stringsAsFactors=F, na.strings=c("","NA"))
thca_paper.info <- thca_paper.info[, -c(11)]
thca_paper.info$sample <- substring(thca_paper.info$sample, 1, 15)
metaData.THCA$samplesClinical.info <- merge(metaData.THCA$samplesClinical.info, thca_paper.info, by.x = "submitter_id", by.y = "sample")


rownames(metaData.THCA$samplesClinical.info) <-metaData.THCA$samplesClinical.info$submitter_id

metaData.THCA$samplesClinical.info<-metaData.THCA$samplesClinical.info[,c("submitter_id","sample_type_id","gender","race","ethnicity","age",
	"age_at_diagnosis","tumor_stage","tss","portion","plate","histological_type")]


any(is.na(metaData.THCA$samplesClinical.info))
#[1] TRUE

metaData.THCA$samplesClinical.info[is.na(metaData.THCA$samplesClinical.info[,"race"]),"race"] <- "unknown"
metaData.THCA$samplesClinical.info[is.na(metaData.THCA$samplesClinical.info[,"tumor_stage"]),"tumor_stage"] <- "unknown"
metaData.THCA$samplesClinical.info[is.na(metaData.THCA$samplesClinical.info[,"ethnicity"]),"ethnicity"] <- "unknown"
metaData.THCA$samplesClinical.info[is.na(metaData.THCA$samplesClinical.info[,"histological_type"]),"histological_type"] <- "unknown"
metaData.THCA$samplesClinical.info[is.na(metaData.THCA$samplesClinical.info[,"age"]),"age"] <- as.character(round(mean(as.numeric(metaData.THCA$samplesClinical.info[!is.na(metaData.THCA$samplesClinical.info[,"age"]),"age"]))))
metaData.THCA$samplesClinical.info$age<-as.numeric(metaData.THCA$samplesClinical.info$age)
metaData.THCA$samplesClinical.info$age_at_diagnosis<-as.numeric(metaData.THCA$samplesClinical.info$age_at_diagnosis)
metaData.THCA$samplesClinical.info$sample_type_id<-as.character(metaData.THCA$samplesClinical.info$sample_type_id)
metaData.THCA$samplesClinical.info$race<-as.character(metaData.THCA$samplesClinical.info$race)
metaData.THCA$samplesClinical.info$ethnicity<-as.character(metaData.THCA$samplesClinical.info$ethnicity)
metaData.THCA$samplesClinical.info$tumor_stage<-as.character(metaData.THCA$samplesClinical.info$tumor_stage)
metaData.THCA$samplesClinical.info$gender<-as.character(metaData.THCA$samplesClinical.info$gender)
metaData.THCA$samplesClinical.info$tss<-as.character(metaData.THCA$samplesClinical.info$tss)
metaData.THCA$samplesClinical.info$portion<-as.character(metaData.THCA$samplesClinical.info$portion)
metaData.THCA$samplesClinical.info$plate<-as.character(metaData.THCA$samplesClinical.info$plate)
metaData.THCA$samplesClinical.info$histological_type<-as.character(metaData.THCA$samplesClinical.info$histological_type)
any(is.na(metaData.THCA$samplesClinical.info))
#[1] FALSE


# filtered genes, exons and junxtion annotatios by genes present in diferential expression only

DEGS.THCA<-read.table("genes_dges_THCA.txt",header=F,sep="\t")
head(DEGS.THCA)
ALL.genes<-read.table("/home/mferreira/Abel_project/splicing_data/all_genes.txt",header=T,sep="\t")
head(ALL.genes)

DEGS.annot<-ALL.genes[ALL.genes$geneID %in% DEGS.THCA[,1],]
head(DEGS.annot)

exonIDs.annot$geneName<-strsplit(exonIDs.annot$gene,"\\|.*")
head(exonIDs.annot)

exonIDs.annot_new<-exonIDs.annot[exonIDs.annot$geneName %in% DEGS.annot$geneName,]
exonIDs.annot_new$geneName<-as.character(exonIDs.annot_new$geneName)
exonIDs_new<-exonIDs[exonIDs %in% exonIDs.annot_new$exon]

geneIDs.annot$geneName<-strsplit(geneIDs.annot$gene,"\\|.*")
geneIDs.annot_new<-geneIDs.annot[geneIDs.annot$geneName %in% DEGS.annot$geneName,]
geneIDs_new<-geneIDs[geneIDs %in% geneIDs.annot_new$gene]


junctionIDs.annot_new<-junctionIDs.annot[junctionIDs.annot$gene %in% geneIDs.annot_new$gene,]

junctionIDs_new<-junctionIDs[junctionIDs %in% junctionIDs.annot_new$junction]
length(junctionIDs_new)

#read expression data
expressionData.THCA <- retrieve.expressionData(metaData.THCA$exprFiles.info, files.dir.THCA, geneIDs, exonIDs, junctionIDs)

expressionData.THCA_new<-expressionData.THCA
expressionData.THCA_new$gene.counts <- expressionData.THCA_new$gene.counts[rownames(expressionData.THCA_new$gene.counts) %in% geneIDs_new,]
expressionData.THCA_new$gene.norm <- expressionData.THCA_new$gene.norm[rownames(expressionData.THCA_new$gene.norm) %in% geneIDs_new,]
expressionData.THCA_new$exon.counts <- expressionData.THCA_new$exon.counts[rownames(expressionData.THCA_new$exon.counts) %in% exonIDs_new,]
expressionData.THCA_new$exon.rpkm <- expressionData.THCA_new$gene.counts[rownames(expressionData.THCA_new$exon.rpkm) %in% exonIDs_new,]
expressionData.THCA_new$junction.counts  <- expressionData.THCA_new$junction.counts[rownames(expressionData.THCA_new$junction.counts) %in% junctionIDs_new,]
dim(expressionData.THCA$gene.counts)
dim(expressionData.THCA$gene.norm)
dim(expressionData.THCA$exon.counts)
dim(expressionData.THCA$exon.rpkm )
dim(expressionData.THCA$junction.counts)
dim(expressionData.THCA_new$gene.counts)
dim(expressionData.THCA_new$gene.norm )
dim(expressionData.THCA_new$exon.counts)
dim(expressionData.THCA_new$exon.rpkm )
dim(expressionData.THCA_new$junction.counts)



#filter out low expressed genes
#always confirm how many genes are removed before proceed with the differential expression analyses
 

keep.geneids <- filter.feature(expressionData.THCA_new$gene.counts, metaData.THCA$samplesClinical.info, geneIDs_new, 1, 0.1)
gene.counts <- expressionData.THCA_new$gene.counts[keep.geneids, ]
gene.norm <- expressionData.THCA_new$gene.norm[keep.geneids, ]


#filter out low expressed exons

#always confirm how many exons are removed before proceed with the differential splicing analyses
keep.exonids <- filter.feature(expressionData.THCA_new$exon.counts, metaData.THCA$samplesClinical.info, exonIDs_new, 1, 0.1)
#keep only the exons with gene information

keep.exonids <- intersect(keep.exonids, exonIDs.annot_new$exon)
exon.counts <- expressionData.THCA_new$exon.counts[keep.exonids, ]
exon.rpkm <- expressionData.THCA_new$exon.rpkm[keep.exonids, ]
#length(unique(exonIDs.annot[exonIDs.annot$exon %in% keep.exonids, "gene"]))


#prepare data for differential expression/splicing analysis
allData.THCA <- split.data(gene.counts, gene.norm, exon.counts, exon.rpkm, expressionData.THCA_new$junction.counts, metaData.THCA$samplesClinical.info)

rm(gene.counts, gene.norm)
rm(exon.counts, exon.rpkm)




#====================================
#differential exon splicing/inclusion
#====================================

#Tumor: Male vs Female



design.tumour.MvsF  <- model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + histological_type + tss + portion + plate + gender, allData.THCA$exon$counts$tumour$clin)
THCA.tumour.MvsF.limma <- run.limma.diff.spl(allData.THCA$exon$counts$tumour$expr,design.tumour.MvsF ,allData.THCA$exon$counts$tumour$clin$gender, exonIDs.annot_new, "THCA.Tumor.MvsF")
THCA.tumour.MvsF.limma.degs.fil <- THCA.tumour.MvsF.limma$diff.ex.usage[[2]][THCA.tumour.MvsF.limma$diff.ex.usage[[2]]$FDR<0.05, ]
dim(THCA.tumour.MvsF.limma.degs.fil )
#[1] 1 6
THCA.tumour.MvsF.limma.degs.fil
              gene geneName chrom NExons      P.Value          FDR
175305 EIF4A1|1973   EIF4A1 chr17     19 9.412904e-11 1.011699e-06


png("THCA.tumour.MvsF.limma.EIF4A1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.tumour.MvsF.limma$diff.ex.usage[[3]], geneid="EIF4A1", genecol="geneName")
dev.off()


#Normal: Male vs Female

design.normal.MvsF <- model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + portion + gender, allData.THCA$exon$counts$normal$clin)
THCA.normal.MvsF.limma <- run.limma.diff.spl(allData.THCA$exon$counts$normal$expr, design.normal.MvsF ,allData.THCA$exon$counts$normal$clin$gender, exonIDs.annot_new, "THCA.Normal.MvsF")
THCA.normal.MvsF.limma.degs.fil <- THCA.normal.MvsF.limma$diff.ex.usage[[2]][THCA.normal.MvsF.limma$diff.ex.usage[[2]]$FDR<0.05, ]
dim(THCA.normal.MvsF.limma.degs.fil)
#[1] 5 6

THCA.normal.MvsF.limma.degs.fil
#               gene geneName chrom NExons      P.Value        FDR
# 190930 PIAS4|51588    PIAS4 chr19     11 1.010963e-05 0.04012751
# 36610   RQCD1|9125    RQCD1  chr2     10 1.076041e-05 0.04012751
# 1411    LZIC|84328     LZIC  chr1     10 1.484415e-05 0.04012751
# 83133    CHN2|1124     CHN2  chr7     18 1.741125e-05 0.04012751
# 60410  SNX25|83891    SNX25  chr4     21 1.866743e-05 0.04012751



png("THCA.normal.MvsF.limma.PIAS4.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.normal.MvsF.limma$diff.ex.usage[[3]], geneid="PIAS4", genecol="geneName")
dev.off()

png("THCA.normal.MvsF.limma.RQCD1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.normal.MvsF.limma$diff.ex.usage[[3]], geneid="RQCD1", genecol="geneName")
dev.off()

png("THCA.normal.MvsF.limma.LZIC.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.normal.MvsF.limma$diff.ex.usage[[3]], geneid="LZIC", genecol="geneName")
dev.off()

png("THCA.normal.MvsF.limma.CHN2.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.normal.MvsF.limma$diff.ex.usage[[3]], geneid="CHN2", genecol="geneName")
dev.off()

png("THCA.normal.MvsF.limma.SNX25.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.normal.MvsF.limma$diff.ex.usage[[3]], geneid="SNX25", genecol="geneName")
dev.off()

#Tumor vs Normal:Female
exon_counts_tumour_normal_female_clin <- subset(allData.THCA$exon$counts$tumour.normal$clin,allData.THCA$exon$counts$tumour.normal$clin$gender=="female")
exon_counts_tumour_normal_female <- allData.THCA$exon$counts$tumour.normal$expr[,colnames(allData.THCA$exon$counts$tumour.normal$expr) %in% exon_counts_tumour_normal_female_clin$submitter_id]


design.female.TvsN <-model.matrix(~ race + ethnicity + age_at_diagnosis + tumor_stage + tss + portion + plate + sample_type_id, exon_counts_tumour_normal_female_clin)
THCA.female.TvsN.limma <- run.limma.diff.spl(exon_counts_tumour_normal_female, design.female.TvsN, exon_counts_tumour_normal_female_clin$sample_type, exonIDs.annot_new, "THCA.female.TvsN")
THCA.female.TvsN.limma.degs.fil <- THCA.female.TvsN.limma$diff.ex.usage[[2]][THCA.female.TvsN.limma$diff.ex.usage[[2]]$FDR<0.05, ]
dim(THCA.female.TvsN.limma.degs.fil)
#[1] 2 6

THCA.female.TvsN.limma.degs.fil
#               gene geneName chrom NExons      P.Value        FDR
# 150719  SNW1|22938     SNW1 chr14     13 2.393991e-06 0.02040671
# 4083   DHDDS|79947    DHDDS  chr1      9 3.797303e-06 0.02040671


png("THCA.female.TvsN.limma.SNW1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.female.TvsN.limma$diff.ex.usage[[3]], geneid="SNW1", genecol="geneName")
dev.off()

png("THCA.female.TvsN.limma.DHDDS.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.female.TvsN.limma$diff.ex.usage[[3]], geneid="DHDDS", genecol="geneName")
dev.off()


#Tumor vs Normal:Male
exon_counts_tumour_normal_male_clin <- subset(allData.THCA$exon$counts$tumour.normal$clin, allData.THCA$exon$counts$tumour.normal$clin$gender=="male")
exon_counts_tumour_normal_male <- allData.THCA$exon$counts$tumour.normal$expr[,colnames(allData.THCA$exon$counts$tumour.normal$expr) %in% exon_counts_tumour_normal_male_clin$submitter_id]

design.male.TvsN <-model.matrix(~ race + ethnicity + age_at_diagnosis +  tumor_stage + tss + portion + plate + sample_type, exon_counts_tumour_normal_male_clin)
THCA.male.TvsN.limma <- run.limma.diff.spl(exon_counts_tumour_normal_male, design.male.TvsN, exon_counts_tumour_normal_male_clin$sample_type, exonIDs.annot_new, "THCA.male.TvsN")
THCA.male.TvsN.limma.degs.fil <- THCA.male.TvsN.limma$diff.ex.usage[[2]][THCA.male.TvsN.limma$diff.ex.usage[[2]]$FDR<0.05, ]
dim(THCA.male.TvsN.limma.degs.fil)
#[1] 3 6
THCA.male.TvsN.limma.degs.fil
#               gene geneName chrom NExons      P.Value          FDR
# 218148 MAGED1|9500   MAGED1  chrX     15 7.119366e-12 7.651895e-08
# 191475  MLLT1|4298    MLLT1 chr19     12 2.480926e-11 1.333250e-07
# 145086 DZIP1|22873    DZIP1 chr13     24 7.951524e-11 2.848766e-07

png("THCA.male.TvsN.limma.MAGED1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.male.TvsN.limma$diff.ex.usage[[3]], geneid="MAGED1", genecol="geneName")
dev.off()

png("THCA.male.TvsN.limma.MLLT1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.male.TvsN.limma$diff.ex.usage[[3]], geneid="MLLT1", genecol="geneName")
dev.off()


png("THCA.male.TvsN.limma.DZIP1.png", width=800, height =650)
par(mar=c(10,5.1, 4.1, 2.1))
plotSplice(THCA.male.TvsN.limma$diff.ex.usage[[3]], geneid="DZIP1", genecol="geneName")
dev.off()

save(list=ls(), file="THCA.limma.spl.RData")

