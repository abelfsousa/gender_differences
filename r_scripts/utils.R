# Understanding Gender Differential Susceptibility in Cancer

# This script contains some useful functions



wilcoxon_diff_expression <- function(group1, group2, expression.data, gene.annot){

    #differential expression analysis using Wilcoxon signed-rank test
    #compares expression differences between two experimental groups

    # this function expects rpkm values as expression data
    # (at least some kind of library normalization)

    p.val <- c()
    log2FC <- c()

    for(gene in rownames(expression.data)){

        #cat(paste(gene, match(gene, rownames(expression.data)), sep=" "), sep = "\n")

        group1.data <- as.numeric( expression.data[ gene, colnames(expression.data) %in% group1 ] )
        group2.data <- as.numeric( expression.data[ gene, colnames(expression.data) %in% group2 ] )

        pvalue <- wilcox.test(group1.data, group2.data)$p.value
        fc <- log2( (median(group1.data)+0.05) / (median(group2.data)+0.05) )

        p.val <- c(p.val, pvalue)
        log2FC <- c(log2FC, fc)
    }

    FDR <- p.adjust(p.val, method="BH")
    gene.table <- cbind.data.frame(rownames(expression.data), log2FC, p.val, FDR)
    colnames(gene.table)[1] <- "genes"

    gene.table <- merge(gene.table, gene.annot, by.x = "genes", by.y = "geneID")
    gene.table <- gene.table[, c("genes", "geneName", "geneType", "chrom", "log2FC", "p.val", "FDR")]
    gene.table <- gene.table[order(gene.table$FDR,decreasing=F),]

    return(gene.table)
}


edgeR_diff_expression <- function(design, factor, expression.data, gene.annot){

	#differential expression analysis by edgeR

	#adjust for covariates when defining design matrix:
	#model.matrix(~ cov1 + cov2 + cov3 + ..., data)

  # this function expects raw counts values as expression data

	library(edgeR)

	#create a DGEList object
	y <- DGEList( counts=expression.data, group=factor, genes=rownames(expression.data) )

	#normalization TMM
	y <- calcNormFactors(y)
	#estimating NB dispertions
	y <- estimateDisp(y, design)
	#fitting the negative binomial GLM
	fit <- glmFit(y, design)
	#likelihood ratio test
	lrt <- glmLRT(fit)

	degs <- cbind(lrt$genes, lrt$table, p.adjust(lrt$table$PValue, method="BH"))
	colnames(degs)[6] <- "FDR"
	degs <- merge(degs, gene.annot, by.x = "genes", by.y = "geneID")
	degs <- degs[,c("genes","geneName","geneType","chrom","logFC","logCPM","LR","PValue","FDR")]
	degs <- degs[order(degs$FDR,decreasing=F),]

	return(degs)
}


sam_diff_expression <- function(expression.data, factor, gene.annot){

  # this function expects raw counts values as expression data

  library(samr)

  #fit the model
  samfit <- SAMseq(x = expression.data, y = factor, resp.type = "Two class unpaired", geneid = gene.annot$geneID, genenames = gene.annot$geneName)

  return(samfit)
}


noiseq_diff_expression <- function(expression.data, metadata, f, gene.annot){
	library(NOISeq)

  # this function expects raw counts values as expression data

	#create a NOISeq object
	noiseqbioData <- readData(data = expression.data, factors = metadata)

	#run NOISeqBIO
	noiseqbioRun <- noiseqbio(input = noiseqbioData, factor = f, norm = "tmm", filter = 0)

  vals <- levels(as.factor(metadata[, f]))

	degs <- noiseqbioRun@results[[1]]
	degs <- cbind(rownames(degs), degs)
	colnames(degs)[1] <- "genes"
	degs <- merge(degs, gene.annot, by.x = "genes", by.y = "geneID")
	degs <- degs[, c("genes","geneName","geneType","chrom",paste(vals[1], "mean", sep="_"),paste(vals[2], "mean", sep="_"),"theta","prob","log2FC")]
  degs$prob <- (1 - degs$prob)
  colnames(degs)[8] <- "FDR"
	degs <- degs[order(degs$FDR,decreasing=T),]

	return(degs)
}


limma_diff_expression <- function(design, factor, expression.data, gene.annot){

	#differential expression analysis by limma

	#adjust for covariates when defining design matrix:
	#model.matrix(~ cov1 + cov2 + cov3 + ..., data)

  # this function expects raw counts values as expression data

	library(edgeR)
  library(limma)

  # DGEList object using edgeR
  y <- DGEList(counts=expression.data, group=factor)

  # scale normalization using TMM method
  y <- calcNormFactors(y)

  # transform the counts using voom
  y <- voom(y, design)

  # fit the linear model for each gene
  fit <- lmFit(y, design)

  #empirical bayes statistics for differential expression
  fit <- eBayes(fit)

  # extract gene table with statistics
  diff_expr_table <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "p", p.value = 1, lfc = 0, adjust.method="BH")

  diff_expr_table <- diff_expr_table %>%
    mutate(genes = rownames(diff_expr_table)) %>%
    left_join(gene.annot, by=c("genes" = "geneID")) %>%
    dplyr::rename(FDR = adj.P.Val)

	return(diff_expr_table)
}





enrich_test <- function(gene_set, universe, terms1, terms2, p_adj, q_value){

  # functional enrichment function using clusterProfiler
  # hypergeometric test

  library(clusterProfiler)

  enr <- enricher(
    gene = gene_set,
    universe = universe,
    TERM2GENE = terms1,
    TERM2NAME = terms2,
    pvalueCutoff = p_adj,
    qvalueCutoff = q_value,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500)

  enr <- enr@result %>% dplyr::select(Description, ID, Count, p.adjust, GeneRatio, BgRatio, geneID)
  return(enr)
}


gsea_test <- function(data, terms1, terms2, p_adj){

  # functional enrichment function using clusterProfiler
  # gsea test

  library(clusterProfiler)

  data <- as.data.frame(data)
  glist <- data$kme
  names(glist) <- data$geneName

  tryCatch( {
    enr <- GSEA(
    geneList = glist,
    TERM2GENE = terms1,
    TERM2NAME = terms2,
    exponent = 1,
    nPerm = 1000,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = p_adj,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = FALSE,
    by = "fgsea")

    enr <- enr@result

    if( sum(dim(enr) == c(0, 8)) == 2 ){
      enr <- data.frame(ID=NA, Description=NA, setSize=NA, enrichmentScore=NA, NES=NA, pvalue=NA, p.adjust=NA, qvalues=NA, rank=NA, leading_edge=NA, core_enrichment=NA)
    }
    return(enr) },
    error = function(e) return(data.frame(ID=NA, Description=NA, setSize=NA, enrichmentScore=NA, NES=NA, pvalue=NA, p.adjust=NA, qvalues=NA, rank=NA, leading_edge=NA, core_enrichment=NA)) )
}


remove_batch1 <- function(gene, data, covar){

    # remove batch effects by regressing them out


    # gene expression
    y <- data[gene, ] %>% as.numeric()
    names(y) <- data[gene, ] %>% names()

    #data type
    cat <- covar

    # fit the model
    lreg <- lm(y ~ cat)

    return(residuals(lreg))
}


# function to regress-out a set of covariates from a variable fitting a linear model
# the first argument must be a data.frame with the following columns:
# column1 - sample identifiers; column2 - response variable; others - covariates to regress-out
# other covariates can be passed using the covs argument
# covs is a data.frame with one column "sample" containing the same sample identifiers as df
remove_batch2 <- function(df, covs = NULL, ztransf = FALSE){

  # additional covariates to regress-out
  if(!is.null(covs)){
    df <- df %>%
      inner_join(covs, by = "sample")
  }

  # remove NAs from the dataset
  dfM <- na.exclude(df)

  # remove sample column
  samples <- dfM$sample
  dfM <- dfM %>%
    select(-sample)

  # z-score transform numeric variables if specified
  if(ztransf){
    dfM <- dfM %>%
      mutate_if(.predicate = is.numeric, .funs = ~ scale(.x)[,1])
  }

  # define response and explanatory variables
  resp <- colnames(dfM)[1]
  expl <- colnames(dfM)[-1]


  # set up the formula
  f <- as.formula(paste0(resp, "~", paste(expl, collapse = "+")))

  # fit the model and get the residuals
  reg <- lm(f, dfM)
  reg_residual <- residuals(reg)

  res = tibble(sample = samples, residual = reg_residual)

  return(res)
}




get_lower_tri <- function(cormat){

  # get lower triangle of the correlation matrix

  cormat[upper.tri(cormat)] <- NA

  return(cormat)
}



get_upper_tri <- function(cormat){

  # get upper triangle of the correlation matrix

  cormat[lower.tri(cormat)]<- NA

  return(cormat)
}




# WGCNA functions


normalize_countTable <- function(countTable){

    # normalize a count table returned by overlapTable()
    # divide each value by the number of genes in the smaller module

    row_mod <- rowSums(countTable)
    col_mod <- colSums(countTable)

    countTable2 <- countTable

    for( row in 1:nrow(countTable2)){
        for (col in 1:ncol(countTable2)){
            countTable2[row, col] <- countTable2[row, col]/min(row_mod[row], col_mod[col])
        }
    }

    return(countTable2)
}




overlap_table_wgcna <- function(network1, network2, net1.label, net2.label, file.name, width, height){

    library(WGCNA)

    # calculate and plot the overlap between two networks

    network1.colors = labels2colors(network1$colors)
    network2.colors = labels2colors(network2$colors)

    # Calculate the contingency table and p-values
    overlap.table = overlapTable(network1.colors, network2.colors)
    # The numMat will encode color. We use -log of the p value.
    numMat = -log10(overlap.table$pTable)
    numMat[numMat >50] = 50
    # Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of
    # counts and corresponding p-values.
    textMat = paste(overlap.table$countTable, "\n", signif(overlap.table$pTable, 2))
    dim(textMat) = dim(numMat)
    # Additional information for the plot. These will be used shortly.
    yLabels = paste("M", sort(unique(network1.colors)))
    ySymbols = paste(sort(unique(network1.colors)), ": ", table(network1.colors), sep = "")
    xLabels = paste("M", sort(unique(network2.colors)))
    xSymbols = paste(sort(unique(network2.colors)), ": ", table(network2.colors), sep = "")

    # Plot the overlap table
    pdf(file = file.name, w = width, h = height)
    #par(oma=c(0,20,4,0))
    fp = TRUE
    fcex = 1.2
    pcex = 1.0
    fcexl = 1.5
    pcexl = 1.0
    par(mar = c(13, 15, 3, 0.5))
    labeledHeatmap(Matrix = numMat,
        yLabels = yLabels, ySymbols = ySymbols,
        xLabels = xLabels, xSymbols = xSymbols,
        colorLabels = TRUE,
        colors = greenWhiteRed(100)[50:100],
        textMatrix = textMat, cex.text = if (fp) fcex else pcex, setStdMargins = FALSE,
        cex.lab = if (fp) fcexl else pcexl,
        xColorWidth = 0.08,
        main = paste( paste(net1.label, "modules (rows)", sep=" "), paste(net2.label, "modules (columns)", sep=" "), sep=" vs. " ), cex.main = 1.5)
    dev.off()

    return(overlap.table)

}




networks_comparison <- function(network1.fpkm.mat, network2.fpkm.mat, network1, network2, overlap.table, pval = 0.05, kme = 0.7){

  library(WGCNA)

  #build a correlation table filled with zeros;
  #for common modules between datasets, a correlation of kME < 0.7 (correlation between the common genes) will indicate a loss of topology
  #while a correlation of kME >= 0.7 will indicate a preservation of topology

  network1.colors = labels2colors(network1$colors)
  network2.colors = labels2colors(network2$colors)

  color.code = standardColors(n = NULL)
  color.code = c("grey", color.code)

  cor.table.1 = as.data.frame( matrix( data = 0, nrow = length(table(network1.colors))-1, ncol = length(table(network2.colors))-1 ) )
  rownames(cor.table.1) = names(table(network1.colors))[which(names(table(network1.colors)) != "grey")]
  colnames(cor.table.1) = names(table(network2.colors))[which(names(table(network2.colors)) != "grey")]
  cor.table.2 = cor.table.1

  #calculate kME for net1 and net2 modules using signedKME function from WGCNA
  kME.table.network1 = signedKME(network1.fpkm.mat, network1$MEs)
  kME.table.network2 = signedKME(network2.fpkm.mat, network2$MEs)

  #assign names to the vectors of colors in order to be able to identify the genes (ensembl IDs) that intersect in each module overlap
  names(network1.colors) = colnames(network1.fpkm.mat)
  names(network2.colors) = colnames(network2.fpkm.mat)

  #iterate over the net1 modules (rows)
  for(net1.module in rownames(overlap.table$countTable)){
    if (net1.module != "grey"){
      #iterate over the net2 modules (columns)
      for(net2.module in colnames(overlap.table$countTable)){
        if (net2.module != "grey"){
          #calculate if the p-value of the overlap is lower than pval
          if( as.numeric(overlap.table$pTable[net1.module, net2.module]) < 0.05 ) {

            #calculate the percentage of overlap
            ov <- overlap.table$countTable[net1.module, net2.module]/min(as.numeric(table(network1.colors)[net1.module]), as.numeric(table(network2.colors)[net2.module]))
            #intersect the genes in both modules in order to get the common gene IDs
            overlaped.genes = intersect( names( which(network1.colors == net1.module) ), names( which(network2.colors == net2.module) ) )
            #obtain the kMEs (correlation of the gene profiles with the module eigengene) of the intersected genes in net1 and net2 datasets
            #selecting the correct module eigengene, and correlate them
            cor.kME = cor(
              kME.table.network1[overlaped.genes , paste("kME", match(net1.module, color.code)-1, sep="")],
              kME.table.network2[overlaped.genes , paste("kME", match(net2.module, color.code)-1, sep="")] )

            if( ov >= 0.2 & ov < 0.5 ){
                cor.table.1[net1.module, net2.module] = 2
                cor.table.2[net1.module, net2.module] = cor.kME
              } else if( ov >= 0.5 & ov < 0.7 ){
                cor.table.1[net1.module, net2.module] = 3
                cor.table.2[net1.module, net2.module] = cor.kME
              } else if( ov >= 0.7 ){
                cor.table.1[net1.module, net2.module] = 4
                cor.table.2[net1.module, net2.module] = cor.kME
              }
            }
          }
        }
      }
    }

  rows = which(rowSums(cor.table.1) == 0)
  cols = which(colSums(cor.table.1) == 0)
  cor.table.1[rows, ] <- 1
  cor.table.1[, cols] <- 1

  data.list <- list(cor.table.1, cor.table.2)
  names(data.list) <- c("corTable1", "corTable2")

  return(data.list)

}
