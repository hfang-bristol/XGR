#' Function to perform eQTL analysis using MatrixEQTL
#'
#' \code{xME} is supposed to perform eQTL analysis using MatrixEQTL.
#
#' @param expression a data matrix/frame with genes in rows and samples in columns
#' @param genotype a data (sparse) matrix/frame with SNPs in rows and samples in columns
#' @param GR.Gene an GR object storing the genomic regions of genes
#' @param GR.SNP an GR object storing the genomic regions of SNPs
#' @param mode a character specifying the eQTL mode. It can be 'cis' for cis-eQTL analysis, 'both' for both cis-eQTL and trans-eQTL analysis
#' @param distance.cis a distance defining cis-eQTL. By default, it is 1Mb
#' @param pvalue.cis the p-value threshold outputing the cis-eQTL. By default, it is 1e-2
#' @param pvalue.trans the p-value threshold outputing the trans-eQTL. By default, it is 1e-5
#' @param which.PCs a vector specifying which principle components (PCs) are used for expression being regressed out. If NULL (by default), no gression is done
#' @param mydir a temporary directory file
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @return 
#' a "MatrixEQTL" object, a list including componets 'cis' and 'trans' (if mode is 'both'). The 'cis' (or 'trans') contains a data frame 'eqtls' with following columns ("snps","gene","statistic","pvalue","FDR","beta")
#' @note none
#' @export
#' @seealso \code{\link{xRegress}}
#' @include xME.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }

xME <- function(expression, genotype, GR.Gene, GR.SNP, mode=c('cis','both'), distance.cis=1e6, pvalue.cis=1e-2, pvalue.trans=1e-5, which.PCs=NULL, mydir=tempdir(check=T), verbose=T, silent=F)
{
    startT <- Sys.time()
    if(!silent){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }else{
    	verbose <- FALSE
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    mode <- match.arg(mode)
	
	if(any(which.PCs==0)){
		which.PCs <- NULL
	}
	if(!is.null(which.PCs)){
		if(verbose){
			message(sprintf("0. Expression matrix regressed out from '%s' PCs (%s) ...", paste(which.PCs,collapse=","), as.character(Sys.time())), appendLF=TRUE)
		}
		res <- xRegress(expression, center=TRUE, scale=TRUE, which.PCs=which.PCs)
		data_expression <- res$regressed
	}else{
		data_expression <- expression
	}
	#### only genes in common ('data_expression' and 'GR.Gene')
	ind <- match(rownames(data_expression), names(GR.Gene))
	data_expression <- data_expression[!is.na(ind),]
	GR.Gene <- GR.Gene[ind[!is.na(ind)]]
	
	if(class(genotype) %in% c('dgCMatrix','data.frame')){
		data_genotype <- as.matrix(genotype)
	}else{
		data_genotype <- genotype
	}
	#### only snps in common ('data_genotype' and 'GR.SNP')
	ind <- match(rownames(data_genotype), names(GR.SNP))
	data_genotype <- data_genotype[!is.na(ind),]
	GR.SNP <- GR.SNP[ind[!is.na(ind)]]
	
	################
	### only samples in common ('data_expression' and 'data_genotype')
	ind <- match(colnames(data_expression), colnames(data_genotype))
	data_expression <- data_expression[,!is.na(ind)]
	data_genotype <- data_genotype[,ind[!is.na(ind)]]
	################
	
	if(is.null(mydir)){
		mydir <- file.path(tempdir(check=T),)
	}
	# random numbers derived from time stamps
	tempnum <- as.character(as.numeric(Sys.time()))
	
	###########
	# gene file
	###########
	genefile <- file.path(mydir, paste0(tempnum,'.expression.txt'))
	df_output <- data.frame(id=rownames(data_expression), data_expression, stringsAsFactors=F, check.names=F)
	utils::write.table(df_output, file=genefile, row.names=F, col.names=T, quote=F, sep="\t")
    if(verbose){
        message(sprintf("1. Expression matrix of %d rows/genes X %d columns saved into '%s' (%s) ...", nrow(data_expression), ncol(data_expression), genefile, as.character(Sys.time())), appendLF=TRUE)
    }
    
	###########
	# snp file
	###########
	snpfile <- file.path(mydir, paste0(tempnum,'.genotype.txt'))
	df_output <- data.frame(id=rownames(data_genotype), data_genotype, stringsAsFactors=F, check.names=F)
	utils::write.table(df_output, file=snpfile, row.names=F, col.names=T, quote=F, sep="\t")
    if(verbose){
        message(sprintf("2. Genotype matrix of %d rows/SNPs X %d columns saved into '%s' (%s) ...", nrow(data_genotype), ncol(data_genotype), snpfile, as.character(Sys.time())), appendLF=TRUE)
    }

	###########
	# snpspos and genepos
	###########
	df <- GenomicRanges::as.data.frame(GR.Gene)
	genepos <- data.frame(gene=rownames(df), df[,c(1:3)], stringsAsFactors=F)
	df <- GenomicRanges::as.data.frame(GR.SNP)
	snpspos <- data.frame(snp=rownames(df), df[,c(1:2)], stringsAsFactors=F)
	
	###########################################
	# library(MatrixEQTL)
	###########################################
	# snps
	snps <- MatrixEQTL::SlicedData$new()
	snps$fileDelimiter <- "\t"
	snps$fileOmitCharacters <- "NA"
	snps$fileSkipRows <- 1
	snps$fileSkipColumns <- 1
	snps$fileSliceSize <- 2000
	suppressMessages(snps$LoadFile(snpfile))
	# gene
	gene <- MatrixEQTL::SlicedData$new()
	gene$fileDelimiter <- "\t"
	gene$fileOmitCharacters <- "NA"
	gene$fileSkipRows <- 1
	gene$fileSkipColumns <- 1
	gene$fileSliceSize <- 2000
	suppressMessages(gene$LoadFile(genefile))

	if(mode=='cis'){
		if(verbose){
			message(sprintf("3. cis-eQTL analysis: within %d distance window and pvalue < %.1e reported (%s) ...", distance.cis, pvalue.cis, as.character(Sys.time())), appendLF=TRUE)
		}
		if(silent){
			suppressMessages(me <- MatrixEQTL::Matrix_eQTL_main(snps=snps, gene=gene, useModel=MatrixEQTL::modelLINEAR, pvOutputThreshold=0, output_file_name.cis=NULL, pvOutputThreshold.cis=pvalue.cis, snpspos=snpspos, genepos=genepos, cisDist=distance.cis, min.pv.by.genesnp=F, noFDRsaveMemory=F, verbose=F))
		}else{
			me <- MatrixEQTL::Matrix_eQTL_main(snps=snps, gene=gene, useModel=MatrixEQTL::modelLINEAR, pvOutputThreshold=0, output_file_name.cis=NULL, pvOutputThreshold.cis=pvalue.cis, snpspos=snpspos, genepos=genepos, cisDist=distance.cis, min.pv.by.genesnp=F, noFDRsaveMemory=F, verbose=F)
		}
		
	}else{
		if(verbose){
			message(sprintf("3. eQTL analysis: cis (within %d distance window and pvalue < %.1e reported) and trans (pvalue < %.1e reported) (%s) ...", distance.cis, pvalue.cis, pvalue.trans, as.character(Sys.time())), appendLF=TRUE)
		}
		if(silent){
			suppressMessages(me <- MatrixEQTL::Matrix_eQTL_main(snps=snps, gene=gene, useModel=MatrixEQTL::modelLINEAR, output_file_name=NULL, pvOutputThreshold=pvalue.trans, output_file_name.cis=NULL, pvOutputThreshold.cis=pvalue.cis, snpspos=snpspos, genepos=genepos, cisDist=distance.cis, min.pv.by.genesnp=F, noFDRsaveMemory=F, verbose=F))
		}else{
			me <- MatrixEQTL::Matrix_eQTL_main(snps=snps, gene=gene, useModel=MatrixEQTL::modelLINEAR, output_file_name=NULL, pvOutputThreshold=pvalue.trans, output_file_name.cis=NULL, pvOutputThreshold.cis=pvalue.cis, snpspos=snpspos, genepos=genepos, cisDist=distance.cis, min.pv.by.genesnp=F, noFDRsaveMemory=F, verbose=F)
		}
	}

    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xEnricherGenes): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }

    invisible(me)
}
