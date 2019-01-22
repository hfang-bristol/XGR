#' Function to perform eQTL analysis using MatrixEQTL
#'
#' \code{xME} is supposed to perform eQTL analysis using MatrixEQTL.
#
#' @param expression a data matrix/frame with genes in rows and samples in columns
#' @param genotype a data (sparse) matrix/frame with SNPs in rows and samples in columns
#' @param GR.Gene an GR object storing the genomic regions of genes
#' @param GR.SNP an GR object storing the genomic regions of SNPs
#' @param mode a character specifying the eQTL mode. It can be 'cis' for cis-eQTL analysis, 'both' for both cis-eQTL and trans-eQTL analysis, and 'none' (no eQTL analysis only returning variables actually used for eQTL analysis: expression, genotype, GR.Gene and GR.SNP)
#' @param distance.cis a distance defining cis-eQTL. By default, it is 1Mb
#' @param pvalue.cis the p-value threshold outputing the cis-eQTL. By default, it is 1e-2
#' @param pvalue.trans the p-value threshold outputing the trans-eQTL. By default, it is 1e-5
#' @param num.PCs an integer specifying the number of principle components (PCs) used for expression being regressed out. If zero (by default), no gression is done
#' @param mydir a temporary directory file
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @return 
#' a "MatrixEQTL" object, a list including componets 'cis' and 'trans' (if mode is 'both'). The 'cis' (or 'trans') contains a data frame 'eqtls' with following columns ("snps","gene","statistic","pvalue","FDR","beta"). If the mode is 'none', return a list with four components 'expression', 'genotype', 'GR.Gene' and 'GR.SNP'
#' @note none
#' @export
#' @seealso \code{\link{xRegress}}
#' @include xME.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' res <- xME(expression, genotype, GR.Gene, GR.SNP, mode='none')
#' }

xME <- function(expression, genotype, GR.Gene, GR.SNP, mode=c('cis','both','none'), distance.cis=1e6, pvalue.cis=1e-2, pvalue.trans=1e-5, num.PCs=0, mydir=tempdir(check=T), verbose=T, silent=F)
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
	
	if(verbose){
		message(sprintf("0. Do checking (%s) ...", as.character(Sys.time())), appendLF=TRUE)
	}

	### only samples in common ('expression' and 'genotype')
	ind <- match(colnames(expression), colnames(genotype))
	data_expression <- expression[,!is.na(ind)]
	data_genotype <- as.matrix(genotype[,ind[!is.na(ind)]])
	
	### only snps in common ('data_genotype' and 'GR.SNP')
	ind <- match(rownames(data_genotype), names(GR.SNP))
	data_genotype <- data_genotype[!is.na(ind),]
	GR.SNP <- GR.SNP[ind[!is.na(ind)]]
	
	### only genes in common ('data_expression' and 'GR.Gene')
	ind <- match(rownames(data_expression), names(GR.Gene))
	data_expression <- data_expression[!is.na(ind),]
	GR.Gene <- GR.Gene[ind[!is.na(ind)]]
	
	### do regression or not
	num.PCs <- as.integer(num.PCs)
	if(num.PCs>0){
		which.PCs <- 1:num.PCs
		if(verbose){
			message(sprintf("\texpression matrix regressed out from '%d PCs (%s) ...", num.PCs, as.character(Sys.time())), appendLF=TRUE)
		}
		res <- xRegress(data_expression, center=TRUE, scale=TRUE, which.PCs=which.PCs)
		data_expression <- res$regressed
	}
	
	if(mode=='none'){
		me <- list(expression=data_expression, genotype=data_genotype, GR.Gene=GR.Gene, GR.SNP=GR.SNP)
		
		if(verbose){
			message(sprintf("1. Return a list with four components 'expression', 'genotype', 'GR.Gene' and 'GR.SNP' (%s) ...", as.character(Sys.time())), appendLF=TRUE)
			message(sprintf("\txpression matrix of %d rows/genes X %d columns", nrow(data_expression), ncol(data_expression)), appendLF=TRUE)
			message(sprintf("\tgenotype matrix of %d rows/SNPs X %d columns", nrow(data_genotype), ncol(data_genotype)), appendLF=TRUE)
		}
		
	}else{
	
		######################################################################################
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
	
		me <- NULL
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
	
		if(any(class(me) %in% c("MatrixEQTL"))){
			if(verbose){
				message(sprintf("Successfully finish the eQTL analysis (%s)!", as.character(Sys.time())), appendLF=TRUE)
			}
			unlink(genefile)
			unlink(snpfile)
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
