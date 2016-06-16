#' Function to score lead or LD SNPs based on the given significance level
#'
#' \code{xSNPscores} is supposed to score a list of Lead SNPs together with the significance level. It can consider LD SNPs and the given threshold of the significant level.
#'
#' @param data a named input vector containing the sinificance level for nodes (dbSNP). For this named vector, the element names are dbSNP (starting with rs or in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for example, 'chr16:28525386'), the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for dbSNP, 2nd column for the significance level. 
#' @param include.LD additional SNPs in LD with Lead SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, LD SNPs will be included based on one or more of 26 populations and 5 super populations from 1000 Genomics Project data (phase 3). The population can be one of 5 super populations ("AFR", "AMR", "EAS", "EUR", "SAS"), or one of 26 populations ("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"). Explanations for population code can be found at \url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
#' @param LD.customised a user-input matrix or data frame with 3 columns: 1st column for Lead SNPs, 2nd column for LD SNPs, and 3rd for LD r2 value. It is designed to allow the user analysing their precalcuated LD info. This customisation (if provided) has the high priority over built-in LD SNPs
#' @param LD.r2 the LD r2 value. By default, it is 0.8, meaning that SNPs in LD (r2>=0.8) with input SNPs will be considered as LD SNPs. It can be any value from 0.8 to 1
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of SNPs into scores. If given, those SNPs below this are considered significant and thus scored positively. Instead, those above this are considered insigificant and thus receive no score
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a data frame with following columns:
#' \itemize{
#'  \item{\code{SNP}: Lead and/or LD SNPs}
#'  \item{\code{Score}: the scores for SNPs calculated based on p-values taking into account the given threshold of the significant level}
#'  \item{\code{Pval}: the input p-values for Lead SNPs or R2-adjusted p-values for LD SNPs}
#' }
#' @note None
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xSNPscores.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location="~/Sites/SVN/github/RDataCentre/Portal"
#'
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' data <- GenomicRanges::mcols(gr)[,c(1,3)]
#'
#' # b) calculate SNP scores (considering significant cutoff 5e-5)
#' ## without inclusion of LD SNPs
#' df_SNP <- xSNPscores(data=data, significance.threshold=5e-5, RData.location=RData.location)
#' ## include LD SNPs (calculated based on European populations)
#' df_SNP <- xSNPscores(data=data, significance.threshold=5e-5, include.LD="EUR", RData.location=RData.location)
#' }

xSNPscores <- function(data, include.LD=NA, LD.customised=NULL, LD.r2=0.8, significance.threshold=5e-5, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/Portal")
{

    if(is.null(data)){
        stop("The input data must be not NULL.\n")
    }else{
    	
    	if(class(data)=='DataFrame'){
    		data <- S4Vectors::as.matrix(data)
    	}
    
		if (is.vector(data)){
			if(length(data)>1){
				# assume a vector
				if(is.null(names(data))){
					stop("The input data must have names with attached dbSNP ID.\n")
				}
			}else{
				# assume a file
				data <- utils::read.delim(file=data, header=F, row.names=NULL, stringsAsFactors=F)
			}
		}
		
		if (is.vector(data)){
			pval <- data
		}else if(is.matrix(data) | is.data.frame(data)){
			data <- as.matrix(data)
			data_list <- split(x=data[,2], f=as.character(data[,1]))
			res_list <- lapply(data_list, function(x){
				x <- as.numeric(x)
				x <- x[!is.na(x)]
				if(length(x)>0){
					min(x)
				}else{
					NULL
				}
			})
			pval <- unlist(res_list)
		}
		
		# force those zeros to be miminum of non-zeros
		#tmp <- as.numeric(format(.Machine)['double.xmin'])
		tmp <- min(pval[pval!=0])
		pval[pval < tmp] <- tmp
	}
	
	## replace '_' with ':'
	tmp <- names(pval)
	tmp <- gsub("_", ":", tmp, perl=T)
	## replace 'imm:' with 'chr'
	names(pval) <- gsub("imm:", "chr", tmp, perl=T)
	
	Lead_Sig <- data.frame(SNP=names(pval), Sig=pval, row.names=NULL, stringsAsFactors=F)
	leads <- Lead_Sig[,1]
	sigs <- Lead_Sig[,2]

	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d Lead SNPs are input", length(leads)), appendLF=T)
	}

	###########################
	## include additional SNPs that are in LD with input SNPs
	if(LD.r2>=0.8 & LD.r2<=1){
		default.include.LD <- c("ACB","AFR","AMR","ASW","BEB","CDX","CEU","CHB","CHS","CLM","EAS","ESN","EUR","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","SAS","STU","TSI","YRI")
		ind <- match(default.include.LD, include.LD)
		include.LD <- default.include.LD[!is.na(ind)]
	}else{
		include.LD <- NULL
	}
	
	LLR <- NULL
	if(length(include.LD) > 0 & is.null(LD.customised)){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Inclusion of LD SNPs is based on population (%s) with R2 >= %f", paste(include.LD, collapse=','), LD.r2), appendLF=T)
		}
	
		GWAS_LD <- xRDataLoader(RData.customised='GWAS_LD', RData.location=RData.location, verbose=verbose)
		res_list <- lapply(include.LD, function(x){
			data_ld <- ''
			eval(parse(text=paste("data_ld <- GWAS_LD$", x, sep="")))
			ind <- match(rownames(data_ld), leads)
			ind_lead <- which(!is.na(ind))
			
			if(length(ind_lead) > 2){
				ind_ld <- which(Matrix::colSums(data_ld[ind_lead,]>=LD.r2)>0)
				sLL <- data_ld[ind_lead, ind_ld]
				summ <- summary(sLL)
				res <- data.frame(Lead=rownames(sLL)[summ$i], LD=colnames(sLL)[summ$j], R2=summ$x, stringsAsFactors=F)
			}else if(length(ind_lead) == 1){
				ind_ld <- which(data_ld[ind_lead,]>=LD.r2)
				sLL <- data_ld[ind_lead, ind_ld]
				res <- data.frame(Lead=rep(rownames(data_ld)[ind_lead],length(sLL)), LD=names(sLL), R2=sLL, stringsAsFactors=F)
			}else{
				NULL
			}
		})
		## get data frame (Lead LD R2)
		LLR <- do.call(rbind, res_list)
		
		###########################
		## also based on ImmunoBase
		if(1){
			ImmunoBase_LD <- xRDataLoader(RData.customised='ImmunoBase_LD', RData.location=RData.location, verbose=verbose)
			res_list <- lapply(include.LD, function(x){
				data_ld <- ''
				eval(parse(text=paste("data_ld <- ImmunoBase_LD$", x, sep="")))
				ind <- match(rownames(data_ld), leads)
				ind_lead <- which(!is.na(ind))
				
				if(length(ind_lead) > 2){
					ind_ld <- which(Matrix::colSums(data_ld[ind_lead,]>=LD.r2)>0)
					sLL <- data_ld[ind_lead, ind_ld]
					summ <- summary(sLL)
					res <- data.frame(Lead=rownames(sLL)[summ$i], LD=colnames(sLL)[summ$j], R2=summ$x, stringsAsFactors=F)
				}else if(length(ind_lead) == 1){
					ind_ld <- which(data_ld[ind_lead,]>=LD.r2)
					sLL <- data_ld[ind_lead, ind_ld]
					res <- data.frame(Lead=rep(rownames(data_ld)[ind_lead],length(sLL)), LD=names(sLL), R2=sLL, stringsAsFactors=F)
				}else{
					NULL
				}
				
			})
			## get data frame (Lead LD R2)
			LLR_tmp <- do.call(rbind, res_list)
			LLR <- rbind(LLR, LLR_tmp)
		}
		
	###########################	
	}else if(!is.null(LD.customised)){
		if (is.vector(LD.customised)){
			# assume a file
			LLR <- utils::read.delim(file=LD.customised, header=F, row.names=NULL, stringsAsFactors=F)
		}else if(is.matrix(LD.customised) | is.data.frame(LD.customised)){
			LLR <- LD.customised
		}
		
		if(!is.null(LLR)){
			flag <- LLR[,3]>=LD.r2
			if(sum(flag)>0){
				LLR <- LLR[LLR[,3]>=LD.r2,]
				colnames(LLR) <- c("Lead", "LD", "R2")
			
				if(verbose){
					now <- Sys.time()
					message(sprintf("Inclusion of LD SNPs is based on customised data (%d Lead SNPs and %d LD SNPs) with R2>=%f", length(unique(LLR[,1])), length(unique(LLR[,2])), LD.r2), appendLF=T)
				}
			}else{
				LLR <- NULL
			}
		}
		
	}
	
	if(!is.null(LLR)){
		## get data frame (LD Sig)
		ld_list <- split(x=LLR[,-2], f=LLR[,2])
		res_list <- lapply(ld_list, function(x){
			ind <- match(x$Lead, leads)
			## power transformation of p-values X R2, then keep the min (the most significant)
			min(sigs[ind] ^ x$R2)
		})
		vec <- unlist(res_list)
		LD_Sig <- data.frame(SNP=names(vec), Sig=vec, row.names=NULL, stringsAsFactors=F)

		## merge Lead and LD
		df <- rbind(Lead_Sig, as.matrix(LD_Sig))
		res_list <- split(x=df$Sig, f=df$SNP)
		res <- lapply(res_list, function(x){
			min(x)
		})
		vec <- unlist(res)
		SNP_Sig <- data.frame(SNP=names(vec), FDR=vec, row.names=NULL, stringsAsFactors=F)
	}else{
		if(verbose){
			now <- Sys.time()
			message(sprintf("Do not include any LD SNPs"), appendLF=T)
		}
	
		SNP_Sig <- Lead_Sig
	}
	###########################
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d Lead/LD SNPs are considered", nrow(SNP_Sig)), appendLF=T)
	}
	
	pval <- as.numeric(SNP_Sig[,2])
	names(pval) <- SNP_Sig[,1]
	
	# transformed into scores according to log-likelihood ratio between the true positives and the false positivies
    ## also take into account the given threshold of the significant level
    ## SNPs with p-value below this are considered significant and thus scored positively
    ## Instead, SNPs with p-values fdr above this are considered insigificant and thus scored negatively (zero-out)
	
	if(is.null(significance.threshold)){
        scores <- log10((1-pval)/pval)
        #scores <- log10(1/pval)
    }else{
		scores <- log10((1-pval)/pval) - log10((1-significance.threshold)/significance.threshold)
		#scores <- log10(1/pval) - log10(1/significance.threshold)
	}
    ## replace those infinite values with the next finite ones
    tmp_max <- max(scores[!is.infinite(scores)])
    tmp_min <- min(scores[!is.infinite(scores)])
    scores[scores>tmp_max] <- tmp_max
    scores[scores<tmp_min] <- tmp_min
	## zero-out SNPs with negative scores
	scores[scores<0] <- 0
	seeds.snps <- scores
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d Lead/LD SNPs are scored positively", sum(seeds.snps>0)), appendLF=T)
	}
    
    #########
    df_SNP <- data.frame(SNP=names(pval), Score=seeds.snps, Pval=pval, row.names=NULL, stringsAsFactors=F)
    #########
    
    invisible(df_SNP)
}
