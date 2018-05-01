#' Function to score genomic regions based on the given significance level
#'
#' \code{xGRscores} is supposed to score a list of genomic regions together with the significance level.
#'
#' @param data a named input vector containing the sinificance level for genomic regions (GR). For this named vector, the element names are GR, in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. The element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for GR, 2nd column for the significance level. 
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of GR into scores. If given, those GR below this are considered significant and thus scored positively. Instead, those above this are considered insigificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to 10. If NULL, no capping is applied
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a data frame with following columns:
#' \itemize{
#'  \item{\code{GR}: genomic regions}
#'  \item{\code{Score}: the scores for GR calculated based on p-values taking into account the given threshold of the significant level}
#'  \item{\code{Pval}: the input p-values for GR}
#' }
#' @note None
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xGRscores.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#'
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' df <- as.data.frame(gr, row.names=NULL)
#' chr <- df$seqnames
#' start <- df$start
#' end <- df$end
#' sig <- df$Pvalue
#' GR <- paste(chr,':',start,'-',end, sep='')
#' data <- cbind(GR=GR, Sig=sig)
#'
#' # b) calculate GR scores (considering significant cutoff 5e-5)
#' df_GR <- xGRscores(data=data, significance.threshold=5e-5, RData.location=RData.location)
#' }

xGRscores <- function(data, significance.threshold=5e-2, score.cap=10, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
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

	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d GR are input", length(pval)), appendLF=T)
	}
	
	# transformed into scores according to log-likelihood ratio between the true positives and the false positivies
    ## also take into account the given threshold of the significant level
    ## GR with p-value below this are considered significant and thus scored positively
    ## Instead, GR with p-values above this are considered insigificant and thus scored negatively (zero-out)
	
	if(is.null(significance.threshold)){
        scores <- log10((1-pval)/pval)
    }else{
		scores <- log10((1-pval)/pval) - log10((1-significance.threshold)/significance.threshold)
	}
    ## replace those infinite values with the next finite ones
    tmp_max <- max(scores[!is.infinite(scores)])
    tmp_min <- min(scores[!is.infinite(scores)])
    scores[scores>tmp_max] <- tmp_max
    scores[scores<tmp_min] <- tmp_min
	## zero-out SNPs with negative scores
	ind_remained <- which(scores>0)
	seeds.snps <- scores[ind_remained]
	pval <- pval[ind_remained]
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d GR are scored positively", sum(seeds.snps>0)), appendLF=T)
	}
    
    #########
    
    df_GR <- data.frame(GR=names(pval), Score=seeds.snps, Pval=pval, row.names=NULL, stringsAsFactors=F)
    
    
    ##############################
    ## cap the maximum score
    if(!is.null(score.cap)){
    	score.cap <- as.numeric(score.cap)
    	if(score.cap <= max(df_GR$Score)){
    		df_GR$Score[df_GR$Score>=score.cap] <- score.cap
    		
			if(verbose){
				now <- Sys.time()
				message(sprintf("GR score capped to the maximum score %d.", score.cap), appendLF=T)
			}
    	}
    }
    ##############################
    
    df_GR <- df_GR[order(df_GR$Score,-df_GR$Pval,df_GR$GR,decreasing=TRUE),]
    #########
    
    invisible(df_GR)
}
