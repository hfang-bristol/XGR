#' Function to define genes crosslinking to an input list of genomic regions and the customised crosslink info
#'
#' \code{xGR2xGenes} is supposed to define genes crosslinking to an input list of genomic regions (GR). Also required is the crosslink info with a score crosslinking a GR to a gene. 
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink.customised the link info with a score linking a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if nor provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param cdf.function a character specifying how to transform the input crosslink score. It can be one of 'original' (no such transformation), and 'empirical'  for looking at empirical Cumulative Distribution Function (cdf; as such it is converted into pvalue-like values [0,1])
#' @param scoring logical to indicate whether gene-level scoring will be further calculated. By default, it sets to false
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param scoring.rescale logical to indicate whether gene scores will be further rescaled into the [0,1] range. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' If scoring sets to false, a data frame with following columns:
#' \itemize{
#'  \item{\code{GR}: genomic regions}
#'  \item{\code{Gene}: crosslinked genes}
#'  \item{\code{Score}: the score between the gene and the GR}
#'  \item{\code{Context}: the context}
#' }
#' If scoring sets to true, a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: crosslinked genes}
#'  \item{\code{Score}: gene score summarised over its list of crosslinked GR}
#'  \item{\code{Context}: the context}
#' }
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xGR}}
#' @include xGR2xGenes.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#'
#' # a) provide the genomic regions
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' names(gr) <- NULL
#' dGR <- xGR(gr, format="GRanges")
#'
#' # b) provide crosslink.customised
#' ## illustration purpose only (see the content of 'crosslink.customised')
#' df <- xGR2nGenes(dGR, format="GRanges", RData.location=RData.location)
#' crosslink.customised <- data.frame(GR=df$GR, Gene=df$Gene, Score=df$Weight, Context=rep('C',nrow(df)), stringsAsFactors=F)
#' #crosslink.customised <- data.frame(GR=df$GR, Gene=df$Gene, Score=df$Weight, stringsAsFactors=F)
#' 
#' # c) define crosslinking genes
#' # without gene scoring
#' df_xGenes <- xGR2xGenes(dGR, format="GRanges", crosslink.customised=crosslink.customised, RData.location=RData.location)
#' # with their scores
#' df_xGenes <- xGR2xGenes(dGR, format="GRanges", crosslink.customised=crosslink.customised, scoring=T, scoring.scheme="max", RData.location=RData.location)
#' }

xGR2xGenes <- function(data, format=c("chr:start-end","data.frame","bed","GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring=F, scoring.scheme=c("max","sum","sequential"), scoring.rescale=F, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
	
	dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
	
	###########################	
	# customised df_SGS
	###########################
	df_SGS_customised <- NULL
	if(!is.null(crosslink.customised)){

		if(is.vector(crosslink.customised)){
			# assume a file
			df <- utils::read.delim(file=crosslink.customised, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
		}else if(is.matrix(crosslink.customised) | is.data.frame(crosslink.customised)){
			df <- crosslink.customised
		}
		
		if(!is.null(df) & (ncol(df)==4 | ncol(df)==3)){
		
			if(ncol(df)==4){
				SGS_customised <- df
			}else{
				SGS_customised <- df
				SGS_customised$Context <- 'C'
			}
			colnames(SGS_customised) <- c("GR", "Gene", "Score", "Context")

			############################
			# remove Gene if NA
			# remove GR if NA
			df_SGS_customised <- SGS_customised[!is.na(SGS_customised[,1]) & !is.na(SGS_customised[,2]),]
			############################
			
			if(verbose){
				message(sprintf("%d genes are customised for contexts (n=%d) (%s) ...", length(unique(df_SGS_customised[,2])), length(unique(df_SGS_customised[,4])), as.character(Sys.time())), appendLF=TRUE)
			}
		}
	}
	
	#########################################
	if(!is.null(df_SGS_customised)){
		############################
		# remove Gene if ''
		# remove GR if ''
		df_SGS_customised <- df_SGS_customised[df_SGS_customised[,1]!='' & df_SGS_customised[,2]!='',]
		############################
	}
	
	if(is.null(df_SGS_customised)){
		res_df <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No xGenes are defined"), appendLF=TRUE)
		}
		
	}else{
		
		Gene <- Weight <- Score <- NULL
		
		ls_df_SGS <- split(x=df_SGS_customised, f=df_SGS_customised$Context)
		ls_res_df <- lapply(1:length(ls_df_SGS), function(j){
			df_SGS <- ls_df_SGS[[j]]
			
			if(cdf.function=="empirical"){
				## Compute an empirical cumulative distribution function
				my.CDF <- stats::ecdf(df_SGS$Score)
				df_SGS$Weight <- my.CDF(df_SGS$Score)
				
			}else{
				df_SGS$Weight <- df_SGS$Score
			}
			
			gr <- xGR(df_SGS$GR, format="chr:start-end", verbose=verbose, RData.location=RData.location)
			
			q2r <- as.data.frame(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=0, minoverlap=1L, type="any", select="all", ignore.strand=TRUE))
			q2r$gr <- names(gr[q2r[,2]])
			q2r$dgr <- names(dGR[q2r[,1]])
			
			ls_dgr <- split(x=q2r$gr, f=q2r$dgr)
			ls_df <- lapply(1:length(ls_dgr), function(i){
				ind <- match(ls_dgr[[i]], df_SGS$GR)
				df <- df_SGS[ind, ]
				
				#################################
				## keep maximum weight if there are many overlaps
				#################################
				df <- as.data.frame(df %>% dplyr::group_by(Gene) %>% dplyr::summarize(Score=max(Weight)))
				data.frame(GR=rep(names(ls_dgr)[i],nrow(df)), df, stringsAsFactors=F)
			})
			df_xGenes <- do.call(rbind, ls_df)
	
			if(verbose){
				message(sprintf("\t%d xGenes are defined for the context '%s' (%s)", length(unique(df_xGenes$Gene)), names(ls_df_SGS)[j], as.character(Sys.time())), appendLF=T)
			}
		
			############################################
			## whether gene scoring
			if(scoring){
			
				## calculate genetic influence score under a set of SNPs for each seed gene
				if(scoring.scheme=="max"){
					summaryFun <- max
				}else if(scoring.scheme=="sum"){
					summaryFun <- sum
				}else if(scoring.scheme=="sequential"){
					summaryFun <- function(x){
						base::sum(x / base::rank(-x,ties.method="min"))
					}
				}
				
				df_xGenes <- as.data.frame(df_xGenes %>% dplyr::group_by(Gene) %>% dplyr::summarise(Score=summaryFun(Score)))

				if(verbose){
					now <- Sys.time()
					message(sprintf("\t%d xGenes are scored using '%s' scoring scheme (%s)", length(unique(df_xGenes$Gene)), scoring.scheme, as.character(now)), appendLF=T)
				}

				if(scoring.rescale){
					if(verbose){
						now <- Sys.time()
						message(sprintf("\talso rescale score into the [0,1] range (%s)", as.character(now)), appendLF=T)
					}
					# rescale to [0 1]
					rescaleFun <- function(x){
						(x - min(x))/(max(x) - min(x))
					}
					
					df_xGenes$Score <- rescaleFun(df_xGenes$Score)
				}
		
			}
	
			data.frame(df_xGenes, Context=rep(names(ls_df_SGS)[j],nrow(df_xGenes)), stringsAsFactors=F)
		})
		res_df <- do.call(rbind, ls_res_df)
	}
	
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
  	
	invisible(res_df)
}
