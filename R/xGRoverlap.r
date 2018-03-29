#' Function to extract overlap-based scores given a list of genomic regions
#'
#' \code{xGRoverlap} is supposed to extract overlap-based scores given a list of genomic regions. Scores are extracted for overlapped sub-regions only; otherwise NA. It returns a GR object.
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param GR.annotation the genomic regions of annotation data. By default, it is 'NA' to disable this option. Pre-built genomic annotation data: 'RecombinationRate' (recombintion rate, \url{http://www.ncbi.nlm.nih.gov/pubmed/17943122})). Beyond pre-built annotation data, the user can specify the customised input: load your customised GR object directly
#' @param scoring.scheme the method used to calculate scores spanning a set of GR. It can be one of "mean", "median", "max", "min" and "sum"
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a GenomicRanges object, appended with a meta-column 'GScore'
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xGRoverlap.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#'
#' # a) provide the genomic regions
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS
#' data <- ImmunoBase$AS$variant
#'
#' # b) extract recombination rate
#' gr <- xGRoverlap(data=data, format="GRanges", GR.annotation="RecombinationRate", scoring.scheme="mean", RData.location=RData.location)
#' }

xGRoverlap <- function(data, format=c("chr:start-end","data.frame","bed","GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), GR.annotation=c(NA, "RecombinationRate"), scoring.scheme=c("mean","median","max","min","sum"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata_dev")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    #GS.annotation <- match.arg(GS.annotation)
    scoring.scheme <- match.arg(scoring.scheme)
	
	##########################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("First, prepare a GR object from the input file formatted as '%s' (%s) ...", format, as.character(now)), appendLF=T)
	}
	dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
	
	#####################################
	## A function to return an GR object storing overlapped regions (ie only overlapped regions!)
	mergeOverlaps <- function(qGR, sGR, maxgap=-1L, minoverlap=0L){
		hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query=qGR, subject=sGR, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
		qhits <- qGR[hits[,1]]
		shits <- sGR[hits[,2]]

		oGR <- IRanges::pintersect(qhits, shits, ignore.strand=T)
		#IRanges::reduce(oGR)
	}
	
	##########################################
	
	if(class(GR.annotation) == "GRanges"){
		if(verbose){
			now <- Sys.time()
			message(sprintf("Second, load the customised genomic annotation (%s) ...",  as.character(now)), appendLF=T)
		}
			
		###################################
		## now GR.annotation can be directly provided as a GR object
		###################################
		qGR <- GR.annotation
	}else{
		if(is.na(GR.annotation)){
			stop("Please specify annotation RData!\n")
		}else{
		
			if(verbose){
				now <- Sys.time()
				message(sprintf("Second, load the genomic annotation '%s' (%s) ...", GS.annotation[1], as.character(now)), appendLF=T)
			}
		
			GR.annotation <- GR.annotation[1]
			qGR <- xRDataLoader(GS.annotation, verbose=verbose, RData.location=RData.location)
		}
	}
	
	if(GS.annotation=='RecombinationRate'){
		qGR$Value <- qGR$Rate
	}else if(!is.null(GenomicRanges::mcols(qGR))){
		## only the first column
		qGR$Value <- GenomicRanges::mcols(qGR)[,1]
	}else{
		## otherwise 1
		qGR$Value <- 1
	}
	
	oGR <- mergeOverlaps(qGR=qGR, sGR=dGR, maxgap=-1L, minoverlap=0L)
    
	##########################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("Last, calculate the '%s' scores for %d genomic regions (%s) ...", scoring.scheme, length(dGR), as.character(now)), appendLF=T)
	}
	if(scoring.scheme=="mean"){
		summaryFun <- mean
	}else if(scoring.scheme=="median"){
		summaryFun <- stats::median
	}else if(scoring.scheme=="max"){
		summaryFun <- max
	}else if(scoring.scheme=="min"){
		summaryFun <- min
	}else if(scoring.scheme=="sum"){
		summaryFun <- sum
	}
	
    hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query=dGR, subject=oGR, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
    ls_vec <- split(x=hits[,2], f=hits[,1])
    ls_res <- lapply(ls_vec, function(x){
    	do.call(summaryFun, list(x))
    })
    vec_res <- unlist(ls_res)
    dGR$GScore <- NA
    dGR$GScore[as.numeric(names(vec_res))] <- vec_res

	return(dGR)
}
