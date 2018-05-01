#' Function to extract overlap-based scores given a list of genomic regions
#'
#' \code{xGRoverlap} is supposed to extract overlap-based scores given a list of genomic regions. Scores are extracted for overlapped sub-regions only, valued at the mean per base; otherwise NA. It returns a GR object.
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param GR.score the genomic regions together with score data. By default, it is 'NA' to disable this option. Pre-built genomic score data: 'RecombinationRate' (recombintion rate, \url{http://www.ncbi.nlm.nih.gov/pubmed/17943122})), 'phastCons100way', 'phyloP100way', 'GERP'. Beyond pre-built data, the user can specify the customised input: load your customised GR object directly (with the first meta column for scores; if not provided, it will be valued at 1)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a GenomicRanges object, appended with a meta-column 'GScore'. If input data contains only a genomic region, then outputs are all overlapped regions from GR.score; otherwise all overlapped regions from input data will be output.
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xGRoverlap.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#'
#' # a) provide the genomic regions
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS
#' data <- ImmunoBase$AS$variant
#'
#' # b) extract recombination rate
#' gr <- xGRoverlap(data=data, format="GRanges", GR.score="RecombinationRate", RData.location=RData.location)
#'
#' ############################################
#' # gene-centric genomic score (per base)
#' gr_Gene <- xRDataLoader('UCSC_knownGene', RData.location=RData.location)
#' ## recombination rate
#' gr_rr <- xGRoverlap(data=gr_Gene, format="GRanges", GR.score="RecombinationRate", RData.location=RData.location)
#' ## phastCons100way
#' gr_phast <- xGRoverlap(data=gr_Gene, format="GRanges", GR.score="phastCons100way", RData.location=RData.location)
#' ## phyloP100way
#' gr_phylo <- xGRoverlap(data=gr_Gene, format="GRanges", GR.score="phyloP100way", RData.location=RData.location)
#' }

xGRoverlap <- function(data, format=c("chr:start-end","data.frame","bed","GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), GR.score=c(NA,"RecombinationRate","phastCons100way","phyloP100way","GERP"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    #GR.score <- match.arg(GR.score)
	
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
	
	if(class(GR.score) == "GRanges"){
		if(verbose){
			now <- Sys.time()
			message(sprintf("Second, load the customised genomic annotation (%s) ...",  as.character(now)), appendLF=T)
		}
		
		###################################
		## now GR.score can be directly provided as a GR object
		###################################
		qGR <- GR.score
	}else{
		if(is.na(GR.score)){
			warning("Please specify GR.score")
			return(NULL)
		}else{
		
			if(verbose){
				now <- Sys.time()
				message(sprintf("Second, load the genomic annotation '%s' (%s) ...", GR.score[1], as.character(now)), appendLF=T)
			}
		
			qGR <- xRDataLoader(GR.score, verbose=verbose, RData.location=RData.location)
		}
	}
	
	if(!is.null(GenomicRanges::mcols(qGR))){
		## only the first column
		qGR$Value <- GenomicRanges::mcols(qGR)[,1]
	}else{
		## otherwise 1
		qGR$Value <- 1
	}
	
	oGR <- mergeOverlaps(qGR=qGR, sGR=dGR, maxgap=-1L, minoverlap=0L)

    if(length(dGR)==1){
    	## if dGR is 1 in length, it will return overlapped part from GR.score
		oGR$GScore <- oGR$Value
		oGR$hit <- NULL
		return(oGR)
		
	}else{
		## otherwise, it will return overlapped part from dGR
		
		##########################################
		if(verbose){
			now <- Sys.time()
			message(sprintf("Last, calculate the '%s' scores for %d genomic regions (%s) ...", GR.score, length(dGR), as.character(now)), appendLF=T)
		}
	
		hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query=dGR, subject=oGR, maxgap=-1L, minoverlap=0L, type="any", select="all", ignore.strand=T)))
		ls_vec <- split(x=hits[,2], f=hits[,1])
		width_oGR <- as.numeric(IRanges::width(oGR))
		value_oGR <- oGR$Value
		ls_res <- lapply(ls_vec, function(x){
			#y <- oGR[x]
			#z <- as.numeric(IRanges::width(y))
			z <- width_oGR[x]
			sum(value_oGR[x] * z) / sum(z)
		})
		vec_res <- unlist(ls_res)
		dGR$GScore <- NA
		dGR$GScore[as.numeric(names(vec_res))] <- vec_res
	}
	
	return(dGR)
}
