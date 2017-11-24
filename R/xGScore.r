#' Function to extract genomic scores given a list of genomic regions
#'
#' \code{xGScore} is supposed to extract genomic scores given a list of genomic regions. Scores for genomic regions/variants can be constraint/conservation or impact/pathogenicity. It returns a GR object.
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param GS.annotation which genomic scores (GS) annotaions used. It can be 'fitCons' (the probability of fitness consequences for point mutations; \url{http://www.ncbi.nlm.nih.gov/pubmed/25599402}), 'phastCons' (the probability that each nucleotide belongs to a conserved element/negative selection [0,1]: \url{http://www.ncbi.nlm.nih.gov/pubmed/25599402}), 'phyloP' (conservation at individual sites representing -log p-values under a null hypothesis of neutral evolution, positive scores for conservation and negative scores for acceleration), 'mcap' (eliminating a majority of variants with uncertain significance in clinical exomes at high sensitivity: \url{http://www.ncbi.nlm.nih.gov/pubmed/27776117}), and 'cadd' (combined annotation dependent depletion for estimating relative levels of pathogenicity of potential human variants: \url{http://www.ncbi.nlm.nih.gov/pubmed/24487276})
#' @param scoring.scheme the method used to calculate scores spanning a set of GR. It can be one of "mean", "median", "max", "min" and "sum"
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a GenomicRanges object 
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xGScore.r
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
#' # b) extract fitness consequence score
#' gr <- xGScore(data=data, format="GRanges", GS.annotation="fitCons", scoring.scheme="mean", RData.location=RData.location)
#' }

xGScore <- function(data, format=c("chr:start-end","data.frame","bed","GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), GS.annotation=c("fitCons","phastCons","phyloP","mcap","cadd"), scoring.scheme=c("mean","median","max","min","sum"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    GS.annotation <- match.arg(GS.annotation)
    scoring.scheme <- match.arg(scoring.scheme)
	
	##########################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("First, prepare a GR object from the input file formatted as '%s' (%s) ...", format, as.character(now)), appendLF=T)
	}
	dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
	
	##########################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("Second, load the genomic scores '%s' (%s) ...", GS.annotation, as.character(now)), appendLF=T)
	}
	if(GS.annotation=='fitCons'){
		gsco <- xRDataLoader('hg19.fitCons.UCSC', RData.location=RData.location)
	}else if(GS.annotation=='phastCons'){
		gsco <- xRDataLoader('hg19.phastCons100way.UCSC', RData.location=RData.location)
	}else if(GS.annotation=='phyloP'){
		gsco <- xRDataLoader('hg19.phyloP100way.UCSC', RData.location=RData.location)
	}else if(GS.annotation=='mcap'){
		gsco <- xRDataLoader('mcap.v1.0.hg19', RData.location=RData.location)
	}else if(GS.annotation=='cadd'){
		gsco <- xRDataLoader('cadd.v1.3.hg19', RData.location=RData.location)
	}
	
	#############
	## replace '/var/www/bigdata_dev' with RData.location
	gsco@data_dirpath <- gsub('/var/www/bigdata_dev', RData.location, gsco@data_dirpath)
	#############
	
	if(verbose){
		message(sprintf("\tIntended directory: '%s' (%s) ...", gsco@data_dirpath, as.character(Sys.time())), appendLF=T)
	}
	
	## check file exists
	url <- gsco@data_dirpath
	if(!file.exists(url)){
		getLinks <- function() {
			links <- character(0)
			list(a = function(node, ...) {
					   links <<- c(links, XML::xmlGetAttr(node, "href"))
					   node
					 },
				 links = function() links)
		}
		h1 <- getLinks()
		XML::htmlTreeParse(url, handlers=h1)
		res <- h1$links()
		vec_files <- res[!(res %in% c("?C=N;O=D", "?C=M;O=A", "?C=S;O=A", "?C=D;O=A", "/bigdata_dev/", "/bigdata/"))]
		
		## create a new directory to hold the downloads
		my_dir <- file.path(getwd(), GS.annotation)
		if(verbose){
			message(sprintf("\tActual directory: '%s' (%s) ...", my_dir, as.character(Sys.time())), appendLF=T)
		}
		if(!file.exists(my_dir)){
		
			if(verbose){
				message(sprintf("\tDownloading files (once and only once) into '%s' (%s) ...", my_dir, as.character(Sys.time())), appendLF=T)
			}
		
			dir.create(my_dir)
			#my_files <- file.path(my_dir, vec_files)
			#file.create(my_files)
			## download all files
			ls_tmp <- lapply(file.path(url,vec_files), function(x){
				source <- x
				target <- file.path(my_dir, basename(x))
				utils::download.file(source, target, quiet=T)
			})
		}

		gsco@data_dirpath <- my_dir
	}
    
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
	
	dGR$GScore <- suppressWarnings(GenomicScores::scores(gsco, dGR, scores.only=TRUE, summaryFun=summaryFun))

	return(dGR)
}
