#' Function to create a GRanges object given a list of genomic regions
#'
#' \code{xGR} is supposed to create a GRanges object given a list of genomic regions. 
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "chr:start-end", "data.frame", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param add.name logical to add names. By default, it sets to true
#' @param remove.mcol logical to remove meta-columns. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return a GenomicRanges object 
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xGR.r
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
#' df <- as.data.frame(gr, row.names=NULL)
#' chr <- df$seqnames
#' start <- df$start
#' end <- df$end
#' data <- paste(chr,':',start,'-',end, sep='')
#'
#' # b) create a GRanges object
#' GR <- xGR(data=data, format="chr:start-end", RData.location=RData.location)
#' }

xGR <- function(data, format=c("chr:start-end","data.frame","bed","GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), add.name=T, remove.mcol=F, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
	
	###################
	if(is.null(data)){
		return(NULL)
	}
	###################
		
    ## import data
    if(is.matrix(data) | is.data.frame(data) | class(data)=="GRanges"){
        data <- data
    }else if(!is.null(data) & any(!is.na(data))){
    	if(length(data)==1){
    		if(file.exists(data)){
    			data <- utils::read.delim(file=data, header=F, row.names=NULL, stringsAsFactors=F)
    			data <- unique(data[,1])
    		}else{
				data <- data
			}
		}else{
			data <- data
		}
    }else{
		warning("The file 'data' must be provided!\n")
		return(NULL)
    }
	
    ## construct GR
	if(format=="data.frame"){
		## construct data GR
		if(ncol(data)>=3){
			data <- data
		}else if(ncol(data)==2){
			data <- cbind(data, data[,2])
		}else{
			warning("Your input 'data.file' is not as expected!\n")
			return(NULL)
		}
		## make sure positions are numeric
		ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
		dGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(data[ind,1]),
			ranges = IRanges::IRanges(start=as.numeric(data[ind,2]), end=as.numeric(data[ind,3])),
			strand = S4Vectors::Rle(rep('*',length(ind)))
		)
		if(add.name){
			names(dGR) <- paste(data[ind,1], ':', data[ind,2], '-', data[ind,3], sep='')
		}
		
	}else if(format=="chr:start-end"){
		data <- unique(data[!is.na(data)])
		input <- do.call(rbind, strsplit(data, ":|-|,"))
		if(ncol(input)>=3){
			data <- matrix(input[,1:3], nrow=nrow(input))
		}else if(ncol(input)==2){
			data <- matrix(input[,c(1,2,2)], nrow=nrow(input))
		}else{
			warning("Your input 'data' does not meet the format 'chr:start-end'!\n")
			return(NULL)
		}
		## make sure positions are numeric
		ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
		dGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(data[ind,1]),
			ranges = IRanges::IRanges(start=as.numeric(data[ind,2]), end=as.numeric(data[ind,3])),
			strand = S4Vectors::Rle(rep('*',length(ind)))
		)
		if(add.name){
			names(dGR) <- paste(data[ind,1], ':', data[ind,2], '-', data[ind,3], sep='')
		}
		
	}else if(format=="bed"){
		## construct data GR
		## make sure positions are numeric
		ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
		dGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(data[ind,1]),
			ranges = IRanges::IRanges(start=as.numeric(data[ind,2])+1, end=as.numeric(data[ind,3])),
			strand = S4Vectors::Rle(rep('*',length(ind)))
		)
		if(add.name){
			names(dGR) <- paste(data[ind,1], ':', data[ind,2]+1, '-', data[ind,3], sep='')
		}
		
	}else if(format=="GRanges"){
		dGR <- data
		
		if(remove.mcol){
			GenomicRanges::mcols(dGR) <- NULL
		}
		
		if(is.null(names(dGR)) & add.name){
			df <- as.data.frame(dGR, row.names=NULL)
			names(dGR) <- paste(df$seqnames,':',df$start,'-',df$end, sep='')
		}
	}

	# lift over
	if(!is.na(build.conversion)){
		if(verbose){
			message(sprintf("\tdata genomic regions: lifted over via genome build conversion `%s`", build.conversion), appendLF=T)
		}
		dGR <- xLiftOver(data.file=dGR, format.file="GRanges", build.conversion=build.conversion, merged=F, verbose=verbose, RData.location=RData.location)
	}
  	#######################################################
  	
  	
  	invisible(dGR)

}
