#' Function to lift genomic intervals from one genome build to another.
#'
#' \code{xLiftOver} is supposed to lift genomic intervals from one genome build to another. Supported are the conversion between genom builds 'hg38' (GRCh38), 'hg19' (GRCh37) and 'h18'. 
#'
#' @param data.file an input data file, containing a list of genomic regions to test. If the input file is formatted as a 'data.frame' (specified by the parameter 'format.file' below), the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. If the format is indicated as "chr:start-end", instead of using the first 3 columns, only the first column will be used and processed. If the file also contains other columns, these additional columns will be ignored. Alternatively, the input file can be the content itself assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param format.file the format for input files. It can be one of "data.frame", "chr:start-end", "bed"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19", "hg19.to.hg38", "hg19.to.hg18", "hg18.to.hg38" and "hg18.to.hg19". By default it is NA, forcing the user to specify the corrent one.
#' @param merged logical to indicate whether multiple ranges should be merged into the one per a range in query. By default, it sets to true
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' an GR oject storing converted genomic intervals.
#' @export
#' @seealso \code{\link{xLiftOver}}
#' @include xLiftOver.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location="~/Sites/SVN/github/bigdata"
#' 
#' # Provide UCSC genes (hg19)
#' UCSC_genes <- xRDataLoader(RData.customised='UCSC_genes', RData.location=RData.location)
#' UCSC_genes
#' 
#' # Lift over to hg38
#' gr <- xLiftOver(UCSC_genes, format.file="GRanges", build.conversion="hg19.to.hg38", RData.location=RData.location)
#' gr
#' }

xLiftOver <- function(data.file, format.file=c("data.frame", "bed", "chr:start-end", "GRanges"), build.conversion=c(NA, "hg38.to.hg19","hg19.to.hg38","hg19.to.hg18","hg18.to.hg38","hg18.to.hg19"), merged=T, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/Portal")
{

    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format.file <- match.arg(format.file)
    build.conversion <- match.arg(build.conversion)
    
    ###################
	if(verbose){
		now <- Sys.time()
		message(sprintf("First, import the files formatted as '%s' (%s) ...", format.file, as.character(now)), appendLF=T)
	}
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\timport the data file (%s) ...", as.character(now)), appendLF=T)
	}
    ## import data file
    if(is.matrix(data.file) | is.data.frame(data.file) | class(data.file)=="GRanges"){
        data <- data.file
    }else if(!is.null(data.file) & !is.na(data.file)){
		data <- utils::read.delim(file=data.file, header=F, row.names=NULL, stringsAsFactors=F)
    }else{
    	stop("The file 'data.file' must be provided!\n")
    }
    
    ###################
	if(verbose){
		now <- Sys.time()
		message(sprintf("Second, construct GenomicRanges object (%s) ...", as.character(now)), appendLF=T)
	}
    
	if(format.file=="data.frame"){
		## construct data GR
		if(ncol(data)>=3){
			data <- data
		}else if(ncol(data)==2){
			data <- cbind(data, data[,2])
		}else{
			stop("Your input 'data.file' is not as expected!\n")
		}
		## make sure positions are numeric
		ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
		data <- data[ind,]
		dGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(data[,1]),
			ranges = IRanges::IRanges(start=as.numeric(data[,2]), end=as.numeric(data[,3])),
			strand = S4Vectors::Rle(rep('*',nrow(data)))
		)
		
	}else if(format.file=="chr:start-end"){
		
		## construct data GR
		input <- do.call(rbind, strsplit(data[,1], ":|-"))
		if(ncol(input)>=3){
			data <- input
		}else if(ncol(input)==2){
			data <- cbind(input, input[,2])
		}else{
			stop("Your input 'data.file' does not meet the format 'chr:start-end'!\n")
		}
		## make sure positions are numeric
		ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
		data <- data[ind,]
		dGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(data[,1]),
			ranges = IRanges::IRanges(start=as.numeric(data[,2]+1), end=as.numeric(data[,3])),
			strand = S4Vectors::Rle(rep('*',nrow(data)))
		)
		
	}else if(format.file=="bed"){
		## construct data GR
		## make sure positions are numeric
		ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
		data <- data[ind,]
		dGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(data[,1]),
			ranges = IRanges::IRanges(start=as.numeric(data[,2])+1, end=as.numeric(data[,3])),
			strand = S4Vectors::Rle(rep('*',nrow(data)))
		)

	}else if(format.file=="GRanges"){
		## construct data GR
		dGR <- data
	}
	
	#####################################
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Third, lift intervals between genome builds '%s' (%s) ...", build.conversion, as.character(now)), appendLF=T)
	}
    
    chains <- xRDataLoader(RData.customised='chain', RData.location=RData.location, verbose=verbose)
	chain <- ''
	eval(parse(text=paste("chain <- chains$", build.conversion, sep="")))
	suppressMessages(res_GRL <- rtracklayer::liftOver(dGR, chain))
	res_GR <- GenomicRanges::unlist(res_GRL)
	
	if(merged){	
		mcols_data <- GenomicRanges::mcols(dGR)
		if(is.null(names(dGR))){
			names(dGR) <- 1:length(dGR)
		}
		names_data <- names(dGR)
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Finally, keep the first range if multiple found (%s) ...", as.character(now)), appendLF=T)
		}
	
		## keep only the first range (if multiple)
		res_df <- GenomicRanges::as.data.frame(res_GR, row.names=NULL)
		uid <- names(res_GR)
		res_ls <- split(x=res_df[,c(1:3,5)], f=uid)
		ls_df <- lapply(res_ls, function(x){
			c(as.character(unique(x$seqnames))[1],min(x$start), max(x$end), as.character(unique(x$strand))[1])
		})
		df <- do.call(rbind, ls_df)
	
		## construct GR object
		gr <- GenomicRanges::GRanges(
				seqnames=S4Vectors::Rle(df[,1]),
				ranges = IRanges::IRanges(start=as.numeric(df[,2]), end=as.numeric(df[,3])),
				strand = S4Vectors::Rle(df[,4])
			)
	
		## append back meta data
		ind <- as.numeric(rownames(df))
		names(gr) <- names_data[ind]
		GenomicRanges::mcols(gr) <- mcols_data[ind,]
		
		res_GR <- gr
	}
	
####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
	
	invisible(res_GR)
}
