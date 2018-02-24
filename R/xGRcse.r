#' Function to create a vector for genomic regions
#'
#' \code{xGRcse} is supposed to create genomic regions in the format of 'chr:start-end'. 
#'
#' @param data input genomic regions (GR). If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "bed" or "GRanges"
#' @return a vector for genomic regions the format of 'chrN:start-end'  
#' @export
#' @seealso \code{\link{xGRcse}}
#' @include xGRcse.r
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
#'
#' # b) create a GRanges object
#' cse <- xGRcse(gr)
#' }

xGRcse <- function(data, format=c("GRanges","data.frame","bed"))
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
	
	dGR <- xGR(data, format=format)
	df <- as.data.frame(dGR, row.names=NULL)
	vec_cse <- paste(df$seqnames,':',df$start,'-',df$end, sep='')
	
  	invisible(vec_cse)

}
