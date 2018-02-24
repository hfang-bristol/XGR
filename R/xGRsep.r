#' Function to obtain separator index.
#'
#' \code{xGRsep} is supposed to obtain separator index. 
#'
#' @param data input genomic regions (GR). GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'
#' @return a vector for separator index
#' @export
#' @seealso \code{\link{xGRsep}}
#' @include xGRsep.r
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
#' cse <- xGRcse(gr)
#'
#' # b) sort index
#' ind <- xGRsort(cse)
#' data <- cse[ind]
#'
#' # c) get separator index
#' vec_sep <- xGRsep(data)
#' }

xGRsep <- function(data)
{
	
	x <- gsub(":.*|chr","",data)
	x[x=='X'] <- 23
	x[x=='Y'] <- 24
	x <- as.numeric(x)
	vec_sep <- cumsum(table(x))
	vec_sep <- vec_sep[-length(vec_sep)]
	
  	invisible(vec_sep)

}
