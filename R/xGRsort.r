#' Function to sort by chromosomes/seqnames, start and end coordinates of the intervals.
#'
#' \code{xGRsort} is supposed to sort by chromosomes/seqnames, start and end coordinates of the intervals. 
#'
#' @param data input genomic regions (GR). GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'
#' @return index
#' @export
#' @seealso \code{\link{xGRsort}}
#' @include xGRsort.r
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
#' ind <- xGRsort(data)
#' }

xGRsort <- function(data)
{
	data <- gsub(',.*','',data)
	## gr has unique regions
	gr <- xGR(data, format='chr:start-end')
	gr <- GenomeInfoDb::sortSeqlevels(gr)
	gr <- sort(gr)
	
	############################
	## very important (because of non-redundant)
	df <- data.frame(sid=1:length(gr), gr=names(gr), stringsAsFactors=F)
	ind <- match(data, df$gr)
	df <- df[ind, ]
	df$oid <- 1:length(data)
	sid <- NULL
	df <- df %>% dplyr::arrange(sid)
	############################
	
  	invisible(df$oid)

}
