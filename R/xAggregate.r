#' Function to aggregate data respecting number of features
#'
#' \code{xAggregate} is supposed to aggregate data respecting number of features. Per row, the aggregated is the sum of two items: the number of features, and the sum of all but scaled into [0,0.9999999]
#'
#' @param data a data frame. The aggregation is done across columns per row
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a data frame with an appended column called 'Aggregate'
#' @note None
#' @export
#' @seealso \code{\link{xAggregate}}
#' @include xAggregate.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' # HiC-gene pairs per cell types/states
#' g <- xRDataLoader(RData.customised='ig.PCHiC', RData.location=RData.location)
#' df <- do.call(cbind, igraph::edge_attr(g))
#'
#' # aggregate over cell types/states
#' data <- df
#' data[data<5] <- NA
#' res <- xAggregate(data)
#' }

xAggregate <- function(data, verbose=T)
{
	
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
	
	if(is.matrix(data)){
		data <- as.data.frame(data)
	}
	
	num_nonzero <- apply(data, 1, function(x) sum(!is.na(x)))
	sum_nonzero <- apply(data, 1, function(x) sum(x,na.rm=T))
	scale_sum <- (sum_nonzero - min(sum_nonzero)) / (max(sum_nonzero) - min(sum_nonzero)) * 0.9999999
	data$Aggregate <- num_nonzero + scale_sum
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(data)
}
