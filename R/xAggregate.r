#' Function to aggregate data respecting number of features
#'
#' \code{xAggregate} is supposed to aggregate data respecting number of features. Per row, the aggregated is the sum of two items: the number of features, and the sum of all but scaled into [0,0.9999999]. Also supported is the rank-transformation of the input data per column, binned into the predefined number of discrete bins.
#'
#' @param data a data frame. The aggregation is done across columns per row. Each cell should contain positive values or NA; if infinite, it will be replaced with the maximum finite value
#' @param bin logical to indicate whether the input data per column is rank-transformed into the predefined number of discrete bins. By default, it sets to false
#' @param nbin the number of discrete bins. By default, it sets to 10 (only works when bin is true)
#' @param scale.log logical to indicate whether the per-row sum is log-scaled. By default, it sets to true
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a data frame with an appended column 'Aggregate'
#' @note None
#' @export
#' @seealso \code{\link{xAggregate}}
#' @include xAggregate.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' # HiC-gene pairs per cell types/states
#' g <- xRDataLoader(RData.customised='ig.PCHiC', RData.location=RData.location)
#' df <- do.call(cbind, igraph::edge_attr(g))
#' # aggregate over cell types/states
#' data <- df
#' data[data<5] <- NA
#' res <- xAggregate(data)
#' }

xAggregate <- function(data, bin=F, nbin=10, scale.log=T, verbose=T)
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
	
	data_input <- data
	
	## whether to bin first
	if(bin){
		if(verbose){
			message(sprintf("The input data of %d X %d is first binned into %d ranks per column (%s) ...", nrow(data), ncol(data), nbin, as.character(Sys.time())), appendLF=T)
		}
	
		for(j in 1:ncol(data)){
			x <- data[,j]
			y <- x[!is.na(x)]
			cut_index <- as.numeric(cut(y, breaks=min(y)+(max(y)-min(y))*seq(0, 1, length.out=nbin+1)))
			cut_index[is.na(cut_index)] <- 1
			x[!is.na(x)] <- cut_index
			data[,j] <- x
		}
	
	}else{
		if(verbose){
			message(sprintf("The input data of %d X %d is directly used for aggregation (%s) ...", nrow(data), ncol(data), as.character(Sys.time())), appendLF=T)
		}
		
	}
	
    ####################################################################################
	num_nonna <- apply(data, 1, function(x) sum(!is.na(x)))
	sum_nonna <- apply(data, 1, function(x) sum(x,na.rm=T))
	
	### avoid Inf
	tmp <- max(sum_nonna[!is.infinite(sum_nonna)])
	sum_nonna[is.infinite(sum_nonna)] <- tmp
	###
	
	if(scale.log){
		if(verbose){
			message(sprintf("The per-row sum is log-scaled"), appendLF=T)
		}
		sum_nonna[sum_nonna==0] <- NA
		sum_nonna <- log(sum_nonna)
	}else{
		if(verbose){
			message(sprintf("The per-row sum is WITHOUT log-scaled"), appendLF=T)
		}
	}
	
	scale_sum <- (sum_nonna - min(sum_nonna,na.rm=T)) / (max(sum_nonna,na.rm=T) - min(sum_nonna,na.rm=T)) * 0.9999999
	data_input$Aggregate <- num_nonna + scale_sum
    
  ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(data_input)
}
