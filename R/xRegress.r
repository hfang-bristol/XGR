#' Function to regress data according to principle components (PCs)
#'
#' \code{xRegress} is supposed to regress data according to principle components (PCs). 
#
#' @param data a data matrix/frame with, for exampe, genes in rows and samples in columns
#' @param center logical to indicate whether the input data columns should be shifted to be zero centered when calculating PCs
#' @param scale logical to indicate whether the input data columns should have unit variance when calculating PCs
#' @param which.PCs a vector specifying which PCs are used for being regressed out. If NULL (by default), no gression is done
#' @return 
#' a list with two componets:
#' \itemize{
#'  \item{\code{regressed}: the regressed data with the same dimension as the input data}
#'  \item{\code{PCs}: a data matrix of PCs X samples}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xRegress}}
#' @include xRegress.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' \dontrun{
#' data(Fang)
#' ls_res <- xRegress(Fang, which.PCs=1)
#' gp <- xHeatmap(ls_res$PCs)
#' gp
#' }

xRegress <- function(data, center=TRUE, scale=TRUE, which.PCs=NULL)
{

	data <- t(data)
	
	# Identify PCs
	res <- stats::prcomp(as.matrix(data), center=center, scale.=scale)
	PCs <- t(res$x)
	
	if(length(which.PCs)==0){
		regressed.data <- t(data)
		
	}else{
		# which PCs to regress out
		if(length(which.PCs)==1){
			expression.pc <- matrix(PCs[which.PCs,], ncol=1)
			rownames(expression.pc) <- colnames(PCs)
		}else{
			expression.pc <- t(PCs[which.PCs,])
		}
	
		# regress out
		regressed.data <- apply(data, 2, function(x){
			fit <- stats::lm(as.matrix(x) ~ as.matrix(expression.pc))
			if(ncol(expression.pc)==1){
				x - expression.pc * fit$coefficients[-1]
			}else{
				x - expression.pc %*% fit$coefficients[-1]
			}
		})
		rownames(regressed.data)<-rownames(data)
		colnames(regressed.data)<-colnames(data)
		regressed.data <- t(regressed.data)
	}
    
    ls_res <- list(regressed = regressed.data,
    			   PCs = PCs
                 )
    
    invisible(ls_res)
}
