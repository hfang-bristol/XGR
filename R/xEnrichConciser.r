#' Function to make enrichment results conciser by removing redundant terms
#'
#' \code{xEnrichConciser} is supposed to make enrichment results conciser by removing redundant terms. 
#'
#' @param eTerm an object of class "eTerm"
#' @param cutoff the cutoff used to remove redunant terms. For a term (less significant), if there is a more significant term having members that cover at least the cutoff percentage (0.95 by default) of its members, then this term will be removed
#' @return
#' an object of class "eTerm", after redundant terms being removed.
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}
#' @include xEnrichConciser.r
#' @examples
#' \dontrun{
#' eTerm_concise <- xEnrichConciser(eTerm)
#' }

xEnrichConciser <- function(eTerm, cutoff=0.95) 
{
    
    if(is.logical(eTerm)){
        stop("There is no enrichment in the 'eTerm' object.\n")
    }
    
    if(class(eTerm) == "eTerm" ){
		
		cross <- eTerm$cross
		
		df <- xEnrichViewer(eTerm, top_num='all', sortBy="pvalue")
		ind <- match(rownames(df), rownames(cross))
		cross <- cross[ind, ind]
		
		nRedundant <- rep(0, ncol(cross))
		for(j in seq(2, ncol(cross))){
			i <- seq(1, j-1)
			nRedundant[j] <- sum(cross[i, j] >= cross[j,j]*cutoff)
		}
		names(nRedundant) <- colnames(cross)
		
		ind <- match(names(eTerm$adjp), names(nRedundant))
		nRedundant <- nRedundant[ind]
		
		## update eTerm by removing redundant terms
		flag <- nRedundant == 0
		eTerm$term_info <- eTerm$term_info[flag, ]
		eTerm$annotation <- eTerm$annotation[flag]
		eTerm$overlap <- eTerm$overlap[flag]
		eTerm$fc <- eTerm$fc[flag]
		eTerm$zscore <- eTerm$zscore[flag]
		eTerm$pvalue <- eTerm$pvalue[flag]
		eTerm$adjp <- eTerm$adjp[flag]
		eTerm$cross <- eTerm$cross[flag, flag]
		
		res <- eTerm
		
	}else{
		res <- NULL
	}
    
    res
}
