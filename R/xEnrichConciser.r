#' Function to make enrichment results conciser by removing redundant terms
#'
#' \code{xEnrichConciser} is supposed to make enrichment results conciser by removing redundant terms. A redundant term (called 'B') is claimed if its overlapped part (A&B) with a more significant term (called 'A') meets both criteria: 1) |A&B| > 0.9*|B|; and 2) |A&B| > 0.5*|A|.
#'
#' @param eTerm an object of class "eTerm"
#' @param cutoff a cutoff vector used to remove redundant terms. By default, it has the first element 0.9 and the second element 0.5. It means, for a term (less significant; called 'B'), if there is a more significant term (called 'A'), their overlapped members cover at least 90% of the B's members, and also at least 50% of the A's members, then this term B will be defined as redundant and thus being removed
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
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

xEnrichConciser <- function(eTerm, cutoff=c(0.9,0.5), verbose=T) 
{
    
    if(is.null(eTerm)){
        warnings("There is no enrichment in the 'eTerm' object.\n")
        return(NULL)
    }
    
    if(class(eTerm) == "eTerm" ){
		
		cross <- eTerm$cross
		
		df <- xEnrichViewer(eTerm, top_num='all', sortBy="pvalue")
		ind <- match(rownames(df), rownames(cross))
		cross <- cross[ind, ind]
		
		nRedundant_1 <- matrix(0, nrow=ncol(cross), ncol=ncol(cross))
		nRedundant_2 <- matrix(0, nrow=ncol(cross), ncol=ncol(cross))
		for(j in seq(2, ncol(cross))){
			for(i in seq(1, j-1)){
				nRedundant_1[i,j] <- cross[i,j] >= cross[j,j]*cutoff[1]
				nRedundant_2[i,j] <- cross[i,j] >= cross[i,i]*cutoff[2]
			}
		}
		nRedundant <- apply(nRedundant_1 & nRedundant_2, 2, sum)
		names(nRedundant) <- colnames(cross)
		
		ind <- match(names(eTerm$adjp), names(nRedundant))
		nRedundant <- nRedundant[ind]
		
		## update eTerm by removing redundant terms
		flag <- nRedundant == 0
		
        if(verbose){
            now <- Sys.time()
            message(sprintf("\tAmong %d terms, there are %d non-redundant terms", length(flag), sum(flag)), appendLF=T)
        }
		
		eTerm$term_info <- eTerm$term_info[flag, ]
		eTerm$annotation <- eTerm$annotation[flag]
		eTerm$overlap <- eTerm$overlap[flag]
		eTerm$fc <- eTerm$fc[flag]
		eTerm$zscore <- eTerm$zscore[flag]
		eTerm$pvalue <- eTerm$pvalue[flag]
		eTerm$adjp <- eTerm$adjp[flag]
		eTerm$or <- eTerm$or[flag]
		eTerm$CIl <- eTerm$CIl[flag]
		eTerm$CIu <- eTerm$CIu[flag]
		eTerm$cross <- eTerm$cross[flag, flag]
		
		res <- eTerm
		
	}else{
		res <- NULL
		message("The provided object cannot be used for removing redundant terms.\n")
	}
    
    res
}
