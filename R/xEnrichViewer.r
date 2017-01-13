#' Function to view enrichment results
#'
#' \code{xEnrichViewer} is supposed to view results of enrichment analysis. 
#'
#' @param eTerm an object of class "eTerm"
#' @param top_num the number of the top terms (sorted according to 'sortBy' below) will be viewed
#' @param sortBy which statistics will be used for sorting and viewing gene sets (terms). It can be "adjp" or "fdr" for adjusted p value (FDR), "pvalue" for p value, "zscore" for enrichment z-score, "fc" for enrichment fold change, "nAnno" for the number of sets (terms), "nOverlap" for the number in overlaps, and "none" for ordering according to ID of terms
#' @param decreasing logical to indicate whether to sort in a decreasing order. If it is null, it would be true for "zscore", "nAnno" or "nOverlap"; otherwise it would be false
#' @param details logical to indicate whether the detailed information of gene sets (terms) is also viewed. By default, it sets to false for no inclusion
#' @return
#' a data frame with following components:
#' \itemize{
#'  \item{\code{id}: term ID; as rownames}
#'  \item{\code{name}: term name}
#'  \item{\code{nAnno}: number in members annotated by a term}
#'  \item{\code{nOverlap}: number in overlaps}
#'  \item{\code{fc}: enrichment fold changes}
#'  \item{\code{zscore}: enrichment z-score}
#'  \item{\code{pvalue}: nominal p value}
#'  \item{\code{adjp}: adjusted p value (FDR)}
#'  \item{\code{distance}: term distance; optional, it is only appended when "details" is true}
#'  \item{\code{members}: members (represented as Gene Symbols) in overlaps; optional, it is only appended when "details" is true}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}
#' @include xEnrichViewer.r
#' @examples
#' \dontrun{
#' xEnrichViewer(eTerm)
#' }

xEnrichViewer <- function(eTerm, top_num=10, sortBy=c("adjp","fdr","pvalue","zscore","fc","nAnno","nOverlap","none"), decreasing=NULL, details=F) 
{
    
    sortBy <- match.arg(sortBy)
    
    if(sortBy=='fdr'){
    	sortBy <- 'adjp'
    }
    
    if(is.logical(eTerm)){
        stop("There is no enrichment in the 'eTerm' object.\n")
    }
    
    if(class(eTerm) == "data.frame" ){
        warnings("The function apply to a 'data.frame' object.\n")
        
		if( is.null(top_num) ){
			top_num <- nrow(eTerm)
		}
		if ( top_num > nrow(eTerm) ){
			top_num <- nrow(eTerm)
		}
		top_num <- as.integer(top_num)
        
        tab <- eTerm
        
    }else if(class(eTerm) == "eTerm" ){
	
		if( is.null(top_num) ){
			top_num <- length(eTerm$term_info$id)
		}
		if ( top_num > length(eTerm$term_info$id) ){
			top_num <- length(eTerm$term_info$id)
		}
		top_num <- as.integer(top_num)
	
		if(dim(eTerm$term_info)[1]==1){
			tab <- data.frame( name         = as.character(eTerm$term_info$name),
							   nAnno        = as.numeric(sapply(eTerm$annotation,length)),
							   nOverlap     = as.numeric(sapply(eTerm$overlap,length)),
							   fc         	= as.numeric(eTerm$fc),
							   zscore       = as.numeric(eTerm$zscore),
							   pvalue       = as.numeric(eTerm$pvalue),
							   adjp         = as.numeric(eTerm$adjp),
							   distance     = eTerm$term_info$distance,
							   members      = sapply(eTerm$overlap, function(x) paste(x,collapse=', ')),
							   stringsAsFactors=F
							  )
		}else{
			tab <- data.frame( name         = as.character(eTerm$term_info$name),
							   nAnno        = as.numeric(sapply(eTerm$annotation,length)),
							   nOverlap     = as.numeric(sapply(eTerm$overlap,length)),
							   fc         	= as.numeric(eTerm$fc),
							   zscore       = as.numeric(eTerm$zscore),
							   pvalue       = as.numeric(eTerm$pvalue),
							   adjp         = as.numeric(eTerm$adjp),
							   distance     = eTerm$term_info$distance,
							   members      = sapply(eTerm$overlap, function(x) paste(x,collapse=', ')),
							   stringsAsFactors=F
							  )
		}
	
		rownames(tab) <- eTerm$term_info$id
	}
	    
    if(details == T){
        res <- tab[,c(1:9)]
    }else{
        res <- tab[,c(1:7)]
    }
    
    if(is.null(decreasing)){
        if(sortBy=="zscore" | sortBy=="nAnno" | sortBy=="nOverlap"){
            decreasing <- T
        }else{
            decreasing <- F
        }
    }
    
    switch(sortBy, 
    	adjp={res <- res[with(res,order(adjp,-zscore))[1:top_num],]},
    	pvalue={res <- res[with(res,order(pvalue,-zscore))[1:top_num],]},
    	zscore={res <- res[with(res,order(-zscore,adjp))[1:top_num],]},
    	fc={res <- res[with(res,order(-fc,adjp))[1:top_num],]},
    	nAnno={res <- res[with(res,order(-nAnno,adjp,-zscore))[1:top_num],]},
    	nOverlap={res <- res[with(res,order(-nOverlap,adjp,-zscore))[1:top_num],]},
        none={res <- res[order(rownames(res), decreasing=decreasing)[1:top_num],]}
    )
    
    if(sortBy=='none'){
    	suppressWarnings(flag <- all(!is.na(as.numeric(rownames(res)))))
    	if(flag){
    		res <- res[order(as.numeric(rownames(res)), decreasing=decreasing)[1:top_num],]
    	}
    }
    
    res
}
