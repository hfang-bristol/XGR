#' Function to convert gene symbols to entrez geneid
#'
#' \code{xGeneID2Symbol} is supposed to convert gene symbols to entrez geneid.
#'
#' @param data an input vector containing gene symbols
#' @param org a character specifying an organism. Currently supported organisms are 'human' and 'mouse'. It can be an object 'EG'
#' @param details logical to indicate whether to result in a data frame (in great details). By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return a vector containing symbol with 'NA' for the unmatched if (details set to false); otherwise, a data frame is returned
#' @note none.
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xSocialiserGenes}}
#' @include xGeneID2Symbol.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' 
#' # a) provide the input Genes of interest (eg 100 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg')
#' GeneID <- sample(org.Hs.eg$gene_info$GeneID, 100)
#' GeneID
#' 
#' # b) convert into GeneID
#' Symbol <- xGeneID2Symbol(GeneID)
#' 
#' # c) convert into a data frame
#' df <- xGeneID2Symbol(GeneID, details=TRUE)
#' 
#' 
#' # advanced use
#' df <- xGeneID2Symbol(GeneID, org=org.Hs.eg, details=TRUE)
#' }

xGeneID2Symbol <- function(data, org=c("human","mouse"), details=F, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    
    if (!is.vector(data)){
        stop("The input data must be a vector.\n")
    }
    GeneID <- as.numeric(data)
    
    if(verbose){
        message(sprintf("%d GeneID of input data (%s)", length(data), as.character(Sys.time())), appendLF=T)
    }
    
    ## load Enterz Gene information
	if(class(org) == "EG"){
		df_eg <- org$gene_info
		if(verbose){
			message(sprintf("Customised organism (%s)", as.character(Sys.time())), appendLF=T)
		}
	}else{
		org <- org[1]
		if(org=='human'){
			df_eg <- xRDataLoader(RData.customised='org.Hs.eg', RData.location=RData.location, verbose=verbose)$gene_info
		}else if(org=='mouse'){
			df_eg <- xRDataLoader(RData.customised='org.Mm.eg', RData.location=RData.location, verbose=verbose)$gene_info
		}
		if(verbose){
			message(sprintf("%s organism (%s)", org, as.character(Sys.time())), appendLF=T)
		}
	}
	
	ind <- match(data, df_eg$GeneID)
	df_res <- df_eg[ind, ]
    
    if(verbose){
        message(sprintf("%d mappable but %d left unmappable (%s)", sum(!is.na(ind)), sum(is.na(ind)), as.character(Sys.time())), appendLF=T)
    }
	
	if(details){
		df_res <- data.frame(Input=GeneID, df_res, stringsAsFactors=F)
		return(df_res)
	}else{
		return(df_res$Symbol)
	}
}
