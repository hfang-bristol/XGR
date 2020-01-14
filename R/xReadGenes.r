#' Function to prepare a gene_info RData file
#'
#' \code{xReadGenes} is supposed to prepare a gene_info RData file.
#'
#' @param gene_info.file an input file for gene_info
#' @param output.prefix a prefix for the output file
#' @param output.dir the output directory
#' @return an EG object
#' @note None
#' @export
#' @seealso \code{\link{xReadGenes}}
#' @include xReadGenes.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' 
#' EG <- xReadGenes(gene_info.file="GeneInfo.9606", output.prefix="org.Hs.eg")
#' }

xReadGenes <- function(gene_info.file=NULL, output.prefix=NULL, output.dir="./")
{
    
    GeneID <- NULL
    
    if(!is.null(gene_info.file) & !is.null(output.prefix)){
    	output.file <- paste0(output.dir, output.prefix, ".RData")
    	
    	gene_info <- readr::read_delim(gene_info.file, delim="\t")
    	gene_info %>% arrange(GeneID) -> gene_info
    	
    	EG <- list(gene_info=gene_info)
   	 	class(EG) <- "EG"
   	 	
   	 	do.call(assign, list(output.prefix, EG))
    	save(list=output.prefix, file=output.file, compress="xz")
    	
    	message(sprintf("Saved into '%s' (%s)!", output.file, as.character(Sys.time())), appendLF=T)
    	invisible(EG)
    	
    }
}
