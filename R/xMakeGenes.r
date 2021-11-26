#' Function to prepare a gene_info RDS file
#'
#' \code{xMakeGenes} is supposed to prepare a gene_info RDS file.
#'
#' @param gene_info.file an input file for gene_info
#' @param output.prefix a prefix for the output file
#' @param output.dir the output directory
#' @param stamp the stamp associated with this RDS file. By default it is the date when the file created
#' @return an EG object
#' @note None
#' @export
#' @seealso \code{\link{xMakeGenes}}
#' @include xMakeGenes.r
#' @examples
#' \dontrun{
#' EG <- xMakeGenes(gene_info.file="GeneInfo.9606", output.prefix="org.Hs.eg")
#' }

xMakeGenes <- function(gene_info.file=NULL, output.prefix=NULL, output.dir="./", stamp=as.Date(Sys.time()))
{
    
    GeneID <- NULL
    
    if(!is.null(gene_info.file) & !is.null(output.prefix)){
    	output.file <- paste0(output.dir, output.prefix, ".RDS")
    	
    	gene_info <- readr::read_delim(gene_info.file, delim="\t")
    	gene_info %>% arrange(GeneID) -> gene_info
    	
    	EG <- list(gene_info=gene_info, stamp=stamp)
   	 	class(EG) <- "EG"
    	saveRDS(EG, file=output.file, compress="gzip")
    	#readr::write_rds(EG, path=output.file, compress="gz")
    	
    	message(sprintf("Saved into '%s' (%s)!", output.file, as.character(Sys.time())), appendLF=T)
    	invisible(EG)
    	
    }
}
