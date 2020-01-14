#' Function to prepare a gene_info RData file
#'
#' \code{xReadGenesets} is supposed to prepare a gene_info RData file.
#'
#' @param association.file an input file for association
#' @param set_info.file an input file for set_info
#' @param output.prefix a prefix for the output file
#' @param output.dir the output directory
#' @return a GS object
#' @note None
#' @export
#' @seealso \code{\link{xReadGenesets}}
#' @include xReadGenesets.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' 
#' GS <- xReadGenesets(association.file="GENE2GOMF.txt", set_info.file="GO.txt", output.prefix="org.Hs.egGOMF")
#' }

xReadGenesets <- function(association.file=NULL, set_info.file=NULL, output.prefix=NULL, output.dir="./")
{
    
    setID <- NULL
    
    if(!is.null(association.file) & !is.null(set_info.file) & !is.null(output.prefix)){
    	output.file <- paste0(output.dir, output.prefix, ".RData")
    	
    	## import set info
    	set_info <- readr::read_delim(set_info.file, delim="\t")
    	colnames(set_info) <- c("setID","name","namespace","distance")

   	 	## import genesets
   	 	#association <- readr::read_delim(association.file, col_names=F, delim="\t", skip=1)
   	 	association <- readr::read_delim(association.file, delim="\t")
		colnames(association) <- c("GeneID","setID")
		
		## focus those gene sets in common
		set_info %>% dplyr::semi_join(association, by="setID") %>% dplyr::arrange(setID) -> set_info
    	## define genesets
    	association %>% dplyr::semi_join(set_info, by="setID") -> df
    	gs <- split(x=df$GeneID, f=df$setID)
   	
    	## GS
    	GS <- list(gs=gs, set_info=set_info)
   	 	class(GS) <- "GS"
   	 	
   	 	do.call(assign, list(output.prefix, GS))
    	save(list=output.prefix, file=output.file, compress="xz")
    	
    	message(sprintf("Saved into '%s' (%s)!", output.file, as.character(Sys.time())), appendLF=T)
    	invisible(GS)
    	
    }
}
