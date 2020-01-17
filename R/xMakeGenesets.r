#' Function to prepare a gene_info RDS file
#'
#' \code{xMakeGenesets} is supposed to prepare a gene_info RDS file.
#'
#' @param association.file an input file for association
#' @param set_info.file an input file for set_info
#' @param output.prefix a prefix for the output file
#' @param output.dir the output directory
#' @param stamp the stamp associated with this RDS file. By default it is the date when the file created
#' @return a GS object
#' @note None
#' @export
#' @seealso \code{\link{xMakeGenesets}}
#' @include xMakeGenesets.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' 
#' GS <- xMakeGenesets(association.file="GENE2GOMF.txt", set_info.file="GO.txt", output.prefix="org.Hs.egGOMF")
#' }

xMakeGenesets <- function(association.file=NULL, set_info.file=NULL, output.prefix=NULL, output.dir="./", stamp=as.Date(Sys.time()))
{
    
    setID <- NULL
    
    if(!is.null(association.file) & !is.null(set_info.file) & !is.null(output.prefix)){
    	output.file <- paste0(output.dir, output.prefix, ".RDS")
    	
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
    	gs %>% tibble::enframe(name='setID', value='GeneID') -> gs

    	## GS
    	GS <- list(gs=gs, set_info=set_info, stamp=stamp)
   	 	class(GS) <- "GS"
    	saveRDS(GS, file=output.file, compress="gzip")
    	#readr::write_rds(GS, path=output.file, compress="gz")
    	
    	message(sprintf("Saved into '%s' (%s)!", output.file, as.character(Sys.time())), appendLF=T)
    	invisible(GS)
    	
    }
}
