#' Function to convert an igraph into a tibble for nodes or edges
#'
#' \code{xIG2TB} is supposed to convert an igraph into a tibble for nodes or edges.
#'
#' @param ig an "igraph" object
#' @param what what to extract. It can be "edges" for edges and "nodes" for nodes
#' @return
#' a tibble object
#' @note none
#' @export
#' @seealso \code{\link{xTB2IG}}
#' @include xIG2TB.r
#' @examples	
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' ig <- xDefineNet(network="KEGG", RData.location=RData.location)
#' ig %>% xIG2TB('edges')
#' ig %>% xIG2TB('nodes')
#' }

xIG2TB <- function(ig, what=c('edges','nodes'))
{
    what <- match.arg(what)
    
   	if(any(class(ig) %in% c("igraph"))){
   		
   		if(what=='edges'){
   			edges <- igraph::as_data_frame(ig, what="edges") %>% tibble::as_tibble()
   			return(edges)
   		}else if(what=='nodes'){
   			nodes <- igraph::as_data_frame(ig, what="vertices") %>% tibble::as_tibble()
   			return(nodes)
   		}
	}else{
		return(NULL)
	}
}
