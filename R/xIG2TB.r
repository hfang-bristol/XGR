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
#' @seealso \code{\link{xIG2TB}}
#' @include xIG2TB.r
#' @examples
#' set.seed(825)
#' ig <- sample_pa(20)
#' V(ig)$name <- seq(1,vcount(ig))
#' ig %>% xIG2TB('edges')
#' ig %>% xIG2TB('nodes')

xIG2TB <- function(ig, what=c('edges','nodes'))
{
    what <- match.arg(what)
    
   	if(is(ig,"igraph")){
   		
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
