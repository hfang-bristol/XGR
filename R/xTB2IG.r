#' Function to convert an igraph from one or two tibbles
#'
#' \code{xTB2IG} is supposed to convert an igraph from one or two tibbles.
#'
#' @param edges a tibble or data frame for edge attributes
#' @param nodes a tibble or data frame for node attributes. It can be NULL
#' @param directed a logic specifying whether to create a directed graph. By default it is false
#' @param intersected a logic or NULL specifying whether to restrict to only the intersected vertices (from input nodes and edges). By default it is NULL; it will be forced true (if the input nodes do not contain all vertices in the input edges); however, the input nodes can contain vertices not in the input edges. It only works when the input nodes are not NULL
#' @return
#' an igraph object
#' @note none
#' @export
#' @seealso \code{\link{xIG2TB}}
#' @include xTB2IG.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' ig <- xDefineNet(network="KEGG", RData.location=RData.location)
#' ig %>% xIG2TB('edges') -> edges
#' ig %>% xIG2TB('nodes') -> nodes
#'
#' xTB2IG(edges, nodes) -> ig
#' }

xTB2IG <- function(edges, nodes=NULL, directed=F, intersected=NULL)
{
    if(any(class(edges) %in% c("tbl"))){
    	edges <- edges %>% as.data.frame()
    }
    
    if(any(class(nodes) %in% c("tbl"))){
    	nodes <- nodes %>% as.data.frame()
    }
    
    if(!is.null(nodes)){
    	colnames(edges)[1:2] <- c("from","to")
    	colnames(nodes)[1] <- "name"
    	
    	from <- to <- name <- NULL
    	
    	# if intersected is NULL, it will be forced to true if the input nodes do not contain all vertices in the input edges
    	if(is.null(intersected)){
    		dplyr::bind_rows(edges %>% dplyr::transmute(name=from), edges %>% dplyr::transmute(name=to)) %>% dplyr::n_distinct() -> num_nodes_in_edges
    		nodes %>% dplyr::select(name) %>% dplyr::n_distinct() -> num_nodes
    		if(num_nodes < num_nodes_in_edges){
    			# force to the restriction
    			intersected <- T
    		}
    	}
    	
    	if(!is.null(intersected) && intersected){
			# edges and nodes restricted to intersected nodes
			edges %>% dplyr::semi_join(nodes, by=c("from"="name")) %>% 
	dplyr::semi_join(nodes, by=c("to"="name")) -> edges
	
			nodes %>% dplyr::semi_join(dplyr::bind_rows(edges %>% dplyr::transmute(name=from), edges %>% dplyr::transmute(name=to)) %>% dplyr::distinct(name), by="name") -> nodes    	
    	}
    }
    
    ig <- igraph::graph_from_data_frame(d=edges, directed=directed, vertices=nodes)
    return(ig)
}
