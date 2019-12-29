#' Function to mark a network within another network
#'
#' \code{xMarkNet} is supposed to mark a network within another network.
#'
#' @param ig1 a "igraph" object within which the mark happens
#' @param ig2 a "igraph" object to be marked. Only overlapped nodes and edges are marked
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' an object of class "igraph" appended with a node attribute 'mark' (0 for the background, 1 for the marked) and an edge attribute 'mark' (0 for the background, 1 for the marked)
#' @note none
#' @export
#' @seealso \code{\link{xMarkNet}}
#' @include xMarkNet.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' ig1 <- xDefineNet(network="KEGG_environmental", RData.location=RData.location)
#' ig2 <- xDefineNet(network="KEGG_organismal", RData.location=RData.location)
#' ig <- xMarkNet(ig1, ig2)
#' 
#' E(ig)$color <- ifelse(E(ig)$mark==0, 'lightblue1', 'darkgreen')
#' E(ig)$color.alpha <- ifelse(E(ig)$mark==0, 0.3, 0.7)
#' ig <- ig %>% xLayout("gplot.layout.fruchtermanreingold")
#' gp <- xGGnetwork(ig, node.xcoord='xcoord', node.ycoord='ycoord', node.color="mark", colormap="orange-darkgreen", node.color.alpha=0.7, edge.color="color", edge.color.alpha="color.alpha", edge.arrow.gap=0) + theme(legend.position='none')
#' }

xMarkNet <- function(ig1, ig2, verbose=TRUE)
{

   	if(!(any(class(ig1) %in% "igraph")) & !(any(class(ig2) %in% "igraph"))){
		return(NULL)
	}
	
	flag_direct <- 'undirect'
	vec <- sapply(list(ig1,ig2), igraph::is_directed)
	if(sum(vec)==length(vec)){
		flag_direct <- 'direct'
	}
	
	######################
	mark <- name <- NULL
	## ig1
    node1 <- ig1 %>% igraph::as_data_frame("vertices")
	edge1 <- ig1 %>% igraph::as_data_frame("edges")
	## ig2
    node2 <- ig2 %>% igraph::as_data_frame("vertices") %>% dplyr::mutate(mark=1) %>% dplyr::select(name,mark)
	edge2 <- ig2 %>% igraph::as_data_frame("edges") %>% dplyr::mutate(mark=1)
	
	## ig1 marked by ig2
	node1 <- node1 %>% dplyr::left_join(node2, by=c("name")) %>% dplyr::mutate(mark=ifelse(is.na(mark),0,mark))
	edge1 <- edge1 %>% dplyr::left_join(edge2, by=c("from","to")) %>% dplyr::mutate(mark=ifelse(is.na(mark),0,mark))
	if(flag_direct=='direct'){
		ig <- igraph::graph_from_data_frame(d=edge1, directed=TRUE, vertices=node1)
	}else{
		ig <- igraph::graph_from_data_frame(d=edge1, directed=FALSE, vertices=node1)
	}

	if(verbose){
		message(sprintf("The '%s' network (%d nodes and %d edges) has %d nodes and %d edges marked", flag_direct, vcount(ig), ecount(ig), sum(V(ig)$mark), sum(E(ig)$mark), appendLF=TRUE))
	}

    invisible(ig)
}
