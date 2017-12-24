#' Function to split a graph according to a node attribute
#'
#' \code{xGraphSplit} is supposed to split a graph according to a node attribute such as community or comp.
#'
#' @param g an object of class "igraph" (or "graphNEL") for a graph with such as a 'community' node attribute
#' @param node.attr a charatter specifying a node attribute. If NULL or no match, it returns NULL
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a list of igraph objects.
#' @export
#' @seealso \code{\link{xA2Net}}
#' @include xGraphSplit.r
#' @examples
#' # 1) generate a random bipartite graph
#' set.seed(123)
#' g <- sample_bipartite(100, 50, p=0.1)
#' V(g)$name <- V(g)
#' 
#' \dontrun{
#' # 2) obtain and append the community
#' cs <- igraph::cluster_louvain(g)
#' V(g)$community <- cs$membership
#' ls_ig <- xGraphSplit(g, node.attr="community")
#' }

xGraphSplit <- function(g, node.attr=NULL, verbose=TRUE)
{
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
	
	node.attrs <- igraph::vertex_attr_names(ig)
	flag <- F
	if(is.null(node.attr)){
		flag <- T
	}else if(!(node.attr %in% node.attrs)){
		flag <- T
	}
	if(flag){
		message("The igraph object must have a vertex attribute of interest.\n")
		return(NULL)
	}
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("The igraph object has %d nodes and %d edges (%s) ...", vcount(ig), ecount(ig), as.character(now)), appendLF=T)
	}
	
	df_nodes <- igraph::get.data.frame(ig,what="vertices")
	df_nodes$node.attr <- df_nodes[,node.attr]
	
	ls_tmp <- split(x=df_nodes$name, f=df_nodes$node.attr)
	ls_subg <- lapply(ls_tmp, function(x){
		g <- dnet::dNetInduce(ig, nodes_query=x, knn=0, largest.comp=TRUE)
	})
    
    invisible(ls_subg)
}


