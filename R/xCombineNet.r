#' Function to combine networks from a list of igraph objects
#'
#' \code{xCombineNet} is supposed to combine networks from a list of igraph objects.
#'
#' @param list_ig a list of "igraph" objects or a "igraph" object
#' @param combineBy how to resolve edges from a list of "igraph" objects. It can be "intersect" for intersecting edges and "union" for unionising edges (by default)
#' @param attrBy the method used to extract node attributes. It can be "intersect" for intersecting node attributes (by default) and "union" for unionising node attributes
#' @param keep.all.vertices logical to indicate whether all nodes are kept when intersecting edges. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' an object of class "igraph"
#' @note none
#' @export
#' @seealso \code{\link{xDefineNet}}
#' @include xCombineNet.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' g1 <- xDefineNet(network="KEGG_environmental", RData.location=RData.location)
#' g2 <- xDefineNet(network="KEGG_organismal", RData.location=RData.location)
#' ls_ig <- list(g1, g2)
#' ig <- xCombineNet(ls_ig, combineBy='union', attrBy="intersect", verbose=TRUE)
#' }

xCombineNet <- function(list_ig, combineBy=c('union','intersect'), attrBy=c("intersect","union"), keep.all.vertices=FALSE, verbose=TRUE)
{
    combineBy <- match.arg(combineBy)
    attrBy <- match.arg(attrBy) 
    
   	if(any(class(list_ig) %in% c("igraph"))){
		ls_ig <- list(list_ig)
	}else if(class(list_ig)=="list"){
		## Remove null elements in a list
		ls_ig <- base::Filter(base::Negate(is.null), list_ig)
		if(length(ls_ig)==0){
			return(NULL)
		}
	}else{
		stop("The function must apply to 'list' of 'igraph' objects or a 'igraph' object.\n")
	}
	
	flag_direct <- 'undirect'
	vec <- sapply(ls_ig, igraph::is_directed)
	if(sum(vec)==length(vec)){
		flag_direct <- 'direct'
	}
	
	######################
	## get node attributes
	ls_node_attr <- lapply(1:length(ls_ig), function(i){
		ig <- ls_ig[[i]]
		igraph::vertex_attr_names(ig)
	})
	if(attrBy=='intersect'){
		node_attr <- base::Reduce(intersect, ls_node_attr)
	}else if(attrBy=='union'){
		node_attr <- base::Reduce(union, ls_node_attr)
	}

	######################	
	## df_node
	ls_node <- lapply(1:length(ls_ig), function(i){
		ig <- ls_ig[[i]]
		nodes <- igraph::get.data.frame(ig, what="vertices")[,node_attr]
	})
	df_node <- unique(do.call(rbind, ls_node))
	
	## df_edge
	ls_edge <- lapply(1:length(ls_ig), function(i){
		ig <- ls_ig[[i]]
		relations <- igraph::get.data.frame(ig, what="edges")[,c(1,2)]
		if(flag_direct=='undirect'){
			ind <- which(relations[,1] > relations[,2])
			relations[ind, c(1:2)] <- relations[ind, c(2,1)]
			relations <- unique(relations)
		}
		return(relations)
	})
	df_edge <- do.call(rbind, ls_edge)
	if(combineBy=='intersect'){
		from <- to <- NULL
		#ftc <- plyr::ddply(df_edge, .(from,to), nrow)
		ftc <- plyr::ddply(df_edge, c('from','to'), nrow)
		ind <- which(ftc[,3]==length(ls_ig))
		df_edge <- ftc[ind, c(1,2)]
		if(!keep.all.vertices){
			tmp_nodes <- unique(c(df_edge[,1], df_edge[,2]))
			ind <- match(df_node$name, tmp_nodes)
			df_node <- df_node[!is.na(ind),]
		}
	}else if(combineBy=='union'){
		df_edge <- unique(df_edge)
	}

	######################	
	## combined_ig
	if(flag_direct=='direct'){
		combined_ig <- igraph::graph.data.frame(d=df_edge, directed=TRUE, vertices=df_node)
	}else{
		combined_ig <- igraph::graph.data.frame(d=df_edge, directed=FALSE, vertices=df_node)
	}
	
	if(verbose){
		message(sprintf("%d network(s) are combined (via '%s') into a '%s' network (%d nodes and %d edges) with %d node attributes (via '%s')", length(ls_ig), combineBy, flag_direct, vcount(combined_ig), ecount(combined_ig), length(node_attr), attrBy, appendLF=TRUE))
	}

    invisible(combined_ig)
}