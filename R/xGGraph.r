#' Function to visualise an igraph object using ggraph
#'
#' \code{xGGraph} is supposed to visualise an igraph object using ggraph, with nodes/tips labelled (aligned to left-right or top-bottom edges).
#'
#' @param ig an object of class "igraph" with node attribute 'name'. It could be a 'phylo' object converted to. Note: the node/leave labels would be the node attribute 'name' unless the node attribute 'label' is explicitely provided
#' @param layout the layout supported in ggraph::create_layout. This can be ggraph layouts 'partition' (by default), 'dendrogram', 'circlepack', 'treemap' (-1,1). This can be also igraph-supported layout ('nicely','fr','kk','sugiyama','randomly','star','circle','gem','dh','graphopt','grid','mds','drl','lgl','sphere')
#' @param circular the logic specifying whether or not circular representations. This will be disabled implicitly if the layout does not support circularity
#' @param leave the logic specifying whether or not only leaves (nodes/labellings) shown. This can be disenabled if the layout does not support tips
#' @param node.label.size the text size of the leave labelings. By default, it is 2. If 0, all labellings will be disabled
#' @param node.label.direction the leave label direction. It can be "none", "leftright" (aligned to the left- and right-most edge) and "topbottom" (aligned to the top- and bottom-most edge)
#' @param node.label.color the color of the leave labelings
#' @param node.label.alpha the alpha of the leave labelings
#' @param node.label.wrap the wrap width of the leave labelings
#' @param node.label.offset the offset of the leave labelings aligned to the edge. It is defined as relative to the range of limits (x-limit for left-right, and y-limit for top-bottom)
#' @param node.size the size of the leave nodes. By default, it is 0
#' @param limit.expansion the x- and y-limit expansion. By default, it is NULL, decided by "node.label.offset"
#' @param edge the edge type. It can be "diagonal" (default) , "link" (straight lines), "arc", "fan" (curves of different curvature), "elbow"
#' @param edge.color the color of edges
#' @param edge.alpha the alpha of edges
#' @param edge.width the width of edges
#' @param ... additional graphic parameters (such as size, color) used in ggrepel::geom_text_repel to control labels
#' @return 
#' a ggplot2 object appended with 'ig' and 'data' which should contain columns 'x','y','name' (the same as V(ig)$name), 'label' (if not given in ig, a 'name' varient). Also contain 'leaf' (T/F), 'depth' (the number of step to the root) for tree-like graph with certain layouts.
#' @note none
#' @export
#' @seealso \code{\link{xGGraph}}
#' @include xGGraph.r
#' @examples
#' \dontrun{
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' 
#' AA.template <- xRDataLoader("AA.template", RData.location=RData.location)
#' # consensus tree
#' ig <- AA.template$consensus$ig
#'
#' # Default: partition-like circular layout
#' # none
#' gp <- xGGraph(ig, node.label.direction="none", node.label.wrap=50)
#' # leftright
#' gp <- xGGraph(ig, node.label.direction="leftright", node.label.wrap=50, node.label.offset=0.5)
#' # topbottom
#' gp <- xGGraph(ig, node.label.direction="topbottom", node.label.wrap=50, node.label.offset=0.5)
#' 
#' # advanced usage
#' ## ggraph layouts
#' gp <- xGGraph(ig, layout='dendrogram', node.label.direction="leftright")
#' gp <- xGGraph(ig, layout='treemap')
#' gp <- xGGraph(ig, layout='circlepack')
#' ## igraph layouts
#' set.seed(825)
#' gp <- xGGraph(ig, layout='nicely', node.label.direction="leftright")
#' gp <- xGGraph(ig, layout='kk')
#' gp <- xGGraph(ig, layout='fr', node.label.direction="leftright")
#' gp <- xGGraph(ig, layout='gem')
#' }

xGGraph <- function(ig, layout='partition', circular=T, leave=T, node.label.size=2, node.label.direction=c('none','leftright','topbottom'), node.label.color="steelblue", node.label.alpha=0.7, node.label.wrap=NULL, node.label.offset=0.5, node.size=2, limit.expansion=NULL, edge=c("diagonal","link","arc","fan","elbow"), edge.color='grey', edge.alpha=0.5, edge.width=0.5, ...)
{

    node.label.direction <- match.arg(node.label.direction)
    edge <- match.arg(edge)
	
	## how to convert a phylo object 'tree' into igraph object 'ig'
	if(0){
		## if internal node labels are duplicated, remove all (and add by as.igraph)
		if(any(duplicated(tree$node.label))) tree$node.label<-NULL
		ig <- as.igraph(tree, directed=T, use.labels=T)
	}
	
    if(!is(ig,"igraph")){
        warnings("The function must apply to the 'igraph' object.\n")
        return(NULL)
    }else{
		if(!all(c("name") %in% igraph::vertex_attr_names(ig))){
			warnings("The igraph object must have vertex attribute 'name'.\n")
			return(NULL)
		}
    }

	##################
	if(!("label" %in% igraph::vertex_attr_names(ig))){
		# append 'label'
		V(ig)$label <- V(ig)$name
	}
	## label wrap
	if(!is.null(node.label.wrap)){
		width <- as.integer(node.label.wrap)
		res_list <- lapply(V(ig)$label, function(x){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=width)
			if(length(y)>1){
				paste(y,collapse='\n')
			}else{
				y
			}
		})
		V(ig)$label <- unlist(res_list)
	}

	x <- y <- leaf <- label <- name <- NULL
	
	############
	#gp <- ggraph::ggraph(ig, layout=layout, circular=circular)
	## disable the circular layout if error
	if(length(suppressWarnings(tryCatch(gp <- ggraph::ggraph(ig, layout=layout, circular=circular), error=function(e) e, warning=function(w) w)))==2){
		warning("The layout does not support circularity!")
		circular <- F
		gp <- ggraph::ggraph(ig, layout=layout, circular=circular)
	}
	############
	
	## nodes/leave
	if(leave & !is.null(gp$data$leaf)){
		gp <- gp + ggraph::geom_node_point(aes(filter=leaf),size=node.size, color=edge.color,alpha=edge.alpha)
	}else{
		gp <- gp + ggraph::geom_node_point(size=node.size, color=edge.color,alpha=edge.alpha)	
	}

	## edges
	if(edge=="diagonal"){
		gp <- gp + ggraph::geom_edge_diagonal(color=edge.color,alpha=edge.alpha,width=edge.width)
	}else if(edge=="link"){
		gp <- gp + ggraph::geom_edge_link(color=edge.color,alpha=edge.alpha,width=edge.width)
	}else if(edge=="arc"){
		gp <- gp + ggraph::geom_edge_arc(color=edge.color,alpha=edge.alpha,width=edge.width)
	}else if(edge=="fan"){
		gp <- gp + ggraph::geom_edge_fan(color=edge.color,alpha=edge.alpha,width=edge.width)
	}else if(edge=="elbow"){
		gp <- gp + ggraph::geom_edge_elbow(color=edge.color,alpha=edge.alpha,width=edge.width)
	}
	
	if(node.label.size>0){
	
		## node/leave labels
		if(leave & !is.null(gp$data$leaf)){
			df <- subset(gp$data, leaf==T)
		}else{
			df <- gp$data
		}
	
		if(node.label.direction=='none'){
			gp <- gp + ggrepel::geom_text_repel(data=df, aes(x=x, y=y, label=label), color=node.label.color, size=node.label.size, alpha=node.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), ...)
		
		}else if(node.label.direction=='leftright'){
			offset <- (range(gp$data$x)[2]-range(gp$data$x)[1]) * node.label.offset
	
			root <- subset(gp$data, name==dnet::dDAGroot(ig))
		
			## left
			df1 <- subset(df, x < root$x)
			df1$nudge_x <- -1 * (df1$x - min(gp$data$x)) - offset
			gp <- gp + ggrepel::geom_text_repel(data=df1, aes(x=x, y=y, label=label), color=node.label.color, size=node.label.size, alpha=node.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), direction="y", hjust=0, nudge_x=df1$nudge_x)
		
			## right
			df2 <- subset(df, x >= root$x)
			df2$nudge_x <- (max(gp$data$x)-df2$x) + offset
			gp <- gp + ggrepel::geom_text_repel(data=df2, aes(x=x, y=y, label=label), color=node.label.color, size=node.label.size, alpha=node.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), direction="y", hjust=1, nudge_x=df2$nudge_x)
		
			if(is.null(limit.expansion)){
				gp <- gp + expand_limits(x=range(gp$data$x)*(1+2*node.label.offset), y=range(gp$data$y))
			}else{
				gp <- gp + expand_limits(x=c(-limit.expansion, limit.expansion), y=c(-limit.expansion, limit.expansion))
			}
	
		}else if(node.label.direction=='topbottom'){
			offset <- (range(gp$data$y)[2]-range(gp$data$y)[1]) * node.label.offset
	
			root <- subset(gp$data, name==dnet::dDAGroot(ig))
		
			## bottom
			df1 <- subset(df, y < root$y)
			df1$nudge_y <- -1 * (df1$y - min(gp$data$y)) - offset
			gp <- gp + ggrepel::geom_text_repel(data=df1, aes(x=x, y=y, label=label), color=node.label.color, size=node.label.size, alpha=node.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), direction="x", hjust=0, nudge_y=df1$nudge_y, angle=90)
		
			## top
			df2 <- subset(df, y >= root$y)
			df2$nudge_y <- (max(gp$data$y)-df2$y) + offset
			gp <- gp + ggrepel::geom_text_repel(data=df2, aes(x=x, y=y, label=label), color=node.label.color, size=node.label.size, alpha=node.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), direction="x", hjust=1, nudge_y=df2$nudge_y, angle=90)
		
			if(is.null(limit.expansion)){
				gp <- gp + expand_limits(x=range(gp$data$x), y=range(gp$data$y)*(1+2*node.label.offset))
			}else{
				gp <- gp + expand_limits(x=c(-limit.expansion, limit.expansion), y=c(-limit.expansion, limit.expansion))
			}
		
		}
	}
		
	gp <- gp + ggraph::theme_graph(base_family="sans",plot_margin=margin(0,0,0,0))
	
	if(0){
		# order by tipid
		tipid <- NULL
		#gp$data <- gp$data %>% dplyr::arrange(tipid)
	}

	gp$ig <- ig
	
	invisible(gp)
}
