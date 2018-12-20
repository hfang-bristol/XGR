#' Function to draw a partition-like circular plot
#'
#' \code{xCpartition} is supposed to draw a partition-like circular plot (partition circular layout), with tips labelled (aligned to left-right or top-bottom edges).
#'
#' @param ig an object of class "igraph" with node attribute 'name'. It could be a 'phylo' object converted to. Note: the leave labels would be the node attribute 'name' unless the node attribute 'label' is explicitely provided
#' @param leave.label.direction the leave label direction. It can be "none", "leftright" (aligned to the left- and right-most edge) and "topbottom" (aligned to the top- and bottom-most edge)
#' @param leave.label.size the text size of the leave labelings. By default, it is 2
#' @param leave.label.color the color of the leave labelings
#' @param leave.label.alpha the alpha of the leave labelings
#' @param leave.label.wrap the wrap width of the leave labelings
#' @param leave.label.offset the offset of the leave labelings aligned to the edge. It is defined as relative to the range of limits (x-limit for left-right, and y-limit for top-bottom)
#' @param leave.size the size of the leave nodes. By default, it is 0
#' @param limit.expansion the x- and y-limit expansion. By default, it is NULL, decided by "leave.label.offset"
#' @param edge.color the color of edges
#' @param edge.alpha the alpha of edges
#' @param edge.width the width of edges
#' @param ... additional graphic parameters (such as size, color) used in ggrepel::geom_text_repel to control labels
#' @return 
#' a ggplot2 object appended with 'ig' and 'data' which should contain columns 'x','y', 'leaf' (T/F), 'depth' (the number of step to the root), 'name' (the same as V(ig)$name), 'tipid' (tip id), 'label' (if not given in ig, a 'name' varient).
#' @note none
#' @export
#' @seealso \code{\link{xCpartition}}
#' @include xCpartition.r
#' @examples
#' \dontrun{
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' 
#' AA.path <- xRDataLoader("AA.path", RData.location=RData.location)
#' # consensus tree
#' ig <- AA.path$consensus$ig
#' 
#' # none
#' gp <- xCpartition(ig, leave.label.direction="none", leave.label.wrap=50, leave.label.offset=0.5, leave.size=2)
#' head(gp$data %>% dplyr::arrange(tipid))
#' 
#' # leftright
#' gp <- xCpartition(ig, leave.label.direction="leftright", leave.label.wrap=50, leave.label.offset=0.5, leave.size=2)
#' 
#' # topbottom
#' gp <- xCpartition(ig, leave.label.direction="topbottom", leave.label.wrap=50, leave.label.offset=0.5, leave.size=2)
#' }

xCpartition <- function(ig, leave.label.direction=c('none','leftright','topbottom'), leave.label.size=2, leave.label.color="darkblue", leave.label.alpha=0.7, leave.label.wrap=NULL, leave.label.offset=0.5, leave.size=0, limit.expansion=NULL, edge.color='grey', edge.alpha=0.5, edge.width=0.5, ...)
{

    leave.label.direction <- match.arg(leave.label.direction)
	
	## how to convert a phylo object 'tree' into igraph object 'ig'
	if(0){
		## if internal node labels are duplicated, remove all (and add by as.igraph)
		if(any(duplicated(tree$node.label))) tree$node.label<-NULL
		ig <- as.igraph(tree, directed=T, use.labels=T)
	}
	
    if(class(ig) != "igraph"){
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
	if(!is.null(leave.label.wrap)){
		width <- as.integer(leave.label.wrap)
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
	##################
	if(!("tipid" %in% igraph::vertex_attr_names(ig))){
		# append 'tipid': NA for internal node, sequential order for tips
		ind <- match(V(ig)$name, dnet::dDAGtip(ig))
		V(ig)$tipid <- NA
		V(ig)$tipid[!is.na(ind)] <- 1:sum(!is.na(ind))
	}
	#####################################
	x <- y <- leaf <- label <- name <- NULL
	
	gp <- ggraph::ggraph(ig, layout='partition', circular=TRUE)
	gp <- gp + ggraph::geom_edge_diagonal(color=edge.color,alpha=edge.alpha,width=edge.width)
	gp <- gp + ggraph::geom_node_point(aes(filter=leaf),size=leave.size, color=edge.color,alpha=edge.alpha)
	
	if(leave.label.direction=='none'){
		df <- subset(gp$data, leaf==T)
		gp <- gp + ggrepel::geom_text_repel(data=df, aes(x=x, y=y, label=label), color=leave.label.color, size=leave.label.size, alpha=leave.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), ...)
		
	}else if(leave.label.direction=='leftright'){
		offset <- (range(gp$data$x)[2]-range(gp$data$x)[1]) * leave.label.offset
	
		root <- subset(gp$data, name==dnet::dDAGroot(ig))
		
		## left
		df1 <- subset(gp$data, leaf==T & x < root$x)
		df1$nudge_x <- -1 * (df1$x - min(gp$data$x)) - offset
		gp <- gp + ggrepel::geom_text_repel(data=df1, aes(x=x, y=y, label=label), color=leave.label.color, size=leave.label.size, alpha=leave.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), direction="y", hjust=0, nudge_x=df1$nudge_x)
		
		## right
		df2 <- subset(gp$data, leaf==T & x >= root$x)
		df2$nudge_x <- (max(gp$data$x)-df2$x) + offset
		gp <- gp + ggrepel::geom_text_repel(data=df2, aes(x=x, y=y, label=label), color=leave.label.color, size=leave.label.size, alpha=leave.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), direction="y", hjust=1, nudge_x=df2$nudge_x)
		
		if(is.null(limit.expansion)){
			gp <- gp + expand_limits(x=range(gp$data$x)*(1+2*leave.label.offset), y=range(gp$data$y))
		}else{
			gp <- gp + expand_limits(x=c(-limit.expansion, limit.expansion), y=c(-limit.expansion, limit.expansion))
		}
	
	}else if(leave.label.direction=='topbottom'){
		offset <- (range(gp$data$y)[2]-range(gp$data$y)[1]) * leave.label.offset
	
		root <- subset(gp$data, name==dnet::dDAGroot(ig))
		
		## left
		df1 <- subset(gp$data, leaf==T & y < root$y)
		df1$nudge_y <- -1 * (df1$y - min(gp$data$y)) - offset
		gp <- gp + ggrepel::geom_text_repel(data=df1, aes(x=x, y=y, label=label), color=leave.label.color, size=leave.label.size, alpha=leave.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), direction="x", hjust=0, nudge_y=df1$nudge_y, angle=90)
		
		## top
		df2 <- subset(gp$data, leaf==T & y >= root$y)
		df2$nudge_y <- (max(gp$data$y)-df2$y) + offset
		gp <- gp + ggrepel::geom_text_repel(data=df2, aes(x=x, y=y, label=label), color=leave.label.color, size=leave.label.size, alpha=leave.label.alpha, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), direction="x", hjust=1, nudge_y=df2$nudge_y, angle=90)
		
		if(is.null(limit.expansion)){
			gp <- gp + expand_limits(x=range(gp$data$x), y=range(gp$data$y)*(1+2*leave.label.offset))
		}else{
			gp <- gp + expand_limits(x=c(-limit.expansion, limit.expansion), y=c(-limit.expansion, limit.expansion))
		}
		
	}
	
	gp <- gp + ggraph::theme_graph()
	
	if(0){
		# order by tipid
		tipid <- NULL
		#gp$data <- gp$data %>% dplyr::arrange(tipid)
	}

	gp$ig <- ig
	
	invisible(gp)
}
