#' Function to draw a tree-like circular plot
#'
#' \code{xCtree} is supposed to draw a tree-like circular plot (dendrogram circular layout), with tips labelled (outwards or inwards). The tree is provided as an object of class "igraph".
#'
#' @param ig an object of class "igraph" with node attribute 'name'. It could be a 'phylo' object converted to. Note: the leave labels would be the node attribute 'name' unless the node attribute 'label' is explicitely provided
#' @param leave.label.orientation the leave label orientation. It can be "outwards" and "inwards"
#' @param leave.label.size the text size of the leave labelings. By default, it is 2
#' @param leave.label.color the color of the leave labelings
#' @param leave.label.alpha the alpha of the leave labelings
#' @param leave.label.wrap the wrap width of the leave labelings
#' @param leave.label.expansion the x- and y-expansion of the leave labelings. The value of 1 for the exact location of the leave, and the outwards (>1; by default 1.05 if NULL) and the inwards (<1; by default 0.98 if NULL)
#' @param leave.size the size of the leave nodes. By default, it is 0
#' @param limit.expansion the x- and y-limit expansion. By default, it is 1.1. Beware the orignial limit is [-1,1]
#' @param edge.color the color of edges
#' @param edge.alpha the alpha of edges
#' @param edge.width the width of edges
#' @return 
#' a ggplot2 object appended with 'ig' and 'data' which should contain columns 'x','y', 'leaf' (T/F), 'name' (the same as V(ig)$name), 'tipid' (tip id), 'label' (if not given in ig, a 'name' varient), 'angle' and 'hjust' (assist in leave label orientation).
#' @note none
#' @export
#' @seealso \code{\link{xCtree}}
#' @include xCtree.r
#' @examples
#' \dontrun{
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' 
#' AA.path <- xRDataLoader("AA.path", RData.location=RData.location)
#' # consensus tree
#' ig <- AA.path$consensus$ig
#' 
#' # outwards
#' gp <- xCtree(ig, leave.label.orientation="outwards", leave.label.wrap=50, limit.expansion=1.5, leave.size=2)
#' head(gp$data %>% dplyr::arrange(tipid))
#' 
#' # inwards
#' gp <- xCtree(ig, leave.label.orientation="inwards", leave.label.wrap=30)
#' 
#' # obtain 'xcoord' and 'ycoord'
#' gp <- ggraph::ggraph(ig, layout='dendrogram', circular=TRUE)
#' data <- gp$data %>% dplyr::arrange(ggraph.orig_index)
#' V(ig)$xcoord <- data[,'x']
#' V(ig)$ycoord <- data[,'y']
#' }

xCtree <- function(ig, leave.label.orientation=c('outwards','inwards'), leave.label.size=2, leave.label.color="darkblue", leave.label.alpha=0.7, leave.label.wrap=NULL, leave.label.expansion=NULL, leave.size=0, limit.expansion=1.1, edge.color='grey', edge.alpha=0.5, edge.width=0.5)
{

    leave.label.orientation <- match.arg(leave.label.orientation)
	
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
	## orientation: angle and hjust			
	V(ig)$angle <- 90 - 360 * V(ig)$tipid / sum(!is.na(V(ig)$tipid))
	if(leave.label.orientation=='outwards'){
		V(ig)$hjust <- ifelse(V(ig)$angle < -90, 1, 0)
		if(is.null(leave.label.expansion)){
			leave.label.expansion <- 1.05
		}
	}else if(leave.label.orientation=='inwards'){
		V(ig)$hjust <- ifelse(V(ig)$angle < -90, 0, 1)
		if(is.null(leave.label.expansion)){
			leave.label.expansion <- 0.98
		}
	}
	V(ig)$angle <- ifelse(V(ig)$angle < -90, V(ig)$angle+180, V(ig)$angle)
	#####################################
	x <- y <- leaf <- label <- angle <- hjust <- NULL
	
	gp <- ggraph::ggraph(ig, layout='dendrogram', circular=TRUE)
	gp <- gp + ggraph::geom_edge_diagonal(color=edge.color,alpha=edge.alpha,width=edge.width)
	gp <- gp + ggraph::geom_node_point(aes(filter=leaf),size=leave.size, color=edge.color,alpha=edge.alpha)
	gp <- gp + ggraph::geom_node_text(aes(x=x*leave.label.expansion, y=y*leave.label.expansion, filter=leaf, label=label, angle=angle, hjust=hjust),show.legend=F, color=leave.label.color, size=leave.label.size, alpha=leave.label.alpha) + expand_limits(x=c(-limit.expansion, limit.expansion), y=c(-limit.expansion, limit.expansion))
	
	gp <- gp + coord_fixed() + theme(legend.position="bottom") + ggraph::theme_graph(base_family="sans")
	
	if(0){
		# order by tipid
		tipid <- NULL
		#gp$data <- gp$data %>% dplyr::arrange(tipid)
	}

	if(1){
		ggraph.orig_index <- NULL
		# append 'xcoord' and 'ycoord'
		data <- gp$data %>% dplyr::arrange(ggraph.orig_index)
		V(ig)$xcoord <- data[,'x']
		V(ig)$ycoord <- data[,'y']
	}

	gp$ig <- ig
	
	invisible(gp)
}
