#' Function to add coordinates into a graph according to a node attribute
#'
#' \code{xAddCoords} is supposed to add coordinates into a graph according to a node attribute such as community or comp.
#'
#' @param g an object of class "igraph" (or "graphNEL") for a graph with such as a 'community' node attribute
#' @param node.attr a charatter specifying a node attribute. If NULL or no match, it returns NULL
#' @param glayout a graph layout function. This function can be one of "layout_nicely" (previously "layout.auto"), "layout_randomly" (previously "layout.random"), "layout_in_circle" (previously "layout.circle"), "layout_on_sphere" (previously "layout.sphere"), "layout_with_fr" (previously "layout.fruchterman.reingold"), "layout_with_kk" (previously "layout.kamada.kawai"), "layout_as_tree" (previously "layout.reingold.tilford"), "layout_with_lgl" (previously "layout.lgl"), "layout_with_graphopt" (previously "layout.graphopt"), "layout_with_sugiyama" (previously "layout.sugiyama"), "layout_with_dh" (previously "layout.davidson.harel"), "layout_with_drl" (previously "layout.drl"), "layout_with_gem" (previously "layout.gem"), "layout_with_mds", and "layout_as_bipartite". A full explanation of these layouts can be found in \url{http://igraph.org/r/doc/layout_nicely.html}
#' @param edge.color.alternative two alternative colors for edges within the community (grey70 by default) and edges between communities (grey95 by default)
#' @param seed an integer specifying the seed
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns an igraph object, appended by node attributes including "xcoord" for x-coordinates, "ycoord" for y-coordiates, and by edge attributes including "color" for between-community edges ('grey95') and within-community edges ('grey70').
#' @export
#' @seealso \code{\link{xGGnetwork}}
#' @include xAddCoords.r
#' @examples
#' # 1) generate a random bipartite graph
#' set.seed(825)
#' g <- sample_bipartite(100, 50, p=0.1)
#' V(g)$name <- V(g)
#' 
#' \dontrun{
#' # 2) obtain and append the community
#' cs <- igraph::cluster_louvain(g)
#' set.seed(825); cs <- igraph::cluster_spinglass(g)
#' V(g)$community <- cs$membership
#' ig <- xAddCoords(g, node.attr="community", edge.color.alternative=c("grey50","grey95"))
#' if(class(V(ig)$community)=='character') V(ig)$community <- as.factor(V(ig)$community)
#' gp <- xGGnetwork(ig, node.label='name', node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0, node.label.force=0.002, node.xcoord='xcoord', node.ycoord='ycoord', node.color='community', node.color.title='Community', colormap='jet.both', ncolors=64, zlim=NULL, edge.color="color",edge.color.alpha=0.5,edge.curve=0,edge.arrow.gap=0)
#' 
#' ## make it discrete for the colorbar
#' gp + scale_colour_gradientn(colors=xColormap('jet')(64),breaks=seq(1,9)) + guides(color=guide_legend(title="Community"))
#' 
#' ## add vertex hull for each community
#' df <- gp$data_nodes
#' ls_res <- lapply(split(x=df,f=df$community), function(z) z[chull(z$x,z$y),])
#' data <- do.call(rbind, ls_res)
#' gp + geom_polygon(data=data, aes(x=x,y=y,group=community), alpha=0.1)
#' gp + geom_polygon(data=data, aes(x=x,y=y,group=community,fill=community), alpha=0.1) + scale_fill_gradientn(colors=xColormap('jet.both')(64)) + guides(fill="none")
#' }

xAddCoords <- function(g, node.attr=NULL, glayout=layout_with_kk, edge.color.alternative=c("grey70","grey95"), seed=825, verbose=TRUE)
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
	
	type <- community <- contribution <- name <- NULL
	df_nodes <- igraph::get.data.frame(ig,what="vertices")
	df_nodes$node.attr <- df_nodes[,node.attr]
	
	if(all(c("type",node.attr,"contribution") %in% colnames(df_nodes))){
		## sorted by type, community (1,2,3,...), -contribution
		df_nodes <- df_nodes %>% dplyr::arrange(type,node.attr,-contribution)
		
	}else{
		## sorted by community (1,2,3,...)
		df_nodes <- df_nodes %>% dplyr::arrange(node.attr)
	}
	
	#############################
	## append node attributes: xcoord, ycoord
	#############################
	if(1){
		ls_tmp <- split(x=df_nodes$name, f=df_nodes$node.attr)
		ls_ig <- lapply(ls_tmp, function(x){
			g <- dnet::dNetInduce(ig, nodes_query=x, knn=0, largest.comp=TRUE)
		})
		
		## collectively
		if(1){
			layouts <- lapply(ls_ig, glayout)
		}else{
			layouts <- lapply(ls_ig, function(g) {
				igraph::layout_as_bipartite(g, types=as.logical(V(g)$type=='xnode'))
				igraph::layout_in_circle(g)
				igraph::layout_with_kk(g)
			})
		}
		set.seed(seed)
		glayout <- igraph::merge_coords(ls_ig, layouts)
		node.xcoord <- glayout[,1]
		node.ycoord <- glayout[,2]
		## scale into [-1,1]
		if(max(node.xcoord) != min(node.xcoord)){
			node.xcoord <- (node.xcoord - min(node.xcoord)) / (max(node.xcoord) - min(node.xcoord)) * 2 - 1
		}
		if(max(node.ycoord) != min(node.ycoord)){
			node.ycoord <- (node.ycoord - min(node.ycoord)) / (max(node.ycoord) - min(node.ycoord)) * 2 - 1
		}
		glayout <- cbind(node.xcoord, node.ycoord)
		## glayout sorted by ig
		g_tmp <- igraph::disjoint_union(ls_ig)
		ind <- match(V(ig)$name, V(g_tmp)$name)
		glayout <- glayout[ind,]
		V(ig)$xcoord <- glayout[,1]
		V(ig)$ycoord <- glayout[,2]

		#############################
		## append edge attribute: color
		#############################
		## define edge.color
		### within and between: grey70 and grey90
		m <- df_nodes$node.attr
		names(m) <- df_nodes$name
		### edges
		df_edges <- igraph::get.data.frame(ig,what="edges")
		### decide on edges within and between community
		ind <- match(df_edges$from, names(m))
		m1 <- m[ind]
		ind <- match(df_edges$to, names(m))
		m2 <- m[ind]
		res <- m1!=m2
		if(!is.null(names(m1))){
			names(res) <- paste(names(m1), names(m2), sep="|")
		}
		
		if(length(edge.color.alternative)==2){
			edge.color.alternative <- edge.color.alternative
		}else{
			edge.color.alternative <- c("grey70","grey95")
		}
		
		edge.color <- edge.color.alternative[res+1]
		E(ig)$color <- edge.color

		#gp <- xGGnetwork(ig, node.label='name', node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0, node.label.force=0.002, node.shape='type', node.shape.title='Type', node.xcoord='xcoord', node.ycoord='ycoord', node.color='community', node.color.title='Community', colormap='jet.both', ncolors=64, zlim=NULL, node.size='contribution', node.size.range=c(1,3), node.size.title='Contribution', slim=c(0,0.1), edge.color="color",edge.color.alpha=0.5,edge.curve=0,edge.arrow.gap=0, title='')
		
    }
    
    invisible(ig)
}


