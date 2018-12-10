#' Function to visualise a graph with communities using hierarchical edge bundling
#'
#' \code{xHEB} is supposed to visualise a graph with communities using hierarchical edge bundling (HEB), an effective way to visualise connections between leaves of a hierarchical/tree graph (representing the community structure). The connections are curved and follow the tree structure (in a circular layout). It returns a ggplot object.
#'
#' @param g an object of class "igraph" with node attributes 'name' and 'community'
#' @param leave.label.size the text size of the leave labelings. By default, it is 3
#' @param leave.label.color the color of the leave labelings. By default, it is 'black'. If NULL, the label will be colored by the community
#' @param leave.label.wrap the wrap width of the leave labelings
#' @param leave.size the size of the leave nodes. By default, it is 3 if there is no node attribute 'size'
#' @param edge.tension the bundling strength of edges. 1 for very tight bundles, 0 for no bundle (straight lines). By defaults it is 0.8
#' @param edge.alpha the alpha of edges
#' @param edge.width the width of edges
#' @param edge.palette the palette defining edge color. It is correponding to the edge attribute 'weight' for the input graph (if any). By default, it is NULL: if the edge attribute 'weight' exists for the input graph, it will be 'RdPu' (RColorBrewer::display.brewer.all()); otherwise 'skyblue'
#' @return
#' a ggplot2 object
#' @export
#' @seealso \code{\link{xHEB}}
#' @include xHEB.r
#' @examples
#' # 1) generate a random bipartite graph
#' set.seed(123)
#' g <- sample_bipartite(50, 20, p=0.1)
#' V(g)$name <- paste0('node_',1:vcount(g))
#' 
#' \dontrun{
#' # 2) obtain its community
#' ig <- xBigraph(g)
#' 
#' # 3) HEB visualisation
#' library(ggraph)
#' E(ig)$weight <- runif(ecount(ig))
#' gp <- xHEB(ig)
#' }

xHEB <- function(g, leave.label.size=3, leave.label.color="black", leave.label.wrap=NULL, leave.size=NULL, edge.tension=0.8, edge.alpha=1, edge.width=0.5, edge.palette=NULL)
{
    
    if (class(g) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }else{
    	ig <- g
    }
	
	##################
	## text wrap
	if(!is.null(leave.label.wrap)){
		width <- as.integer(leave.label.wrap)
		res_list <- lapply(V(ig)$name, function(x){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=width)
			if(length(y)>1){
				paste(y,collapse='\n')
			}else{
				y
			}
		})
		V(ig)$name <- unlist(res_list)
	}
	##################
		
	if(!all(c("community","name") %in% igraph::vertex_attr_names(ig))){
		stop("The igraph object must have vertex attributes 'community'.\n")
		
	}else{
		
		## calculate within-community degree
		if(1){
			# degree within each community
			df_tmp <- igraph::get.data.frame(ig,what="vertices")
			ls_tmp <- split(x=df_tmp$name, f=df_tmp$community)
			ls_degree <- lapply(ls_tmp, function(x){
				ig_tmp <- igraph::induced.subgraph(ig, x)
				within.degrees <- igraph::degree(ig_tmp)
				data.frame(name=names(within.degrees), degree=within.degrees, stringsAsFactors=F)
			})
			df_withindegree <- do.call(rbind, ls_degree)
			ind <- match(V(ig)$name, df_withindegree$name)
			V(ig)$withindegree <- df_withindegree$degree[ind]
		}
		
		if(!all("size" %in% igraph::vertex_attr_names(ig))){
			if(0){
				V(ig)$size <- igraph::degree(ig)
			}else{
				V(ig)$size <- V(ig)$withindegree
			}
		}
		
		if(!all("order" %in% igraph::vertex_attr_names(ig))){
			V(ig)$order <- V(ig)$withindegree
		}
		
		message(sprintf("The igraph object has %d nodes and %d edges (%s) ...", vcount(ig), ecount(ig), as.character(Sys.time())), appendLF=T)
		
		community <- order <- NULL
		df_nodes <- igraph::get.data.frame(ig,what="vertices")
		#df_nodes$community <- factor(df_nodes$community, levels=unique(df_nodes$community))
		df_nodes <- df_nodes %>% dplyr::arrange(community,order)

		## a data structure contains both hierarchical components (i.e., parent-child relations between data items) and non-hierarchical components (representing additional relations between data items). Parent-child relations are called hierarchy relations, whereas additional, non-hierarchical relations are called adjacency relations
		
		##################
		## Construct hierarchy relations
		##################
		df1 <- data.frame(from='ROOT', to=sort(unique(df_nodes$community)), stringsAsFactors=F)
		df2 <- data.frame(from=df_nodes$community, to=df_nodes$name, stringsAsFactors=F)
		df_hierarchy <- rbind(df1, df2)
		
		# df_vertices
		df_vertices <- data.frame(name=unique(c(df_hierarchy$from, df_hierarchy$to)), stringsAsFactors=F)
		## df_vertices$community
		ind <- match(df_vertices$name, df_hierarchy$to)
		df_vertices$community <- df_hierarchy$from[ind]
		## df_vertices$size
		ind <- match(df_vertices$name, df_nodes$name)
		df_vertices$size <- df_nodes$size[ind]
		## orientation
		myleaves <- which(is.na( match(df_vertices$name, df_hierarchy$from) ))
		nleaves <- length(myleaves)
		df_vertices$id <- NA
		df_vertices$id[ myleaves ] <- seq(1:nleaves)
		df_vertices$angle <- 90 - 360 * df_vertices$id / nleaves
		df_vertices$hjust<-ifelse( df_vertices$angle < -90, 1, 0)
		df_vertices$angle<-ifelse(df_vertices$angle < -90, df_vertices$angle+180, df_vertices$angle)

		# Create ig_hierarchy
		ig_hierarchy <- igraph::graph_from_data_frame(df_hierarchy, vertices=df_vertices)
		
		##################
		## adjacency relations
		##################
		df_edges <- igraph::get.data.frame(ig)
		from <- match(df_edges$from, df_vertices$name)
		to <- match(df_edges$to, df_vertices$name)

		##################
		##################
		leaf <- x <- y <- NULL
		weight <- NULL
		size <- name <- angle <- hjust <- NULL
		
		gp <- ggraph::ggraph(ig_hierarchy, layout='dendrogram', circular=TRUE)
		## edge weight
		if(!is.null(df_edges$weight) & length(unique(df_edges$weight))>1){
			gp <- gp + ggraph::geom_conn_bundle(data=ggraph::get_con(from=from, to=to, weight=df_edges$weight), aes(colour=weight), alpha=edge.alpha, width=edge.width, tension=edge.tension)
			if(is.null(edge.palette)){
				#RColorBrewer::display.brewer.all()
				edge.palette <- "RdPu"
			}
			gp <- gp + ggraph::scale_edge_colour_distiller(palette=edge.palette, direction=1)
			
		}else{
			if(is.null(edge.palette)){
				edge.palette <- "skyblue"
			}
			gp <- gp + ggraph::geom_conn_bundle(data=ggraph::get_con(from=from,to=to), colour=edge.palette, alpha=edge.alpha, width=edge.width, tension=edge.tension)
		}
		## leaf size
		if(length(unique(V(ig_hierarchy)$size))>2){
			gp <- gp + ggraph::geom_node_point(aes(filter=leaf, x=x*1.05, y=y*1.05, colour=community, size=size), alpha=0.3) 
			#gp <- gp + scale_size_continuous(range=c(1,7))
		}else{
			if(is.null(leave.size)){
				leave.size <- 3
			}
			gp <- gp + ggraph::geom_node_point(aes(filter=leaf, x=x*1.05, y=y*1.05, colour=community), size=leave.size, alpha=0.3)
		}
		## leaf label
		if(is.null(leave.label.color)){
			gp <- gp + ggraph::geom_node_text(aes(x=x*1.1, y=y*1.1, filter=leaf, label=name, angle=angle, hjust=hjust, color=community), size=leave.label.size, alpha=1)
		}else{
			gp <- gp + ggraph::geom_node_text(aes(x=x*1.1, y=y*1.1, filter=leaf, label=name, angle=angle, hjust=hjust), color=leave.label.color, size=leave.label.size, alpha=1)
		}
		
		## theme
		gp <- gp + theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"cm")) + expand_limits(x=c(-1.2, 1.2), y=c(-1.2, 1.2))
		gp <- gp + theme_void() + coord_fixed()
	}
	
    invisible(gp)
}


