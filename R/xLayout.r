#' Function to define graph node coordinates according to igraph- or sna-style layout
#'
#' \code{xLayout} is supposed to define graph node coordinates according to igraph- or sna-style layout.
#'
#' @param g an object of class "igraph" (or "graphNEL") for a graph with such as a 'community' node attribute
#' @param layout a character specifying graph layout function. This character can be used to indicate igraph-style layout ("layout_nicely","layout_randomly","layout_in_circle","layout_on_sphere","layout_with_fr","layout_with_kk","layout_as_tree","layout_with_lgl","layout_with_graphopt","layout_with_sugiyama","layout_with_dh","layout_with_drl","layout_with_gem","layout_with_mds","layout_as_bipartite"), or sna-style layout ("gplot.layout.adj","gplot.layout.circle","gplot.layout.circrand","gplot.layout.eigen","gplot.layout.fruchtermanreingold","gplot.layout.geodist","gplot.layout.hall","gplot.layout.kamadakawai","gplot.layout.mds","gplot.layout.princoord","gplot.layout.random","gplot.layout.rmds","gplot.layout.segeo","gplot.layout.seham","gplot.layout.spring","gplot.layout.springrepulse","gplot.layout.target")
#' @param seed an integer specifying the seed
#' @return It returns an igraph object, appended by node attributes including "xcoord" for x-coordinates, "ycoord" for y-coordiates.
#' @export
#' @seealso \code{\link{xGGnetwork}}
#' @include xLayout.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' # load REACTOME
#' # restricted to Immune System ('R-HSA-168256') or Signal Transduction ('R-HSA-162582')
#' g <- xRDataLoader('ig.REACTOME', RData.location=RData.location)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="R-HSA-168256", mode="out")
#' nodeInduced <- V(g)[unique(unlist(neighs.out))]$name
#' ig <- igraph::induced.subgraph(g, vids=nodeInduced)
#' 
#' # compare Fruchterman and Reingold force-directed placement algorithm
#' ## based on igraph layout
#' ig1 <- xLayout(ig, layout="layout_with_fr")
#' gp1 <- xGGnetwork(ig1, node.xcoord="xcoord", node.ycoord="ycoord")
#' ## based on sna layout
#' ig2 <- xLayout(ig, layout="gplot.layout.fruchtermanreingold")
#' gp2 <- xGGnetwork(ig2, node.xcoord="xcoord", node.ycoord="ycoord")
#' 
#' # compare Kamada-Kawai force-directed placement algorithm
#' ## based on igraph layout
#' ig1 <- xLayout(ig, layout="layout_with_kk")
#' gp1 <- xGGnetwork(ig1, node.xcoord="xcoord", node.ycoord="ycoord")
#' ## based on sna layout
#' ig2 <- xLayout(ig, layout="gplot.layout.kamadakawai")
#' gp2 <- xGGnetwork(ig2, node.xcoord="xcoord", node.ycoord="ycoord")
#' }

xLayout <- function(g, layout=c("layout_nicely","layout_randomly","layout_in_circle","layout_on_sphere","layout_with_fr","layout_with_kk","layout_as_tree","layout_with_lgl","layout_with_graphopt","layout_with_sugiyama","layout_with_dh","layout_with_drl","layout_with_gem","layout_with_mds","layout_as_bipartite", "gplot.layout.adj","gplot.layout.circle","gplot.layout.circrand","gplot.layout.eigen","gplot.layout.fruchtermanreingold","gplot.layout.geodist","gplot.layout.hall","gplot.layout.kamadakawai","gplot.layout.mds","gplot.layout.princoord","gplot.layout.random","gplot.layout.rmds","gplot.layout.segeo","gplot.layout.seham","gplot.layout.spring","gplot.layout.springrepulse","gplot.layout.target"), seed=825)
{
    
    layout <- match.arg(layout)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if(class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
	
	glayout <- NULL
	
	if(grepl("gplot",layout)){
		## based on sna
		m <- as.matrix(xConverter(ig, from='igraph', to='dgCMatrix', verbose=F))
		set.seed(seed)
		eval(parse(text=paste0('glayout <- sna::',layout,'(m, NULL)')))
	}else{
		## based on igraph
		set.seed(seed)
		eval(parse(text=paste0('glayout <- ',layout,'(ig)')))
	}
	
	if(!is.null(glayout)){
		## scale into [-1,1]
		node.xcoord <- glayout[,1]
		node.ycoord <- glayout[,2]
		if(max(node.xcoord) != min(node.xcoord)){
			node.xcoord <- (node.xcoord - min(node.xcoord)) / (max(node.xcoord) - min(node.xcoord)) * 2 - 1
		}
		if(max(node.ycoord) != min(node.ycoord)){
			node.ycoord <- (node.ycoord - min(node.ycoord)) / (max(node.ycoord) - min(node.ycoord)) * 2 - 1
		}
		glayout <- cbind(node.xcoord, node.ycoord)
		
		## append 'xcoord' and 'ycoord' into ig
		V(ig)$xcoord <- glayout[,1]
		V(ig)$ycoord <- glayout[,2]
	}
    
    invisible(ig)
}


