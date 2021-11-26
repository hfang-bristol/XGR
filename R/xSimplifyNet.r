#' Function to simplify networks from an igraph object
#'
#' \code{xSimplifyNet} is supposed to simplify networks from an igraph object by keeping root-tip shortest paths only.
#'
#' @param g an "igraph" object
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' an object of class "igraph"
#' @note none
#' @export
#' @seealso \code{\link{xSimplifyNet}}
#' @include xSimplifyNet.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' g <- xRDataLoader(RData.customised='ig.DO', RData.location=RData.location)
#' ig <- xSimplifyNet(g)
#' }

xSimplifyNet <- function(g, verbose=TRUE)
{

    if(is(g,"igraph")){
    	ig <- g
    
		if(verbose){
			now <- Sys.time()
			message(sprintf("The input graph has %d nodes and %d edges (%s) ...", vcount(ig),ecount(ig),as.character(now)), appendLF=TRUE)
		}
		
	}else{
		stop("The function must apply to the 'igraph' object.\n")
	}
	
	if(!igraph::is_directed(ig)){
		stop("The function must apply to the direct graph.\n")
	}
	root <- dnet::dDAGroot(ig)
	if(verbose){
		message(sprintf("\tthe number of roots: %d", length(root)), appendLF=TRUE)
	}
	tips <- dnet::dDAGtip(ig)
	if(verbose){
		message(sprintf("\tthe number of tips: %d", length(tips)), appendLF=TRUE)
	}
	epaths <- igraph::get.shortest.paths(ig, from=root, to=tips, output="epath")
	subg <- igraph::subgraph.edges(ig, eids=unlist(epaths$epath), delete.vertices=TRUE)
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("The simplified graph has %d nodes and %d edges (%s) ...", vcount(subg),ecount(subg),as.character(now)), appendLF=TRUE)
	}
    
    invisible(subg)
}
