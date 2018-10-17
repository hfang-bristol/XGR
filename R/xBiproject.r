#' Function to obtain a projected graph from a bipartitle graph
#'
#' \code{xBiproject} is supposed to obtain a projected graph from a bipartitle graph.
#'
#' @param g an object of class "igraph" (or "graphNEL") for a bipartitel graph with a 'type' node attribute
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns an igraph object.
#' @note The input graph will has an equal weight if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{xBiproject}}
#' @include xBiproject.r
#' @examples
#' # 1) generate a random bipartite graph
#' set.seed(123)
#' g <- sample_bipartite(100, 50, p=0.05)
#' V(g)$name <- V(g)
#' 
#' \dontrun{
#' # 2) obtain its projected graph
#' ig <- xBiproject(g)
#' 
#' # 3) estimate pairwise affinity between nodes
#' mat <- xRWkernel(ig)
#' }

xBiproject <- function(g, verbose=TRUE)
{
    
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
	
	if(is.null(V(ig)$type)){
		stop("The igraph object must have a 'type' vertex attribute.\n")
	}
	
	vec_type <- sort(table(V(ig)$type))
	if(verbose){
		message(sprintf("The igraph object has %d nodes (%d xnode '%s' type and %d ynode '%s' type) and %d edges (%s) ...", vcount(ig), vec_type[1],names(vec_type)[1], vec_type[2],names(vec_type)[2], ecount(ig), as.character(Sys.time())), appendLF=T)
	}
	
	if(0){
	## only keep the largest component
	ig <- dnet::dNetInduce(ig, nodes_query=V(ig)$name, knn=0, largest.comp=TRUE)
	if(verbose){
		message(sprintf("\tThe largest component has %d nodes and %d edges", vcount(ig), ecount(ig)), appendLF=T)
	}
	}
	
	if(vcount(ig)==0){
		return(NULL)
	}
	
	## append the weight edge attribute if not found
	if(is.null(E(ig)$weight)){
		if(verbose){
			message(sprintf("\tWithout a 'weight' vertex attributes"), appendLF=T)
		}
		E(ig)$weight <- rep(1, ecount(ig))
	}
	
	df_edges <- get.data.frame(ig, what="edges")[,c("from","to","weight")]
	
	## ynode versus xnode
	if(length(unique(df_edges$from)) < length(unique(df_edges$to))){
		elist <- df_edges[, c("to","from","weight")]
	}else{
		elist <- df_edges
	}
	colnames(elist) <- c("ynode","xnode","weight")
	
	## Sparse matrix: ynode versus xnode
	### ynode indices
	ynodes <- as.integer(factor(elist[,1]))
	ynode.names <- levels(factor(elist[,1]))
	### xnode indices (Genes)
	xnodes <- as.integer(factor(elist[,2]))
	xnode.names <- levels(factor(elist[,2]))
	### weights
	weights <- elist[,3]
	### ynodes versus xnodes
	sM <- Matrix::sparseMatrix(i=ynodes, j=xnodes, x=weights, dims=c(length(unique(ynodes)),length(unique(xnodes))))
	
	## Projected into xnode space with projected adjacency matrix
	xM <- Matrix::t(sM) %*% sM
	colnames(xM) <- rownames(xM) <- xnode.names
	
	## ig_x: igraph from projected adjacency matrix on xnodes
	ig_x <- igraph::graph.adjacency(xM, mode="undirected", weighted=TRUE, diag=FALSE)
	### remove loops and multiple edges
	#ig_x <- dnet::dNetInduce(ig_x, nodes_query=V(ig_x)$name, knn=0, largest.comp=TRUE)

	if(verbose){
		message(sprintf("Projected onto '%s' node type creates a graph of %d nodes and %d edges (%s) ...", names(vec_type)[1], vcount(ig_x), ecount(ig_x), as.character(Sys.time())), appendLF=T)
	}

 ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)

    invisible(ig_x)
}


