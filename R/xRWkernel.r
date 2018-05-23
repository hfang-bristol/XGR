#' Function to calculate random walk kernel on the input graph
#'
#' \code{xRWkernel} is supposed to calculate a weighted random walk kernel (a predefined number of steps) for estimating pairwise connectivity between nodes.
#'
#' @param g an object of class "igraph" or "graphNEL". It will be a weighted graph if having an edge attribute 'weight'. The edge directions are ignored for directed graphs
#' @param steps an integer specifying the number of steps that random walk performs. By default, it is 4
#' @param a an integer specifying the chance of remaining at the same vertex. By default, it is 2, the higher the higher change
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a connectivity matrix for pairwise affinity between nodes
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{xRWkernel}}
#' @include xRWkernel.r
#' @examples
#' # 1) generate a random graph according to the ER model
#' set.seed(123)
#' g <- erdos.renyi.game(10, 3/10)
#' V(g)$name <- paste0('n',1:vcount(g))
#'
#' \dontrun{
#' # 2) obtain the pre-computated affinity matrix
#' mat <- xRWkernel(g)
#' # visualise affinity matrix
#' visHeatmapAdv(as.matrix(mat), colormap="wyr", KeyValueName="Affinity")
#' }

xRWkernel <- function(g, steps=4, a=2, verbose=TRUE)
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
	
	if(igraph::is_directed(ig)){
		ig <- igraph::as.undirected(ig, mode="collapse", edge.attr.comb="max")
	}
	
    if(verbose){
        now <- Sys.time()
        message(sprintf("First, get the adjacency matrix of the input graph (%s) ...", as.character(now)), appendLF=TRUE)
    }
    
    if ("weight" %in% list.edge.attributes(ig)){
        adjM <- get.adjacency(ig, type="both", attr="weight", edges=FALSE, names=TRUE, sparse=getIgraphOpt("sparsematrices"))
        if(verbose){
            message(sprintf("\tNotes: using weighted graph!"), appendLF=TRUE)
        }
    }else{
        adjM <- get.adjacency(ig, type="both", attr=NULL, edges=FALSE, names=TRUE, sparse=getIgraphOpt("sparsematrices"))
        if(verbose){
            message(sprintf("\tNotes: using unweighted graph!"), appendLF=TRUE)
        }
    }
    
    if(verbose){
        message(sprintf("Then, laplacian normalisation of the adjacency matrix (%s) ...", as.character(Sys.time())), appendLF=TRUE)
    }
    
    A <- adjM!=0
    ## library(Matrix) 
    ## vignette("Intro2Matrix")
    D <- Matrix::Diagonal(x=(Matrix::colSums(A))^(-0.5))
    nadjM <- D %*% adjM %*% D
    #nadjM <- as.matrix(nadjM)
    
    steps <- as.integer(steps)
    if(verbose){
        message(sprintf("Last, %d-step random walk kernel (%s) ...", steps, as.character(Sys.time())), appendLF=TRUE)
    }
    
    if(verbose){
        message(sprintf("\tstep 1 (%s) ...", as.character(Sys.time())), appendLF=TRUE)
    }
	## one-step random walk kernel
	I <- Matrix::Matrix(diag(x=a-1,nrow=vcount(g)), sparse=TRUE)
	RW <- I + nadjM
    res <- RW
    
    ## p-step random walk kernel
    if(steps >=2){
		for (i in 2:steps){
			if(verbose){
				message(sprintf("\tstep %d (%s) ...", i, as.character(Sys.time())), appendLF=TRUE)
			}
		  	res <- res %*% RW
		}
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)

    invisible(as.matrix(res))
}


