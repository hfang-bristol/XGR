#' Function to implement Random Walk with Restart (RWR) on the input graph solved iteratively
#'
#' \code{xGraphRWR} is supposed to implement Random Walk with Restart (RWR) on the input graph solved iteratively. If the seeds (i.e. a set of starting nodes) are given, it intends to calculate the affinity score of all nodes in the graph to the seeds. If the seeds are not given, it will pre-compute affinity matrix for nodes in the input graph with respect to each starting node (as a seed) by looping over every node in the graph. Parallel computing is also supported.
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param setSeeds an input matrix used to define sets of starting seeds. One column corresponds to one set of seeds that a walker starts with. The input matrix must have row names, coming from node names of input graph, i.e. V(g)$name, since there is a mapping operation. The non-zero entries mean that the corresonding rows (i.e. the gene/row names) are used as the seeds, and non-zero values can be viewed as how to weight the relative importance of seeds. By default, this option sets to "NULL", suggesting each node in the graph will be used as a set of the seed to pre-compute affinity matrix for the input graph. This default does not scale for large input graphs since it will loop over every node in the graph; however, the pre-computed affinity matrix can be extensively reused for obtaining affinity scores between any combinations of nodes/seeds, allows for some flexibility in the downstream use, in particular when sampling a large number of random node combinations for statistical testing
#' @param restart the restart probability used for RWR. The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param max.step an integer specifying the maximum number of steps that random walk performs. By default, it is 50
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a sparse matrix, called 'PTmatrix':
#' \itemize{
#'  \item{When the seeds are NOT given: a pre-computated affinity matrix with the dimension of n X n, where n is the number of nodes in the input graph. Columns stand for starting nodes walking from, and rows for ending nodes walking to. Therefore, a column for a starting node represents a steady-state affinity vector that the starting node will visit all the ending nodes in the graph}
#'  \item{When the seeds are given: an affinity matrix with the dimension of n X nset, where n is the number of nodes in the input graph, and nset for the number of the sets of seeds (i.e. the number of columns in setSeeds). Each column stands for the steady probability vector, storing the affinity score of all nodes in the graph to the starting nodes/seeds. This steady probability vector can be viewed as the "influential impact" over the graph imposed by the starting nodes/seeds.}
#' }
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{xGraphRWR}}
#' @include xGraphRWR.r
#' @examples
#' g <- make_graph("Zachary")
#' gp <- xGGnetwork(g,node.label='name',node.label.size=3, node.label.force=0.01)
#'
#' \dontrun{
#' # obtain the pre-computated affinity matrix
#' PTmatrix <- xGraphRWR(g, normalise="laplacian", restart=0.75)
#' # visualise affinity matrix
#' visHeatmapAdv(as.matrix(PTmatrix), Rowv=FALSE, Colv=FALSE, colormap="yr", KeyValueName="Affinity")
#' }

xGraphRWR <- function(g, normalise=c("laplacian","row","column","none"), setSeeds=NULL, restart=0.75, max.step=50, verbose=TRUE)
{
    
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    normalise <- match.arg(normalise)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }

    if(verbose){
        now <- Sys.time()
        message(sprintf("First, get the adjacency matrix of the input graph (%s) ...", as.character(now)), appendLF=TRUE)
    }
    
    if ("weight" %in% igraph::list.edge.attributes(ig)){
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
        now <- Sys.time()
        message(sprintf("Then, normalise the adjacency matrix using %s normalisation (%s) ...", normalise, as.character(now)), appendLF=TRUE)
    }
    
    A <- adjM!=0
    if(normalise == "row"){
        D <- Matrix::Diagonal(x=(Matrix::rowSums(A))^(-1))
        nadjM <- adjM %*% D
    }else if(normalise == "column"){
        D <- Matrix::Diagonal(x=(Matrix::colSums(A))^(-1))
        nadjM <- D %*% adjM
    }else if(normalise == "laplacian"){
        D <- Matrix::Diagonal(x=(Matrix::colSums(A))^(-0.5))
        nadjM <- D %*% adjM %*% D
    }else{
        nadjM <- adjM
    }

    
    ####################################################
    # A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=FALSE){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d seed sets (%s)", i, B, as.character(Sys.time())), appendLF=TRUE)
            }
        }
    }
    
    ## A function to make sure the sum of elements in each steady probability vector is one
    sum2one <- function(PTmatrix){
        col_sum <- apply(PTmatrix, 2, sum)
        col_sum_matrix <- matrix(rep(col_sum, nrow(PTmatrix)), ncol=ncol(PTmatrix), nrow=nrow(PTmatrix), byrow=TRUE)
        res <- as.matrix(PTmatrix)/col_sum_matrix
        res[is.na(res)] <- 0
        return(res)
    }
    
    ####################################################
    
    ################## RWR    
    if(is.null(setSeeds)){
        P0matrix <- Matrix::Matrix(diag(vcount(ig)), sparse=TRUE)
        rownames(P0matrix) <- V(ig)$name
        colnames(P0matrix) <- V(ig)$name
        
    }else{
        ## check input data
        if(is.matrix(setSeeds) | is.data.frame(setSeeds)){
            data <- as.matrix(setSeeds)
        }else if(is.vector(setSeeds)){
            data <- as.matrix(setSeeds, ncol=1)
        }

        if(is.null(rownames(data))) {
            stop("The function must require the row names of the input setSeeds.\n")
        }else if(any(is.na(rownames(data)))){
            warning("setSeeds with NA as row names will be removed")
            data <- data[!is.na(rownames(data)),]
        }
        cnames <- colnames(data)
        if(is.null(cnames)){
            cnames <- seq(1,ncol(data))
        }
    
        ## check mapping between input data and graph
        ind <- match(rownames(data), V(ig)$name)
        nodes_mapped <- V(ig)$name[ind[!is.na(ind)]]
        if(length(nodes_mapped)!=vcount(ig)){
            warning("The row names of input setSeeds do not contain all those in the input graph.\n")
        }
        P0matrix <- matrix(0,nrow=nrow(nadjM),ncol=ncol(data))
        P0matrix[ind[!is.na(ind)],] <- as.matrix(data[!is.na(ind),])
        
        ## make sure the sum of elements in each steady probability vector is one
        P0matrix <- sum2one(P0matrix)
        
        ## convert to sparse matrix
        P0matrix <- Matrix::Matrix(P0matrix, sparse=TRUE)
    }
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Third, RWR of %d sets of seeds using %1.1e restart probability (with maximum step) (%s) ...", ncol(P0matrix), restart, as.character(now)), appendLF=TRUE)
    }
    
    if(restart==1){
        ## just seeds themselves
        PTmatrix <- P0matrix
    }else{
    
		## stopping critera
		stop_delta <- 1e-6   # L1 norm of successive estimates 'PT' below the threshold 'stop_delta'
		stop_step <- min(50, as.integer(max.step)) # maximum steps of iterations

		ls_res <- lapply(1:ncol(P0matrix), function(j){
			progress_indicate(j, ncol(P0matrix), 10, flag=TRUE)
			
           	P0 <- P0matrix[,j]
        
            ## Initializing variables
            step <- 0
            delta <- 1
        
            PT <- P0
            ## Iterative update till convergence (delta<=1e-10)
            while (delta>stop_delta && step<=stop_step){
                PX <- (1-restart) * nadjM %*% PT + restart * P0

                # p-norm of v: sum((abs(v).p)^(1/p))
                delta <- sum(abs(PX-PT))
    
                PT <- PX
                step <- step+1
            }
            PT
		})
		PTmatrix <- do.call(cbind, ls_res)
    }
    
    ## make sure the sum of elements in each steady probability vector is one
    if(verbose){
        now <- Sys.time()
        message(sprintf("Fourth, rescale steady probability vector (%s) ...", as.character(now)), appendLF=TRUE)
    }
    PTmatrix <- sum2one(PTmatrix) # input/output: full matrix
    #PTmatrix <- Matrix::Matrix(PTmatrix, sparse=TRUE)
    rownames(PTmatrix) <- rownames(P0matrix)
    colnames(PTmatrix) <- colnames(P0matrix)
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)

    invisible(PTmatrix)
}


