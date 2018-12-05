#' Function to calculate node similarity based on the network structure
#'
#' \code{xGraphSim} is supposed to calculate node similarity based on the network structure.
#'
#' @param g an object of class "igraph" (or "graphNEL") for a graph
#' @param measure the similarity metrics used to measure similarity between nodes. It can be node-dependent similarity, path-dependent similarity, and random walk-based similarity
#' @param type the type defining nodes to be calcuated. It can be 'full' for all possible pairs of nodes (used for link prediction) and 'edge' for paris of nodes in edges (used for link recommendation)
#' @param edge.fast logical to indicate whether a vectorised fast computation is used. By default, it sets to true. It is always advisable to use this vectorised fast computation; though the conventional computation is just used for script understanding
#' @param measure.para a list of measure-specific parameters.
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a symmetric sparse matrix involving all network nodes.
#' @note The similarity measures are detailed as follows.\cr
#' 1. Node-dependent similarity
#' \itemize{
#' \item{"Common Neighbors (CN)": the number of common neighbors.}
#' \item{"Salton Index (SI)": the cosine similarity, that is, the number of common neighbors normalised by square root of degree product of two nodes.}
#' \item{"Jaccard Index (JI)": the number of common neighbors normalised by the union of neighbors.}
#' \item{"Dice Index (DI)": mainly for ecological community data, that is, the number of common neighbors normalised by degree average of two nodes.}
#' \item{"Hub Promoted Index (HPI)": the topological overlap, that is, the number of common neighbors normalised by minimum degree of two nodes (with hubs receiving higher score).}
#' \item{"Hub Depressed Index (HPI)": the number of common neighbors normalised by maximum degree of two nodes (with hubs receiving lower score).}
#' \item{"Leicht-Holme-Newman Index (LHN)": inversely proportional to the expected number of common neighbors in the configuration model, that is, the number of common neighbors normalised by degree product of two nodes.}
#' \item{"Adamic-Adar Index (AA)": the number of common neighbors with a logarithmic weighting, that is, weighted by 1/log(degree).}
#' \item{"Resource Allocation Index (RA)": the number of common neighbors with an inverse weighting, that is, weighted by 1/degree.}
#' }
#' 2. Path-dependent similarity
#' \itemize{
#' \item{"Local Path Index (LP)": the number of two-step paths and three-step paths (weighted by the eps).}
#' \item{"Katz Index (KI)": sums over the collection of paths exponentially damped by length to give the shorter paths more weights.}
#' \item{"Leicht-Holme-Newman Index global similarity (LHNg)": two nodes are similar if their immediate neighbors are similar by themselves.}
#' \item{"Shortest Path Index (SP)": the length of the shortest paths.}
#' }
#' 3. Random walk-based similarity
#' \itemize{
#' \item{"Average Commute Time (ACT)": the average number of steps required by a random walker starting from a node to reach another node.}
#' \item{"Cosine based on Laplacian matrix (CL)": the cosine of the angle between columns of the pseudoinverse of the Laplacian matrix.}
#' \item{"Matrix Forest Index (MFI)": the ratio of the number of spanning rooted forests belonging to the same tree to all spanning rooted forests of the network.}
#' \item{"Random Walk with Restart (RWR)": the probability at the steady state, solved analytically (not iteratively).}
#' \item{"Short Random Walk (SRW)": p-step random walk, solved analytically (not iteratively).}
#' }
#' @export
#' @seealso \code{\link{xGraphSim}}
#' @include xGraphSim.r
#' @examples
#' g <- make_graph("Zachary")
#' gp <- xGGnetwork(g,node.label='name',node.label.size=3, node.label.force=0.01)
#' 
#' \dontrun{
#' # Common Neighbors (CN)
#' res <- xGraphSim(g, measure="CN")
#' # Resource Allocation Index (RA)
#' res <- xGraphSim(g, measure="RA")
#' # Local Path Index (LP)
#' res <- xGraphSim(g, measure="LP", measure.para=list(LP.eps=0.01))
#' # Katz Index (KI)
#' res <- xGraphSim(g, measure="KI", measure.para=list(KI.beta=0.001))
#' # Random Walk with Restart (RWR)
#' res <- xGraphSim(g, measure="RWR", measure.para=list(RWR.r=0.5))
#' # Short Random Walk (SRW)
#' res <- xGraphSim(g, measure="SRW", measure.para=list(SRW.p=4))
#' }

xGraphSim <- function(g, measure=c("CN","SI","JI","DI","HPI","HDI","LHN","AA","RA","LP","KI","LHNg","SP","ACT","CL","MFI","RWR","SRW"), type=c('full','edge'), edge.fast=F, measure.para=list(LP.eps=0.01,KI.beta=0.001,LHNg.theta=0.5,RWR.r=0.5,SRW.p=4), verbose=TRUE)
{
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    measure <- match.arg(measure)
    type <- match.arg(type)

    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    if(type=='edge' & !edge.fast){
		if(is.null(V(ig)$iid)){
			V(ig)$iid <- 1:vcount(ig)
		}
	
		vec_degree <- igraph::degree(ig)
		names(vec_degree) <- V(ig)$iid
	
		# a node: neighbors
		ls_neighbors <- lapply(V(ig), function(v){
			res <- igraph::neighbors(ig, v)
		})
		names(ls_neighbors) <- V(ig)$iid
	
		df_edges <- get.data.frame(ig, what='edge')
		## add degree
		ind <- match(df_edges$from, names(vec_degree))
		df_edges$degree_from <- vec_degree[ind]
		ind <- match(df_edges$to, names(vec_degree))
		df_edges$degree_to <- vec_degree[ind]
	
		# two-end nodes: common neighbors
		ls_common <- lapply(1:nrow(df_edges), function(i){
			from <- ls_neighbors[[df_edges$from[i]]]
			to <- ls_neighbors[[df_edges$to[i]]]
			intersection(from, to)
		})
    }
	
	n <- vcount(ig)
	m <- igraph::ecount(ig)
	
	###############
	if(measure=='CN'){
		# Common Neighbors (CN)
		if(type=='full' | (type=='edge' & edge.fast)){
			score <- igraph::cocitation(ig)
			score <- Matrix::Matrix(score, sparse=TRUE)
			if(type=='edge'){
				score  <- igraph::get.adjacency(ig) * score
			}
		}else{
			df_edges$sim <- sapply(1:length(ls_common), function(i){
				length(ls_common[[i]])
			})
			mat_triu <- as.matrix(xSparseMatrix(df_edges[,c('from','to','sim')], rows=1:n, columns=1:n, verbose=F))
			score <- Matrix::Matrix(mat_triu + t(mat_triu), sparse=TRUE)
		}
		
	}else if(measure=='SI'){
		# Salton Index (SI): cosine similarity
		if(type=='full' | (type=='edge' & edge.fast)){
			deg <- igraph::degree(ig)
			score <- igraph::cocitation(ig) / outer(deg, deg, function(x,y) sqrt(x*y))
			score <- Matrix::Matrix(score, sparse=TRUE)
			if(type=='edge'){
				score  <- igraph::get.adjacency(ig) * score
			}
		}else{
			df_edges$sim <- sapply(1:length(ls_common), function(i){
				length(ls_common[[i]]) / (sqrt(df_edges$degree_from[i] * df_edges$degree_to[i]))
			})
			mat_triu <- as.matrix(xSparseMatrix(df_edges[,c('from','to','sim')], rows=1:n, columns=1:n, verbose=F))
			score <- Matrix::Matrix(mat_triu + t(mat_triu), sparse=TRUE)
		}
		
	}else if(measure=='JI'){
		# Jaccard Index (JI)
		if(type=='full' | (type=='edge' & edge.fast)){
			score <- igraph::similarity.jaccard(ig)
			score <- Matrix::Matrix(score, sparse=TRUE)
			if(type=='edge'){
				score  <- igraph::get.adjacency(ig) * score
			}
		}else{
			ls_union <- lapply(1:nrow(df_edges), function(i){
				from <- ls_neighbors[[df_edges$from[i]]]
				to <- ls_neighbors[[df_edges$to[i]]]
				union(from, to)
			})
			df_edges$sim <- sapply(1:length(ls_common), function(i){
				length(ls_common[[i]]) / length(ls_union[[i]])
			})
			mat_triu <- as.matrix(xSparseMatrix(df_edges[,c('from','to','sim')], rows=1:n, columns=1:n, verbose=F))
			score <- Matrix::Matrix(mat_triu + t(mat_triu), sparse=TRUE)
		}
		
	}else if(measure=='DI'){
		# Dice Index (DI): mainly for ecological community data
		if(type=='full' | (type=='edge' & edge.fast)){
			score <- igraph::similarity.dice(ig)
			score <- Matrix::Matrix(score, sparse=TRUE)
			if(type=='edge'){
				score  <- igraph::get.adjacency(ig) * score
			}
		}else{
			df_edges$sim <- sapply(1:length(ls_common), function(i){
				2 * length(ls_common[[i]]) / (df_edges$degree_from[i] + df_edges$degree_to[i])
			})
			mat_triu <- as.matrix(xSparseMatrix(df_edges[,c('from','to','sim')], rows=1:n, columns=1:n, verbose=F))
			score <- Matrix::Matrix(mat_triu + t(mat_triu), sparse=TRUE)
		}
		
	}else if(measure=='HPI'){
		# Hub Promoted Index (HPI): quantifying the topological overlap of pairs of substrates in metabolic networks, with links adjacent to hubs receiving higher scores
		if(type=='full' | (type=='edge' & edge.fast)){
			deg <- igraph::degree(ig)
			score <- igraph::cocitation(ig) / outer(deg, deg, base::pmin)
			score <- Matrix::Matrix(score, sparse=TRUE)
			if(type=='edge'){
				score  <- igraph::get.adjacency(ig) * score
			}
		}else{
			df_edges$sim <- sapply(1:length(ls_common), function(i){
				length(ls_common[[i]]) / min(df_edges$degree_from[i], df_edges$degree_to[i])
			})
			mat_triu <- as.matrix(xSparseMatrix(df_edges[,c('from','to','sim')], rows=1:n, columns=1:n, verbose=F))
			score <- Matrix::Matrix(mat_triu + t(mat_triu), sparse=TRUE)
		}
		
	}else if(measure=='HDI'){
		# Hub Depressed Index (HPI): links adjacent to hubs receiving lower scores
		if(type=='full' | (type=='edge' & edge.fast)){
			deg <- igraph::degree(ig)
			score <- igraph::cocitation(ig) / outer(deg, deg, base::pmax)
			score <- Matrix::Matrix(score, sparse=TRUE)
			if(type=='edge'){
				score  <- igraph::get.adjacency(ig) * score
			}
		}else{
			df_edges$sim <- sapply(1:length(ls_common), function(i){
				length(ls_common[[i]]) / max(df_edges$degree_from[i], df_edges$degree_to[i])
			})
			mat_triu <- as.matrix(xSparseMatrix(df_edges[,c('from','to','sim')], rows=1:n, columns=1:n, verbose=F))
			score <- Matrix::Matrix(mat_triu + t(mat_triu), sparse=TRUE)
		}
		
	}else if(measure=='LHN'){
		# Leicht-Holme-Newman Index (LHN): inversely proportional to the expected number of common neighbors in the configuration model
		if(type=='full' | (type=='edge' & edge.fast)){
			deg <- igraph::degree(ig)
			score <- igraph::cocitation(ig) / outer(deg,deg)
			score <- Matrix::Matrix(score, sparse=TRUE)
			if(type=='edge'){
				score  <- igraph::get.adjacency(ig) * score
			}
		}else{
			df_edges$sim <- sapply(1:length(ls_common), function(i){
				length(ls_common[[i]]) / (df_edges$degree_from[i] * df_edges$degree_to[i])
			})
			mat_triu <- as.matrix(xSparseMatrix(df_edges[,c('from','to','sim')], rows=1:n, columns=1:n, verbose=F))
			score <- Matrix::Matrix(mat_triu + t(mat_triu), sparse=TRUE)
		}
		
	}else if(measure=='AA'){
		# Adamic-Adar Index (AA): penalising the high-degree common neighbors
		if(type=='full' | (type=='edge' & edge.fast)){
			score <- igraph::similarity.invlogweighted(ig)
			score <- Matrix::Matrix(score, sparse=TRUE)
			if(type=='edge'){
				score  <- igraph::get.adjacency(ig) * score
			}
		}else{
			df_edges$sim <- sapply(1:length(ls_common), function(i){
				x <- ls_common[[i]]
				y <- vec_degree[x]
				sum(1/log(y))
			})
			mat_triu <- as.matrix(xSparseMatrix(df_edges[,c('from','to','sim')], rows=1:n, columns=1:n, verbose=F))
			score <- Matrix::Matrix(mat_triu + t(mat_triu), sparse=TRUE)
		}
		
	}else if(measure=='RA'){
		# Resource Allocation Index (RA): penalising even more the high-degree common neighbors. Motivated by the resource allocation dynamics on complex networks, assuming that each common neighbor (transmitter) has a unit of resource being distributed equally to all its neighbors. The similarity defined as the amount of resource transmitted betweeen two
		if(type=='edge'){ 
			df_edges$sim <- sapply(1:length(ls_common), function(i){
				x <- ls_common[[i]]
				y <- vec_degree[x]
				sum(1/y)
			})
			mat_triu <- as.matrix(xSparseMatrix(df_edges[,c('from','to','sim')], rows=1:n, columns=1:n, verbose=F))
			score <- Matrix::Matrix(mat_triu + t(mat_triu), sparse=TRUE)
		}else if(type=='full'){
			score <- matrix(0, nrow=n, ncol=n)
			for(i in 1:(n-1)){
				from <- ls_neighbors[[i]]
				for(j in (i+1):n){
					to <- ls_neighbors[[j]]
					x <- intersection(from, to)
					y <- vec_degree[x]
					score[i,j] <- score [j,i] <- sum(1/y)
				}
			}
			score <- Matrix::Matrix(score, sparse=TRUE)
		}
		
	}else if(measure=='LP'){
		# Local Path Index (LP): the number of two-step paths and three-step paths (weighted by the eps), providing a tradeoff of accuracy and computational complexity
		eps <- 0.01
		if(names(measure.para) %in% "LP.eps"){
			if(!is.null(measure.para$LP.eps)){
				eps <- measure.para$LP.eps
			}
		}
		A <- igraph::get.adjacency(ig)
		A2 <- A %*% A
		A3 <- A2 %*% A
		score <- A2 + A3 * eps
		if(type=='edge'){
			score  <- igraph::get.adjacency(ig) * score
		}
		
	}else if(measure=='KI'){
		# Katz Index (KI): based on the ensemble of all paths, which directly sums over the collection of paths and is exponentially damped by length to give the shorter paths more weights
		beta <- 0.001
		if(names(measure.para) %in% "KI.beta"){
			if(!is.null(measure.para$KI.beta)){
				beta <- measure.para$KI.beta
			}
		}
		A  <- igraph::get.adjacency(ig)
    	I <- diag(n)
    	score <- solve(I - beta * A) - I
		score <- Matrix::Matrix(score, sparse=TRUE)
		if(type=='edge'){
			score  <- igraph::get.adjacency(ig) * score
		}
    	
	}else if(measure=='LHNg'){
		# Leicht-Holme-Newman Index global similarity (LHNg): two nodes are similar if their immediate neighbors are similar by themselves
		theta <- 0.5
		if(names(measure.para) %in% "LHNg.theta"){
			if(!is.null(measure.para$LHNg.theta)){
				theta <- measure.para$LHNg.theta
			}
		}
		lambda <- igraph::graph.eigen(ig)$value
		deg <- igraph::degree(ig)
		A  <- igraph::get.adjacency(ig)
    	I <- diag(n)
    	score <- solve(I - (theta / lambda) * A) / outer(deg,deg)
		score <- Matrix::Matrix(score, sparse=TRUE)
		if(type=='edge'){
			score  <- igraph::get.adjacency(ig) * score
		}
    	
	}else if(measure=='SP'){
		# Shortest Path Index (SP): based on the length of the shortest paths
		score <- igraph::shortest.paths(ig)
    	score[is.infinite(score)] <- igraph::vcount(ig) + 1
    	score <- 1 / score
    	diag(score) <- 0
		score <- Matrix::Matrix(score, sparse=TRUE)
		if(type=='edge'){
			score  <- igraph::get.adjacency(ig) * score
		}
    	
	}else if(measure=='ACT'){
		# Average Commute Time (ACT): the average number of steps required by a random walker starting from a node to reach another node
  		L <- igraph::graph.laplacian(ig)
  		L_psinv <- solve(L - 1/n) + 1/n
  		score <- 1 / (2 * m * (diag(L_psinv) %*% t(rep(1,n)) + rep(1,n) %*% t(diag(L_psinv)) - 2 * L_psinv))
    	diag(score) <- 0
		score <- Matrix::Matrix(score, sparse=TRUE)
		if(type=='edge'){
			score  <- igraph::get.adjacency(ig) * score
		}
    	
	}else if(measure=='CL'){
		# Cosine based on Laplacian matrix (CL): the cosine of the angle between columns of the pseudoinverse of the Laplacian matrix
  		L <- igraph::graph.laplacian(ig)
  		L_psinv <- solve(L - 1/n) + 1/n
  		dL <- diag(L_psinv)
  		score <- L_psinv / outer(dL,dL,function(x,y) sqrt(x*y))
  		diag(score) <- 0
		score <- Matrix::Matrix(score, sparse=TRUE)
		if(type=='edge'){
			score  <- igraph::get.adjacency(ig) * score
		}
    	
	}else if(measure=='MFI'){
		# Matrix Forest Index (MFI): the ratio of the number of spanning rooted forests belonging to the same tree to all spanning rooted forests of the network
		L <- igraph::graph.laplacian(ig)
		I <- diag(n)
    	score <- solve(I + L)
		score <- Matrix::Matrix(score, sparse=TRUE)
		if(type=='edge'){
			score  <- igraph::get.adjacency(ig) * score
		}
    	
	}else if(measure=='RWR'){
		# Random Walk with Restart (RWR)
		# solved analytically (not iteratively)
		r <- 0.5
		if(names(measure.para) %in% "RWR.r"){
			if(!is.null(measure.para$RWR.r)){
				r <- measure.para$RWR.r
			}
		}
  		if(1){
			if ("weight" %in% igraph::list.edge.attributes(ig)){
				adjM <- igraph::get.adjacency(ig, type="both", attr="weight", edges=FALSE, names=TRUE, sparse=getIgraphOpt("sparsematrices"))
				if(verbose){
					message(sprintf("\tNotes: using weighted graph!"), appendLF=TRUE)
				}
			}else{
				adjM <- igraph::get.adjacency(ig, type="both", attr=NULL, edges=FALSE, names=TRUE, sparse=getIgraphOpt("sparsematrices"))
				if(verbose){
					message(sprintf("\tNotes: using unweighted graph!"), appendLF=TRUE)
				}
			}
			## D is the degree matrix of the graph (^-1/2)
			A <- adjM!=0
			D <- Matrix::Diagonal(x=(Matrix::colSums(A))^(-0.5))
			nadjM <- as.matrix(D %*% adjM %*% D)			
			## steady probability
			score <- solve((diag(nrow(nadjM)) - (1-r) * t(nadjM)) / r)
  		}else{
  			P <- as.matrix(igraph::get.stochastic(ig)) # the transition matrix
			score <- solve((diag(nrow(P)) - (1 - r) * t(P)) / r)
			score <- score + t(score)
  		}
		score <- Matrix::Matrix(score, sparse=TRUE)
		if(type=='edge'){
			score  <- igraph::get.adjacency(ig) * score
		}
    	
	}else if(measure=='SRW'){
		# Short Random Walk (SRW)
		# solved analytically
		p <- 4
		if(names(measure.para) %in% "SRW.p"){
			if(!is.null(measure.para$SRW.p)){
				p <- measure.para$SRW.p
			}
		}
		score <- xRWkernel(ig, steps=p, verbose=F)
		if(type=='edge'){
			score  <- igraph::get.adjacency(ig) * score
		}
    	
	}
	
	colnames(score) <- rownames(score) <- 1:n
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)

    invisible(score)
}


