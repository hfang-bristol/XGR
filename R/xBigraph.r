#' Function to obtain communities from a bipartitle graph
#'
#' \code{xBigraph} is supposed to obtain communities from a bipartitle graph.
#'
#' @param g an object of class "igraph" (or "graphNEL") for a bipartitel graph with a 'type' node attribute
#' @param algorithm the algorithm to initialise community from a projected graph. It can be 'louvain' using multi-level optimization, 'leading_eigen' using leading eigenvector, 'fast_greedy' using greedy optimization, and 'walktrap' via short random walks
#' @param seed an integer specifying the seed
#' @param glayout a graph layout function. This function can be one of "layout_nicely" (previously "layout.auto"), "layout_randomly" (previously "layout.random"), "layout_in_circle" (previously "layout.circle"), "layout_on_sphere" (previously "layout.sphere"), "layout_with_fr" (previously "layout.fruchterman.reingold"), "layout_with_kk" (previously "layout.kamada.kawai"), "layout_as_tree" (previously "layout.reingold.tilford"), "layout_with_lgl" (previously "layout.lgl"), "layout_with_graphopt" (previously "layout.graphopt"), "layout_with_sugiyama" (previously "layout.kamada.kawai"), "layout_with_dh" (previously "layout.davidson.harel"), "layout_with_drl" (previously "layout.drl"), "layout_with_gem" (previously "layout.gem"), "layout_with_mds", and "layout_as_bipartite". A full explanation of these layouts can be found in \url{http://igraph.org/r/doc/layout_nicely.html}
#' @param edge.color.alternative two alternative colors for edges within the community (grey70 by default) and edges between communities (grey95 by default)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns an igraph object, appended by node attributes including "community" for communit memberships, "contribution" for contribution to its community, "xcoord" for x-coordinates, "ycoord" for y-coordiates, and by edge attributes including "color" for between-community edges ('grey95') and within-community edges ('grey70').
#' @note The input graph will has an equal weight if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{xGGnetwork}}
#' @include xBigraph.r
#' @examples
#' # 1) generate a random bipartite graph
#' set.seed(123)
#' g <- sample_bipartite(100, 50, p=0.1)
#' V(g)$name <- V(g)
#' 
#' \dontrun{
#' # 2) obtain its community
#' ig <- xBigraph(g)
#' gp <- xGGnetwork(ig, node.label='name', node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0, node.label.force=0.002, node.shape='type', node.shape.title='Type', node.xcoord='xcoord', node.ycoord='ycoord', node.color='community', node.color.title='Community', colormap='jet.both', ncolors=64, zlim=NULL, node.size='contribution', node.size.range=c(1,4), node.size.title='Contribution', slim=NULL, edge.color="color",edge.color.alpha=0.5,edge.curve=0,edge.arrow.gap=0)
#' 
#' ## make it discrete for the colorbar
#' gp + scale_colour_gradientn(colors=xColormap('jet.both')(64),breaks=seq(1,max(V(ig)$community))) + guides(color=guide_legend(title="Community"))
#' }

xBigraph <- function(g, algorithm=c("louvain","leading_eigen","fast_greedy","walktrap"), seed=825, glayout=layout_with_kk, edge.color.alternative=c("grey70","grey95"), verbose=TRUE)
{
    
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    algorithm <- match.arg(algorithm)
    
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
	ls_res <- igraph::bipartite_mapping(ig)
	if(!(ls_res$res)){
		stop("The igraph object is not bipartitle.\n")
	}
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("The igraph object has %d nodes and %d edges (%s) ...", vcount(ig), ecount(ig), as.character(now)), appendLF=T)
	}
	
	## only keep the largest component
	ig <- dnet::dNetInduce(ig, nodes_query=V(ig)$name, knn=0, largest.comp=TRUE)
	
	if(verbose){
		message(sprintf("\tThe largest component has %d nodes and %d edges", vcount(ig), ecount(ig)), appendLF=T)
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
	xM <- t(sM) %*% sM
	colnames(xM) <- rownames(xM) <- xnode.names
	
	## ig_x: igraph from projected adjacency matrix on xnodes
	ig_x <- igraph::graph.adjacency(xM, mode="undirected", weighted=TRUE, diag=FALSE)
	### remove loops and multiple edges
	ig_x <- dnet::dNetInduce(ig_x, nodes_query=V(ig_x)$name, knn=0, largest.comp=TRUE)

	# find community for ig_x
	if(algorithm=="louvain"){
		## multi-level optimization of modularity
		cs0 <- igraph::cluster_louvain(ig_x, weights=E(ig_x)$weight)
	}else if(algorithm=="leading_eigen"){
		## leading eigenvector of the community matrix
		cs0 <- igraph::cluster_leading_eigen(ig_x, weights=E(ig_x)$weight)
	}else if(algorithm=="fast_greedy"){
		## greedy optimization of modularity
		cs0 <- igraph::cluster_fast_greedy(ig_x, weights=E(ig_x)$weight)
	}else if(algorithm=="walktrap"){
		## via short random walks
		cs0 <- igraph::cluster_walktrap(ig_x, weights=E(ig_x)$weight)
	}
	
	## initial condition for ig_x community membership
	T0 <- data.frame(as.integer(factor(xnode.names)), as.vector(membership(cs0)))
	## bipartite modularity maximization using initial assignments T0
	if(1){
		#define a local function
		machine.which.max = function(X){
			thresh <- .Machine$double.eps
			X[abs(X) <= thresh] <- 0
			true_max <- which.max(X)
			return(true_max)
		}
		machine.max = function(X){
			thresh <- .Machine$double.eps
			X[abs(X) <= thresh] <- 0
			true_max <- max(X)
			return(true_max)
		}

		## number of edges
		elist.index <- cbind(ynodes, xnodes)
		m <- nrow(elist)
		deltaQmin <- min(10^-4, 1/m)

		A <- as.matrix(sM)
		rownames(A) <- ynode.names
		colnames(A) <- xnode.names

		p <- nrow(A)
		q <- ncol(A)
		N <- p+q

		#initialize community assignments for ynodes
		T0[,1] <- as.integer(factor(T0[,1]))
		T0 <- T0[order(T0[,1]),]
		Tind <- data.matrix(T0)
		Rind <- data.matrix(cbind(1:p,rep(0,length=p)))
		cs <- sort(unique(Tind[,2]))
	
		#variables to track modularity changes after each loop
		Qhist <- vector()
		Qnow  <- 0
		deltaQ <- 1
	
		#Sparse matrix with the upper right block of the true Adjacency matrix
		ki <- rowSums(A)
		dj <- colSums(A)
		m <- sum(ki)

		#Make B tilde, Btilde = Aij - (ki*dj)/m
		Btilde <- -(ki %o% dj)/m
		Btilde[elist.index] <- weights + Btilde[elist.index]
		
		#initialize Tm
		#Tm = sparseMatrix(i=Tind[,1],j=Tind[,2],x=1,index1=TRUE)
		Tm <- matrix(0,nrow=q,ncol=max(Tind[,2]))
		Tm[Tind] <- 1 
		
		###########################################################
		iter <- 1
		while(deltaQ > deltaQmin){
	
			### Step 1, assign red nodes
			Ttilde <- Btilde %*% Tm
	
			#Find first max for all rows, update R
			Rind[,2] <- unname(apply(Ttilde, 1, machine.which.max))
	
			#Special condition if all nodes are stuck in one large community,
			#if TRUE, randomly assign two nodes to new communities.
			if(iter > 1 && length(unique(c(Rind[,2],Tind[,2]))) == 1){
				random_nodes <- sample(1:nrow(Rind),2)
				Rind[random_nodes,2] <- max(c(Rind[,2],Tind[,2])) + 1:2
			}
	
			#Check to see if new communities should be made
			negative_contribution <- unname(apply(Ttilde, 1,machine.max)) < 0
			if( any( negative_contribution )){
				#add new communities
				num_new_com <- sum(negative_contribution)
				cs <- length(unique(c(Tind[,2],Rind[,2])))
				Rind[negative_contribution,2] <- (cs + 1):(cs + num_new_com)
			}
	
			#Rm = sparseMatrix(i=Rind[,1],j=Rind[,2],x=1,index1=TRUE)
			Rm <- matrix(0,nrow=p,ncol=max(Rind[,2]))
			Rm[data.matrix(Rind)] <- 1
			
			### Step 2, assign blue nodes
			# R tilde = transpose(Btilde) %*% R
			Rtilde <- crossprod(Btilde,Rm)
	
			#Find first max for all rows, update T
			Tind[,2] <- unname(apply(Rtilde, 1, machine.which.max))
	
			#Check to see if new communities should be made
			negative_contribution <- unname(apply(Rtilde, 1,machine.max)) < 0
			if( any( negative_contribution )){
				#add new communities
				num_new_com <- sum(negative_contribution)
				cs <- length(unique(c(Tind[,2],Rind[,2])))
				Tind[negative_contribution,2] <- (cs + 1):(cs + num_new_com)
			}
			#Tm dimensions, note the extra empty community.
			#Tm = sparseMatrix(i=Tind[,1],j=Tind[,2],x=1,index1=TRUE) 
			Tm <- matrix(0,nrow=q,ncol=max(Tind[,2]))
			Tm[data.matrix(Tind)] <- 1 
	
			Qthen <- Qnow
			#replace with diag(crossprod(T,BTR))/m
			Qcom <- diag(crossprod(Rm,Btilde %*% Tm))/m
			Qnow <- sum(Qcom)
			if(abs(Qnow) < .Machine$double.eps){Qnow <- 0}
			Qhist <- c(Qhist,Qnow)

			#print(paste("Q =",Qnow,sep=" "))
			if(Qnow != 0){
				deltaQ <- Qnow - Qthen
			}
			iter <- iter+1
		}
		###########################################################
	
		cs <- 1:ncol(Tm)
		if(any(sort(unique(Tind[,2])) != sort(unique(Rind[,2])))){
			stop("number of ynode and xnode communities unequal")
		}
	
		#drop empty communities
		qcom_temp <- cbind(Qcom,cs)
		qcom_out <- qcom_temp[abs(Qcom) > .Machine$double.eps,]
	
		#if communities were dropped, relabel so community labels can function
		#as row/column indices in the modularity matrix, B_ij.
		if(nrow(qcom_out) < nrow(qcom_temp)){
			qcom_out[,2] <- as.integer(factor(qcom_out[,2]))
			Rind[,2] <- as.integer(factor(Rind[,2]))
			Tind[,2] <- as.integer(factor(Tind[,2]))
		}
		colnames(qcom_out) <- c("Qcom","community")
		
		Qcoms <- as.data.frame(qcom_out)
		modularity <- Qhist
		ynode.memb	<- data.frame(ynode.names=ynode.names[Rind[,1]], com=Rind[,2])
		xnode.memb	<- data.frame(xnode.names=xnode.names[Tind[,1]], com=Tind[,2])
		
		#######################################################
		
		## contribution to its modularity
		R1 <- ynode.memb
		T1 <- xnode.memb
		r1 <- cbind(as.numeric(factor(R1[,1])), R1[,2])
		t1 <- cbind(as.numeric(factor(T1[,1])), T1[,2]) 
		Rtrans <- Matrix::sparseMatrix(i=r1[,2], j=r1[,1], x=1, dims=c(max(r1[,2]),length(unique(r1[,1]))))  
		T2 <- Matrix::sparseMatrix(i=t1[,1], j=t1[,2], x=1, dims=c(max(t1[,1]),max(t1[,2])))

		Qjk <- vector(length=q)
		for(j in 1:max(t1[,1])){
			if(j %% 1000 == 0){print(paste(j,t1[j,]))}
			Bj <- A[,j] - (ki*dj[j])/m;
			Qjk[j] = ((Rtrans[t1[j,2],] %*% Bj)/(2*m))*(1/Qcoms[t1[j,2],1])
		}
		Qik <- vector(length=p)
		for(i in 1:max(r1[,1])){
			if(i %% 1000 == 0){print(i)}
			Bi <- A[i,] - (ki[i]*dj)/m;
			Qik[i] <- ((Bi %*% T2[,r1[i,2]])/(2*m))*(1/Qcoms[r1[i,2],1])  
		}
		
		xnode.memb <- data.frame(T1, Qjk)
		ynode.memb <- data.frame(R1, Qik)
		
		out <- list(Qcoms		= Qcoms,
					ynode.memb 	= ynode.memb,
					xnode.memb  = xnode.memb
					)
		
		## rename the modules
		df <- out$Qcoms
		df$xnode <- table(out$xnode.memb$com)
		df$ynode <- table(out$ynode.memb$com)
		## sorted by ynode, Qcom
		Qcom <- ynode <- xnode <- NULL
		df <- df %>% dplyr::arrange(-ynode,-Qcom) %>% dplyr::mutate(new_community=1:nrow(out$Qcoms))
		### rename
		ind <- match(out$ynode.memb$com, df$com)
		out$ynode.memb$com <- df$new_community[ind]
		ind <- match(out$xnode.memb$com, df$com)
		out$xnode.memb$com <- df$new_community[ind]
	}
	
	############
	############
	ig_converted <- igraph::graph.data.frame(elist, directed=FALSE)
	ind_x <- match(V(ig_converted)$name, out$xnode.memb$xnode.names)
	ind_y <- match(V(ig_converted)$name, out$ynode.memb$ynode.names)
	## type
	V(ig_converted)$type <- NA
	V(ig_converted)$type[!is.na(ind_x)] <- 'xnode'
	V(ig_converted)$type[!is.na(ind_y)] <- 'ynode'
	## community
	V(ig_converted)$community <- 0
	V(ig_converted)$community[!is.na(ind_x)] <- out$xnode.memb$com[ind_x[!is.na(ind_x)]]
	V(ig_converted)$community[!is.na(ind_y)] <- out$ynode.memb$com[ind_y[!is.na(ind_y)]]
	## contribution	
	V(ig_converted)$contribution <- 0
	V(ig_converted)$contribution[!is.na(ind_x)] <- out$xnode.memb$Qjk[ind_x[!is.na(ind_x)]]
	V(ig_converted)$contribution[!is.na(ind_y)] <- out$ynode.memb$Qik[ind_y[!is.na(ind_y)]]
	############
	############
	
	#############################
	## append node attributes: xcord, ycord
	## append edge attribute: color
	#############################
	ig_converted <- xAddCoords(ig_converted, node.attr="community", seed=seed, glayout=glayout, edge.color.alternative=edge.color.alternative, verbose=verbose)
	
 ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)

    invisible(ig_converted)
}


