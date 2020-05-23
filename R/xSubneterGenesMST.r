#' Function to identify a minimum spanning tree for subnetworks
#'
#' \code{xSubneterGenesAdv} is supposed to identify a minimum spanning tree for subnetworks. It returns an object of class "igraph". 
#'
#' @param isubg an "iSubg" object resulting from \code{\link{xSubneterGenesAdv}}
#' @param metric the distance metric for subnetworks. It can be either "max" for the maximum distance between any two nodes (one from a subnetwork, and other from another subnetwork) based on the whole network, or "jaccard"  for jaccard distance between two subnetworks (nodes overlapped), or "hybrid" (that is, "max" multiplied by "jaccard")
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a minimum spanning tree, an object of class "igraph". It has graph attributes ('summary', 'detail' and 'matrix'), and node attributes ('xcoord', 'ycoord', 'num_nodes', 'num_edges', 'weight' and 'weight_scaled' [1,5] for visualisation).
#' @export
#' @seealso \code{\link{xSubneterGenesAdv}}
#' @include xSubneterGenesMST.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#'
#' # a) provide the input nodes/genes with the significance info
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' sig <- rbeta(500, shape1=0.5, shape2=1)
#' data <- data.frame(symbols=org.Hs.eg$gene_info$Symbol[1:500], sig)
#' 
#' # b) find a series of maximum-scoring subnets with the desired node number=50
#' isubg <- xSubneterGenesAdv(data=data, network="STRING_high", subnet.size=50, RData.location=RData.location)
#'
#' # c) represent a series of subnets as a minimum spanning tree
#' mst <- xSubneterGenesMST(isubg)
#' mst$summary
#' head(mst$detail)
#' head(mst$matrix)
#' gp_mst <- xGGnetwork(mst, node.label='name', node.label.size=3, node.label.force=1, node.xcoord='xcoord', node.ycoord='ycoord', edge.size='weight_scaled', node.size='num_edges', node.size.title="Num of \nedges", node.size.range=c(1,4))
#' }

xSubneterGenesMST <- function(isubg, metric=c("hybrid","max","jaccard"), verbose=T)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    metric <- match.arg(metric)
    
   	if(is(isubg,"iSubg")){
   		ls_subg <- isubg$ls_subg
   		g <- isubg$g
		## Remove null elements in a list
		ls_subg <- base::Filter(base::Negate(is.null), ls_subg)
		if(length(ls_subg)<=1){
			return(NULL)
		}
	}else{
		stop("The function must apply to a 'iSubg' objects.\n")
	}
    
    ## df_summary
	ls_df <- lapply(1:length(ls_subg), function(i){
		x <- ls_subg[[i]]
		y <- range(as.numeric(V(x)$significance))
		data.frame(iteration=i, num.node=vcount(x), num.edge=ecount(x), num.linker=sum(V(x)$type=='linker'), threshold=x$threshold, sig.min=y[1], sig.max=y[2], stringsAsFactors=F)
	})
	df_summary <- do.call(rbind, ls_df)
	
	## df_detail
	ls_df <- lapply(1:length(ls_subg), function(i){
		x <- ls_subg[[i]]
		y <- as_data_frame(x, what='vertices')
		data.frame(iteration=i, y, stringsAsFactors=F)
	})
	df_detail <- do.call(rbind, ls_df)
    
	## mat
	name <- iteration <- score <- NULL
	mat <- df_detail %>% dplyr::select(name,iteration,score) %>% tidyr::spread(key=iteration, value=score)
	rownames(mat) <- mat[,1]
	mat <- mat[,-1]
	mat[is.na(mat)] <- 0
	mat[mat<0] <- -1
	mat[mat>0] <- 1
    mat <- data.matrix(mat)
    
    ## mst
    ### calculate distance between subg
    n <- ncol(mat)
	mat_max <- matrix(0, nrow=n, ncol=n)
	mat_jaccard <- matrix(0, nrow=n, ncol=n)
	for(i in 1:(n-1)){
		name_i <- names(which(mat[,i]!=0))
		v <- match(name_i, V(g)$name)
		for(j in (i+1):n){
			name_j <- names(which(mat[,j]!=0))
			to <- match(name_j, V(g)$name)
			dist_tmp <- igraph::distances(g, v=v, to=to, mode='all')
			mat_max[i,j] <- max(dist_tmp)
			mat_jaccard[i,j] <- 1 - length(intersect(v, to)) / length(union(v,to))
		}
	}
	mat_max <- mat_max + t(mat_max)
	colnames(mat_max) <- rownames(mat_max) <- colnames(mat)
	mat_jaccard <- mat_jaccard + t(mat_jaccard)
	colnames(mat_jaccard) <- rownames(mat_jaccard) <- colnames(mat)
	if(metric=='max'){
		mat_dist <- mat_max
	}else if(metric=='jaccard'){
		mat_dist <- mat_jaccard
	}else if(metric=='hybrid'){
		mat_dist <- mat_max * mat_jaccard
	}
	### obtain mst
	#### first, convert distnace matrix into a graph
	g_subg <- xConverter(mat_dist, from="dgCMatrix", to="igraph", verbose=verbose)
	#### second, get mst (connected graph with distance miminised)
	mst <- igraph::minimum.spanning.tree(g_subg, weights=E(g_subg)$weight)
	#### third, do visualisation
	##### rescale the weight [1,5]
	w <- E(mst)$weight
	E(mst)$weight_scaled <- 1 + 4*(w - min(w))/(max(w) - min(w))
	##### coords determined considering edge weight 
	##### after transformation, nodes with less weights are closer to each other
	set.seed(825)
	coords <- igraph::layout_with_fr(mst, weights=5-(E(mst)$weight_scaled-1))
	V(mst)$xcoord <- coords[,1]
	V(mst)$ycoord <- coords[,2]
	##### append node attributes ('num_nodes' and 'num_edges')
	ind <- match(V(mst)$name, df_summary$iteration)
	V(mst)$num_nodes <- df_summary$num.node[ind]
	V(mst)$num_edges <- df_summary$num.edge[ind]
	##### append graph attributes ('summary', 'detail' and 'matrix')
	mst$summary <- df_summary
	mst$detail <- df_detail
	mst$matrix <- mat
	
	#gp_mst <- xGGnetwork(mst, node.label='name', node.label.size=3, node.label.force=1, node.xcoord='xcoord', node.ycoord='ycoord', edge.size='weight_scaled', node.size='num_edges', node.size.title="Num of \nedges", node.size.range=c(1,4))
    
    return(mst)
}
