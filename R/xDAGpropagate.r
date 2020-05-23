#' Function to generate a subgraph of a direct acyclic graph (DAG) propagaged by the input annotation data
#'
#' \code{xDAGpropagate} is supposed to produce a subgraph induced by the input annotation data, given a direct acyclic graph (DAG; an ontology). The input is a graph of "igraph", a list of the vertices containing annotation data, and the mode defining the paths to the root of DAG. The annotations are propagated to the ontology root (eg, retaining the minmum pvalue). The propagaged subgraph contains vertices (with annotation data) and their ancestors along with the defined paths to the root of DAG. The annotations at these vertices (including their ancestors) can also be updated according to the true-path rule: those annotated to a term should also be annotated by its all ancestor terms.
#'
#' @param g an object of class "igraph" to represent DAG
#' @param annotation the vertices/nodes for which annotation data are provided. It can be a sparse Matrix of class "dgCMatrix" (with variants/genes as rows and terms as columns), or a data frame with three columns: 1st column for variants/genes, 2nd column for terms, and 3rd column for values, or a list (with the name for terms) each element storing a named vactor (that is, value for the content and variants/genes as names)
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param propagation how to propagate the score. It can be "max" for retaining the maximum value from its children, "min" for retaining the minimum value from its children, and 'all' for retaining all from its children (by default)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' \itemize{
#'  \item{\code{subg}: an induced/propagated subgraph, an object of class "igraph". In addition to the original attributes to nodes and edges, the return subgraph is also appended by two node attributes: 1) "anno" containing a list of variants/genes (with numeric values as elements); 2) "IC" standing for information content defined as negative 10-based log-transformed frequency of variants/genes annotated to that term.}
#' }
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xDAGpropagate.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' # 1) EF ontology
#' # ig.EF (an object of class "igraph" storing as a directed graph)
#' ig <- xRDataLoader('ig.EF', RData.location=RData.location)
#' ## optional: extract the disease part (EFO:0000408)
#' neighs.out <- igraph::neighborhood(ig, order=vcount(ig), nodes='EFO:0000408', mode="out")
#' vids <- V(ig)[unique(unlist(neighs.out))]$name
#' g <- igraph::induced.subgraph(ig, vids=vids)
#' 
#' ################
#' # 2a) load GWAS SNPs annotated by EF (an object of class "dgCMatrix" storing a sparse matrix)
#' annotation <- xRDataLoader(RData='GWAS2EF', RData.location=RData.location)
#' ## only significant
#' annotation[as.matrix(annotation>5e-8)] <- 0
#'
#' # 2b) propagation based on shortest paths (ie the most concise subgraph)
#' dag <- xDAGpropagate(g, annotation, path.mode="shortest_paths", propagation="min")
#'
#' # 2c) color-code nodes/terms according to the number of annotations
#' data <- sapply(V(dag)$anno, length)
#' names(data) <- V(dag)$name
#' ## only those GWAS>=100
#' nodes <- V(dag)$name[data>=100]
#' dagg <- igraph::induced.subgraph(dag, vids=nodes)
#' ### DAG plot
#' dnet::visDAG(dagg, data, node.info="both")
#' ### Net plot
#' set.seed(825); glayout <- layout_with_kk(dagg)
#' xVisNet(dagg, pattern=data, colormap="yr", glayout=glayout, vertex.label=V(dagg)$term_name, vertex.shape="sphere", vertex.label.font=2, vertex.label.dist=0.2, vertex.label.cex=0.5, zlim=c(100,300))
#' ### interpolation plot
#' set.seed(825); glayout <- layout_with_kk(dagg)
#' pattern <- sapply(V(dagg)$anno, length)
#' ls_xyz <- data.frame(x=glayout[,1], y=glayout[,2], z=log10(pattern))
#' xVisInterp(ls_xyz, nD="auto", image=TRUE)
#' 
#' ################
#' 3a) load ChEMBL targets annotated by EF (an object of class "dgCMatrix" storing a sparse matrix)
#' annotation <- xRDataLoader(RData='Target2EF', RData.location=RData.location)
#' ## only approved (phase=4)
#' annotation[as.matrix(annotation<4)] <- 0
#'
#' 3b) propagation based on all paths
#' dag <- xDAGpropagate(g, annotation, path.mode="all_paths", propagation="max")
#'
#' 3c) color-code nodes/terms according to the number of annotations
#' data <- sapply(V(dag)$anno, length)
#' names(data) <- V(dag)$name
#' ## only those Targets>=50
#' nodes <- V(dag)$name[data>=50]
#' dagg <- igraph::induced.subgraph(dag, vids=nodes)
#' ### DAG plot
#' dnet::visDAG(dagg, data, node.info="both")
#' ### Net plot
#' set.seed(825); glayout <- layout_with_kk(dagg)
#' xVisNet(dagg, pattern=data, colormap="yr", glayout=glayout, vertex.label=V(dagg)$term_name, vertex.shape="sphere", vertex.label.font=2, vertex.label.dist=0.2, vertex.label.cex=0.5, zlim=c(50,300))
#' ### interpolation plot
#' set.seed(825); glayout <- layout_with_kk(dagg)
#' pattern <- sapply(V(dagg)$anno, length)
#' ls_xyz <- data.frame(x=glayout[,1], y=glayout[,2], z=log10(pattern))
#' xVisInterp(ls_xyz, nD="3D", contour=TRUE)
#' }

xDAGpropagate <- function (g, annotation, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), propagation=c("all","min","max"), verbose=TRUE)
{
    
    path.mode <- match.arg(path.mode)
    propagation <- match.arg(propagation)
    
    ig <- g
    if (!is(ig,"igraph")){
        stop("The function must apply to the 'igraph' object.\n")
    }
    
    if(is(annotation,"list")){
    	oAnnos <- annotation
    }else{
		if(is(annotation,"data.frame")){
			annotation <- xSparseMatrix(annotation, verbose=FALSE)
		}
	
		if(is(annotation,"dgCMatrix")){
			D <- annotation
			oAnnos <- lapply(1:ncol(D), function(j){
				ind <- which(D[,j]!=0)
				x <- D[ind,j]
				names(x) <- names(ind)
				x
			})
			names(oAnnos) <- colnames(annotation)
	
		}else{
			stop("The input annotation must be either 'GS' or 'list' or 'dgCMatrix' object.\n")
		}
	}    
    
    ## check nodes in annotation
    if (is.list(oAnnos)){
        originNodes <- names(oAnnos)
        
        ind <- match(originNodes, V(ig)$name)
        nodes_mapped <- originNodes[!is.na(ind)]
        if(length(nodes_mapped)==0){
            stop("The input annotation data do not contain terms matched to the nodes/terms in the input graph.\n")
        }
    }
    
	#######################
    ## generate a subgraph of a direct acyclic graph (DAG) induced by terms from input annotations
    dag <- dnet::dDAGinduce(ig, originNodes, path.mode=path.mode)
    allNodes <- V(dag)$name
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Do propagation via '%s' operation (%s) ...", propagation, as.character(now)), appendLF=T)
    }
    
    ## Propagate using list
    if(1){
        ## initialise annotations    
        pAnnos <- oAnnos[allNodes]
        names(pAnnos) <- allNodes

        ## get the levels list
        level2node <- dnet::dDAGlevel(dag, level.mode="longest_path", return.mode="level2node")
        nLevels <- length(level2node)
        for(i in nLevels:2) {
            currNodes <- level2node[[i]]

            ## get the incoming neighbors (excluding self) that are reachable (i.e. nodes from i-1 level)
            adjNodesList <- lapply(currNodes, function(node){
                neighs.in <- igraph::neighborhood(dag, order=1, nodes=node, mode="in")
                setdiff(V(dag)[unlist(neighs.in)]$name, node)
            })
            names(adjNodesList) <- currNodes

            ## inherit the annotations from level i to level i - 1
            for(k in 1:length(currNodes)){
                node <- currNodes[k]
                ## get the domain annotations from this current node
                nowDomain <- pAnnos[[node]]
                nowDomain_mat <- cbind(names(nowDomain), as.numeric(nowDomain))
            
                ## assigin inherit annotations to all its adjacent nodes
                adjNodes <- adjNodesList[[node]]
                res <- lapply(adjNodes, function(adjNode){
                    ## get the domain annotations from this adjacent node
                    adjDomain <- pAnnos[[adjNode]]
                    adjDomain_mat <- cbind(names(adjDomain), as.numeric(adjDomain))
            
                    ### update: keep the largest score if overlap
                    all_mat <- rbind(nowDomain_mat, adjDomain_mat)
                    all_list <- base::split(x=as.numeric(all_mat[,2]), f=all_mat[,1])
                    output_list <- lapply(1:length(all_list), function(i){
                    	x <- all_list[[i]]
                        if(propagation=='max'){
                            y <- max(x)
                            #data.frame(name=names(all_list)[i], value=y, stringsAsFactors=FALSE)
                        }else if(propagation=='min'){
                            y <- min(x)
                            #data.frame(name=names(all_list)[i], value=y, stringsAsFactors=FALSE)
                        }else{
                        	y <- unique(x)
                        }
                        data.frame(name=rep(names(all_list)[i],length(y)), value=y, stringsAsFactors=FALSE)
                    })
                    
                    if(0){
						x_mat <- base::do.call(base::rbind, output_list)
						output <- as.vector(x_mat)
						names(output) <- rownames(x_mat)
                	}else{
						x_mat <- base::do.call(base::rbind, output_list)
						output <- x_mat$value
						names(output) <- x_mat$name
                	}
                	
                    return(output)
                    
                })
                pAnnos[adjNodes] <- res
            }
        
            if(verbose){
                message(sprintf("\tAt level %d, there are %d nodes, and %d incoming neighbors (%s).", i, length(currNodes), length(unique(unlist(adjNodesList))), as.character(Sys.time())), appendLF=T)
            }
        
        }
        
    }
    
    ## append 'anno' attributes to the graph
    ind <- match(V(dag)$name, names(pAnnos))
    V(dag)$anno <- pAnnos[ind]
    
    ## append 'IC' attributes to the graph
    counts <- sapply(V(dag)$anno, length)
    IC <- -1*log10(counts/max(counts))
    ### force those 'Inf' to be 'zero'
    if(1){
    	IC[is.infinite(IC)] <- 0
    }
    V(dag)$IC <- IC
    
    return(dag)
}
