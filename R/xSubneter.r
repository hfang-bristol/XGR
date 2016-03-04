#' Function to identify a subnetwork from an input network and the signficance level imposed on its nodes
#'
#' \code{xSubneter} is supposed to identify maximum-scoring subnetwork from an input graph with the node information on the significance (measured as p-values or fdr). It returns an object of class "igraph". 
#'
#' @param data a named input vector containing the sinificance level for nodes (gene symbols). For this named vector, the element names are gene symbols, the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for gene symbols, 2nd column for the significance level
#' @param network the ontology supported currently. It can be "STRING" for human functional protein association network from the STRING database (version 10), "PCOMMONS_UN" for human protein interaction undirected network from Pathway Commons (version 7), "PCOMMONS_DN" for human protein interaction directed network from Pathway Commons (version 7). In order to provide source-specific directed network for "PCOMMONS_DN", it can be "PCOMMONS_DN_Reactome" for those from Reactome, "PCOMMONS_DN_KEGG" for those from KEGG, "PCOMMONS_DN_HumanCyc" for those from HumanCyc, "PCOMMONS_DN_PID" for those froom PID, "PCOMMONS_DN_PANTHER" for those from PANTHER, "PCOMMONS_DN_ReconX" for those from ReconX, "PCOMMONS_DN_TRANSFAC" for those from TRANSFAC, "PCOMMONS_DN_PhosphoSite" for those from PhosphoSite, and "PCOMMONS_DN_CTD" for those from CTD
#' @param network.confidence the parameter used to control the level of confidence of network edges that are to be further used. When "network" is "STRING", network edges with highest confidence are those with confidence scores=900, high confidence for confidence scores=700, medium confidence for confidence scores=400. When "network" is either "PCOMMONS_UN" or "PCOMMONS_DN", network edges with highest confidence are those supported with the PubMed references plus at least 3 different data sources, high confidence for those supported with the PubMed references plus at least 2 different data sources, medium confidence for those supported with at least 2 different data sources
#' @param subnet.significance the given significance threshold. By default, it is set to NULL, meaning there is no constraint on nodes/genes. If given, those nodes/genes with p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively
#' @param subnet.size the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a subgraph with a maximum score, an object of class "igraph"
#' @note The algorithm identifying a subnetwork is implemented in the dnet package (http://genomemedicine.biomedcentral.com/articles/10.1186/s13073-014-0064-8). In brief, from an input network with input node/gene information (the significant level; p-values or FDR), the way of searching for a maximum-scoring subnetwork is done as follows. Given the threshold of tolerable p-value, it gives positive scores for nodes with p-values below the threshold (nodes of interest), and negative scores for nodes with threshold-above p-values (intolerable). After score transformation, the search for a maximum scoring subnetwork is deduced to find the connected subnetwork that is enriched with positive-score nodes, allowing for a few negative-score nodes as linkers. This objective is met through minimum spanning tree finding and post-processing, previously used as a heuristic solver of prize-collecting Steiner tree problem. The solver is deterministic, only determined by the given tolerable p-value threshold. For identification of the subnetwork with a desired number of nodes, an iterative procedure is also developed to fine-tune tolerable thresholds. This explicit control over the node size may be necessary for guiding follow-up experiments.
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xSubneter.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' library(dnet)
#'
#' # a) provide the input nodes/genes with the significance info
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg')
#' sig <- rbeta(500, shape1=0.5, shape2=1)
#' data <- data.frame(symbols=org.Hs.eg$gene_info$Symbol[1:500], sig)
#' 
#' # b) perform network analysis
#' # b1) find maximum-scoring subnet based on the given significance threshold
#' subnet <- xSubneter(data=data, network="STRING", subnet.significance=0.01)
#' # b2) find maximum-scoring subnet with the desired node number=50
#' subnet <- xSubneter(data=data, network="STRING", subnet.size=50)
#'
#' # c) save subnet results to the files called 'subnet_edges.txt' and 'subnet_nodes.txt'
#' output <- igraph::get.data.frame(subnet, what="edges")
#' utils::write.table(output, file="subnet_edges.txt", sep="\t", row.names=FALSE)
#' output <- igraph::get.data.frame(subnet, what="vertices")
#' utils::write.table(output, file="subnet_nodes.txt", sep="\t", row.names=FALSE)
#'
#' # d) visualise the identified subnet
#' ## do visualisation
#' xVisNet(g=subnet, pattern=V(subnet)$score, vertex.shape="sphere")
#' 
#' # e) visualise the identified subnet as a circos plot
#' library(RCircos)
#' library(GenomicRanges)
#' xCircos(g=subnet, entity="Gene")
#' }

xSubneter <- function(data, network=c("STRING","PCOMMONS_UN","PCOMMONS_DN","PCOMMONS_DN_Reactome","PCOMMONS_DN_Reactome","PCOMMONS_DN_KEGG","PCOMMONS_DN_HumanCyc","PCOMMONS_DN_PID","PCOMMONS_DN_PANTHER","PCOMMONS_DN_ReconX","PCOMMONS_DN_TRANSFAC","PCOMMONS_DN_PhosphoSite","PCOMMONS_DN_CTD"), network.confidence=c("highest","high","medium"), subnet.significance=0.01, subnet.size=NULL, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/0.99.0")
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    network.confidence <- match.arg(network.confidence)
    
    if(is.null(data)){
        stop("The input data must be not NULL.\n")
    }
    if (is.vector(data)){
    	if(length(data)>1){
    		# assume a vector
			if(is.null(names(data))){
				stop("The input data must have names with attached gene symbols.\n")
			}
		}else{
			# assume a file
			data <- utils::read.delim(file=data, header=F, row.names=NULL, stringsAsFactors=F)
		}
    }
    if (is.vector(data)){
    	pval <- data
    }else if(is.matrix(data) | is.data.frame(data)){
        data <- as.matrix(data)
		data_list <- split(x=data[,2], f=as.character(data[,1]))
		res_list <- lapply(data_list, function(x){
			x <- as.numeric(x)
			x <- x[!is.na(x)]
			if(length(x)>0){
				min(x)
			}else{
				NULL
			}
		})
		pval <- unlist(res_list)
    }
    
    if(!is.na(network)){
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the network %s (%s) ...", network, as.character(now)), appendLF=T)
		}
        
        if(network=='STRING'){
			g <- xRDataLoader(RData.customised='org.Hs.string', RData.location=RData.location, verbose=verbose)
			## restrict to those edges with given confidence
			if(network.confidence=='highest'){
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=900])"))
			}else if(network.confidence=='high'){
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=700])"))
			}else if(network.confidence=='medium'){
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=400])"))
			}
			
        }else if(network=='PCOMMONS_UN'){
			g <- xRDataLoader(RData.customised='org.Hs.PCommons_UN', RData.location=RData.location, verbose=verbose)
			if(network.confidence=='highest'){
				# restrict to those edges with physical interactions and with highest confidence (score>=103)
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[in_complex_with>=103 | interacts_with>=103])"))
			}else if(network.confidence=='high'){
				# restrict to those edges with physical interactions and with high confidence score>=102
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[in_complex_with>=102 | interacts_with>=102])"))
			}else if(network.confidence=='medium'){
				# restrict to those edges with physical interactions and with median confidence score>=2
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[in_complex_with>=2 | interacts_with>=2])"))
			}
			
        }else if(network=='PCOMMONS_DN'){
			g <- xRDataLoader(RData.customised='org.Hs.PCommons_DN', RData.location=RData.location, verbose=verbose)
			if(network.confidence=='highest'){
				# restrict to those edges with highest confidence (score>=103)
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[catalysis_precedes>=103 | controls_expression_of>=103 | controls_phosphorylation_of>=103 | controls_state_change_of>=103 | controls_transport_of>=103])"))
			}else if(network.confidence=='high'){
				# restrict to those edges with high confidence score>=102
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[catalysis_precedes>=102 | controls_expression_of>=102 | controls_phosphorylation_of>=102 | controls_state_change_of>=102 | controls_transport_of>=102])"))
			}else if(network.confidence=='medium'){
				# restrict to those edges with median confidence score>=2
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[catalysis_precedes>=2 | controls_expression_of>=2 | controls_phosphorylation_of>=2 | controls_state_change_of>=2 | controls_transport_of>=2])"))
			}
        }else{
        	g <- xRDataLoader(RData.customised='org.Hs.PCommons_DN.source', RData.location=RData.location, verbose=verbose)
        	g <- g[[ unlist(strsplit(network,"_"))[3] ]]
			if(network.confidence=='highest'){
				# restrict to those edges with highest confidence (score>=102)
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[catalysis_precedes>=102 | controls_expression_of>=102 | controls_phosphorylation_of>=102 | controls_state_change_of>=102 | controls_transport_of>=102])"))
			}else if(network.confidence=='high'){
				# restrict to those edges with high confidence score>=101
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[catalysis_precedes>=101 | controls_expression_of>=101 | controls_phosphorylation_of>=101 | controls_state_change_of>=101 | controls_transport_of>=101])"))
			}else if(network.confidence=='medium'){
				# restrict to those edges with median confidence score>=1
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[catalysis_precedes>=1 | controls_expression_of>=1 | controls_phosphorylation_of>=1 | controls_state_change_of>=1 | controls_transport_of>=1])"))
			}
        }
	
	}else{
		stop("There is no input for the ontology.\n")
	}
    
	## further restrict the network to only nodes/genes with pval values
	ind <- match(V(g)$symbol, names(pval))
	## for extracted graph
	nodes_mapped <- V(g)$name[!is.na(ind)]
	g <- dnet::dNetInduce(g=g, nodes_query=nodes_mapped, knn=0, remove.loops=F, largest.comp=T)
	V(g)$name <- V(g)$symbol
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("start to identify subnet (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
    subnet <- dNetPipeline(g=g, pval=pval, method="customised", significance.threshold=subnet.significance, nsize=subnet.size, plot=F, verbose=verbose)

	# extract relevant info
	if(1 & class(subnet)=="igraph"){
		relations <- igraph::get.data.frame(subnet, what="edges")[,c(1,2)]
		nodes <- igraph::get.data.frame(subnet, what="vertices")
		nodes <- cbind(symbol=nodes$symbol, description=nodes$description, significance=pval[rownames(nodes)], score=nodes$score)
		if(is.directed(subnet)){
			subg <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
		}else{
			subg <- igraph::graph.data.frame(d=relations, directed=F, vertices=nodes)
		}
	}

	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("subnet has been identified (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    return(subg)
}
