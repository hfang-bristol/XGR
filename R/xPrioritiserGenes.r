#' Function to priorise genes from an input network and the weight info imposed on its nodes
#'
#' \code{xPrioritiserGenes} is supposed to prioritise genes given an input graph and a list of seed nodes. It implements Random Walk with Restart (RWR) and calculates the affinity score of all nodes in the graph to the seeds. The priority score is the affinity score. Parallel computing is also supported for Linux or Mac operating systems. It returns an object of class "pNode". 
#'
#' @param data a named input vector containing a list of seed nodes (ie gene symbols). For this named vector, the element names are seed/node names (e.g. gene symbols), the element (non-zero) values used to weight the relative importance of seeds. Alternatively, it can be a matrix or data frame with two columns: 1st column for seed/node names, 2nd column for the weight values
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathways Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), and "STRING_medium" for interactions with medium confidence (confidence scores>=400). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addtion to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for Random Walk with Restart (RWR). The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "pNode", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a matrix of nNode X 4 containing node priority information, where nNode is the number of nodes in the input graph, and the 4 columns are "name" (node names), "seed" (1 for seeds, 0 for non-seeds), "weight" (weight values),  "priority" (the priority scores that are rescaled to the range [0,1]), "rank" (ranks of the priority scores)}
#'  \item{\code{g}: an input "igraph" object}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xPrioritiserSNPs}}, \code{\link{xPrioritiser}}, \code{\link{xPrioritiserPathways}}
#' @include xPrioritiserGenes.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' library(dnet)
#' library(GenomicRanges)
#'
#' # a) provide the seed nodes/genes with the weight info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' ## get genes within 500kb away from AS GWAS lead SNPs
#' seeds.genes <- ImmunoBase$AS$genes_variants
#' ## seeds weighted according to distance away from lead SNPs
#' data <- 1- seeds.genes/500000
#'
#' # b) perform priority analysis
#' pNode <- xPrioritiserGenes(data=data, network="PCommonsDN_medium",restart=0.7)
#'
#' # c) save to the file called 'Genes_priority.txt'
#' write.table(pNode$priority, file="Genes_priority.txt", sep="\t", row.names=FALSE)
#' }

xPrioritiserGenes <- function(data, network=c("STRING_highest","STRING_high","STRING_medium","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD"), network.customised=NULL, normalise=c("laplacian","row","column","none"), restart=0.75, normalise.affinity.matrix=c("none","quantile"), parallel=TRUE, multicores=NULL, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    if(!is.null(network.customised) && class(network.customised)=="igraph"){
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the customised network (%s) ...", as.character(now)), appendLF=T)
		}
		g <- network.customised
		
	}else{
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the network %s (%s) ...", network, as.character(now)), appendLF=T)
		}
		
		if(length(grep('STRING',network,perl=T)) > 0){
			g <- xRDataLoader(RData.customised='org.Hs.string', RData.location=RData.location, verbose=verbose)
			
			## restrict to those edges with given confidence
			flag <- unlist(strsplit(network,"_"))[2]
			if(flag=='highest'){
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=900])"))
			}else if(flag=='high'){
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=700])"))
			}else if(flag=='medium'){
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=400])"))
			}
			
			V(g)$name <- V(g)$symbol
			if(0){
				relations <- igraph::get.data.frame(g, what="edges")[, c(1,2,10)]
				colnames(relations) <- c("from","to","weight")
			}else{
				relations <- igraph::get.data.frame(g, what="edges")[, c(1,2)]
				colnames(relations) <- c("from","to")
			}
			g <- igraph::graph.data.frame(d=relations, directed=F)
			
        }else if(length(grep('PCommonsUN',network,perl=T)) > 0){
			g <- xRDataLoader(RData.customised='org.Hs.PCommons_UN', RData.location=RData.location, verbose=verbose)
			
			flag <- unlist(strsplit(network,"_"))[2]
			if(flag=='high'){
				# restrict to those edges with physical interactions and with score>=102
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[in_complex_with>=102 | interacts_with>=102])"))
			}else if(flag=='medium'){
				# restrict to those edges with physical interactions and with score>=101
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[in_complex_with>=101 | interacts_with>=101])"))
			}
			
			relations <- igraph::get.data.frame(g, what="edges")[, c(1,2)]
			colnames(relations) <- c("from","to")
			g <- igraph::graph.data.frame(d=relations, directed=F)
			
        }else if(length(grep('PCommonsDN',network,perl=T)) > 0){
			flag <- unlist(strsplit(network,"_"))[2]
			if(flag=='high'){
				g <- xRDataLoader(RData.customised='org.Hs.PCommons_DN', RData.location=RData.location, verbose=verbose)
				# restrict to those edges with high confidence score>=102
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[catalysis_precedes>=102 | controls_expression_of>=102 | controls_phosphorylation_of>=102 | controls_state_change_of>=102 | controls_transport_of>=102])"))
			}else if(flag=='medium'){
				g <- xRDataLoader(RData.customised='org.Hs.PCommons_DN', RData.location=RData.location, verbose=verbose)
				# restrict to those edges with median confidence score>=101
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[catalysis_precedes>=101 | controls_expression_of>=101 | controls_phosphorylation_of>=101 | controls_state_change_of>=101 | controls_transport_of>=101])"))
			}else{
				g <- xRDataLoader(RData.customised='org.Hs.PCommons_DN.source', RData.location=RData.location, verbose=verbose)
				g <- g[[ flag ]]
				# restrict to those edges with high confidence score>=101
				eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[catalysis_precedes>=101 | controls_expression_of>=101 | controls_phosphorylation_of>=101 | controls_state_change_of>=101 | controls_transport_of>=101])"))
			}
			
			relations <- igraph::get.data.frame(g, what="edges")[, c(1,2)]
			colnames(relations) <- c("from","to")
			g <- igraph::graph.data.frame(d=relations, directed=F)
        }
	
	}
    ######################################################################################
    
    pNode <- xPrioritiser(seeds=data, g=g, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose)
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(pNode)
}
