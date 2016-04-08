#' Function to identify a gene network from top prioritised genes
#'
#' \code{xPrioritiserSubnet} is supposed to identify maximum-scoring gene subnetwork from a graph with the node information on priority scores, both are part of an object of class "pNode". It returns an object of class "igraph". 
#'
#' @param pNode an object of class "pNode"
#' @param priority.quantite the quantite of the top priority genes. By default, 10% of top prioritised genes will be used for network analysis. If NULL or NA, all prioritised genes will be used
#' @param network the built-in network. If NULL, the network used for prioritisation will be used, which is part of the object of class "pNode". Otherwise, choose the other network of interest. Currently two sources of network information are supported: the STRING database (version 10) and the Pathways Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), and "STRING_medium" for interactions with medium confidence (confidence scores>=400). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addtion to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param subnet.significance the given significance threshold. By default, it is set to NULL, meaning there is no constraint on nodes/genes. If given, those nodes/genes with p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively
#' @param subnet.size the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a subgraph with a maximum score, an object of class "igraph". It has ndoe attributes: signficance, score, priority (part of the "pNode" object)
#' @note The priority score will be first scaled to the range x=[0 100] and then is converted to pvalue-like significant level: 10^(-x). Next, \code{\link{xSubneterGenes}} is used to identify a maximum-scoring gene subnetwork that contains as many highly prioritised genes as possible but a few lowly prioritised genes as linkers.
#' @export
#' @seealso \code{\link{xSubneterGenes}}
#' @include xPrioritiserSubnet.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' library(dnet)
#' library(ggbio)
#'
#' # a) provide the SNPs with the significance info
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' AS <- read.delim(file.path(path.package("XGR"),"AS.txt"), stringsAsFactors=FALSE)
#'
#' # b) perform priority analysis
#' pNode <- xPrioritiserSNPs(data=AS, network="PCommonsUN_medium",restart=0.7)
#' 
#' # c) perform network analysis
#' # find maximum-scoring subnet with the desired node number=50
#' subnet <- xPrioritiserSubnet(pNode, priority.quantite=0.1, subnet.size=50)
#'
#' # d) save subnet results to the files called 'subnet_edges.txt' and 'subnet_nodes.txt'
#' output <- igraph::get.data.frame(subnet, what="edges")
#' utils::write.table(output, file="subnet_edges.txt", sep="\t", row.names=FALSE)
#' output <- igraph::get.data.frame(subnet, what="vertices")
#' utils::write.table(output, file="subnet_nodes.txt", sep="\t", row.names=FALSE)
#'
#' # e) visualise the identified subnet
#' ## do visualisation with nodes colored according to the priority
#' xVisNet(g=subnet, pattern=V(subnet)$priority, vertex.shape="sphere")
#' ## do visualisation with nodes colored according to pvalue-like signficance
#' xVisNet(g=subnet, pattern=-log10(as.numeric(V(subnet)$significance)), vertex.shape="sphere", colormap="wyr")
#' 
#' # f) visualise the identified subnet as a circos plot
#' library(RCircos)
#' xCircos(g=subnet, entity="Gene")
#' }

xPrioritiserSubnet <- function(pNode, priority.quantite=0.1, network=c(NULL,"STRING_highest","STRING_high","STRING_medium","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD"), network.customised=NULL, subnet.significance=0.01, subnet.size=NULL, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    
    if(is.null(network)){
    	network.customised <- pNode$g
    }
    
    if (class(pNode) != "pNode" ){
        stop("The function must apply to a 'pNode' object.\n")
    }
    
    priority <- pNode$priority$priority
    names(priority) <- rownames(pNode$priority)
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("The 'pNode' object contains %d prioritised genes", length(pval)), appendLF=T)
	}
    
    ## scale to the range [0 100] and then convert to pvalue-like signficant level
	x <- priority
	y <- (x - min(x,na.rm=T)) / (max(x,na.rm=T) - min(x,na.rm=T))
	pval <- 10^(-100*y)
	
	## only keep the top priority (qunatile)
	## priority quantite
	priority.quantite <- as.numeric(priority.quantite)
	if(length(priority.quantite>0 & priority.quantite<1) & !is.na(priority.quantite)){
		cf <- quantile(pval, priority.quantite, na.rm=T)
		ind <- which(pval<cf)
		pval <- pval[ind]
		
		priority <- priority[ind]
	}
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Among, %d genes are used for network analysis", length(pval)), appendLF=T)
	}
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\t maximum priority: %1.2e; minimum priority: %1.2e", max(priority), min(priority)), appendLF=T)
		message(sprintf("\t minimum p-value: %1.2e; maximum p-value: %1.2e", min(pval), max(pval)), appendLF=T)
	}
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("xSubneterGenes is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
    subg <- xSubneterGenes(data=pval, network=network, network.customised=network.customised, seed.genes=T, subnet.significance=subnet.significance, subnet.size=subnet.size, verbose=verbose, RData.location=RData.location)
	
	# extract relevant info
	if(ecount(subg)>0 && class(subg)=="igraph"){
		relations <- igraph::get.data.frame(subg, what="edges")[,c(1,2)]
		nodes <- igraph::get.data.frame(subg, what="vertices")
		nodes <- cbind(symbol=nodes$name, description=nodes$description, significance=nodes$significance, score=nodes$score, priority=priority[rownames(nodes)])
		if(is.directed(subg)){
			subg <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
		}else{
			subg <- igraph::graph.data.frame(d=relations, directed=F, vertices=nodes)
		}
	}else{
		subg <- NULL
	}
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("xSubneterGenes has finished (%s)!", as.character(now)), appendLF=T)
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
