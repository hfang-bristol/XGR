#' Function to iteratively identify subnetworks from an input network and the signficance level imposed on its nodes
#'
#' \code{xSubneterGenesAdv} is supposed to iteratively identify subnetworks from an input network and the signficance level imposed on its nodes. It is an advanced version of the function \code{xSubneterGenes}. It returns a "iSubg" object
#'
#' @param data a named input vector containing the significance level for nodes (gene symbols). For this named vector, the element names are gene symbols, the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for gene symbols, 2nd column for the significance level
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways. 'REACTOME' for protein-protein interactions derived from Reactome pathways
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param seed.genes logical to indicate whether the identified network is restricted to seed genes (ie input genes with the signficant level). By default, it sets to true
#' @param subnet.size the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option
#' @param test.permutation logical to indicate whether the permutation test is perform to estimate the significance of identified network with the same number of nodes. By default, it sets to false
#' @param num.permutation the number of permutations generating the null distribution of the identified network
#' @param respect how to respect nodes to be sampled. It can be one of 'none' (randomly sampling) and 'degree' (degree-preserving sampling)
#' @param aggregateBy the aggregate method used to aggregate edge confidence p-values. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here
#' @param num.subnets the number of subnets to be iteratively identified. If NULL, all subnets will be identified until subnet.significance is no less than 0.05
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an "iSubg" object, with two components ('g' and 'ls_subg'). The 'g', a "igraph" objects for the whole network. The 'ls_subg', a list of "igraph" objects, with each element for a subgraph with a maximum score, having node attributes (significance, score, type) and a graph attribute (threshold; determined when scanning 'subnet.size'). If permutation test is enabled, it also has a graph attribute (combinedP) and an edge attribute (edgeConfidence).
#' @note The algorithm uses \code{\link{xSubneterGenes}} for identifying the first gene subnetwork from the input whole network. The second subnetwork is identified from the whole subnetwork subtracted by all edges in the first identified subnetwork. And do so till subnet.significance is no less than 0.05
#' @export
#' @seealso \code{\link{xSubneterGenes}}
#' @include xSubneterGenesAdv.r
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
#' }

xSubneterGenesAdv <- function(data, network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], network.customised=NULL, seed.genes=T, subnet.size=50, test.permutation=F, num.permutation=100, respect=c("none","degree"), aggregateBy=c("Ztransform","fishers","logistic","orderStatistic"), num.subnets=NULL, verbose=T, silent=F, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    respect <- match.arg(respect)
    aggregateBy <- match.arg(aggregateBy)
    
    subnet.size <- as.integer(subnet.size)
    if(subnet.size<=0){
    	return(NULL)
    }
    
	########################################################################
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
        g <- xDefineNet(network=network, STRING.only=STRING.only, weighted=FALSE, verbose=FALSE, RData.location=RData.location)
	}
    if(verbose){
        message(sprintf("The network you choose has %d nodes and %d edges", vcount(g),ecount(g)), appendLF=T)
    }
	########################################################################
	
	## maximum of num.subnets is 100
    if(is.null(num.subnets)){
    	num.subnets <- 100
    }
	num.subnets <- as.integer(num.subnets)
	
	ig_tmp <- g
	ls_subg <- list()
    k <- 1
    flag <- T
	while(flag && k <= num.subnets){
		subg_tmp <- xSubneterGenes(data=data, network.customised=ig_tmp, seed.genes=seed.genes, subnet.significance=NULL, subnet.size=subnet.size, test.permutation=test.permutation, num.permutation=num.permutation, respect=respect, aggregateBy=aggregateBy, verbose=verbose, silent=T, RData.location=RData.location)
		if(verbose){
			message(sprintf("Iteration %d: thresholded at %1.2e (%s) ...", k, subg_tmp$threshold, as.character(Sys.time())), appendLF=T)
		}
		## make sure threshold is less than 5e-2
		if(subg_tmp$threshold < 5e-2){
			ls_subg[[k]] <- subg_tmp
			k <- k+1
			ig_tmp <- difference(ig_tmp,subg_tmp)		
		}else{
			flag <- F
		}
	}
    
	isubg <- list(g=g, ls_subg=ls_subg, call=match.call())
    class(isubg) <- "iSubg"
    
    invisible(isubg)
}
