#' Function to define a gene network
#'
#' \code{xDefineNet} is supposed to define a gene network sourced from the STRING database or the Pathway Commons database. It returns an object of class "igraph". 
#'
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param weighted logical to indicate whether edge weights should be considered. By default, it sets to false. If true, it only works for the network from the STRING database 
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "igraph"
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xDefineNet.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' g <- xDefineNet(network="PCommonsDN_medium", RData.location=RData.location)
#' }

xDefineNet <- function(network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], weighted=FALSE, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load the network %s (%s) ...", network, as.character(now)), appendLF=TRUE)
	}
		
	if(length(grep('STRING',network,perl=TRUE)) > 0){
		g <- xRDataLoader(RData.customised='org.Hs.string', RData.location=RData.location, verbose=verbose)
		
		## restrict to those edges with given confidence
		flag <- unlist(strsplit(network,"_"))[2]
		if(flag=='highest'){
			eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=900])"))
		}else if(flag=='high'){
			eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=700])"))
		}else if(flag=='medium'){
			eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=400])"))
		}else if(flag=='low'){
			eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=150])"))
		}
		
		## further restricted by the evidence type
		default.STRING.only <- c("neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")
		ind <- match(default.STRING.only, STRING.only)
		STRING.only <- default.STRING.only[!is.na(ind)]
		if(length(STRING.only)>0){
			x <- sapply(STRING.only, function(x) paste0(x,'>0'))
			x <- paste0(x, collapse=' | ')
			x <- paste0("g <- igraph::subgraph.edges(g, eids=E(g)[",x,"])")
			eval(parse(text=x))
		}
		
		########################
		# because of the way storing the network from the STRING database
		## extract relations (by symbol)
		V(g)$name <- V(g)$symbol
		if(weighted){
			relations <- igraph::get.data.frame(g, what="edges")[, c(1,2,10)]
			colnames(relations) <- c("from","to","weight")
		}else{
			relations <- igraph::get.data.frame(g, what="edges")[, c(1,2)]
			colnames(relations) <- c("from","to")
			relations$weight <- rep(1, nrow(relations))
		}
		## do removal for node extraction (without 'name'; otherwise failed to do so using the function 'igraph::get.data.frame')
		g <- igraph::delete_vertex_attr(g, "name")
		g <- igraph::delete_vertex_attr(g, "seqid")
		g <- igraph::delete_vertex_attr(g, "geneid")
		nodes <- igraph::get.data.frame(g, what="vertices")
		### remove the duplicated
		nodes <- nodes[!duplicated(nodes), ]			
		########################
		
		g <- igraph::graph.data.frame(d=relations, directed=FALSE, vertices=nodes)
			
    }else if(length(grep('PCommonsUN',network,perl=TRUE)) > 0){
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
		nodes <- igraph::get.data.frame(g, what="vertices")[, c(3,4)]
		g <- igraph::graph.data.frame(d=relations, directed=FALSE, vertices=nodes)
		
    }else if(length(grep('PCommonsDN',network,perl=TRUE)) > 0){
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
		relations$weight <- rep(1, nrow(relations))
		nodes <- igraph::get.data.frame(g, what="vertices")[, c(3,4)]
		g <- igraph::graph.data.frame(d=relations, directed=FALSE, vertices=nodes)
    }
    
    invisible(g)
}
