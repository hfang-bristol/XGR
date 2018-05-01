#' Function to identify a subnetwork from an input network and the signficance level imposed on its nodes
#'
#' \code{xSubneterGenes} is supposed to identify maximum-scoring subnetwork from an input graph with the node information on the significance (measured as p-values or fdr). It returns an object of class "igraph". 
#'
#' @param data a named input vector containing the significance level for nodes (gene symbols). For this named vector, the element names are gene symbols, the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for gene symbols, 2nd column for the significance level
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways. 'REACTOME' for protein-protein interactions derived from Reactome pathways
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param seed.genes logical to indicate whether the identified network is restricted to seed genes (ie input genes with the signficant level). By default, it sets to true
#' @param subnet.significance the given significance threshold. By default, it is set to NULL, meaning there is no constraint on nodes/genes. If given, those nodes/genes with p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively
#' @param subnet.size the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a subgraph with a maximum score, an object of class "igraph". It has ndoe attributes: significance, score, type
#' @note The algorithm identifying a subnetwork is implemented in the dnet package (http://genomemedicine.biomedcentral.com/articles/10.1186/s13073-014-0064-8). In brief, from an input network with input node/gene information (the significant level; p-values or FDR), the way of searching for a maximum-scoring subnetwork is done as follows. Given the threshold of tolerable p-value, it gives positive scores for nodes with p-values below the threshold (nodes of interest), and negative scores for nodes with threshold-above p-values (intolerable). After score transformation, the search for a maximum scoring subnetwork is deduced to find the connected subnetwork that is enriched with positive-score nodes, allowing for a few negative-score nodes as linkers. This objective is met through minimum spanning tree finding and post-processing, previously used as a heuristic solver of prize-collecting Steiner tree problem. The solver is deterministic, only determined by the given tolerable p-value threshold. For identification of the subnetwork with a desired number of nodes, an iterative procedure is also developed to fine-tune tolerable thresholds. This explicit control over the node size may be necessary for guiding follow-up experiments.
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xDefineNet}}
#' @include xSubneterGenes.r
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
#' # b) perform network analysis
#' # b1) find maximum-scoring subnet based on the given significance threshold
#' subnet <- xSubneterGenes(data=data, network="STRING_high", subnet.significance=0.01, RData.location=RData.location)
#' # b2) find maximum-scoring subnet with the desired node number=50
#' subnet <- xSubneterGenes(data=data, network="STRING_high", subnet.size=50, RData.location=RData.location)
#'
#' # c) save subnet results to the files called 'subnet_edges.txt' and 'subnet_nodes.txt'
#' output <- igraph::get.data.frame(subnet, what="edges")
#' utils::write.table(output, file="subnet_edges.txt", sep="\t", row.names=FALSE)
#' output <- igraph::get.data.frame(subnet, what="vertices")
#' utils::write.table(output, file="subnet_nodes.txt", sep="\t", row.names=FALSE)
#'
#' # d) visualise the identified subnet
#' ## do visualisation with nodes colored according to the significance (you provide)
#' xVisNet(g=subnet, pattern=-log10(as.numeric(V(subnet)$significance)), vertex.shape="sphere", colormap="wyr")
#' ## do visualisation with nodes colored according to transformed scores
#' xVisNet(g=subnet, pattern=as.numeric(V(subnet)$score), vertex.shape="sphere")
#' 
#' # e) visualise the identified subnet as a circos plot
#' library(RCircos)
#' xCircos(g=subnet, entity="Gene", colormap="white-gray", RData.location=RData.location)
#' 
#' # g) visualise the subnet using the same layout_with_kk
#' df_tmp <- df[match(V(subnet)$name,df$Symbol),]
#' vec_tmp <- colnames(df_tmp)
#' names(vec_tmp) <- vec_tmp
#' glayout <- igraph::layout_with_kk(subnet)
#' V(subnet)$xcoord <- glayout[,1]
#' V(subnet)$xcoord <- glayout[,2]
#' # g1) colored according to FDR
#' ls_ig <- lapply(vec_tmp, function(x){
#' 	ig <- subnet
#' 	V(ig)$fdr <- -log10(as.numeric(df_tmp[,x]))
#' 	ig
#' })
#' gp_FDR <- xA2Net(g=ls_g, node.label='name', node.label.size=2, node.label.color='blue', node.label.alpha=0.8, node.label.padding=0.25, node.label.arrow=0, node.label.force=0.1, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='fdr', node.color.title=expression(-log[10]('FDR')), colormap='grey-yellow-orange', ncolors=64, zlim=c(0,3), node.size.range=4, edge.color="black",edge.color.alpha=0.3,edge.curve=0.1,edge.arrow.gap=0.025)
#' # g2) colored according to FC
#' ls_ig <- lapply(vec_tmp, function(x){
#' 	ig <- subnet
#' 	V(ig)$lfc <- as.numeric(df_tmp[,x])
#' 	ig
#' })
#' gp_FC <- xA2Net(g=ls_g, node.label='name', node.label.size=2, node.label.color='blue', node.label.alpha=0.8, node.label.padding=0.25, node.label.arrow=0, node.label.force=0.1, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='lfc', node.color.title=expression(log[2]('FC')), colormap='cyan1-grey-pink1', ncolors=64, zlim=c(-3,3),  node.size.range=4, edge.color="black",edge.color.alpha=0.3,edge.curve=0.1,edge.arrow.gap=0.025)
#' # g3) colored according to FC
#' gridExtra::grid.arrange(grobs=list(gp_FDR, gp_FC), ncol=2, as.table=TRUE)
#' }

xSubneterGenes <- function(data, network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], network.customised=NULL, seed.genes=T, subnet.significance=0.01, subnet.size=NULL, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    
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
		## keep the miminum p-values per gene
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
	
	if(seed.genes){
		## further restrict the network to only nodes/genes with pval values
		ind <- match(V(g)$name, names(pval))
		## for extracted graph
		nodes_mapped <- V(g)$name[!is.na(ind)]
		g <- dnet::dNetInduce(g=g, nodes_query=nodes_mapped, knn=0, remove.loops=F, largest.comp=T)
	}else{
		ind <- match(V(g)$name, names(pval))
		nodes_not_mapped <- V(g)$name[is.na(ind)]
		pval_not_mapped <- rep(1, length(nodes_not_mapped))
		names(pval_not_mapped) <- nodes_not_mapped
		pval <- c(pval, pval_not_mapped)
	}
	
	    
    if(verbose){
        message(sprintf("Restricted to data/nodes of interest, the network (with the largest interconnected component) has %d nodes and %d edges", vcount(g),ecount(g)), appendLF=T)
    }
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("Start to identify a subnetwork (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
    subnet <- dNetPipeline(g=g, pval=pval, method="customised", significance.threshold=subnet.significance, nsize=subnet.size, plot=F, verbose=verbose)

	# extract relevant info
	if(ecount(subnet)>0 && class(subnet)=="igraph"){
		relations <- igraph::get.data.frame(subnet, what="edges")[,c(1,2)]
		nodes <- igraph::get.data.frame(subnet, what="vertices")
		nodes <- cbind(name=nodes$name, description=nodes$description, significance=pval[rownames(nodes)], score=nodes$score, type=nodes$type)
		#nodes <- cbind(name=nodes$name, significance=pval[rownames(nodes)], score=nodes$score)
		if(is.directed(subnet)){
			subg <- igraph::graph.data.frame(d=relations, directed=TRUE, vertices=nodes)
		}else{
			subg <- igraph::graph.data.frame(d=relations, directed=FALSE, vertices=nodes)
		}
	}else{
		subg <- NULL
	}

	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("The subnetwork has been identified (%s)!", as.character(now)), appendLF=T)
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
