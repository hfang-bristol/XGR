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
#' @param test.permutation logical to indicate whether the permutation test is perform to estimate the significance of identified network with the same number of nodes. By default, it sets to false
#' @param num.permutation the number of permutations generating the null distribution of the identified network
#' @param respect how to respect nodes to be sampled. It can be one of 'none' (randomly sampling) and 'degree' (degree-preserving sampling)
#' @param aggregateBy the aggregate method used to aggregate edge confidence p-values. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return
#' a subgraph with a maximum score, an object of class "igraph". It has node attributes (significance, score, type) and a graph attribute (threshold; determined when scanning 'subnet.size'). If permutation test is enabled, it also has a graph attribute (combinedP) and an edge attribute (edgeConfidence).
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

xSubneterGenes <- function(data, network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], network.customised=NULL, seed.genes=T, subnet.significance=0.01, subnet.size=NULL, test.permutation=F, num.permutation=100, respect=c("none","degree"), aggregateBy=c("Ztransform","fishers","logistic","orderStatistic"), verbose=T, silent=F, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{

    startT <- Sys.time()
    if(!silent){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }else{
    	verbose <- FALSE
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    respect <- match.arg(respect)
    aggregateBy <- match.arg(aggregateBy)
    
    if(is.null(data)){
        #stop("The input data must be not NULL.\n")
        return(NULL)
    }
    if (is.vector(data)){
    	if(length(data)>1){
    		# assume a vector
			if(is.null(names(data))){
				#stop("The input data must have names with attached gene symbols.\n")
				return(NULL)
			}
		}else{
			# assume a file
			data <- utils::read.delim(file=data, header=F, row.names=NULL, stringsAsFactors=F)
		}
    }
    if (is.vector(data)){
    	pval <- data[!is.na(data)]
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
        
        g <- xDefineNet(network=network, STRING.only=STRING.only, weighted=FALSE, verbose=FALSE, RData.location=RData.location, guid=guid)

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
    
	if(class(suppressWarnings(try(subnet <- dnet::dNetPipeline(g=g, pval=pval, method="customised", significance.threshold=subnet.significance, nsize=subnet.size, plot=F, verbose=verbose), T)))=="try-error"){
		subnet <- NULL
		
	}else{
		
		if(test.permutation & !is.null(subnet.size)){
		
			if(verbose){
				message(sprintf("Estimate the significance of the identified network (%d nodes) based on %d permutation test respecting '%s' (%s) ...", vcount(subnet), num.permutation, respect, as.character(Sys.time())), appendLF=T)
			}
			
			####################
			if(respect=='degree'){
				## equally binned
				nbin <- 10
				vec_degree <- degree(g)
				breaks <- unique(stats::quantile(vec_degree, seq(0, 1, 1/nbin)))
				cut_index <- as.numeric(cut(vec_degree, breaks=breaks))
				cut_index[is.na(cut_index)] <- 1
				names(cut_index) <- V(g)$name
				# update pval and generate pval_degree
				ind <- match(names(pval), names(cut_index))
				pval <- pval[!is.na(ind)]
				pval_degree <- cut_index[ind[!is.na(ind)]]
				
				## df_ind_B
				B <- num.permutation
				set.seed(825)
				### per node
				ls_df <- lapply(1:length(pval), function(i){
					#message(sprintf("%d (%s) ...", i, as.character(Sys.time())), appendLF=T)
					x <- pval[i]
					## all_to_sample:
					ind <- match(names(x), names(pval_degree))
					all_to_sample <- which(pval_degree == pval_degree[ind])
					## ind_sampled
					ind_sampled <- base::sample(all_to_sample, B, replace=T)
					res <- data.frame(name=names(pval[i]), pval=pval[ind_sampled], B=1:B, stringsAsFactors=F)
				})
				df_ind_B <- do.call(rbind, ls_df)
			
				# Estimate the significance of the identified subnetwork based on data permutation test
				ls_index <- split(x=df_ind_B[,c("name","pval")], f=df_ind_B$B)
				ls_subnet_permutated <- lapply(1:length(ls_index), function(j){
					if(verbose & j%%10==0){
						message(sprintf("\t%d (out of %d) (%s) ...", j, B, as.character(Sys.time())), appendLF=T)
					}
					pval_permutated <- ls_index[[j]]$pval
					names(pval_permutated) <- ls_index[[j]]$name
					# For permutated pval
					#subnet_permutated <- dnet::dNetPipeline(g=g, pval=pval_permutated, method="customised", significance.threshold=subnet.significance, nsize=subnet.size, plot=F, verbose=F)
					if(class(suppressWarnings(try(subnet_permutated <- suppressMessages(dnet::dNetPipeline(g=g, pval=pval_permutated, method="customised", significance.threshold=subnet.significance, nsize=subnet.size, plot=F, verbose=F)), T)))=="try-error"){
						return(NULL)
					}else{
						return(subnet_permutated)
					}
				})
				## Remove null elements in a list
				ls_subnet_permutated <- base::Filter(base::Negate(is.null), ls_subnet_permutated)
			
			}else{
			
				# Estimate the significance of the identified subnetwork based on data permutation test
				B <- num.permutation
				ls_subnet_permutated <- list()
				set.seed(825)
				for(j in 1:B){
					if(verbose & j%%10!=0){
						message(sprintf("\t%d (out of %d) (%s) ...", j, B, as.character(Sys.time())), appendLF=T)
					}
					ind <- base::sample(1:length(pval), replace=T)
					pval_permutated <- pval[ind]
					names(pval_permutated) <- names(pval)
					# For permutated pval
					#ls_subnet_permutated[[j]] <- dnet::dNetPipeline(g=g, pval=pval_permutated, method="customised", significance.threshold=subnet.significance, nsize=subnet.size, plot=F, verbose=F)
					if(class(suppressWarnings(try(ls_subnet_permutated[[j]] <- suppressMessages(dnet::dNetPipeline(g=g, pval=pval_permutated, method="customised", significance.threshold=subnet.significance, nsize=subnet.size, plot=F, verbose=F)), T)))=="try-error"){
						ls_subnet_permutated[[j]] <- NULL
					}
				}
				## Remove null elements in a list
				ls_subnet_permutated <- base::Filter(base::Negate(is.null), ls_subnet_permutated)
				
			}
			
			# append the confidence information from the source graphs into the target graph
			subnet <- dNetConfidence(target=subnet, sources=ls_subnet_permutated, plot=F)
			E(subnet)$edgeConfidence <- E(subnet)$edgeConfidence/100
			
			# combined P-values for aggregated/global p-value
			p_combined <- dnet::dPvalAggregate(pmatrix=matrix(E(subnet)$edgeConfidence,nrow=1), method=aggregateBy)
			subnet$combinedP <- p_combined
			
			if(verbose){
				message(sprintf("\t'%s' combined p-value (%1.2e) of the identified network (%d nodes) based on %d permutation test respecting '%s' (%s)", aggregateBy, subnet$combinedP, vcount(subnet), num.permutation, respect, as.character(Sys.time())), appendLF=T)
			}
			
		}
	}
    #subnet <- dnet::dNetPipeline(g=g, pval=pval, method="customised", significance.threshold=subnet.significance, nsize=subnet.size, plot=F, verbose=verbose)

	# extract relevant info
	#if(igraph::ecount(subnet)>0 && class(subnet)=="igraph"){
	if(class(subnet)=="igraph"){
		relations <- igraph::get.data.frame(subnet, what="edges")[,c(1,2)]
		if(!is.null(subnet$combinedP)){
			relations$edgeConfidence <- igraph::get.data.frame(subnet, what="edges")[,"edgeConfidence"]
		}
		nodes <- igraph::get.data.frame(subnet, what="vertices")
		nodes <- cbind(name=nodes$name, description=nodes$description, significance=pval[rownames(nodes)], score=nodes$score, type=nodes$type)
		#nodes <- cbind(name=nodes$name, significance=pval[rownames(nodes)], score=nodes$score)
		if(igraph::is.directed(subnet)){
			subg <- igraph::graph.data.frame(d=relations, directed=TRUE, vertices=nodes)
		}else{
			subg <- igraph::graph.data.frame(d=relations, directed=FALSE, vertices=nodes)
		}
		if(!is.null(subnet$combinedP)){
			subg$combinedP <- subnet$combinedP
		}
		subg$threshold <- subnet$threshold
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
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xSubneterGenes): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    return(subg)
}
