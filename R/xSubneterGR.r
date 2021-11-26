#' Function to identify a gene network from an input network given a list of genomic regions together with the significance level
#'
#' \code{xSubneterGR} is supposed to identify maximum-scoring gene subnetwork from an input graph with the node information on the significance (measured as p-values or fdr). To do so, it defines seed genes and their scores that take into account the distance to and the significance of input genomic regions (GR). It returns an object of class "igraph". 
#'
#' @param data a named input vector containing the sinificance level for genomic regions (GR). For this named vector, the element names are GR, in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. The element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for GR, 2nd column for the significance level. 
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of GR into scores. If given, those GR below this are considered significant and thus scored positively. Instead, those above this are considered insigificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to 10. If NULL, no capping is applied
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways. 'REACTOME' for protein-protein interactions derived from Reactome pathways
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param seed.genes logical to indicate whether the identified network is restricted to seed genes (ie nearby genes that are located within defined distance window centred on lead or LD SNPs). By default, it sets to true
#' @param subnet.significance the given significance threshold. By default, it is set to NULL, meaning there is no constraint on nodes/genes. If given, those nodes/genes with p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively
#' @param subnet.size the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option
#' @param test.permutation logical to indicate whether the permutation test is perform to estimate the significance of identified network with the same number of nodes. By default, it sets to false
#' @param num.permutation the number of permutations generating the null distribution of the identified network
#' @param respect how to respect nodes to be sampled. It can be one of 'none' (randomly sampling) and 'degree' (degree-preserving sampling)
#' @param aggregateBy the aggregate method used to aggregate edge confidence p-values. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here
#' @param num.subnets the number of subnets to be iteratively identified. If NA, this functionality is disabled. If NULL, all subnets will be identified until subnet.significance is no less than 0.05
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return
#' IF num.subnets is NA, a subgraph with a maximum score, an object of class "igraph". It has ndoe attributes: significance, score, type. If permutation test is enabled, it also has a graph attribute (combinedP) and an edge attribute (edgeConfidence)
#' IF num.subnets is not NA (also subnet.size is not NULL), an "iSubg" object, with two components ('g' and 'ls_subg'). The 'g', a "igraph" objects for the whole network. The 'ls_subg', a list of "igraph" objects, with each element for a subgraph with a maximum score, having node attributes (significance, score, type) and a graph attribute (threshold; determined when scanning 'subnet.size'). If permutation test is enabled, it also has a graph attribute (combinedP) and an edge attribute (edgeConfidence).
#' @note The algorithm identifying a gene subnetwork that is likely modulated by input genomic regions (GR) includes two major steps. The first step is to use \code{\link{xGR2GeneScores}} for defining and scoring nearby genes that are located within distance window of input GR. The second step is to use \code{\link{xSubneterGenes}} for identifying a maximum-scoring gene subnetwork that contains as many highly scored genes as possible but a few less scored genes as linkers. Also supported at the second step is to identify multiple subnetwoks using \code{\link{xSubneterGenesAdv}} when num.subnets is not NA.
#' @export
#' @seealso \code{\link{xGR2GeneScores}}, \code{\link{xSubneterGenes}}, \code{\link{xSubneterGenesAdv}}
#' @include xSubneterGR.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#'
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' df <- as.data.frame(gr, row.names=NULL)
#' chr <- df$seqnames
#' start <- df$start
#' end <- df$end
#' sig <- df$Pvalue
#' GR <- paste(chr,':',start,'-',end, sep='')
#' data <- cbind(GR=GR, Sig=sig)
#' 
#' # b) perform network analysis
#' # b1) find maximum-scoring subnet based on the given significance threshold
#' subnet <- xSubneterGR(data=data, network="STRING_high", seed.genes=F, subnet.significance=0.01, RData.location=RData.location)
#' # b2) find maximum-scoring subnet with the desired node number=30
#' subnet <- xSubneterGR(data=data, network="STRING_high", seed.genes=F, subnet.size=30, RData.location=RData.location)
#'
#' # c) save subnet results to the files called 'subnet_edges.txt' and 'subnet_nodes.txt'
#' output <- igraph::get.data.frame(subnet, what="edges")
#' utils::write.table(output, file="subnet_edges.txt", sep="\t", row.names=FALSE)
#' output <- igraph::get.data.frame(subnet, what="vertices")
#' utils::write.table(output, file="subnet_nodes.txt", sep="\t", row.names=FALSE)
#'
#' # d) visualise the identified subnet
#' ## do visualisation with nodes colored according to the significance
#' xVisNet(g=subnet, pattern=-log10(as.numeric(V(subnet)$significance)), vertex.shape="sphere", colormap="wyr")
#' ## do visualisation with nodes colored according to transformed scores
#' xVisNet(g=subnet, pattern=as.numeric(V(subnet)$score), vertex.shape="sphere")
#' 
#' # e) visualise the identified subnet as a circos plot
#' library(RCircos)
#' xCircos(g=subnet, entity="Gene", colormap="white-gray", RData.location=RData.location)
#' }

xSubneterGR <- function(data, significance.threshold=5e-5, score.cap=10, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), distance.max=50000, decay.kernel=c("slow","linear","rapid","constant"), decay.exponent=2, GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), 
scoring.scheme=c("max","sum","sequential"), network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], network.customised=NULL, seed.genes=T, subnet.significance=5e-5, subnet.size=NULL, test.permutation=F, num.permutation=100, respect=c("none","degree"), aggregateBy=c("Ztransform","fishers","logistic","orderStatistic"), num.subnets=NA, verbose=T, silent=F, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
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
    build.conversion <- match.arg(build.conversion)
    decay.kernel <- match.arg(decay.kernel)
    scoring.scheme <- match.arg(scoring.scheme)
    network <- match.arg(network)
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xGR2GeneScores' is being called to score seed genes (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
    
    mSeed <- xGR2GeneScores(data=data, significance.threshold=significance.threshold, score.cap=score.cap, build.conversion=build.conversion, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.Gene=GR.Gene, scoring.scheme=scoring.scheme, verbose=verbose, RData.location=RData.location, guid=guid)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xGR2GeneScores' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
	df_Gene <- mSeed$Gene
	pval <- df_Gene$Pval
	names(pval) <- df_Gene$Gene
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\t\t minimum p-value: %1.2e; maximum p-value: %1.2e", min(pval), max(pval)), appendLF=T)
	}
    
    #############################################################################################
    flag <- F
    if(!is.null(num.subnets)){
    	if(is.na(num.subnets)){
    		flag <- T
    	}
    }
    
    if(flag){
		if(verbose){
			now <- Sys.time()
			message(sprintf("\n#######################################################", appendLF=T))
			message(sprintf("xSubneterGenes is being called (%s):", as.character(now)), appendLF=T)
			message(sprintf("#######################################################", appendLF=T))
		}
		subg <- xSubneterGenes(data=pval, network=network, STRING.only=STRING.only, network.customised=network.customised, seed.genes=seed.genes, subnet.significance=subnet.significance, subnet.size=subnet.size, test.permutation=test.permutation, num.permutation=num.permutation, respect=respect, aggregateBy=aggregateBy, verbose=verbose, silent=!verbose, RData.location=RData.location, guid=guid)
		if(verbose){
			now <- Sys.time()
			message(sprintf("#######################################################", appendLF=T))
			message(sprintf("xSubneterGenes has finished (%s)!", as.character(now)), appendLF=T)
			message(sprintf("#######################################################\n", appendLF=T))
		}
    }else{
    	if(!is.null(subnet.size)){
			if(verbose){
				now <- Sys.time()
				message(sprintf("\n#######################################################", appendLF=T))
				message(sprintf("xSubneterGenesAdv is being called (%s):", as.character(now)), appendLF=T)
				message(sprintf("#######################################################", appendLF=T))
			}
			subg <- xSubneterGenesAdv(data=pval, network=network, STRING.only=STRING.only, network.customised=network.customised, seed.genes=seed.genes, subnet.size=subnet.size, test.permutation=test.permutation, num.permutation=num.permutation, respect=respect, aggregateBy=aggregateBy, num.subnets=num.subnets, verbose=verbose, silent=!verbose, RData.location=RData.location, guid=guid)
			if(verbose){
				now <- Sys.time()
				message(sprintf("#######################################################", appendLF=T))
				message(sprintf("xSubneterGenesAdv has finished (%s)!", as.character(now)), appendLF=T)
				message(sprintf("#######################################################\n", appendLF=T))
			}
    	}else{
    		stop("The function requires subnet.size is NULL when num.subnets is not NA.\n")
    	}
    }
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xSubneterGR): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    return(subg)
}
