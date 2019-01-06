#' Function to identify a gene network from an input network given a list of genomic regions
#'
#' \code{xGR2xNet} is supposed to identify maximum-scoring gene subnetwork from an input graph with the node information (genomic regions with or without the significance). To do so, it defines seed genes and their scores that take into account the distance to and the significance of input genomic regions (GR). It returns an object of class "igraph". 
#'
#' @param data a named input vector containing the significance level for genomic regions (GR). For this named vector, the element names are GR, in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20', the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for GR, 2nd column for the significance level. Also supported is the input with GR only (without the significance level)
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of GR into scores. If given, those GR below this are considered significant and thus scored positively. Instead, those above this are considered insignificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to NULL, meaning that no capping is applied
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. See \code{\link{xGR2xGenes}} for details
#' @param crosslink.customised the crosslink info with a score quantifying the link of a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if nor provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param cdf.function a character specifying how to transform the input crosslink score. It can be one of 'original' (no such transformation), and 'empirical'  for looking at empirical Cumulative Distribution Function (cdf; as such it is converted into pvalue-like values [0,1])
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR (also over Contexts if many). It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param nearby.distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param nearby.decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param nearby.decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways. 'REACTOME' for protein-protein interactions derived from Reactome pathways
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param seed.genes logical to indicate whether the identified network is restricted to seed genes (ie nearby genes that are located within defined distance window centred on lead or LD SNPs). By default, it sets to true
#' @param subnet.significance the given significance threshold. By default, it is set to NULL, meaning there is no constraint on nodes/genes. If given, those nodes/genes with p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively
#' @param subnet.size the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option
#' @param test.permutation logical to indicate whether the permutation test is perform to estimate the significance of identified network with the same number of nodes. By default, it sets to false
#' @param num.permutation the number of permutations generating the null distribution of the identified network
#' @param respect how to respect nodes to be sampled. It can be one of 'none' (randomly sampling) and 'degree' (degree-preserving sampling)
#' @param aggregateBy the aggregate method used to aggregate edge confidence p-values. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a subgraph with a maximum score, an object of class "igraph". It has graph attributes (evidence, gp_evidence) and node attributes (significance, score, type). If permutation test is enabled, it also has a graph attribute (combinedP) and an edge attribute (edgeConfidence).
#' @note The algorithm identifying a gene subnetwork that is likely modulated by input genomic regions (GR) includes two major steps. The first step is to use \code{\link{xGR2xGeneScores}} for defining and scoring nearby genes that are located within distance window of input GR. The second step is to use \code{\link{xSubneterGenes}} for identifying a maximum-scoring gene subnetwork that contains as many highly scored genes as possible but a few less scored genes as linkers.
#' @export
#' @seealso \code{\link{xGR2xGeneScores}}, \code{\link{xSubneterGenes}}
#' @include xGR2xNet.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' # a) provide the seed SNPs with the significance info
#' data(ImmunoBase)
#' ## only AS GWAS SNPs and their significance info (p-values)
#' df <- as.data.frame(ImmunoBase$AS$variant, row.names=NULL)
#' GR <- paste0(df$seqnames,':',df$start,'-',df$end)
#' data <- cbind(GR=GR, Sig=df$Pvalue)
#' 
#' # b) perform network analysis
#' # b1) find maximum-scoring subnet based on the given significance threshold
#' subnet <- xGR2xNet(data=data, crosslink="genehancer", network="STRING_high", seed.genes=F, subnet.significance=0.01, RData.location=RData.location)
#' # b2) find maximum-scoring subnet with the desired node number=30
#' subnet <- xGR2xNet(data=data, crosslink="genehancer", network="STRING_high", seed.genes=F, subnet.size=30, RData.location=RData.location)
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
#' xCircos(g=subnet, entity="Gene", colormap="orange-darkred", ideogram=F, entity.label.side="out", chr.exclude=NULL, RData.location=RData.location)
#' }

xGR2xNet <- function(data, significance.threshold=NULL, score.cap=NULL, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("genehancer","PCHiC_combined","GTEx_V6p_combined","nearby"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring.scheme=c("max","sum","sequential"), nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), network.customised=NULL, seed.genes=T, subnet.significance=5e-5, subnet.size=NULL, test.permutation=F, num.permutation=100, respect=c("none","degree"), aggregateBy=c("Ztransform","fishers","logistic","orderStatistic"), verbose=T, silent=F, RData.location="http://galahad.well.ox.ac.uk/bigdata")
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
    #crosslink <- match.arg(crosslink)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
    network <- match.arg(network)
    
	if(is.vector(data) && length(data)>1 && is.null(names(data)) && is.character(data)){
		if(verbose){
			message(sprintf("The input has GRs only (without the significance level)"), appendLF=TRUE)
		}
		
		####################################################################################
		if(verbose){
			now <- Sys.time()
			message(sprintf("\n#######################################################", appendLF=T))
			message(sprintf("'xGR2xGenes' is being called to score seed genes (%s):", as.character(now)), appendLF=T)
			message(sprintf("#######################################################", appendLF=T))
		}
		df_xGenes <- xGR2xGenes(data=data, format="chr:start-end", build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function=cdf.function, scoring=T, scoring.scheme=scoring.scheme, scoring.rescale=F, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, verbose=verbose, silent=!verbose, RData.location=RData.location)
		#length(unique(df_xGenes$Gene))
		if(verbose){
			now <- Sys.time()
			message(sprintf("#######################################################", appendLF=T))
			message(sprintf("'xGR2xGenes' has been finished (%s)!", as.character(now)), appendLF=T)
			message(sprintf("#######################################################\n", appendLF=T))
		}
		
		pval <- df_xGenes$Pval
		names(pval) <- df_xGenes$Gene
		
		#########
		df_evidence <- xGR2xGenes(data=data, format="chr:start-end", build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function=cdf.function, scoring=F, scoring.scheme=scoring.scheme, scoring.rescale=F, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, verbose=verbose, silent=!verbose, RData.location=RData.location)
		#length(unique(df_evidence$Gene))
		#########
		
	}else{
	
		####################################################################################
		if(verbose){
			now <- Sys.time()
			message(sprintf("\n#######################################################", appendLF=T))
			message(sprintf("'xGR2xGeneScores' is being called to score seed genes (%s):", as.character(now)), appendLF=T)
			message(sprintf("#######################################################", appendLF=T))
		}
		mSeed <- xGR2xGeneScores(data=data, significance.threshold=significance.threshold, score.cap=score.cap, build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function=cdf.function, scoring.scheme=scoring.scheme, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, verbose=verbose, RData.location=RData.location)
		if(verbose){
			now <- Sys.time()
			message(sprintf("#######################################################", appendLF=T))
			message(sprintf("'xGR2xGeneScores' has been finished (%s)!", as.character(now)), appendLF=T)
			message(sprintf("#######################################################\n", appendLF=T))
		}
	
		df_Gene <- mSeed$Gene
		pval <- df_Gene$Pval
		names(pval) <- df_Gene$Gene
	
		#########
		df_evidence <- mSeed$Link
		#########
	}
    
	####################################################################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("\t\t minimum p-value: %1.2e; maximum p-value: %1.2e", min(pval), max(pval)), appendLF=T)
	}
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("xSubneterGenes is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
    subg <- xSubneterGenes(data=pval, network=network, network.customised=network.customised, seed.genes=seed.genes, subnet.significance=subnet.significance, subnet.size=subnet.size, test.permutation=test.permutation, num.permutation=num.permutation, respect=respect, aggregateBy=aggregateBy, verbose=verbose, silent=!verbose, RData.location=RData.location)
	
	######
	## append graph attribute 'evidence'
	ind <- match(df_evidence$Gene, V(subg)$name)
	evidence <- df_evidence[!is.na(ind), c('GR','Gene','Score')]
	subg$evidence <- evidence
	## append graph attribute 'gp_evidence'
	Gene <- Score <- NULL
	mat_evidence <- tidyr::spread(evidence, key=Gene, value=Score)
	mat <- mat_evidence[,-1]
	rownames(mat) <- mat_evidence[,1]
	#### sort by chromosome, start and end
	ind <- xGRsort(rownames(mat))
	mat <- mat[ind,]
	####
	if(ncol(mat)>=0){
		reorder <- "none"
	}else{
		reorder <- "col"
	}
	gp_evidence <- xHeatmap(mat, reorder=reorder, colormap="spectral", ncolors=64, barwidth=0.4, x.rotate=90, shape=19, size=2, x.text.size=6,y.text.size=6, na.color='transparent')
	gp_evidence <- gp_evidence + theme(legend.title=element_text(size=8), legend.position="left") + scale_y_discrete(position="right")
	subg$gp_evidence <- gp_evidence
	######
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("xSubneterGenes has finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xGR2xNet): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    return(subg)
}
