#' Function to conduct region-based enrichment analysis via crosslinked genes
#'
#' \code{xGR2xGeneAnno} is supposed to conduct region-based enrichment analysis for the input genomic region data (genome build h19), using crosslinked gene annotations. To do so, crosslinked genes are first defined. Currently supported built-in crosslink info is enhancer genes, eQTL genes, conformation genes and nearby genes (purely), though the user can customise it via 'crosslink.customised'; if so, it has priority over the built-in data. Enrichment analysis is then based on either Fisher's exact test or Hypergeometric test for estimating the significance of overlapped crosslinked genes. Test background can be provided; by default, the annotatable genes will be used. 
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param background an input background containing a list of genomic regions as the test background. The file format is the same as 'data' above. By default, it is NULL meaning all annotatable genes are used as background
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. See \code{\link{xGR2xGenes}} for details
#' @param crosslink.customised the crosslink info with a score quantifying the link of a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if not provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param crosslink.top the number of the top genes defined by 'data' will be used for test. By default, it is NULL
#' @param nearby.distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param nearby.decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param nearby.decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param ontology the ontology supported currently. By default, it is 'NA' to disable this option. Pre-built ontology and annotation data are detailed in \code{\link{xDefineOntology}}.
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 2000
#' @param min.overlap the minimum number of overlaps. Only those terms with members that overlap with input data at least min.overlap (3 by default) will be processed
#' @param which.distance which terms with the distance away from the ontology root (if any) is used to restrict terms in consideration. By default, it sets to 'NULL' to consider all distances
#' @param test the test statistic used. It can be "fisher" for using fisher's exact test, "hypergeo" for using hypergeometric test, or "binomial" for using binomial test. Fisher's exact test is to test the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement). Hypergeometric test is to sample at random (without replacement) from the background containing annotated and non-annotated genes, and thus compare sampling to background. Unlike hypergeometric test, binomial test is to sample at random (with replacement) from the background with the constant probability. In terms of the ease of finding the significance, they are in order: hypergeometric test > fisher's exact test > binomial test. In other words, in terms of the calculated p-value, hypergeometric test < fisher's exact test < binomial test
#' @param background.annotatable.only logical to indicate whether the background is further restricted to the annotatable. By default, it is NULL: if ontology.algorithm is not 'none', it is always TRUE; otherwise, it depends on the background (if not provided, it will be TRUE; otherwise FALSE). Surely, it can be explicitly stated
#' @param p.tail the tail used to calculate p-values. It can be either "two-tails" for the significance based on two-tails (ie both over- and under-overrepresentation)  or "one-tail" (by default) for the significance based on one tail (ie only over-representation)
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param ontology.algorithm the algorithm used to account for the hierarchy of the ontology. It can be one of "none", "pc", "elim" and "lea". For details, please see 'Note' below
#' @param elim.pvalue the parameter only used when "ontology.algorithm" is "elim". It is used to control how to declare a signficantly enriched term (and subsequently all genes in this term are eliminated from all its ancestors)
#' @param lea.depth the parameter only used when "ontology.algorithm" is "lea". It is used to control how many maximum depth is used to consider the children of a term (and subsequently all genes in these children term are eliminated from the use for the recalculation of the signifance at this term)
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' an object of class "eTerm", a list with following components:
#' \itemize{
#'  \item{\code{term_info}: a matrix of nTerm X 4 containing snp/gene set information, where nTerm is the number of terms, and the 4 columns are "id" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
#'  \item{\code{annotation}: a list of terms containing annotations, each term storing its annotations. Always, terms are identified by "id"}
#'  \item{\code{g}: an igraph object to represent DAG}
#'  \item{\code{data}: a vector containing input data in consideration. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{background}: a vector containing the background data. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{overlap}: a list of overlapped snp/gene sets, each storing snps overlapped between a snp/gene set and the given input data (i.e. the snps of interest). Always, gene sets are identified by "id"}
#'  \item{\code{fc}: a vector containing fold changes}
#'  \item{\code{zscore}: a vector containing z-scores}
#'  \item{\code{pvalue}: a vector containing p-values}
#'  \item{\code{adjp}: a vector containing adjusted p-values. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{or}: a vector containing odds ratio}
#'  \item{\code{CIl}: a vector containing lower bound confidence interval for the odds ratio}
#'  \item{\code{CIu}: a vector containing upper bound confidence interval for the odds ratio}
#'  \item{\code{cross}: a matrix of nTerm X nTerm, with an on-diagnal cell for the overlapped-members observed in an individaul term, and off-diagnal cell for the overlapped-members shared betwene two terms}
#'  \item{\code{call}: the call that produced this result}
#'  \item{\code{evidence}: a data frame with 3 columns ('GR' for genomic regions, 'Gene' for crosslinked genes, and 'Score' for the score between the gene and the GR)}
#'  \item{\code{gp_evidence}: a ggplot object for evidence data}
#' }
#' @note The interpretation of the algorithms used to account for the hierarchy of the ontology is:
#' \itemize{
#' \item{"none": does not consider the ontology hierarchy at all.}
#' \item{"lea": computers the significance of a term in terms of the significance of its children at the maximum depth (e.g. 2). Precisely, once snps are already annotated to any children terms with a more signficance than itself, then all these snps are eliminated from the use for the recalculation of the signifance at that term. The final p-values takes the maximum of the original p-value and the recalculated p-value.}
#' \item{"elim": computers the significance of a term in terms of the significance of its all children. Precisely, once snps are already annotated to a signficantly enriched term under the cutoff of e.g. pvalue<1e-2, all these snps are eliminated from the ancestors of that term).}
#' \item{"pc": requires the significance of a term not only using the whole snps as background but also using snps annotated to all its direct parents/ancestors as background. The final p-value takes the maximum of both p-values in these two calculations.}
#' \item{"Notes": the order of the number of significant terms is: "none" > "lea" > "elim" > "pc".}
#' }
#' @export
#' @seealso \code{\link{xGR}}, \code{\link{xGR2xGenes}}, \code{\link{xEnricherGenes}}
#' @include xGR2xGeneAnno.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' 
#' # 1) provide the genomic regions
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' names(gr) <- NULL
#' dGR <- xGR(gr, format="GRanges")
#' 
#' ## b) perform DO enrichment analysis
#' ## enhancer genes
#' eTerm <- xGR2xGeneAnno(data=dGR, format="GRanges", crosslink="genehancer", ontology="DO", RData.location=RData.location)
#' ## nearby genes (50kb, decaying rapidly)
#' eTerm <- xGR2xGeneAnno(data=dGR, format="GRanges", crosslink="nearby", ontology="DO", nearby.distance.max=50000, nearby.decay.kernel="rapid", RData.location=RData.location)
#'
#' ## c) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' ## d) save enrichment results to the file called 'Regions2genes_enrichments.txt'
#' output <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' utils::write.table(output, file="Regions2genes_enrichments.txt", sep="\t", row.names=FALSE)
#' 
#' ## e) barplot of significant enrichment results
#' bp <- xEnrichBarplot(eTerm, top_num=10, displayBy="fc")
#' print(bp)
#' 
#' ## f) forest of significant enrichment results
#' gp <- xEnrichForest(eTerm, top_num=10)
#' }

xGR2xGeneAnno <- function(data, background=NULL, format=c("data.frame", "bed", "chr:start-end", "GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("genehancer","PCHiC_combined","GTEx_V6p_combined","nearby"), crosslink.customised=NULL, crosslink.top=NULL, nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, ontology=NA, size.range=c(10,2000), min.overlap=3, which.distance=NULL, test=c("hypergeo","fisher","binomial"), background.annotatable.only=NULL, p.tail=c("one-tail","two-tails"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","pc","elim","lea"), elim.pvalue=1e-2, lea.depth=2, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=F, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    #crosslink <- match.arg(crosslink)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
    test <- match.arg(test)
    p.tail <- match.arg(p.tail)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
    path.mode <- match.arg(path.mode)
    
    ###################
	if(verbose){
		now <- Sys.time()
		message(sprintf("First, import the data/background formatted as '%s' (%s) ...", format, as.character(now)), appendLF=T)
	}
    
	dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
	bGR <- xGR(data=background, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
	
	#####################################
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Second, define crosslinked genes based on '%s' (%s) ...", crosslink, as.character(now)), appendLF=T)
	}
    
    df_xGenes_data <- xGR2xGenes(data=dGR, format="GRanges", crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function="original", scoring=TRUE, scoring.scheme="max", scoring.rescale=F, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, verbose=verbose, RData.location=RData.location)
    df_xGenes_background <- xGR2xGenes(data=bGR, format="GRanges", crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function="original", scoring=TRUE, scoring.scheme="max", scoring.rescale=F, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, verbose=verbose, RData.location=RData.location)
	
	##############################
   	Score <- Gene <- NULL
    
    ## dGR_genes
    df_xGenes_data <- df_xGenes_data %>% dplyr::arrange(-Score)
	if(is.null(crosslink.top)){
		crosslink.top <- nrow(df_xGenes_data)
	}
	if(crosslink.top > nrow(df_xGenes_data)){
		crosslink.top <- nrow(df_xGenes_data)
	}
	crosslink.top <- as.integer(crosslink.top)
	crosslink.cutoff <- df_xGenes_data[crosslink.top,'Score']
    dGR_genes <- df_xGenes_data$Gene[df_xGenes_data$Score >= crosslink.cutoff]

    ## bGR_genes
	if(!is.null(df_xGenes_background)){
		bGR_genes <- (df_xGenes_background %>% dplyr::arrange(-Score))$Gene
	}else{
		bGR_genes <- NULL
	}
	
	if(verbose){
		if(is.null(bGR_genes)){
			message(sprintf("\t%d (out of %d crosslinked genes) are used.", length(dGR_genes), nrow(df_xGenes_data), as.character(Sys.time())), appendLF=T)
		}else{
			message(sprintf("\t%d (out of %d crosslinked genes) and %d background genes are used.", length(dGR_genes), nrow(df_xGenes_data), length(bGR_genes), as.character(Sys.time())), appendLF=T)
		}
	}
	
	#######################################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xEnricherGenes' is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
	eTerm <- xEnricherGenes(data=dGR_genes, background=bGR_genes, ontology=ontology, size.range=size.range, min.overlap=min.overlap, which.distance=which.distance, test=test, background.annotatable.only=background.annotatable.only, p.tail=p.tail, p.adjust.method=p.adjust.method, ontology.algorithm=ontology.algorithm, elim.pvalue=elim.pvalue, lea.depth=lea.depth, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose, RData.location=RData.location)
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xEnricherGenes' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    if(!is.null(eTerm)){
		######
		## append 'evidence'
		df_evidence <- xGR2xGenes(data=dGR, format="GRanges", crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function="original", scoring=FALSE, scoring.scheme="max", scoring.rescale=F, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, verbose=verbose, RData.location=RData.location)
		ind <- match(df_evidence$Gene, dGR_genes)
		evidence <- df_evidence[!is.na(ind), c('GR','Gene','Score')]
		eTerm$evidence <- evidence
		## append 'gp_evidence'
		Gene <- Score <- NULL
		mat_evidence <- tidyr::spread(evidence, key=Gene, value=Score)
		mat <- mat_evidence[,-1]
		rownames(mat) <- mat_evidence[,1]
		#### sort by chromosome, start and end
		ind <- xGRsort(rownames(mat))
		mat <- mat[ind,]
		
		################
		## obtain rowsep
		rowsep <- xGRsep(rownames(mat))
		rowsep <- nrow(mat) - rowsep
		################
				
		####
		if(ncol(mat)>=0){
			reorder <- "none"
		}else{
			reorder <- "col"
		}
		gp_evidence <- xHeatmap(mat, reorder=reorder, colormap="spectral", ncolors=64, barwidth=0.4, x.rotate=90, shape=19, size=2, x.text.size=6,y.text.size=6, na.color='transparent')
		gp_evidence <- gp_evidence + theme(legend.title=element_text(size=8), legend.position="left") + scale_y_discrete(position="right")
		gp_evidence <- gp_evidence + geom_hline(yintercept=rowsep+0.5,color="grey90",size=0.5)
		eTerm$gp_evidence <- gp_evidence
		###### 
    }
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(eTerm)
}
