#' Function to conduct enrichment analysis given a list of gene sets and a list of ontologies
#'
#' \code{xEnricherGenesAdv} is supposed to conduct enrichment analysis given a list of gene sets and a list of ontologies. It is an advanced version of \code{xEnricherGenes}, returning an object of the class 'ls_eTerm'.
#'
#' @param list_vec an input vector containing gene symbols. Alternatively it can be a list of vectors, representing multiple groups of genes
#' @param background a background vector containing gene symbols as the test background. If NULL, by default all annotatable are used as background
#' @param check.symbol.identity logical to indicate whether to match the input data/background via Synonyms for those unmatchable by official gene symbols. By default, it sets to false
#' @param ontologies the ontologies supported currently. Pre-built ontology and annotation data are detailed in \code{\link{xDefineOntology}}.
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
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @param plot logical to indicate whether heatmap plot is drawn
#' @param fdr.cutoff fdr cutoff used to declare the significant terms. By default, it is set to 0.05. This option only works when setting plot (see above) is TRUE
#' @param displayBy which statistics will be used for drawing heatmap. It can be "fc" for enrichment fold change, "fdr" for adjusted p value (or FDR), "pvalue" for p value, "zscore" for enrichment z-score (by default), "or" for odds ratio. This option only works when setting plot (see above) is TRUE
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return 
#' an object of class "ls_eTerm", a list with following components:
#' \itemize{
#'  \item{\code{df}: a data frame of n x 12, where the 12 columns are "group" (the input group names), "ontology" (input ontologies), "id" (term ID), "name" (term name), "nAnno" (number in members annotated by a term), "nOverlap" (number in overlaps), "fc" (enrichment fold changes), "zscore" (enrichment z-score), "pvalue" (nominal p value), "adjp" (adjusted p value (FDR)), "or" (odds ratio), "CIl" (lower bound confidence interval for the odds ratio), "CIu" (upper bound confidence interval for the odds ratio), "distance" (term distance or other information), "members" (members (represented as Gene Symbols) in overlaps)}
#'  \item{\code{mat}: NULL if the plot is not drawn; otherwise, a matrix of term names X groups with numeric values for the signficant enrichment, NA for the insignificant ones}
#'  \item{\code{gp}: NULL if the plot is not drawn; otherwise, a 'ggplot' object}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xEnricherGenes}}, \code{\link{xEnrichViewer}}, \code{\link{xHeatmap}}
#' @include xEnricherGenesAdv.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' # Gene-based enrichment analysis using ontologies (REACTOME and GOMF)
#' # a) provide the input Genes of interest (eg 100 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' set.seed(825)
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 100))
#' data
#' 
#' # optionally, provide the test background (if not provided, all human genes)
#' #background <- as.character(org.Hs.eg$gene_info$Symbol)
#' 
#' # b) perform enrichment analysis
#' ls_eTerm <- xEnricherGenesAdv(data, ontologies=c("REACTOME","GOMF"), RData.location=RData.location)
#' ls_eTerm
#' ## forest plot of enrichment results
#' gp <- xEnrichForest(ls_eTerm, top_num=10)
#' ## heatmap plot of enrichment results
#' gp <- xEnrichHeatmap(ls_eTerm, fdr.cutoff=0.1, displayBy="or")
#' }

xEnricherGenesAdv <- function(list_vec, background=NULL, check.symbol.identity=F, ontologies=NA, size.range=c(10,2000), min.overlap=5, which.distance=NULL, test=c("fisher","hypergeo","binomial"), background.annotatable.only=NULL, p.tail=c("one-tail","two-tails"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","pc","elim","lea"), elim.pvalue=1e-2, lea.depth=2, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=F, verbose=F, silent=FALSE, plot=TRUE, fdr.cutoff=0.05, displayBy=c("zscore","fdr","pvalue","fc","or"), RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
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
    test <- match.arg(test)
    p.tail <- match.arg(p.tail)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
    path.mode <- match.arg(path.mode)
    p.tail <- match.arg(p.tail)
    displayBy <- match.arg(displayBy)
    
    ################################################
    ## support enrichment analysis for modular genes from the cModule object
    if(is(list_vec,'cModule')){
    	list_vec <- split(x=list_vec$mem$nodes, f=list_vec$mem$modules)
    }
    ################################################
        
    ############
    if(length(list_vec)==0){
    	return(NULL)
    }
    ############
    if(is.vector(list_vec) & !is(list_vec,"list")){
    	list_vec <- list(list_vec)
	}else if(is(list_vec,"list")){
		## Remove null elements in a list
		list_vec <- base::Filter(base::Negate(is.null), list_vec)
		if(length(list_vec)==0){
			return(NULL)
		}
    }else{
        stop("The input data must be a vector or a list of vectors.\n")
    }
    
	list_names <- names(list_vec)
	if(is.null(list_names)){
		list_names <- paste0('G', 1:length(list_vec))
		names(list_vec) <- list_names
	}
    
    ls_df <- lapply(1:length(list_vec), function(i){
		
		if(verbose){
			message(sprintf("Analysing group %d ('%s') (%s) ...", i, names(list_vec)[i], as.character(Sys.time())), appendLF=T)
		}
		data <- list_vec[[i]]
    	
    	ls_df <- lapply(1:length(ontologies), function(j){
			if(verbose){
				message(sprintf("\tontology %d ('%s') (%s) ...", j, ontologies[j], as.character(Sys.time())), appendLF=T)
			}
			ontology <- ontologies[j]
			
    		eTerm <- xEnricherGenes(data=data, background=background, check.symbol.identity=check.symbol.identity, ontology=ontology, size.range=size.range, min.overlap=min.overlap, which.distance=which.distance, test=test, background.annotatable.only=background.annotatable.only, p.tail=p.tail, p.adjust.method=p.adjust.method, ontology.algorithm=ontology.algorithm, elim.pvalue=elim.pvalue, lea.depth=lea.depth, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose, silent=!verbose, RData.location=RData.location, guid=guid)
			df <- xEnrichViewer(eTerm, top_num="all", sortBy="or", details=TRUE)
			
			if(is.null(df)){
				return(NULL)
			}else{
				cbind(group=rep(names(list_vec)[i],nrow(df)), ontology=rep(ontology,nrow(df)), id=rownames(df), df, stringsAsFactors=F)
			}
		})
		df <- do.call(rbind, ls_df)
	})
    df_all <- do.call(rbind, ls_df)
    ## group ordered by the input data
    df_all$group <- factor(df_all$group, levels=names(list_vec))
    
    ## heatmap view
    if(plot & !is.null(df_all)){
		
		gp <- xEnrichHeatmap(list_eTerm=df_all, fdr.cutoff=fdr.cutoff, displayBy=displayBy, colormap=NULL, zlim=NULL, reorder="none")
		mat <- gp$mat
		gp$mat <- NULL
		
    }else{
    	mat <- NULL
    	gp <- NULL
    }
    
    ls_eTerm <- list(df = df_all,
    			   mat = mat,
    			   gp = gp
                 )
    class(ls_eTerm) <- "ls_eTerm"
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xEnricherGenesAdv): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(ls_eTerm)
}
