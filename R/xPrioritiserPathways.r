#' Function to prioritise pathways based on GSEA analysis of prioritised genes
#'
#' \code{xPrioritiserPathways} is supposed to prioritise pathways given prioritised genes and the ontology in query. It returns an object of class "eTerm". It is done via gene set enrichment analysis (GSEA). 
#'
#' @param pNode an object of class "pNode"
#' @param ontology the ontology supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "PS" for phylostratific age information, "PS2" for the collapsed PS version (inferred ancestors being collapsed into one with the known taxonomy information), "SF" for domain superfamily assignments, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPCM" for Human Phenotype Clinical Modifier, "HPMA" for Human Phenotype Mortality Aging, "MP" for Mammalian Phenotype, and Drug-Gene Interaction database (DGIdb) for drugable categories, and the molecular signatures database (Msigdb, including "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7")
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 2000
#' @param which.distance which terms with the distance away from the ontology root (if any) is used to restrict terms in consideration. By default, it sets to 'NULL' to consider all distances
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to false
#' @param nperm the number of random permutations. For each permutation, gene-score associations will be permutated so that permutation of gene-term associations is realised
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a data frame of nTerm x 10, where nTerm is the number of terms, and the 10 columns:
#' \itemize{
#'  \item{\code{setID}: Term ID}
#'  \item{\code{setSize}: Number of members in Term}
#'  \item{\code{name}: Term Name}
#'  \item{\code{es}: enrichment score}
#'  \item{\code{nes}: normalised enrichment score. It is the version of enrichment score but after being normalised by gene set size}
#'  \item{\code{pvalue}: nominal p value}
#'  \item{\code{adjp}: adjusted p value. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{fdr}: false discovery rate (FDR). It is the estimated probability that the normalised enrichment score represents a false positive finding}
#'  \item{\code{qvalue}: q value. It is the monotunically increasing FDR}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xPrioritiserSNPs}}, \code{\link{xPrioritiserGenes}}
#' @include xPrioritiserPathways.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' library(dnet)
#' library(GenomicRanges)
#'
#' RData.location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/0.99.0"
#' # a) provide the seed nodes/genes with the weight info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get genes within 500kb away from AS GWAS lead SNPs
#' seeds.genes <- ImmunoBase$AS$genes_variants
#' ## seeds weighted according to distance away from lead SNPs
#' data <- 1- seeds.genes/500000
#'
#' # b) perform priority analysis
#' pNode <- xPrioritiserGenes(data=data, network="PCommonsDN_medium",restart=0.7, RData.location=RData.location)
#' 
#' # c) derive pathway-level priority
#' res <- xPrioritiserPathways(pNode=pNode, ontology="MsigdbC2REACTOME", RData.location=RData.location)
#'
#' # d) save enrichment results to the file called 'Pathways_priority.txt'
#' write.table(res, file="Pathways_priority.txt", sep="\t", row.names=FALSE)
#' }

xPrioritiserPathways <- function(pNode, ontology=c("GOBP","GOMF","GOCC","PS","PS2","SF","DO","HPPA","HPMI","HPCM","HPMA","MP", "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7", "DGIdb"), size.range=c(10,2000), which.distance=NULL, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=F, nperm=1000, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/0.99.0")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    ontology <- match.arg(ontology)
    path.mode <- match.arg(path.mode)
    
    if (class(pNode) != "pNode" ){
        stop("The function must apply to a 'pNode' object.\n")
    }
    
    ## load Enterz Gene information
	EG <- xRDataLoader(RData=paste('org.Hs.eg', sep=''), RData.location=RData.location, verbose=verbose)	
	allGeneID <- EG$gene_info$GeneID
	allSymbol <- as.vector(EG$gene_info$Symbol)
	name <- rownames(pNode$priority)
	ind <- match(name, allSymbol)
	name_GeneID <- allGeneID[ind]
    
    ## remove those zero priority
    score <- pNode$priority[,'priority']
    ### restrict to those seeds
    seed <- pNode$priority[,'seed']
    #ind <- which(score>0 & seed==1)
    ind <- 1:max(which(seed==1))
    data <- data.frame(priority=score[ind], row.names=name_GeneID[ind])
    
    #############################################################################################

    if(!is.na(ontology)){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the ontology %s and its gene annotations (%s) ...", ontology, as.character(now)), appendLF=T)
		}
                                  
		#########
		## load GS information
		## flag the simplified version of PS
		flag_PS2 <- FALSE
		if(ontology=="PS2"){
			flag_PS2 <- TRUE
			ontology <- "PS"
		}
		GS <- xRDataLoader(RData=paste('org.Hs.eg', ontology, sep=''), RData.location=RData.location, verbose=verbose)
		
		################
		if(flag_PS2){
			tmp <- as.character(unique(GS$set_info$name))
			inds <- sapply(tmp,function(x) which(GS$set_info$name==x))
		
			## new set_info
			set_info <- data.frame()
			for(i in 1:length(inds)){
				set_info<- rbind(set_info,as.matrix(GS$set_info[max(inds[[i]]),]))
			}
			## new gs
			gs <- list()
			for(i in 1:length(inds)){
				gs[[i]] <- unlist(GS$gs[inds[[i]]], use.names=F)
			}
			names(gs) <- rownames(set_info)
		
			## new GS
			GS$set_info <- set_info
			GS$gs <- gs
		}
		################
		
		#########
		## get annotation information
		anno <- GS$gs

		#########
		## get ontology information
		## check the eligibility for the ontology
		all.ontologies <- c("GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPCM","HPMA","MP")
		flag_ontology <- ontology %in% all.ontologies
    	
    	if(flag_ontology){
			g <- xRDataLoader(RData=paste('ig.', ontology, sep=''), RData.location=RData.location, verbose=verbose)
		}else{
		
			nodes <- data.frame(name=as.character(GS$set_info$setID), term_id=as.character(GS$set_info$setID), term_name=as.character(GS$set_info$name), term_distance=as.character(GS$set_info$distance), stringsAsFactors=F)
			nodes <- rbind(nodes, c('root','root','root','root'))
			relations <- data.frame(from='root', to=nodes$name)
			g <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
		}
	
	}else{
		stop("There is no input for the ontology.\n")
	}
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("Enrichment analysis is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
    eTerm <- dGSEA(data=data, identity="entrez", check.symbol.identity=F, ontology="Customised", customised.genesets=anno, sizeRange=size.range, which_distance=which.distance, weight=1, nperm=nperm, fast=T, sigTail="two-tails", p.adjust.method="BH", verbose=verbose, RData.location=RData.location)
	
	res <- dGSEAview(eTerm, which_sample=1, top_num=NULL, sortBy="nES", decreasing=T, details=TRUE)
	term_name <- V(g)$term_name
	ind <- match(rownames(res), V(g)$name)
	res$name <- term_name[ind]
	## only those with nES>0
	res <- res[res$nES>0, c(1,11,10,2,3,4,5,8,9)]
	rownames(res) <- NULL

	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("Enrichment analysis has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(res)
}
