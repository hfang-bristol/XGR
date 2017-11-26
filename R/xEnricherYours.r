#' Function to conduct enrichment analysis given YOUR own input data
#'
#' \code{xEnricherYours} is supposed to conduct enrichment analysis given the input data and the ontology in query. It returns an object of class "eTerm". Enrichment analysis is based on either Fisher's exact test or Hypergeometric test.
#'
#' @param data.file an input data file, containing a list of entities (e.g. genes or SNPs) to test. The entities can be anything, for example, in this file \url{http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt}, the entities are InterPro domains (InterPro). As seen in this example, entries in the first column must be domains. If the file also contains other columns, these additional columns will be ignored. Alternatively, the data.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param annotation.file an input annotation file containing annotations between entities and ontology terms. For example, a file containing annotations between InterPro domains and GO Molecular Function (GOMF) terms can be found in \url{http://dcgor.r-forge.r-project.org/data/InterPro/Domain2GOMF.txt}. As seen in this example, the input file must contain two columns: 1st column for domains, 2nd column for ontology terms. If there are additional columns, these columns will be ignored. Alternatively, the annotation.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param background.file an input background file containing a list of entities as the test background. The file format is the same as 'data.file'. By default, it is NULL meaning all annotatable entities (i.g. those entities in 'annotation.file') are used as background
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 2000
#' @param min.overlap the minimum number of overlaps. Only those terms with members that overlap with input data at least min.overlap (3 by default) will be processed
#' @param test the test statistic used. It can be "fisher" for using fisher's exact test, "hypergeo" for using hypergeometric test, or "binomial" for using binomial test. Fisher's exact test is to test the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement). Hypergeometric test is to sample at random (without replacement) from the background containing annotated and non-annotated genes, and thus compare sampling to background. Unlike hypergeometric test, binomial test is to sample at random (with replacement) from the background with the constant probability. In terms of the ease of finding the significance, they are in order: hypergeometric test > fisher's exact test > binomial test. In other words, in terms of the calculated p-value, hypergeometric test < fisher's exact test < binomial test
#' @param background.annotatable.only logical to indicate whether the background is further restricted to annotatable genes (covered by 'annotation.file'). In other words, if the background is provided, the background genes are those after being overlapped with annotatable genes. Notably, if only one annotation is provided in 'annotation.file', it should be false
#' @param p.tail the tail used to calculate p-values. It can be either "two-tails" for the significance based on two-tails (ie both over- and under-overrepresentation)  or "one-tail" (by default) for the significance based on one tail (ie only over-representation)
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
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
#' }
#' @note None
#' @export
#' @seealso \code{\link{xEnricher}}
#' @include xEnricherYours.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' 
#' # Enrichment analysis using your own data
#' # a) provide your own data (eg InterPro domains and their annotations by GO terms)
#' ## All InterPro domains
#' input.file <- "http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt"
#' data <- utils::read.delim(input.file, header=F, row.names=NULL, stringsAsFactors=F)[,1]
#' ## provide the input domains of interest (eg 100 randomly chosen domains)
#' data.file <- sample(data, 100)
#' ## InterPro domains annotated by GO Molecular Function (GOMF) terms
#' annotation.file <- "http://dcgor.r-forge.r-project.org/data/InterPro/Domain2GOMF.txt"
#' 
#' # b) perform enrichment analysis
#' eTerm <- xEnricherYours(data.file=data.file, annotation.file=annotation.file)
#'
#' # c) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' # d) save enrichment results to the file called 'Yours_enrichments.txt'
#' output <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' utils::write.table(output, file="Yours_enrichments.txt", sep="\t", row.names=FALSE)
#' 
#' # e) barplot of significant enrichment results
#' bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="adjp")
#' print(bp)
#' 
#' # Using ImmunoBase SNPs and associations/annotations with disease traits
#' ## get ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' ## get disease associated variants/SNPs
#' variants_list <- lapply(ImmunoBase, function(x) cbind(SNP=names(x$variants), Disease=rep(x$disease,length(x$variants))))
#' ## extract annotations as a data frame: Variant Disease_Name 
#' annotation.file <- do.call(rbind, variants_list)
#' head(annotation.file)
#' ## provide the input SNPs of interest
#' ## for example, cis-eQTLs induced by interferon gamma
#' cis <- xRDataLoader(RData.customised='JKscience_TS2A')
#' data.file <- matrix(cis[which(cis$IFN_t>0),c('variant')], ncol=1)
#' # perform enrichment analysis
#' eTerm <- xEnricherYours(data.file=data.file, annotation.file=annotation.file)
#' # view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#' # barplot of significant enrichment results
#' bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="adjp")
#' print(bp)
#' }

xEnricherYours <- function(data.file, annotation.file, background.file=NULL, size.range=c(10,2000), min.overlap=3, test=c("hypergeo","fisher","binomial"), background.annotatable.only=T, p.tail=c("one-tail","two-tails"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), verbose=T, silent=FALSE)
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
    p.adjust.method <- match.arg(p.adjust.method)
    p.tail <- match.arg(p.tail)
    
    ############
    if(length(data.file)==0){
    	return(NULL)
    }
    ############
    
    ###################
    ## import data file
    if(is.matrix(data.file) | is.data.frame(data.file)){
        data <- unique(data.file[,1])
    }else if(!is.null(data.file) & any(!is.na(data.file))){
    	if(length(data.file)==1){
    	
        	if(class(suppressWarnings(try(data <- utils::read.delim(file=data.file, header=F, row.names=NULL, stringsAsFactors=F), T)))=="try-error"){
        		data <- data.file
        	}else{
        		data <- unique(data[,1])
        	}
    	
			#data <- utils::read.delim(file=data.file, header=F, row.names=NULL, stringsAsFactors=F)
			#data <- unique(data[,1])
		}else{
			data <- data.file
		}
    }else{
    	stop("The file 'data.file' must be provided!\n")
    }
    
	## define background information
    if(is.matrix(background.file) | is.data.frame(background.file)){
        background <- unique(background.file[,1])
    }else if(!is.null(background.file)){
    
    	if(length(background.file)==1){
			background <- utils::read.delim(file=background.file, header=F, row.names=NULL, stringsAsFactors=F)
			background <- unique(background[,1])
		}else{
			background <- background.file
		}
    }else{
    	background <- background.file
    }
    
    ## import annotation file
    if(is.matrix(annotation.file) | is.data.frame(annotation.file)){
        input <- cbind(annotation.file[,1], annotation.file[,2])
    }else if(!is.null(annotation.file) & any(!is.na(annotation.file))){
		input <- utils::read.delim(file=annotation.file, header=F, row.names=NULL, stringsAsFactors=F)
    }else{
    	stop("The file 'annotation.file' must be provided!\n")
    }
    
    ###################################
    # only works when "background.file" is not NULL and "background.annotatable.only" is false
	if(background.annotatable.only==FALSE){
		if(!is.null(background.file)){
			tmp <- cbind(background.file, rep('BG',length(background.file)))
			input <- rbind(input, tmp)
		}
	}
    ###################################    
    ## define annotation information
	anno <- split(x=input[,1], f=input[,2])
    
	## define ontology information (artificially)
	terms <- names(anno)
	nodes <- data.frame(name=terms, term_id=terms, term_name=terms, term_distance=rep(1,length(terms)), stringsAsFactors=F)
	root <- c('Root', 'Root', 'Root', 0)
	nodes <- rbind(nodes, root)
	relations <- data.frame(from='Root', to=nodes$name)
	g <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)


    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xEnricher' is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    eTerm <- xEnricher(data=data, annotation=anno, g=g, background=background, size.range=size.range, min.overlap=min.overlap, test=test, p.tail=p.tail, p.adjust.method=p.adjust.method, ontology.algorithm="none",true.path.rule=F, verbose=verbose)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xEnricher' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(eTerm)
}
