#' Function to identify likely modulated seed genes given a list of genomic regions together with the significance level
#'
#' \code{xGR2GeneScores} is supposed to identify likely modulated seed genes from a list of genomic regions (GR) together with the significance level (measured as p-values or fdr). To do so, it defines seed genes and their scores that take into account the distance to and the significance of input SNPs. It returns an object of class "mSeed". 
#'
#' @param data a named input vector containing the sinificance level for genomic regions (GR). For this named vector, the element names are GR, in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. The element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for GR, 2nd column for the significance level. 
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of GR into scores. If given, those GR below this are considered significant and thus scored positively. Instead, those above this are considered insignificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to 10. If NULL, no capping is applied
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "mSeed", a list with following components:
#' \itemize{
#'  \item{\code{GR}: a matrix of nGR X 3 containing GR information, where nGR is the number of GR, and the 3 columns are "GR" (genomic regions), "Score" (the scores for GR calculated based on p-values taking into account the given threshold of the significant level), "Pval" (the input p-values for GR)}
#'  \item{\code{Gene}: a matrix of nGene X 3 containing Gene information, where nGene is the number of seed genes, and the 3 columns are "Gene" (gene symbol), "Score" (the scores for seed genes), "Pval" (pvalue-like significance level transformed from gene scores)}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note This function uses \code{\link{xGRscores}} and \code{\link{xGR2nGenes}} to define and score nearby genes that are located within distance window of input genomic regions.
#' @export
#' @seealso \code{\link{xGRscores}}, \code{\link{xGR2nGenes}}, \code{\link{xSparseMatrix}}
#' @include xGR2GeneScores.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
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
#' # b) define and score seed geens
#' mSeed <- xGR2GeneScores(data=data, RData.location=RData.location)
#' }

xGR2GeneScores <- function(data, significance.threshold=5e-5, score.cap=10, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), distance.max=50000, decay.kernel=c("slow","linear","rapid","constant"), decay.exponent=2, GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), scoring.scheme=c("max","sum","sequential"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    build.conversion <- match.arg(build.conversion)
    decay.kernel <- match.arg(decay.kernel)
    scoring.scheme <- match.arg(scoring.scheme)
    
  ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xGRscores' is being called to score GR (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
	df_GR <- xGRscores(data=data, significance.threshold=significance.threshold, score.cap=score.cap, verbose=verbose, RData.location=RData.location)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xGRscores' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xGR2nGenes' is being called to define nearby genes (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
	df_nGenes <- xGR2nGenes(data=df_GR$GR, format="chr:start-end", build.conversion=build.conversion, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.Gene=GR.Gene, scoring=F, verbose=verbose, RData.location=RData.location)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xGR2nGenes' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    
    if(1){
    	ind <- match(df_nGenes$GR, df_GR$GR)
    	## Gene2GR
    	score <- df_nGenes$Weight * df_GR$Score[ind]
    	Gene2GR <- data.frame(Gene=df_nGenes$Gene, GR=df_nGenes$GR, Score=score, stringsAsFactors=FALSE)
    	Gene2GR <- Gene2GR[order(Gene2GR$Gene,-Gene2GR$Score,decreasing=FALSE),]
    	
    	ls_gene <- split(x=Gene2GR$Score, f=Gene2GR$Gene)
		## calculate genetic influence score under a set of SNPs for each seed gene
		if(scoring.scheme=='max'){
			seeds.genes <- sapply(ls_gene, max)
		}else if(scoring.scheme=='sum'){
			seeds.genes <- sapply(ls_gene, sum)
		}else if(scoring.scheme=='sequential'){
			seeds.genes <- sapply(ls_gene, function(x){
				base::sum(x / base::rank(-x,ties.method="min"))
			})
		}
    	
    }else{
		## df_GR df_nGenes
		allGenes <- sort(df_nGenes$Gene)
		allGR <- sort(df_GR$GR)
	
		## sparse matrix of nGenes X GR
		G2S_n <- xSparseMatrix(df_nGenes[,c("Gene","GR","Weight")], rows=allGenes, columns=allGR, verbose=verbose)
		G2S <- G2S_n
	
		## consider GR scores
		ind <- match(colnames(G2S), df_GR$GR)
		########
		df_GR <- df_GR[ind,]
		########
		GR_score <- df_GR$Score
		names(GR_score) <- colnames(G2S)
		## convert into matrix
		mat_GR_score <- matrix(rep(GR_score,each=nrow(G2S)), nrow=nrow(G2S))
	
		## calculate genetic influence score for a gene-GR pair
		G2S_score <- G2S * mat_GR_score
    	
		## calculate genetic influence score under a set of SNPs for each seed gene
		if(scoring.scheme=='max'){
			seeds.genes <- apply(G2S_score, 1, function(x) {
				base::max(x)
			})
		}else if(scoring.scheme=='sum'){
			seeds.genes <- apply(G2S_score, 1, function(x) {
				base::sum(x)
			})
		}else if(scoring.scheme=='sequential'){
			seeds.genes <- apply(G2S_score, 1, function(x) {
				base::sum(base::sort(x, decreasing=T) / (1:length(x)))
			})
		}

    }
	
    ################################
    x <- seeds.genes
	## take back to the p-value format
    pval <- 10^(-1*x)
    ################################
	
	#############
	## for output
	df_Gene <- data.frame(Gene=names(seeds.genes), Score=seeds.genes, Pval=pval, row.names=NULL, stringsAsFactors=F)
	df_Gene <- df_Gene[order(df_Gene$Score,decreasing=TRUE),]
	#############
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("In summary, %d Genes are defined as seeds and scored using '%s' scoring scheme", length(seeds.genes), scoring.scheme, as.character(now)), appendLF=T)
	}
    
    df_GR <- df_GR[order(df_GR$Score,df_GR$GR,decreasing=TRUE),]
    
    mSeed <- list(GR = df_GR,
                  Gene = df_Gene,
                  Call = match.call()
                 )
    class(mSeed) <- "mSeed"   
    
    invisible(mSeed)
}
