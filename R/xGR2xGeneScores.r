#' Function to identify likely modulated seed genes from an input list of genomic regions together with the significance level given the crosslink info
#'
#' \code{xGR2xGeneScores} is supposed to identify likely modulated seed genes from a list of genomic regions (GR) together with the significance level (measured as p-values or fdr). To do so, it defines seed genes and their scores given the crosslink info with a score quantifying the link of a GR to a gene. It returns an object of class "mSeed". 
#'
#' @param data a named input vector containing the sinificance level for genomic regions (GR). For this named vector, the element names are GR, in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. The element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for GR, 2nd column for the significance level. 
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of GR into scores. If given, those GR below this are considered significant and thus scored positively. Instead, those above this are considered insignificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to 10. If NULL, no capping is applied
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. It can be one of 'genehancer' (enhancer genes; PMID:28605766) or 'nearby' (nearby genes; if so, please also specify the relevant parameters 'nearby.distance.max', 'nearby.decay.kernel' and 'nearby.decay.exponent' below)
#' @param crosslink.customised the crosslink info with a score quantifying the link of a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if nor provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param cdf.function a character specifying how to transform the input crosslink score. It can be one of 'original' (no such transformation), and 'empirical'  for looking at empirical Cumulative Distribution Function (cdf; as such it is converted into pvalue-like values [0,1])
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR (also over Contexts if many). It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param nearby.distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param nearby.decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param nearby.decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "mSeed", a list with following components:
#' \itemize{
#'  \item{\code{GR}: a matrix of nGR X 3 containing GR information, where nGR is the number of GR, and the 3 columns are "GR" (genomic regions), "Score" (the scores for GR calculated based on p-values taking into account the given threshold of the significant level), "Pval" (the input p-values for GR)}
#'  \item{\code{Gene}: a matrix of nGene X 3 containing Gene information, where nGene is the number of seed genes, and the 3 columns are "Gene" (gene symbol), "Score" (the scores for seed genes), "Pval" (pvalue-like significance level transformed from gene scores)}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note This function uses \code{\link{xGRscores}} and \code{\link{xGR2xGenes}} to define and score seed genes from input genomic regions.
#' @export
#' @seealso \code{\link{xGRscores}}, \code{\link{xGR2xGenes}}, \code{\link{xSparseMatrix}}
#' @include xGR2xGeneScores.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#'
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' df <- as.data.frame(gr, row.names=NULL)
#' GR <- paste0(df$seqnames,':',df$start,'-',df$end)
#' data <- cbind(GR=GR, Sig=df$Pvalue)
#'
#' # b) define and score seed geens
#' mSeed <- xGR2xGeneScores(data=data, crosslink="genehancer", RData.location=RData.location)
#'
#' # c) extract GR info
#' head(mSeed$GR)
#'
#' # d) extract gene info
#' head(mSeed$Gene)
#' }

xGR2xGeneScores <- function(data, significance.threshold=5e-5, score.cap=10, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("genehancer","nearby"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring.scheme=c("max","sum","sequential"), nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    build.conversion <- match.arg(build.conversion)
    crosslink <- match.arg(crosslink)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
    
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
        message(sprintf("'xGR2xGenes' is being called to define crosslinked genes (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
	df_xGenes <- xGR2xGenes(data=df_GR$GR, format="chr:start-end", build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function=cdf.function, scoring=F, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, verbose=verbose, RData.location=RData.location)
	
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xGR2xGenes' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    
    if(1){
    	ind <- match(df_xGenes$GR, df_GR$GR)
    	## Gene2GR
    	score <- df_xGenes$Score * df_GR$Score[ind]
    	Gene2GR <- data.frame(Gene=df_xGenes$Gene, GR=df_xGenes$GR, Score=score, stringsAsFactors=FALSE)
    	Gene2GR <- Gene2GR[order(Gene2GR$Gene,-Gene2GR$Score,decreasing=FALSE),]
    	
    	Gene <- Score <- NULL
    	
    	## summaryFun is applied over GR (also over Contexts if many) for each seed gene
    	
		## calculate genetic influence score under a set of SNPs for each seed gene
		if(scoring.scheme=="max"){
			summaryFun <- max
		}else if(scoring.scheme=="sum"){
			summaryFun <- sum
		}else if(scoring.scheme=="sequential"){
			summaryFun <- function(x){
				base::sum(x / base::rank(-x,ties.method="min"))
			}
		}
    	
		tmp <- as.data.frame(df_xGenes %>% dplyr::group_by(Gene) %>% dplyr::summarise(Score=summaryFun(Score)))
    	seeds.genes <- tmp$Score
    	names(seeds.genes) <- tmp$Gene
    	
    }
	
	##############
	# rescale to [0 1]
	rescaleFun <- function(x){
		(x - min(x))/(max(x) - min(x))
	}
	x <- rescaleFun(seeds.genes)
	# convert into pvalue by 10^(-x*10)
	pval <- 10^(-x*10)
	##############
	
	#############
	## for output
	df_Gene <- data.frame(Gene=names(seeds.genes), Score=seeds.genes, Pval=pval, row.names=NULL, stringsAsFactors=F)
	df_Gene <- df_Gene[order(df_Gene$Score,decreasing=TRUE),]
	#############
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("In summary, %d Genes are defined as seeds from a list of %d (out of %d) genomic regions and scored using '%s' scoring scheme", length(seeds.genes), length(unique(df_xGenes$GR)), nrow(df_GR), scoring.scheme, as.character(now)), appendLF=T)
	}
    
    df_GR <- df_GR[order(df_GR$Score,df_GR$GR,decreasing=TRUE),]
    
    mSeed <- list(GR = df_GR,
                  Gene = df_Gene,
                  Call = match.call()
                 )
    class(mSeed) <- "mSeed"
    
    invisible(mSeed)
}
