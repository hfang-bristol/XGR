#' Function to identify likely modulated seed genes given a list of SNPs together with the significance level (e.g. GWAS reported p-values)
#'
#' \code{xSNP2GeneScores} is supposed to identify likely modulated seed genes from a list of SNPs together with the significance level (measured as p-values or fdr). To do so, it defines seed genes and their scores that take into account the distance to and the significance of input SNPs. It returns an object of class "mSeed". 
#'
#' @param data a named input vector containing the sinificance level for nodes (dbSNP). For this named vector, the element names are dbSNP ID (or in the format such as 'chr16:28525386'), the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for dbSNP, 2nd column for the significance level
#' @param include.LD additional SNPs in LD with Lead SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, LD SNPs will be included based on one or more of 26 populations and 5 super populations from 1000 Genomics Project data (phase 3). The population can be one of 5 super populations ("AFR", "AMR", "EAS", "EUR", "SAS"), or one of 26 populations ("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"). Explanations for population code can be found at \url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
#' @param LD.customised a user-input matrix or data frame with 3 columns: 1st column for Lead SNPs, 2nd column for LD SNPs, and 3rd for LD r2 value. It is designed to allow the user analysing their precalcuated LD info. This customisation (if provided) has the high priority over built-in LD SNPs
#' @param LD.r2 the LD r2 value. By default, it is 0.8, meaning that SNPs in LD (r2>=0.8) with input SNPs will be considered as LD SNPs. It can be any value from 0.8 to 1
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of SNPs into scores. If given, those SNPs below this are considered significant and thus scored positively. Instead, those above this are considered insigificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to 10. If NULL, no capping is applied
#' @param distance.max the maximum distance between genes and SNPs. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby SNPs per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param include.TAD TAD boundary regions are also included. By default, it is 'NA' to disable this option. Otherwise, inclusion of a TAD dataset to pre-filter SNP-nGene pairs (i.e. only those within a TAD region will be kept). TAD datasets can be one of "GM12878"  (lymphoblast), "IMR90" (fibroblast), "MSC" (mesenchymal stem cell) ,"TRO" (trophoblasts-like cell), "H1" (embryonic stem cell), "MES" (mesendoderm) and "NPC" (neural progenitor cell). Explanations can be found at \url{http://dx.doi.org/10.1016/j.celrep.2016.10.061}
#' @param scoring.scheme the method used to calculate seed gene scores under a set of SNPs. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "mSeed", a list with following components:
#' \itemize{
#'  \item{\code{SNP}: a matrix of nSNP X 4 containing SNP information, where nSNP is the number of SNPs, and the 3 columns are "SNP" (Lead and/or LD SNPs), "Score" (the scores for SNPs calculated based on p-values taking into account the given threshold of the significant level), "Pval" (the input p-values for Lead SNPs or R2-adjusted p-values for LD SNPs), "Flag" (indicating  as Lead or LD SNPs)}
#'  \item{\code{Gene}: a matrix of nGene X 3 containing Gene information, where nGene is the number of seed genes, and the 3 columns are "Gene" (gene symbol), "Score" (the scores for seed genes), "Pval" (pvalue-like significance level transformed from gene scores)}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note This function uses \code{\link{xSNPscores}} and \code{\link{xSNP2nGenes}} to define and score nearby genes that are located within distance window of input and/or LD SNPs.
#' @export
#' @seealso \code{\link{xSNPscores}}, \code{\link{xSNP2nGenes}}, \code{\link{xSparseMatrix}}
#' @include xSNP2GeneScores.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#'
#' \dontrun{
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' data <- GenomicRanges::mcols(gr)[,c(1,3)]
#' 
#' # b) define and score seed geens
#' mSeed <- xSNP2GeneScores(data=data, include.TAD="GM12878", RData.location=RData.location)
#'
#' # c) extract SNP info
#' head(mSeed$SNP)
#'
#' # d) extract gene info
#' head(mSeed$Gene)
#' }

xSNP2GeneScores <- function(data, include.LD=NA, LD.customised=NULL, LD.r2=0.8, significance.threshold=5e-5, score.cap=10, distance.max=50000, decay.kernel=c("slow","linear","rapid","constant"), decay.exponent=2, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), include.TAD=c("none","GM12878","IMR90","MSC","TRO","H1","MES","NPC"), scoring.scheme=c("max","sum","sequential"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    decay.kernel <- match.arg(decay.kernel)
    include.TAD <- match.arg(include.TAD)
    scoring.scheme <- match.arg(scoring.scheme)
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xSNPscores' is being called to score SNPs (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
	df_SNP <- xSNPscores(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, verbose=verbose, RData.location=RData.location, guid=guid)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xSNPscores' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xSNP2nGenes' is being called to define nearby genes (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
	df_nGenes <- xSNP2nGenes(data=df_SNP$SNP, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.TAD=include.TAD, verbose=verbose, RData.location=RData.location, guid=guid)
	
	if(include.TAD!='none'){
		TAD <- NULL
		df_nGenes <- base::subset(df_nGenes, TAD!='Excluded')
	}
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xSNP2nGenes' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    
    ## df_SNP df_nGenes
    allGenes <- sort(df_nGenes$Gene)
    allSNPs <- sort(df_SNP$SNP)
    
    ## sparse matrix of nGenes X SNPs
    G2S_n <- xSparseMatrix(df_nGenes[,c("Gene","SNP","Weight")], rows=allGenes, columns=allSNPs, verbose=verbose)
    G2S <- G2S_n
    
    ## consider SNP scores
    ind <- match(colnames(G2S), df_SNP$SNP)
    ########
    df_SNP <- df_SNP[ind,]
    ########
    SNP_score <- df_SNP$Score
    names(SNP_score) <- colnames(G2S)
    ## convert into matrix
    mat_SNP_score <- matrix(rep(SNP_score,each=nrow(G2S)), nrow=nrow(G2S))
    
    ## calculate genetic influence score for a gene-SNP pair
    G2S_score <- G2S * mat_SNP_score
    
    ## Gene2SNP
    Gene2SNP <- xSM2DF(data=G2S_score, verbose=FALSE)
    colnames(Gene2SNP) <- c('Gene','SNP','Score')
	ls_gene <- split(x=Gene2SNP$Score, f=Gene2SNP$Gene)
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
	
	df_SNP <- df_SNP[order(df_SNP$Flag,df_SNP$Score,-df_SNP$Pval,df_SNP$SNP,decreasing=TRUE),]
    
    mSeed <- list(SNP = df_SNP,
                  Gene = df_Gene,
                  Call = match.call()
                 )
    class(mSeed) <- "mSeed"   
    
    invisible(mSeed)
}
