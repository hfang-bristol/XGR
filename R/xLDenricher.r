#' Function to conduct LD-based enrichment analysis using genomic annotations via sampling
#'
#' \code{xLDenricher} is supposed to conduct LD-based enrichment analysis for the input genomic region data (genome build h19), using genomic annotations (eg active chromatin, transcription factor binding sites/motifs, conserved sites). Enrichment analysis is achieved by comparing the observed overlaps against the expected overlaps which are estimated from the null distribution. The null LD block is generated via sampling from the background (for example, all GWAS SNPs or all common SNPs), respecting the maf of the best SNP and/or the distance of the best SNP to the nearest gene, restricting the same chromosome or not. 
#'
#' @param bLD a bLD object, containing a set of blocks based on which to generate a null distribution
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 150) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 150) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised GR object directly
#' @param num.samples the number of samples randomly generated
#' @param respect how to respect the properties of to-be-sampled LD blocks. It can be one of 'maf' (respecting the maf of the best SNP), 'distance' (respecting the distance of the best SNP to the nearest gene), and 'both' (respecting the maf and distance)
#' @param restrict.chr logical to restrict to the same chromosome. By default, it sets to false
#' @param preserve how to preserve the resulting null LD block. It can be one of 'boundary' (preserving the boundary of the LD block), and 'exact' (exactly preserving the relative SNP locations within the LD block). Notably, no huge difference for the boundary preserving when enrichment analysis invovles region-based genomic annotations, but it may make difference when genomic annatations are largely SNP-based (such as eQTLs)
#' @param seed an integer specifying the seed
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param GR.annotation the genomic regions of annotation data. By default, it is 'NA' to disable this option. Pre-built genomic annotation data are detailed in \code{\link{xDefineGenomicAnno}}. Alternatively, the user can also directly provide a customised GR object (or a list of GR objects)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a data frame with 13 columns:
#' \itemize{
#'  \item{\code{name}: the annotation name}
#'  \item{\code{nAnno}: the number of regions from annotation data}
#'  \item{\code{nOverlap}: the observed number of LD blocks overlapped with annotation data}
#'  \item{\code{fc}: fold change}
#'  \item{\code{zscore}: z-score}
#'  \item{\code{pvalue}: p-value}
#'  \item{\code{adjp}: adjusted p-value. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{or}: a vector containing odds ratio}
#'  \item{\code{CIl}: a vector containing lower bound confidence interval for the odds ratio}
#'  \item{\code{CIu}: a vector containing upper bound confidence interval for the odds ratio}
#'  \item{\code{nData}: the number of input LD blocks}
#'  \item{\code{nExpect}: the expected number of LD blocks overlapped with annotation data}
#'  \item{\code{std}: the standard deviation of expected number of LD blocks overlapped with annotation data}
#' }
#' @note Pre-built genomic annotation data are detailed in \code{\link{xDefineGenomicAnno}}.
#' @export
#' @seealso \code{\link{xDefineGenomicAnno}}
#' @include xLDenricher.r
#' @examples
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#'
#' \dontrun{
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' data(ImmunoBase)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' data <- GenomicRanges::mcols(gr)[,c('Variant','Pvalue')]
#'
#' # b) get LD block (EUR population)
#' bLD <- xLDblock(data, include.LD="EUR", LD.r2=0.8, RData.location=RData.location)
#' 
#' ## c) perform enrichment analysis using FANTOM expressed enhancers
#' eTerm <- xLDenricher(bLD, GR.annotation="ReMap_Encode_mergedTFBS", RData.location=RData.location)
#'
#' ## d) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' ## e) barplot of enriched terms
#' bp <- xEnrichBarplot(eTerm, top_num='auto', displayBy="fdr")
#' bp
#'
#' ## f) forest plot of enrichment results
#' gp <- xEnrichForest(eTerm, FDR.cutoff=0.01)
#'
#' ## g) save enrichment results to the file called 'LD_enrichments.txt'
#' output <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' utils::write.table(output, file="LD_enrichments.txt", sep="\t", row.names=FALSE)
#' 
#' ## h) compare boundary and exact
#' GR.SNP <- xRDataLoader("dbSNP_GWAS", RData.location=RData.location)
#' GR.annotation <- xRDataLoader("FANTOM5_CAT_Cell", RData.location=RData.location)
#' eTerm_boundary <- xLDenricher(bLD, GR.SNP=GR.SNP, GR.annotation=GR.annotation, num.samples=20000, preserve="boundary", RData.location=RData.location)
#' eTerm_exact <- xLDenricher(bLD, GR.SNP=GR.SNP, GR.annotation=GR.annotation, num.samples=20000, preserve="exact", RData.location=RData.location)
#' ls_eTerm <- list(boundary=eTerm_boundary, exact=eTerm_exact)
#' ### barplot
#' bp <- xEnrichCompare(ls_eTerm, displayBy="zscore")
#' ### forest plot
#' eTerm_boundary$group <- 'boundary'
#' eTerm_exact$group <- 'exact'
#' df <- rbind(eTerm_boundary, eTerm_exact)
#' gp <- xEnrichForest(df, FDR.cutoff=0.01)
#' }

xLDenricher <- function(bLD, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), num.samples=2000, respect=c("maf","distance","both"), restrict.chr=F, preserve=c("exact","boundary"), seed=825, p.adjust.method=c("BH","BY","bonferroni","holm","hochberg","hommel"), GR.annotation=NA, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    respect <- match.arg(respect)
    preserve <- match.arg(preserve)
    p.adjust.method <- match.arg(p.adjust.method)
    
	#######################################################
	
	if(verbose){
		message(sprintf("First, load GR annotations (%s) ...", as.character(Sys.time())), appendLF=T)
	}
	
	aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location)
	if(is.null(aGRL)){
		return(NULL)
	}
    
	#######################################################
	if(verbose){
		message(sprintf("Second, generate null distribution via doing %d sampling (%s) ...", num.samples, as.character(Sys.time())), appendLF=T)
	}
	
    if(verbose){
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xLDsampling' is being called (%s):", as.character(Sys.time())), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
	grl <- xLDsampling(bLD=bLD, GR.SNP=GR.SNP, num.samples=num.samples, respect=respect, restrict.chr=restrict.chr, preserve=preserve, seed=seed, verbose=verbose, RData.location=RData.location)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xLDsampling' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
	
	#######################################################
	
	if(verbose){
		message(sprintf("Third, perform enrichment analysis (%s) ...", as.character(Sys.time())), appendLF=T)
	}
	
	queryHits <- subjectHits <- B <- n <- best <- NULL
	
	## observed
	### gr_observed
	if(preserve=="boundary"){
		gr_best <- bLD$best
		df_best <- GenomicRanges::as.data.frame(gr_best)
		df_best$start <- df_best$start + df_best$upstream
		df_best$end <- df_best$end + df_best$downstream
		gr_observed <- xGR(df_best, format='data.frame', add.name=F)
		gr_observed$best <- rownames(df_best)
		
	}else if(preserve=="exact"){
		grl_block <- bLD$block
		gr_block <- BiocGenerics::unlist(grl_block,use.names=F)
		gr_observed <- gr_block[,'best']
		
	}
	
	### obs
	system.time({
	q2r <- as.data.frame(GenomicRanges::findOverlaps(gr_observed, aGRL))
	q2r$best <- gr_observed$best[q2r$queryHits]
	df_n_observed <- as.data.frame(q2r %>% dplyr::group_by(subjectHits) %>% dplyr::summarise(n=length(unique(best))))
	obs <- rep(0, length(aGRL))
	obs[df_n_observed$subjectHits] <- df_n_observed$n
	names(obs) <- names(aGRL)
	})
	
	if(verbose){
		message(sprintf("\t%d observed (%s)", length(obs), as.character(Sys.time())), appendLF=T)
	}
	
	## expected
	### gr_expected
	gr_expected <- BiocGenerics::unlist(grl)
	### a2B
	system.time({
	q2r <- as.data.frame(GenomicRanges::findOverlaps(gr_expected, aGRL))
	q2r$best <- gr_expected$best[q2r$queryHits]
	q2r$B <- gr_expected$B[q2r$queryHits]
	})
	system.time({
	df_n_expected <- as.data.frame(q2r %>% dplyr::group_by(subjectHits,B) %>% dplyr::summarise(n=length(unique(best))))
	a2B <- as.matrix(xSparseMatrix(df_n_expected, rows=1:length(aGRL), columns=1:length(grl), verbose=T))
	rownames(a2B) <- names(aGRL)
	})

	if(verbose){
		message(sprintf("\t%d x %d of the expected matrix (%s)", nrow(a2B), ncol(a2B), as.character(now <- Sys.time())), appendLF=T)
	}

	### exp
	exp_mean <- apply(a2B, 1, mean)
	exp_std <- apply(a2B, 1, stats::sd)
	
	## ratio
	ratio <- obs/exp_mean
	### just in case: obs=0 & exp=0
	ratio[obs==0 & exp_mean==0] <- 1
	###
	
    ## zscore
    zscore <- (obs-exp_mean)/exp_std

    ## pvalue
    obs_matrix <- matrix(rep(obs,ncol(a2B)), ncol=ncol(a2B))
    pvalue <- apply((obs_matrix - a2B)<=0, 1, sum) / ncol(a2B)
	####################
	
    zscore[is.na(zscore)] <- 0
    zscore[is.infinite(zscore)] <- max(zscore[!is.infinite(zscore)])
    pvalue[is.na(ratio)] <- 1
    ratio[is.na(ratio)] <- 1
 	
 	####################
 	or <- CIl <- CIu <- NA
 	if(1){
		ls_df <- apply(a2B, 1, function(x){
			y <- stats::t.test(x)
			data.frame(mean=y$estimate, conf_lower=y$conf.int[1], conf_upper=y$conf.int[2], stringsAsFactors=F)
		})
		res_df <- do.call(rbind, ls_df)
		or <- obs / res_df[,1]
		CIu <- obs / res_df[,2]
		CIl <- obs / res_df[,3]
		
		### just in case: obs=0 & exp=0
		ind <- which(obs==0 & res_df[,1]==0)
		or[ind] <- 1
		CIu[ind] <- 1
		CIu[ind] <- 1
		###
 	}
 	####################
 
	enrichment_df <- data.frame(names(aGRL), sapply(aGRL,length), length(unique(gr_observed$best)), obs, exp_mean, exp_std, ratio, zscore, pvalue, or, CIl, CIu, row.names=NULL, stringsAsFactors=F)
	colnames(enrichment_df) <- c("name", "nAnno", "nData", "nOverlap", "nExpect", "std", "fc", "zscore", "pvalue", "or", "CIl", "CIu")

	## Adjust P-values for multiple comparisons
	pvals <- enrichment_df$pvalue
	adjpvals <- stats::p.adjust(pvals, method=p.adjust.method)
	enrichment_df$adjp <- adjpvals
	
	####################################################################################
	
	enrichment_df$zscore <- signif(enrichment_df$zscore, digits=3)
	
	pvals <- enrichment_df$pvalue
	adjpvals <- enrichment_df$adjp
	pvals <- signif(pvals, digits=2)
	adjpvals <- signif(adjpvals, digits=2)
	
	# scientific notations
	pvals  <- base::sapply(pvals, function(x){
		if(x < 0.1 & x!=0){
			as.numeric(format(x,scientific=T))
		}else{
			x
		}
	})
	
	adjpvals <- base::sapply(adjpvals, function(x){
		if(x < 0.1 & x!=0){
			as.numeric(format(x,scientific=T))
		}else{
			x
		}
	})
	
	enrichment_df$pvalue <- pvals
	enrichment_df$adjp <- adjpvals
	
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
	res_df <- enrichment_df[, c("name", "nAnno", "nOverlap", "fc", "zscore", "pvalue", "adjp","or", "CIl", "CIu", "nData", "nExpect","std")]
	
	invisible(res_df)
}
