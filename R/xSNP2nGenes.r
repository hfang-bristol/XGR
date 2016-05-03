#' Function to define nearby genes given a list of SNPs
#'
#' \code{xSNP2nGenes} is supposed to define nearby genes given a list of SNPs within certain distance window. The distance weight is calcualted as a decaying function of the gene-to-SNP distance. 
#'
#' @param data a input vector containing SNPs
#' @param distance.max the maximum distance between genes and SNPs. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby SNPs per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay
#' @param decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_genes', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Even the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: nearby genes}
#'  \item{\code{SNP}: SNPs}
#'  \item{\code{Dist}: the genomic distance between the gene and the SNP}
#'  \item{\code{Weight}: the distance weight based on the gnomic distance}
#' }
#' @note For details on the decay kernels, please refer to \code{\link{xVisKernels}}
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xVisKernels}}
#' @include xSNP2nGenes.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' library(GenomicRanges)
#'
#' # a) provide the SNPs with the significance info
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' AS <- read.delim(file.path(path.package("Pi"),"AS.txt"), stringsAsFactors=FALSE)
#'
#' # b) define nearby genes
#' df_nGenes <- xSNP2nGenes(data=AS$SNP, distance.max=200000, decay.kernel="slow", decay.exponent=2, RData.location=RData.location)
#' }

xSNP2nGenes <- function(data, distance.max=200000, decay.kernel=c("slow","linear","rapid"), decay.exponent=2, GR.SNP="dbSNP_GWAS", GR.Gene="UCSC_genes", verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/Portal")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    decay.kernel <- match.arg(decay.kernel)
	
	data <- unique(data)
	
    ######################################################
    # Link to targets based on genomic distance
    ######################################################
    
  	## load positional information
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for SNPs (%s) ...", as.character(now)), appendLF=T)
	}
    pos_SNP <- xRDataLoader(RData.customised=GR.SNP, verbose=verbose, RData.location=RData.location)
    if(is.null(pos_SNP)){
    	GR.SNP <- "dbSNP_GWAS"
		if(verbose){
			message(sprintf("Instead, %s will be used", GR.SNP), appendLF=T)
		}
    	pos_SNP <- xRDataLoader(RData.customised=GR.SNP, verbose=verbose, RData.location=RData.location)
    }
  	ind <- match(data, names(pos_SNP))
  	ind <- ind[!is.na(ind)]
  	if(length(ind)){
  		gr_SNP <- pos_SNP[ind,]
  		
		if(verbose){
			now <- Sys.time()
			message(sprintf("\tOut of %d input SNPs, %d SNPs have positional info", length(data), length(gr_SNP)), appendLF=T)
		}
  	}
  	
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=T)
	}
    gr_Gene <- xRDataLoader(RData.customised=GR.Gene, verbose=verbose, RData.location=RData.location)
    if(is.null(gr_Gene)){
    	GR.Gene <- "UCSC_genes"
		if(verbose){
			message(sprintf("Instead, %s will be used", GR.Gene), appendLF=T)
		}
    	gr_Gene <- xRDataLoader(RData.customised=GR.Gene, verbose=verbose, RData.location=RData.location)
    }
    
	# genes: get all UCSC genes within defined distance window away from variants
	maxgap <- distance.max
	minoverlap <- 1L # 1b overlaps
	subject <- gr_Gene
	query <- gr_SNP
	q2r <- as.matrix(suppressWarnings(GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
	
	if(length(q2r) > 0){
	
		list_gene <- split(x=q2r[,1], f=q2r[,2])
		ind_gene <- as.numeric(names(list_gene))
		res_list <- lapply(1:length(ind_gene), function(i){
			x <- subject[ind_gene[i],]
			y <- query[list_gene[[i]],]
			dists <- GenomicRanges::distance(x, y, select="all", ignore.strand=T)
			res <- data.frame(Gene=rep(names(x),length(dists)), SNP=names(y), Dist=dists, stringsAsFactors=F)
		})
	
		## weights according to distance away from SNPs
		df_nGenes <- do.call(rbind, res_list)
		x <- df_nGenes$Dist / distance.max
		if(decay.kernel == 'slow'){
			y <- 1-(x)^decay.exponent
		}else if(decay.kernel == 'rapid'){
			y <- (1-x)^decay.exponent
		}else{
			y <- 1-x
		}
		df_nGenes$Weight <- y
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("%d Genes are defined as nearby genes within %d(bp) genomic distance window", length(unique(df_nGenes$Gene)), distance.max), appendLF=T)
		}
	
	}else{
		df_nGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No nearby genes are defined"), appendLF=T)
		}
	}
	
    invisible(df_nGenes)
}
