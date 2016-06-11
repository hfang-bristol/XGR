#' Function to define nearby genes given a list of SNPs
#'
#' \code{xSNP2nGenes} is supposed to define nearby genes given a list of SNPs within certain distance window. The distance weight is calcualted as a decaying function of the gene-to-SNP distance. 
#'
#' @param data a input vector containing SNPs. SNPs should be provided as dbSNP ID (ie starting with rs). Alternatively, they can be in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for example, 'chr16:28525386'
#' @param distance.max the maximum distance between genes and SNPs. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby SNPs per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay
#' @param decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
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
#' library(XGR)
#' RData.location="~/Sites/SVN/github/RDataCentre/Portal"
#'
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' data <- names(gr)
#'
#' # b) define nearby genes
#' df_nGenes <- xSNP2nGenes(data=data, distance.max=200000, decay.kernel="slow", decay.exponent=2, RData.location=RData.location)
#' }

xSNP2nGenes <- function(data, distance.max=200000, decay.kernel=c("rapid","slow","linear"), decay.exponent=2, GR.SNP=c("dbSNP_GWAS","dbSNP_Common"), GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/Portal")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    decay.kernel <- match.arg(decay.kernel)
	
	data <- unique(data)
	
	## replace '_' with ':'
	data <- gsub("_", ":", data, perl=T)
	## replace 'imm:' with 'chr'
	data <- gsub("imm:", "chr", data, perl=T)
	
    ######################################################
    # Link to targets based on genomic distance
    ######################################################
    
  	## load positional information
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for SNPs (%s) ...", as.character(now)), appendLF=T)
	}
    pos_SNP <- xRDataLoader(RData.customised=GR.SNP[1], verbose=verbose, RData.location=RData.location)
    if(is.null(pos_SNP)){
    	GR.SNP <- "dbSNP_GWAS"
		if(verbose){
			message(sprintf("Instead, %s will be used", GR.SNP), appendLF=T)
		}
    	pos_SNP <- xRDataLoader(RData.customised=GR.SNP, verbose=verbose, RData.location=RData.location)
    }
  	ind <- match(data, names(pos_SNP))
  	data_rest <- data[is.na(ind)]
  	ind <- ind[!is.na(ind)]
  	if(length(ind)){
  		gr_SNP <- pos_SNP[ind,]
  		GenomicRanges::mcols(gr_SNP) <- NULL
  	}else{
  		gr_SNP <- NULL
  	}
  	
  	#######################################################
  	### deal with data_rest having the format as: chr\w+:\d+
  	ind <- grep("^chr\\w+:\\d+", data_rest, perl=T)
  	if(length(ind)>0){
		data_rest <- data_rest[ind]
		res_ls <- strsplit(data_rest, ":")
		res_df <- do.call(rbind, res_ls)
		res_gr <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(res_df[,1]),
			ranges = IRanges::IRanges(start=as.numeric(res_df[,2]), end=as.numeric(res_df[,2]), names=data_rest),
			strand = S4Vectors::Rle(rep('*',nrow(res_df)))
		)
	}else{
		res_gr <- NULL
	}
	### combine gr_SNP and res_gr
	if(!is.null(gr_SNP)){
		if(!is.null(res_gr)){
			gr_SNP <- c(gr_SNP, res_gr)
		}
	}else{
		if(!is.null(res_gr)){
			gr_SNP <- res_gr
		} 
	}
  	if(verbose){
		now <- Sys.time()
		message(sprintf("\tOut of %d input SNPs, %d SNPs have positional info", length(data), length(gr_SNP)), appendLF=T)
  	}
  	#######################################################
  	
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=T)
	}
    gr_Gene <- xRDataLoader(RData.customised=GR.Gene[1], verbose=verbose, RData.location=RData.location)
    if(is.null(gr_Gene)){
    	GR.Gene <- "UCSC_knownGene"
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
			message(sprintf("%d Genes are defined as nearby genes within %d(bp) genomic distance window using '%s' decay kernel", length(unique(df_nGenes$Gene)), distance.max, decay.kernel), appendLF=T)
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
