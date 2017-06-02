#' Function to define nearby genes given a list of SNPs
#'
#' \code{xSNP2nGenes} is supposed to define nearby genes given a list of SNPs within certain distance window. The distance weight is calcualted as a decaying function of the gene-to-SNP distance. 
#'
#' @param data a input vector containing SNPs. SNPs should be provided as dbSNP ID (ie starting with rs). Alternatively, they can be in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is genomic positional number; for example, 'chr16:28525386'
#' @param distance.max the maximum distance between genes and SNPs. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby SNPs per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param include.TAD TAD boundary regions are also included. By default, it is 'NA' to disable this option. Otherwise, inclusion of a TAD dataset to pre-filter SNP-nGene pairs (i.e. only those within a TAD region will be kept). TAD datasets can be one of "GM12878"  (lymphoblast), "IMR90" (fibroblast), "MSC" (mesenchymal stem cell) ,"TRO" (trophoblasts-like cell), "H1" (embryonic stem cell), "MES" (mesendoderm) and "NPC" (neural progenitor cell). Explanations can be found at \url{http://dx.doi.org/10.1016/j.celrep.2016.10.061}
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: nearby genes}
#'  \item{\code{SNP}: SNPs}
#'  \item{\code{Dist}: the genomic distance between the gene and the SNP}
#'  \item{\code{Weight}: the distance weight based on the genomic distance}
#'  \item{\code{Gap}: the genomic gap between the gene and the SNP (in the form of 'chr:start-end')}
#'  \item{\code{TAD}: if applied, it can be 'Excluded' or the TAD boundary region (in the form of 'chr:start-end') that the genomic interval falls into}
#' }
#' @note For details on the decay kernels, please refer to \code{\link{xVisKernels}}
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xVisKernels}}
#' @include xSNP2nGenes.r
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
#' data <- names(gr)
#'
#' # b) define nearby genes
#' df_nGenes <- xSNP2nGenes(data=data, distance.max=200000, decay.kernel="slow", decay.exponent=2, RData.location=RData.location)
#'
#' # c) define nearby genes (considering TAD boundary regions in GM12878)
#' df_nGenes <- xSNP2nGenes(data=data, distance.max=200000, decay.kernel="slow", decay.exponent=2, include.TAD='GM12878', RData.location=RData.location)
#' }

xSNP2nGenes <- function(data, distance.max=200000, decay.kernel=c("rapid","slow","linear","constant"), decay.exponent=2, GR.SNP=c("dbSNP_GWAS","dbSNP_Common"), GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), include.TAD=c(NA,"GM12878","IMR90","MSC","TRO","H1","MES","NPC"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    decay.kernel <- match.arg(decay.kernel)
    include.TAD <- match.arg(include.TAD)
	
	## replace '_' with ':'
	data <- gsub("_", ":", data, perl=T)
	## replace 'imm:' with 'chr'
	data <- gsub("imm:", "chr", data, perl=T)
	
	data <- unique(data)
	
    ######################################################
    # Link to targets based on genomic distance
    ######################################################
    
    gr_SNP <- xSNPlocations(data=data, GR.SNP=GR.SNP, verbose=verbose, RData.location=RData.location)

  	#######################################################
  	
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=T)
	}
	if(class(GR.Gene) == "GRanges"){
			gr_Gene <- GR.Gene
	}else{
		gr_Gene <- xRDataLoader(RData.customised=GR.Gene[1], verbose=verbose, RData.location=RData.location)
		if(is.null(gr_Gene)){
			GR.Gene <- "UCSC_knownGene"
			if(verbose){
				message(sprintf("Instead, %s will be used", GR.Gene), appendLF=T)
			}
			gr_Gene <- xRDataLoader(RData.customised=GR.Gene, verbose=verbose, RData.location=RData.location)
		}
    }
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Define nearby genes (%s) ...", as.character(now)), appendLF=T)
	}
    
	# genes: get all UCSC genes within defined distance window away from variants
	maxgap <- distance.max
	minoverlap <- 1L # 1b overlaps
	subject <- gr_Gene
	query <- gr_SNP
	q2r <- as.matrix(as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T))))
	
	if(length(q2r) > 0){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Calculate distance (%s) ...", as.character(now)), appendLF=T)
		}
	
		if(1){
			### very quick
			x <- subject[q2r[,2],]
			y <- query[q2r[,1],]
			dists <- GenomicRanges::distance(x, y, select="all", ignore.strand=T)
			
			###
			df_y <- GenomicRanges::as.data.frame(y, row.names=NULL)
			df_x <- GenomicRanges::as.data.frame(x, row.names=NULL)
			
			if(0){
				## Gap defined as the regions from an SNP to the closest end of a Gene
				df_interval <- data.frame(seqnames=df_y$seqnames, start=df_y$start, end=df_y$end, stringsAsFactors=FALSE)
				ind <- df_y$start < df_x$start
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_y$start[ind], end=df_x$start[ind], stringsAsFactors=FALSE)
				ind <- df_y$start > df_x$end
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_x$end[ind], end=df_y$start[ind], stringsAsFactors=FALSE)
				
			}else{
				## Gap defined as the regions from an SNP to the TSS of a Gene
				
				df_interval <- data.frame(seqnames=df_y$seqnames, start=df_y$start, end=df_y$end, stringsAsFactors=FALSE)
				ind <- df_x$strand=='+'
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_y$start[ind], end=df_x$start[ind], stringsAsFactors=FALSE)
				ind <- df_x$strand=='-'
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_y$start[ind], end=df_x$end[ind], stringsAsFactors=FALSE)
				ind <- df_x$strand=='*'
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_y$start[ind], end=(df_x$start[ind]+df_x$end[ind])/2, stringsAsFactors=FALSE)
				
				## swap the location of start and end
				ind <- df_interval$start > df_interval$end
				df_interval[ind,] <- data.frame(seqnames=df_interval$seqnames[ind], start=df_interval$end[ind], end=df_interval$start[ind], stringsAsFactors=FALSE)
			}
			
			vec_interval <- paste0(df_interval$seqnames,':',df_interval$start,'-',df_interval$end)
			###
			
			df_nGenes <- data.frame(Gene=names(x), SNP=names(y), Dist=dists, Gap=vec_interval, stringsAsFactors=F)
		}
	
		## weights according to distance away from SNPs
		if(distance.max==0){
			x <- df_nGenes$Dist
		}else{
			x <- df_nGenes$Dist / distance.max
		}
		if(decay.kernel == 'slow'){
			y <- 1-(x)^decay.exponent
		}else if(decay.kernel == 'rapid'){
			y <- (1-x)^decay.exponent
		}else if(decay.kernel == 'linear'){
			y <- 1-x
		}else{
			y <- 1
		}
		df_nGenes$Weight <- y
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("%d Genes are defined as nearby genes within %d(bp) genomic distance window using '%s' decay kernel (%s)", length(unique(df_nGenes$Gene)), distance.max, decay.kernel, as.character(now)), appendLF=T)
		}
		
		df_nGenes <- df_nGenes[,c('Gene','SNP','Dist','Weight','Gap')]
		df_nGenes <- df_nGenes[order(df_nGenes$Gene,df_nGenes$Dist,decreasing=FALSE),]
		
	}else{
		df_nGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No nearby genes are defined"), appendLF=T)
		}
	}
	
	###########################
	## include TAD
	default.include.TAD <- c("GM12878","IMR90","MSC","TRO","H1","MES","MES")
	ind <- match(default.include.TAD, include.TAD)
	include.TAD <- default.include.TAD[!is.na(ind)]
	if(length(include.TAD) > 0){
		if(verbose){
			now <- Sys.time()
			message(sprintf("Inclusion of TAD boundary regions is based on '%s'", include.TAD), appendLF=T)
		}
		df_nGenes$TAD <- rep('Excluded', nrow(df_nGenes))
		
		TAD <- xRDataLoader(RData.customised=paste0('TAD.',include.TAD), RData.location=RData.location, verbose=verbose)
		iGR <- xGR(data=df_nGenes$Gap, format="chr:start-end", RData.location=RData.location)
		q2r <- as.matrix(as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=iGR, subject=TAD, type="within", select="all", ignore.strand=T))))
		q2r <- q2r[!duplicated(q2r[,1]), ]
		df_nGenes$TAD[q2r[,1]] <- GenomicRanges::mcols(TAD)[q2r[,2],]
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("\t%d out of %d SNP-nGene pairs are within the same TAD boundary regions", sum(df_nGenes$TAD!='Excluded'), length(iGR)), appendLF=T)
			message(sprintf("\t%d out of %d genes are defined as nearby genes after considering TAD boundary regions", length(unique(df_nGenes[df_nGenes$TAD!='Excluded','Gene'])), length(unique(df_nGenes$Gene))), appendLF=T)
		}
		
	}
	
	####################################
	# only keep those genes with GeneID
	####################################
	if(!is.null(df_nGenes)){
		ind <- xSymbol2GeneID(df_nGenes$Gene, details=FALSE, verbose=verbose, RData.location=RData.location)
		df_nGenes <- df_nGenes[!is.na(ind), ]
		if(nrow(df_nGenes)==0){
			df_nGenes <- NULL
		}
	}
	####################################
	
    invisible(df_nGenes)
}
