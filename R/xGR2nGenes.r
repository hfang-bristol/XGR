#' Function to define nearby genes given a list of genomic regions
#'
#' \code{xGR2nGenes} is supposed to define nearby genes given a list of genomic regions (GR) within certain distance window. The distance weight is calcualted as a decaying function of the gene-to-GR distance. 
#'
#' @param data a input vector containing genomic regions (GR). GR should be provided in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay
#' @param decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: nearby genes}
#'  \item{\code{GR}: genomic regions}
#'  \item{\code{Dist}: the genomic distance between the gene and the GR}
#'  \item{\code{Weight}: the distance weight based on the genomic distance}
#' }
#' @note For details on the decay kernels, please refer to \code{\link{xVisKernels}}
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xVisKernels}}
#' @include xGR2nGenes.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#'
#' # a) provide the genomic regions
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' df <- as.data.frame(gr, row.names=NULL)
#' chr <- df$seqnames
#' start <- df$start
#' end <- df$end
#' data <- paste(chr,':',start,'-',end, sep='')
#'
#' # b) define nearby genes
#' df_nGenes <- xGR2nGenes(data=data, distance.max=10000, decay.kernel="slow", decay.exponent=2, RData.location=RData.location)
#' }

xGR2nGenes <- function(data, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), distance.max=50000, decay.kernel=c("rapid","slow","linear"), decay.exponent=2, GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    build.conversion <- match.arg(build.conversion)
    decay.kernel <- match.arg(decay.kernel)
	
    ######################################################
    # Link to targets based on genomic distance
    ######################################################
    data <- unique(data)
    
	## construct GR
	input <- do.call(rbind, strsplit(data, ":|-"))
	if(ncol(input)>=3){
		data <- input[,1:3]
	}else if(ncol(input)==2){
		data <- input[,c(1,2,2)]
	}else{
		stop("Your input 'data' does not meet the format 'chr:start-end'!\n")
	}
	
	## make sure positions are numeric
	ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
	data <- data[ind,]
	dGR <- GenomicRanges::GRanges(
		seqnames=S4Vectors::Rle(data[,1]),
		ranges = IRanges::IRanges(start=as.numeric(data[,2]), end=as.numeric(data[,3])),
		strand = S4Vectors::Rle(rep('*',nrow(data)))
	)
	names(dGR) <- paste(data[,1], ':', data[,2], '-', data[,3], sep='')
	
	# lift over
	if(!is.na(build.conversion)){
		if(verbose){
			message(sprintf("\tdata genomic regions: lifted over via genome build conversion `%s`", build.conversion), appendLF=T)
		}
		dGR <- xLiftOver(data.file=dGR, format.file="GRanges", build.conversion=build.conversion, merged=F, verbose=verbose, RData.location=RData.location)
	}
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
    
	# genes: get all UCSC genes within defined distance window away from variants
	maxgap <- distance.max
	minoverlap <- 1L # 1b overlaps
	subject <- gr_Gene
	query <- dGR
	q2r <- as.matrix(suppressWarnings(GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
	
	if(length(q2r) > 0){
	
		list_gene <- split(x=q2r[,1], f=q2r[,2])
		ind_gene <- as.numeric(names(list_gene))
		res_list <- lapply(1:length(ind_gene), function(i){
			x <- subject[ind_gene[i],]
			y <- query[list_gene[[i]],]
			dists <- GenomicRanges::distance(x, y, select="all", ignore.strand=T)
			res <- data.frame(Gene=rep(names(x),length(dists)), GR=names(y), Dist=dists, stringsAsFactors=F)
		})
	
		## weights according to distance away from SNPs
		df_nGenes <- do.call(rbind, res_list)
		if(distance.max==0){
			x <- df_nGenes$Dist
		}else{
			x <- df_nGenes$Dist / distance.max
		}
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
		
		df_nGenes <- df_nGenes[order(df_nGenes$Gene,df_nGenes$Dist,decreasing=FALSE),]
		
	}else{
		df_nGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No nearby genes are defined"), appendLF=T)
		}
	}
	
    invisible(df_nGenes)
}
