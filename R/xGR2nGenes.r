#' Function to define nearby genes given a list of genomic regions
#'
#' \code{xGR2nGenes} is supposed to define nearby genes given a list of genomic regions (GR) within certain distance window. The distance weight is calcualted as a decaying function of the gene-to-GR distance. 
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param scoring logical to indicate whether gene-level scoring will be further calculated. By default, it sets to false
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param scoring.rescale logical to indicate whether gene scores will be further rescaled into the [0,1] range. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' If scoring sets to false, a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: nearby genes}
#'  \item{\code{GR}: genomic regions}
#'  \item{\code{Dist}: the genomic distance between the gene and the GR}
#'  \item{\code{Weight}: the distance weight based on the genomic distance}
#' }
#' If scoring sets to true, a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: nearby genes}
#'  \item{\code{Score}: gene score taking into account the distance weight based on the genomic distance}
#' }
#' @note For details on the decay kernels, please refer to \code{\link{xVisKernels}}
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xGR}}
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
#' # b) define nearby genes taking into acount distance weight
#' # without gene scoring
#' df_nGenes <- xGR2nGenes(data=data, format="chr:start-end", distance.max=10000, decay.kernel="slow", decay.exponent=2, RData.location=RData.location)
#' # with their scores
#' df_nGenes <- xGR2nGenes(data=data, format="chr:start-end", distance.max=10000, decay.kernel="slow", decay.exponent=2, scoring=T, scoring.scheme="max", RData.location=RData.location)
#'
#' # c) define nearby genes without taking into acount distance weight
#' # without gene scoring
#' df_nGenes <- xGR2nGenes(data=data, format="chr:start-end", distance.max=10000, decay.kernel="constant", RData.location=RData.location)
#' # with their scores
#' df_nGenes <- xGR2nGenes(data=data, format="chr:start-end", distance.max=10000, decay.kernel="constant", scoring=T, scoring.scheme="max", RData.location=RData.location)
#' }

xGR2nGenes <- function(data, format=c("chr:start-end","data.frame","bed","GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), distance.max=50000, decay.kernel=c("rapid","slow","linear","constant"), decay.exponent=2, GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), scoring=F, scoring.scheme=c("max","sum","sequential"), scoring.rescale=F, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    decay.kernel <- match.arg(decay.kernel)
    scoring.scheme <- match.arg(scoring.scheme)
	
	dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
  	#######################################################
  	
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=T)
	}
	if(class(GR.Gene) == "GRanges"){
		gr_Gene <- xGR(GR.Gene, format="GRanges", build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
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
	#minoverlap <- 1L # 1b overlaps
	minoverlap <- 0L
	subject <- gr_Gene
	query <- dGR
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
			df_nGenes <- data.frame(Gene=names(x), GR=names(y), Dist=dists, stringsAsFactors=F)
		}else{
			### very slow
			list_gene <- split(x=q2r[,1], f=q2r[,2])
			ind_gene <- as.numeric(names(list_gene))
			res_list <- lapply(1:length(ind_gene), function(i){
				x <- subject[ind_gene[i],]
				y <- query[list_gene[[i]],]
				dists <- GenomicRanges::distance(x, y, select="all", ignore.strand=T)
				res <- data.frame(Gene=rep(names(x),length(dists)), GR=names(y), Dist=dists, stringsAsFactors=F)
			})
			df_nGenes <- do.call(rbind, res_list)
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
		
		df_nGenes <- df_nGenes[order(df_nGenes$Gene,df_nGenes$Dist,decreasing=FALSE),]
		
		############################################
		## whether gene scoring
		if(scoring){
			if(0){
			## sparse matrix of nGenes X GR
			G2S_score <- xSparseMatrix(df_nGenes[,c("Gene","GR","Weight")], verbose=verbose)
			}
			
			ls_gene <- split(x=df_nGenes$Weight, f=df_nGenes$Gene)
			
			## calculate genetic influence score under a set of SNPs for each seed gene
			if(scoring.scheme=='max'){
				seeds.genes <- sapply(ls_gene, max)
				
				if(0){
					##########################
					## to be improved (if too many columns)
					#mat <- as.matrix(G2S_score)
					#mat <- G2S_score
					#seeds.genes <- do.call(base::pmax, lapply(1:ncol(mat), function(j) mat[,j]))
					##########################
					ind <- !duplicated(df_nGenes[,1])
					seeds.genes <- df_nGenes[ind, 4]
					names(seeds.genes) <- df_nGenes[ind, 1]					
				
					system.time({
					seeds.genes <- apply(G2S_score, 1, function(x) {
						base::max(x)
					})
					})
				}
				
			}else if(scoring.scheme=='sum'){
				seeds.genes <- sapply(ls_gene, sum)
				
				if(0){
					##########################
					## to be improved (if too many columns)
					seeds.genes <- base::rowSums(as.matrix(G2S_score))
					##########################
					seeds.genes <- apply(G2S_score, 1, function(x) {
						base::sum(x)
					})
				}
				
			}else if(scoring.scheme=='sequential'){
				seeds.genes <- sapply(ls_gene, function(x){
					#base::sum(base::sort(x, decreasing=T) / (1:length(x)))
					base::sum(x / base::rank(-x,ties.method="min"))
				})
				
				if(0){
				seeds.genes <- apply(G2S_score, 1, function(x) {
					base::sum(base::sort(x, decreasing=T) / (1:length(x)))
				})
				}
			}

			if(verbose){
				now <- Sys.time()
				message(sprintf("In summary, %d Genes are defined as seeds and scored using '%s' scoring scheme (%s)", length(seeds.genes), scoring.scheme, as.character(now)), appendLF=T)
			}

			if(scoring.rescale){
				if(verbose){
					now <- Sys.time()
					message(sprintf("Also rescale score into the [0,1] range (%s)", as.character(now)), appendLF=T)
				}
				# rescale to [0 1]
				seeds.genes <- (seeds.genes - min(seeds.genes))/(max(seeds.genes) - min(seeds.genes))
			}

			## for output
			df_Gene <- data.frame(Gene=names(seeds.genes), Score=seeds.genes, stringsAsFactors=F)
			rownames(df_Gene) <- NULL
			df_Gene <- df_Gene[order(df_Gene$Score,decreasing=TRUE),]
			
			invisible(df_Gene)
		
		}else{
			invisible(df_nGenes)
		}
		
	}else{
		df_nGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No nearby genes are defined"), appendLF=T)
		}
		
		invisible(df_nGenes)
	}

}
