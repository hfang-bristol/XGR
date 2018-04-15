#' Function to calculate regulatory potential scores for genomic regions using genomic annotations
#'
#' \code{xRPS} is supposed to calculate regulatory potential scores for genomic regions using genomic annotations. 
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param GR.annotation the genomic regions of annotation data. Pre-built genomic annotation data are detailed in the section 'Note'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a GenomicRanges object appended a metacolumn 'RPS'
#' @note The genomic annotation data are described below according to the data sources and data types.\cr
#' 1. FANTOM5 expressed enhancer atlas
#' \itemize{
#'  \item{\code{FANTOM5_Enhancer_Cell}: a list (71 cell types) of GenomicRanges objects; each is an GR object containing enhancers specifically expressed in a cell type.}
#'  \item{\code{FANTOM5_Enhancer_Tissue}: a list (41 tissues) of GenomicRanges objects; each is an GR object containing enhancers specifically expressed in a tissue.}
#' }
#' 2. FANTOM5 sample-ontology-enriched CAT genes
#' \itemize{
#'  \item{\code{FANTOM5_CAT_Cell}: a list (173 cell types) of GenomicRanges objects; each is an GR object containing CAT genes specifically expressed in a cell type.}
#'  \item{\code{FANTOM5_CAT_Tissue}: a list (174 tissues) of GenomicRanges objects; each is an GR object containing CAT genes specifically expressed in a tissue.}
#' }
#' 3. GWAS Catalog trait-associated SNPs
#' \itemize{
#'  \item{\code{GWAScatalog_alltraits}: a list (390 traits grouped by EFO) of GenomicRanges objects; each is an GR object containing trait-associated SNPs.}
#' }
#' 4. ENCODE DNaseI Hypersensitivity site data
#' \itemize{
#'  \item{\code{ENCODE_DNaseI_ClusteredV3}: a GR object containing clustered peaks, along with a meta-column 'num_cells' telling how many cell types associated with a clustered peak.}
#' }
#' 5. ENCODE Transcription Factor ChIP-seq data
#' \itemize{
#'  \item{\code{ENCODE_TFBS_ClusteredV3}: a list (161 transcription factors) of GenomicRanges objects; each is an GR object containing clustered peaks per transcription factor, along with a meta-column 'cells' telling cell types associated with a clustered peak.}
#' }
#' 6. Roadmap Epigenomics Core 15-state Genome Segmentation data for 127 cell types
#' \itemize{
#' \item{\code{EpigenomeAtlas_15Segments}: a list (127 cell types) of a list (15 categories of segments) of GenomicRanges objects; each is an GR object containing segments per category in the reference epigenome.}
#' }
#' 7. Genomic scores
#' \itemize{
#'  \item{\code{RecombinationRate}: a GR object containing a meta-column for recombination rates.}
#'  \item{\code{phastCons100way}: a GR object containing a meta-column for phastCons100way.}
#'  \item{\code{phyloP100way}: a GR object containing a meta-column for phyloP100way.}
#' }
#' @export
#' @seealso \code{\link{xEnrichViewer}}
#' @include xRPS.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' 
#' # transcribed lncRNAs
#' FANTOM5_CAT <- xRDataLoader('FANTOM5_CAT', RData.location=RData.location)
#' GR_lncRNA <- FANTOM5_CAT[grepl('lncRNA',FANTOM5_CAT$Class)]
#' names(GR_lncRNA) <- NULL
#' data <- GR_lncRNA
#' # RPS calculation
#' dGR <- xRPS(data, format="GRanges", GR.annotation=c("FANTOM5_Enhancer_Cell","FANTOM5_Enhancer_Tissue","FANTOM5_CAT_Cell","FANTOM5_CAT_Tissue","GWAScatalog_alltraits","ENCODE_DNaseI_ClusteredV3","ENCODE_TFBS_ClusteredV3","EpigenomeAtlas_15Segments","RecombinationRate","phastCons100way","phyloP100way"), RData.location=RData.location)
#' }

xRPS <- function(data, format=c("data.frame", "bed", "chr:start-end", "GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), GR.annotation=c("FANTOM5_Enhancer_Cell","FANTOM5_Enhancer_Tissue","FANTOM5_CAT_Cell","FANTOM5_CAT_Tissue","GWAScatalog_alltraits","ENCODE_DNaseI_ClusteredV3","ENCODE_TFBS_ClusteredV3","EpigenomeAtlas_15Segments","RecombinationRate","phastCons100way","phyloP100way"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    
    dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
    
    #########################################################################
    
	default.GR_annotation <- c("FANTOM5_Enhancer_Cell","FANTOM5_Enhancer_Tissue","FANTOM5_CAT_Cell","FANTOM5_CAT_Tissue","GWAScatalog_alltraits","ENCODE_DNaseI_ClusteredV3","ENCODE_TFBS_ClusteredV3","EpigenomeAtlas_15Segments","RecombinationRate","phastCons100way","phyloP100way")
	names(default.GR_annotation) <- c("Fcell","Ftissue","CATcell","CATtissue","Trait","DHS","TF","Epi","RR","phastCons","phyloP")
	ind <- match(default.GR_annotation, GR.annotation)
	GR_annotation <- default.GR_annotation[!is.na(ind)]
	if(length(GR_annotation)==0){
		return(NULL)
	}
    
	if(verbose){
		message(sprintf("First, %d genomic regions are scored using %d annotations '%s' (%s) ...", length(dGR), length(GR_annotation), paste(GR_annotation,collapse=','), as.character(Sys.time())), appendLF=T)
	}
    
	##############
	# Fcell: number of cell types
	##############
	if("FANTOM5_Enhancer_Cell" %in% GR_annotation){
	
		if(verbose){
			message(sprintf("using the annotation '%s' (%s) ...", 'FANTOM5_Enhancer_Cell', as.character(Sys.time())), appendLF=T)
		}
	
		# FANTOM5_Enhancer_Cell: a list (71 cell types) of GenomicRanges objects; each is an GR object containing enhancers specifically expressed in a cell type
		ls_gr <- xRDataLoader('FANTOM5_Enhancer_Cell', verbose=F, RData.location=RData.location)
		## number of cell types
		ls_vec <- lapply(ls_gr, function(gr){
			vec <- rep(0, length(dGR))
			q2r <- as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=1L, type="any", select="all", ignore.strand=T)))
			ind <- unique(q2r[,1])
			vec[ind] <- 1
			return(vec)
		})
		df_res <- do.call(cbind, ls_vec)
		colnames(df_res) <- names(ls_gr)
		vec_res <- apply(df_res, 1, sum)
		vec_res[vec_res==0] <- NA
		dGR$Fcell <- vec_res
    }
    
	##############
	# Ftissue: number of tissues
	##############
	if("FANTOM5_Enhancer_Tissue" %in% GR_annotation){
	
		if(verbose){
			message(sprintf("using the annotation '%s' (%s) ...", 'FANTOM5_Enhancer_Tissue', as.character(Sys.time())), appendLF=T)
		}
	
		# FANTOM5_Enhancer_Tissue: a list (41 tissues) of GenomicRanges objects; each is an GR object containing enhancers specifically expressed in a tissue
		ls_gr <- xRDataLoader('FANTOM5_Enhancer_Tissue', verbose=F, RData.location=RData.location)
		## number of tissues
		ls_vec <- lapply(ls_gr, function(gr){
			vec <- rep(0, length(dGR))
			q2r <- as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=1L, type="any", select="all", ignore.strand=T)))
			ind <- unique(q2r[,1])
			vec[ind] <- 1
			return(vec)
		})
		df_res <- do.call(cbind, ls_vec)
		colnames(df_res) <- names(ls_gr)
		vec_res <- apply(df_res, 1, sum)
		vec_res[vec_res==0] <- NA
		dGR$Ftissue <- vec_res
    }
    
	##############
	# CATcell: number of cell types
	##############
	if("FANTOM5_CAT_Cell" %in% GR_annotation){
	
		if(verbose){
			message(sprintf("using the annotation '%s' (%s) ...", 'FANTOM5_CAT_Cell', as.character(Sys.time())), appendLF=T)
		}
	
		# FANTOM5_CAT_Cell: a list (173 cell types) of GenomicRanges objects; each is an GR object containing CAT genes specifically expressed in a cell type
		ls_gr <- xRDataLoader('FANTOM5_CAT_Cell', verbose=F, RData.location=RData.location)
		## number of cell types
		ls_vec <- lapply(ls_gr, function(gr){
			vec <- rep(0, length(dGR))
			q2r <- as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=1L, type="any", select="all", ignore.strand=T)))
			ind <- unique(q2r[,1])
			vec[ind] <- 1
			return(vec)
		})
		df_res <- do.call(cbind, ls_vec)
		colnames(df_res) <- names(ls_gr)
		vec_res <- apply(df_res, 1, sum)
		vec_res[vec_res==0] <- NA
		dGR$CATcell <- vec_res
    }
    
	##############
	# CATtissue: number of tissues
	##############
	if("FANTOM5_CAT_Tissue" %in% GR_annotation){
	
		if(verbose){
			message(sprintf("using the annotation '%s' (%s) ...", 'FANTOM5_CAT_Tissue', as.character(Sys.time())), appendLF=T)
		}
	
		# FANTOM5_CAT_Tissue: a list (174 tissues) of GenomicRanges objects; each is an GR object containing CAT genes specifically expressed in a tissue
		ls_gr <- xRDataLoader('FANTOM5_CAT_Tissue', verbose=F, RData.location=RData.location)
		## number of tissues
		ls_vec <- lapply(ls_gr, function(gr){
			vec <- rep(0, length(dGR))
			q2r <- as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=1L, type="any", select="all", ignore.strand=T)))
			ind <- unique(q2r[,1])
			vec[ind] <- 1
			return(vec)
		})
		df_res <- do.call(cbind, ls_vec)
		colnames(df_res) <- names(ls_gr)
		vec_res <- apply(df_res, 1, sum)
		vec_res[vec_res==0] <- NA
		dGR$CATtissue <- vec_res
    }
    
	##############
	# Trait: number of traits
	##############
	if("GWAScatalog_alltraits" %in% GR_annotation){
	
		if(verbose){
			message(sprintf("using the annotation '%s' (%s) ...", 'GWAScatalog_alltraits', as.character(Sys.time())), appendLF=T)
		}
	
		# GWAScatalog_alltraits: a list (390 traits grouped by EFO) of GenomicRanges objects; each is an GR object containing trait-associated SNPs
		ls_gr <- xRDataLoader('GWAScatalog_alltraits', verbose=F, RData.location=RData.location)
		## number of tissues
		ls_vec <- lapply(ls_gr, function(gr){
			vec <- rep(0, length(dGR))
			q2r <- as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=1L, type="any", select="all", ignore.strand=T)))
			ind <- unique(q2r[,1])
			vec[ind] <- 1
			return(vec)
		})
		df_res <- do.call(cbind, ls_vec)
		colnames(df_res) <- names(ls_gr)
		vec_res <- apply(df_res, 1, sum)
		vec_res[vec_res==0] <- NA
		dGR$Trait <- vec_res
    }
    
	##############
	# DHS: number of cell types
	##############
	if("ENCODE_DNaseI_ClusteredV3" %in% GR_annotation){
	
		if(verbose){
			message(sprintf("using the annotation '%s' (%s) ...", 'ENCODE_DNaseI_ClusteredV3', as.character(Sys.time())), appendLF=T)
		}
	
		# ENCODE_DNaseI_ClusteredV3: an GR object containing clustered peaks, along with a meta-column 'num_cells' telling how many cell types associated with a clustered peak (125 cell types)
		gr <- xRDataLoader('ENCODE_DNaseI_ClusteredV3', verbose=F, RData.location=RData.location)
		## max value
		q2r <- as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=1L, type="any", select="all", ignore.strand=T)))
		q2r$mdata <- as.data.frame(GenomicRanges::mcols(gr[q2r[,2]]))[,1]
		df_ind <- as.data.frame(q2r %>% dplyr::group_by(queryHits) %>% dplyr::summarise(max(mdata)))
		vec_res <- rep(NA, length(dGR))
		vec_res[df_ind[,1]] <- df_ind[,2]
		dGR$DHS <- vec_res
    }
    
	##############
	# TF: number of cell types * number of TFs
	##############
	if("ENCODE_TFBS_ClusteredV3" %in% GR_annotation){
	
		if(verbose){
			message(sprintf("using the annotation '%s' (%s) ...", 'ENCODE_TFBS_ClusteredV3', as.character(Sys.time())), appendLF=T)
		}
	
		# ENCODE_TFBS_ClusteredV3: a list (161 transcription factors) of GenomicRanges objects
		ls_gr <- xRDataLoader('ENCODE_TFBS_ClusteredV3', verbose=F, RData.location=RData.location)
		## number of cell types * number of TFs
		ls_vec <- lapply(ls_gr, function(gr){
			gr$cells <- sapply(strsplit(gr$cells,','), length)
			q2r <- as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=1L, type="any", select="all", ignore.strand=T)))
			## max value: number of cells
			q2r$mdata <- as.data.frame(GenomicRanges::mcols(gr[q2r[,2]]))[,1]
			df_ind <- as.data.frame(q2r %>% dplyr::group_by(queryHits) %>% dplyr::summarise(max(mdata)))
			## add TF (mcols)
			vec <- rep(0, length(dGR))
			vec[df_ind[,1]] <- df_ind[,2]
			return(vec)
		})
		# entry: number of cells
		df_res <- do.call(cbind, ls_vec)
		colnames(df_res) <- names(ls_gr)
		## add TF (mcols)
		vec_res <- apply(df_res, 1, sum)
		vec_res[vec_res==0] <- NA
		dGR$TF <- vec_res
	}
	
	##############
	# Epi: number of cell types
	##############
	if("EpigenomeAtlas_15Segments" %in% GR_annotation){
	
		if(verbose){
			message(sprintf("using the annotation '%s' (%s) ...", 'EpigenomeAtlas_15Segments', as.character(Sys.time())), appendLF=T)
		}
	
		# Roadmap Epigenomics Core 15-state Genome Segmentation data (127 cell types)
		ls_ls_gr <- xRDataLoader('EpigenomeAtlas_15Segments', verbose=F, RData.location=RData.location)
		## enhancer-like
		ls_vec <- lapply(ls_ls_gr, function(ls_gr){
			## three features containing 'enhancer'
			ind <- grepl('enhancer', names(ls_gr), ignore.case=T)
			grl <- GenomicRanges::GRangesList(ls_gr[ind])
			gr <- IRanges::reduce(BiocGenerics::unlist(grl))
	
			q2r <- as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=1L, type="any", select="all", ignore.strand=T)))
			## a vector
			ind <- unique(q2r[,1])
			vec <- rep(0, length(dGR))
			vec[ind] <- 1
			return(vec)
		})
		df_res <- do.call(cbind, ls_vec)
		colnames(df_res) <- names(ls_ls_gr)
		vec_res <- apply(df_res, 1, sum)
		vec_res[vec_res==0] <- NA
		dGR$Epi <- vec_res
    }
    
	##############
	# RecombinationRate phastCons100way phyloP100way
	##############
	ind <- match(c("RecombinationRate","phastCons100way","phyloP100way"), GR_annotation)
	vec <- GR_annotation[ind[!is.na(ind)]]
	if(length(vec)>0){
		for(j in 1:length(vec)){
		
			if(verbose){
				message(sprintf("using the annotation '%s' (%s) ...", vec[j], as.character(Sys.time())), appendLF=T)
			}
		
			x <- vec[j]
			gr <- xGRoverlap(dGR, format="GRanges", GR.score=x, verbose=F, RData.location=RData.location)
			res <- gr$GScore
			#res[is.na(res)] <- 0
			res[res<=0] <- NA
			GenomicRanges::mcols(dGR)[[names(vec)[j]]] <- res
		}
	}
	
	################################################################
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("Then, scores are combined in a rank aggregation manner (%s) ...", as.character(now)), appendLF=T)
	}
	
	df_mat <- as.data.frame(GenomicRanges::mcols(dGR))
	ind <- match(colnames(df_mat), names(GR_annotation))
	df_mat <- df_mat[,!is.na(ind)]
	
	## Rank aggregation
	res <- xAggregate(df_mat, bin=T, nbin=10, scale.log=T)
	dGR$RPS <- res$Aggregate
	####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
	
	invisible(dGR)
}
