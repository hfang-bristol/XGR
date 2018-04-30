#' Function to calculate per base scores given a list of genomic regions in terms of overlaps with genomic annotations
#'
#' \code{xGScoreAdv} is supposed to calculate per base scores for an input list of genomic regions (genome build 19), using genomic annotations (eg genomic segments, active chromatin, transcription factor binding sites/motifs, conserved sites). The per base scores are calculated for overlaps with each genomic annotation. Scores for genomic regions/variants can be constraint/conservation or impact/pathogenicity.
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param GS.annotation which genomic scores (GS) annotaions used. It can be 'fitCons' (the probability of fitness consequences for point mutations; \url{http://www.ncbi.nlm.nih.gov/pubmed/25599402}), 'phastCons' (the probability that each nucleotide belongs to a conserved element/negative selection [0,1]), 'phyloP' (conservation at individual sites representing -log p-values under a null hypothesis of neutral evolution, positive scores for conservation and negative scores for acceleration), 'mcap' (eliminating a majority of variants with uncertain significance in clinical exomes at high sensitivity: \url{http://www.ncbi.nlm.nih.gov/pubmed/27776117}), and 'cadd' (combined annotation dependent depletion for estimating relative levels of pathogenicity of potential human variants: \url{http://www.ncbi.nlm.nih.gov/pubmed/24487276})
#' @param GR.annotation the genomic regions of annotation data. By default, it is 'NA' to disable this option. Pre-built genomic annotation data are detailed in \code{\link{xDefineGenomicAnno}}. Alternatively, the user can also directly provide a customised GR object (or a list of GR objects)
#' @param details logical to indicate whether the detailed information (ie ratio) is returned. By default, it sets to false for no inclusion
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a data frame with 6 columns:
#' \itemize{
#'  \item{\code{name}: the annotation name}
#'  \item{\code{o_nBase}: the number of bases overlapped between input regions and annotation regions}
#'  \item{\code{o_GS}: the per base genomic scores for overlaps between input regions and annotation regions}
#'  \item{\code{a_nBase}: the number of bases covered by that annotation; optional, it is only appended when "details" is true}
#'  \item{\code{a_GS}: the per base genomic scores for that annotation; optional, it is only appended when "details" is true}
#'  \item{\code{ratio}: ratio of o_GS divided by a_GS; optional, it is only appended when "details" is true}
#' }
#' @note Pre-built genomic annotation data are detailed in \code{\link{xDefineGenomicAnno}}.
#' @export
#' @seealso \code{\link{xGScore}}
#' @include xGScoreAdv.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#'
#' # a) provide the genomic regions
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS
#' data <- ImmunoBase$AS$variant
#'
#' # b) in terms of overlaps with genomic segments (Primary monocytes from peripheral blood)
#' ## fitness consequence score 
#' res_df <- xGScoreAdv(data=data, format="GRanges", GS.annotation="fitCons", GR.annotation="EpigenomeAtlas_15Segments_E029", RData.location=RData.location)
#' ## phastCons conservation score 
#' res_df <- xGScoreAdv(data=data, format="GRanges", GS.annotation="phastCons", GR.annotation="EpigenomeAtlas_15Segments_E029", RData.location=RData.location)
#' 
#' # c) in terms of overlaps with genic annotations
#' ## phyloP conservation score 
#' res_df <- xGScoreAdv(data=data, format="GRanges", GS.annotation="phyloP", GR.annotation="Genic_anno", RData.location=RData.location)
#' }

xGScoreAdv <- function(data, format=c("data.frame", "bed", "chr:start-end", "GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), GS.annotation=c("fitCons","phastCons","phyloP","mcap","cadd"), GR.annotation=NA, details=F, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata_dev")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    GS.annotation <- match.arg(GS.annotation)
    
	if(class(GR.annotation) == "GRanges"){
		###################################
		## now GR.annotation can be directly provided as a GR object
		###################################
		aGR <- GR.annotation
	}else{
		aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location)
		aGR <- lapply(aGRL, function(x) x)
	}
    
	#####################################
	## A function to return an GR object storing overlapped regions (ie only overlapped regions!)
	mergeOverlaps <- function(qGR, sGR, maxgap=-1L, minoverlap=0L){
		hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query=qGR, subject=sGR, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
		qhits <- qGR[hits[,1]]
		shits <- sGR[hits[,2]]

		oGR <- IRanges::pintersect(qhits, shits, ignore.strand=T)
		#IRanges::reduce(oGR)
	}
	#####################################
	
	dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
	
	ls_df <- lapply(1:length(aGR), function(i){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Calculating '%s' for the annotation '%s' (%s) ...", GS.annotation, names(aGR)[i], as.character(now)), appendLF=TRUE)
		}
		gr <- aGR[[i]]
		
		## overlaps
		gr_overlap <- mergeOverlaps(dGR, gr)
		gr_score_overlap <- xGScore(data=gr_overlap, format="GRanges", GS.annotation=GS.annotation, scoring.scheme="sum", RData.location=RData.location, verbose=FALSE)
		ind <- which(!is.na(gr_score_overlap$GScore))
		o_nBase <- sum(as.numeric(IRanges::width(gr_score_overlap[ind])))
		o_GScore <- sum(gr_score_overlap$GScore[ind])
		o_GScore_per_base <- o_GScore / o_nBase
		#o_GScore_per_base[is.na(o_GScore_per_base)] <- 0
		
		if(details){
			## annotation
			gr_score_bg <- xGScore(data=gr, format="GRanges", GS.annotation=GS.annotation, scoring.scheme="sum", RData.location=RData.location, verbose=FALSE)
			ind <- which(!is.na(gr_score_bg$GScore))
			a_nBase <- sum(as.numeric(IRanges::width(gr_score_bg[ind])))
			a_GScore <- sum(gr_score_bg$GScore[ind])
			a_GScore_per_base <- a_GScore / a_nBase
			
			## output
			res <- data.frame(name=names(aGR)[i], o_nBase, o_GS=o_GScore_per_base, a_nBase, a_GS=a_GScore_per_base, ratio=o_GScore_per_base/a_GScore_per_base, stringsAsFactors=F)
			
		}else{
			## output
			res <- data.frame(name=names(aGR)[i], o_nBase, o_GS=o_GScore_per_base, stringsAsFactors=F)
		}
		
		return(res)
	})
	res_df <- do.call(rbind, ls_df)
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
	
	invisible(res_df)
}
