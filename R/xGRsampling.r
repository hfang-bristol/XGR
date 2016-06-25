#' Function to generate random samples for data genomic regions from background genomic regions
#'
#' \code{xGRsampling} is supposed to randomly generate samples for data genomic regions from background genomic regions. To do so, we first identify background islands, that is, non-overlapping regions. Then, we keep only parts of data genomic regions that fall into these background islands. For each kept genomic region, a randomised region of the same length is sampled from the corresponding background islands. If required, the randomised region can be restricted to be no more than (eg 10000bp) away from data genomic regions. 
#'
#' @param GR.data an input data GR object, containing a set of genomic regions based on which to generate a null distribution
#' @param GR.background an input background GR object, containing a set of genomic regions to randomly sample from. It can be a GR list object or a list of GR objects
#' @param num.samples the number of samples randomly generated
#' @param gap.max the maximum distance of background islands to be considered away from data regions. Only background islands no far way from this distance will be considered. For example, if it is 0, meaning that only background islands that overlapp with genomic regions will be considered. By default, it is 50000
#' @param max.distance the maximum distance away from data regions that is allowed when generating random samples. By default, it is NULl meaning no such restriction
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a list of GR ojects, each containing an GR oject storing a sample.
#' @export
#' @seealso \code{\link{xGRsampling}}
#' @include xGRsampling.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location="~/Sites/SVN/github/bigdata"
#' 
#' # Enrichment analysis for GWAS SNPs from ImmunoBase
#' # a) provide input data GR object storing GWAS SNPs
#' dbSNP_GWAS <- xRDataLoader(RData.customised='dbSNP_GWAS', RData.location=RData.location)
#' 
#' # b) provide backgorund data GR object storing FANTOM5 cell-specific enhancers
#' FANTOM5_Enhancer_Cell <- xRDataLoader(RData.customised='FANTOM5_Enhancer_Cell', RData.location=RData.location)
#' 
#' # c) generate random samples as a list of GR objects
#' sGR_List <- xGRsampling(GR.data=dbSNP_GWAS, GR.background=FANTOM5_Enhancer_Cell, num.samples=1000, RData.location=RData.location)
#' }

xGRsampling <- function(GR.data, GR.background, num.samples=100, gap.max=50000, max.distance=NULL, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/Portal")
{

	if(is.null(max.distance)){
		max.distance <- gap.max
	}else if(max.distance > gap.max){
		max.distance <- gap.max
	}
	
  	## Check input GR.data and GR.background
  	if (class(GR.data) != "GRanges") {
    	stop("The function must apply to a 'GRanges' object for input data.\n")
  	}
  	
  	if(class(GR.background) == "GRangesList"){
  		GR.background <- BiocGenerics::unlist(GR.background)
  	}else if(class(GR.background) == "list"){
  		GR.background <- BiocGenerics::unlist(GenomicRanges::GRangesList(GR.background))
  	}
  	if (class(GR.background) != "GRanges") {
    	stop("The function must apply to a 'GRanges' object for input background.\n")
  	}
  	
	if(verbose){
		now <- Sys.time()
		message(sprintf("First, get non-overlapping regions for both input data and background (%s) ...", as.character(now)), appendLF=T)
	}
    ## get reduced ranges (ie non-overlapping regions)
    ### data GR
    dGR_reduced <- IRanges::reduce(GR.data)
    ### background GR
    bGR_reduced <- IRanges::reduce(GR.background)
	if(verbose){
		now <- Sys.time()
		message(sprintf("\tnon-overlapping regions: %d for data, %d for background", length(dGR_reduced), length(bGR_reduced)), appendLF=T)
	}
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Second, keep only data regions that are within background regions (%s) ...", as.character(now)), appendLF=T)
	}
	#####################################
	## A function to return an GR object storing overlapped regions (ie only overlapped regions!)
	mergeOverlaps <- function(qGR, sGR, maxgap=0L, minoverlap=1L){
		hits <- GenomicRanges::findOverlaps(query=qGR, subject=sGR, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)
		qhits <- qGR[S4Vectors::queryHits(hits)]
		shits <- sGR[S4Vectors::subjectHits(hits)]

		oGR <- IRanges::pintersect(qhits, shits)
		IRanges::reduce(oGR)
	}
	#####################################
	## update data GR after considering background
	dGR_reduced <- mergeOverlaps(qGR=dGR_reduced, sGR=bGR_reduced, maxgap=0L, minoverlap=1L)
	if(verbose){
		now <- Sys.time()
		message(sprintf("\t%d within background", length(dGR_reduced)), appendLF=T)
	}

	if(verbose){
		now <- Sys.time()
		message(sprintf("Third, find background islands that contain data regions (%s) ...", as.character(now)), appendLF=T)
	}
	## find islands
	hits <- GenomicRanges::findOverlaps(query=dGR_reduced, subject=GR.background, maxgap=gap.max, minoverlap=1L, type="any", select="all", ignore.strand=T)
	ind_data <- S4Vectors::queryHits(hits)
	ind_background <- S4Vectors::subjectHits(hits)
	dt_ls <- split(x=ind_background, f=ind_data)
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("\t%d background islands", length(unique(ind_background))), appendLF=T)
	}
    
    if(is.null(max.distance)){
		if(verbose){
			now <- Sys.time()
			message(sprintf("Fourth, define the sampling range (%s) ...", as.character(now)), appendLF=T)
		}
    }else{
		max.distance <- as.integer(max.distance)
		if(verbose){
			now <- Sys.time()
			message(sprintf("Fourth, define the sampling range within %d bp distance away from data regions (%s) ...", max.distance, as.character(now)), appendLF=T)
		}
    }
	## convert into data.frame
    df_data <- GenomicRanges::as.data.frame(dGR_reduced, row.names=NULL)
    df_background <- GenomicRanges::as.data.frame(GR.background, row.names=NULL)
    range_ls <- lapply(1:length(dt_ls), function(i){
    	df_dt <- df_data[as.numeric(names(dt_ls)[i]), ]
    	# data width
    	dw <- df_dt$width
    	
    	df_bg <- df_background[dt_ls[[i]], ]
    	res <- lapply(1:nrow(df_bg), function(j){
			# background start and end
			bs <- df_bg$start[j]
			be <- df_bg$end[j]
			
			# range start and end
			rs <- bs
			re <- be-dw
			# within the distance restriction
			if(!is.null(max.distance)){
				if(rs < df_dt$start - max.distance){
					rs <- df_dt$start - max.distance
				}
				if(re > df_dt$end + max.distance){
					re <- df_dt$end + max.distance
				}
			}
			seq(rs, re)
    	})
    	unlist(res)
    })
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Fifth, do '%d' sampling for the starting points (%s) ...", num.samples, as.character(now)), appendLF=T)
	}
    res_ls <- lapply(range_ls, function(x){
    	base::sample(x, num.samples, replace=T)
    })
    df_samples <- do.call(rbind, res_ls)
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Finally, construct GR for each sampling (%s) ...", as.character(now)), appendLF=T)
	}
	## 'df_dt_all' for all data regions   
	df_dt_all <- df_data[as.numeric(names(dt_ls)), ]
    sGR_list <- lapply(1:ncol(df_samples), function(j){
		sGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(df_dt_all$seqnames),
			ranges = IRanges::IRanges(start=df_samples[,j], end=df_samples[,j]+df_dt_all$width-1),
			strand = S4Vectors::Rle(df_dt_all$strand)
		)
    })

	invisible(sGR_list)
}
