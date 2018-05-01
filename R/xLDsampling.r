#' Function to generate randomly sampled LD blocks
#'
#' \code{xLDsampling} is supposed to generate randomly sampled LD blocks. A sample block has the same boundary range as the observed, and can respect maf of the best SNP, and/or distance of the best SNP to the nearest gene. Also it can be restricted to the same chromosome. For each null LD block, it can preserve the boundary only or exactly preserve the relative SNP locations. It returns a GRL object.
#'
#' @param bLD a bLD object, containing a set of blocks based on which to generate a null distribution
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 150) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 150) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised GR object directly
#' @param num.samples the number of samples randomly generated
#' @param respect how to respect the properties of to-be-sampled LD blocks. It can be one of 'maf' (respecting the maf of the best SNP), 'distance' (respecting the distance of the best SNP to the nearest gene), and 'both' (respecting the maf and distance)
#' @param restrict.chr logical to restrict to the same chromosome. By default, it sets to false
#' @param preserve how to preserve the resulting null LD block. It can be one of 'boundary' (preserving the boundary of the LD block), and 'exact' (exactly preserving the relative SNP locations within the LD block). Notably, no huge difference for the boundary preserving when enrichment analysis invovles region-based genomic annotations, but it may make difference when genomic annatations are largely SNP-based (such as eQTLs)
#' @param seed an integer specifying the seed
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a GRL object, each containing an GR oject storing an instance of sampled blocks (with a meta-column 'best' for the identity, and a meta-column 'B' for the instance sequence).
#' @export
#' @seealso \code{\link{xLDsampling}}
#' @include xLDsampling.r
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
#' # c) generate random samples as a GRL object
#' grl <- xLDsampling(bLD, GR.SNP="dbSNP_GWAS", num.samples=2000, RData.location=RData.location)
#' 
#' ##########################
#' ## Advanced use: customised GR.SNP
#' ##########################
#' GR.SNP <- xRDataLoader("dbSNP_GWAS", RData.location=RData.location)
#' grl <- xLDsampling(bLD, GR.SNP=GR.SNP, respect="both", restrict.chr=T, preserve="exact", RData.location=RData.location)
#' }

xLDsampling <- function(bLD, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), num.samples=2000, respect=c("maf","distance","both"), restrict.chr=F, preserve=c("boundary","exact"), seed=825, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    startT <- Sys.time()
    if(verbose){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }
	####################################################################################
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    respect <- match.arg(respect)
    preserve <- match.arg(preserve)
	
	if(class(bLD) == "bLD"){
		gr_best <- bLD$best
		grl_block <- bLD$block
	}else{
		return(NULL)
	}
	
	if(class(gr_best) == "GRanges"){
		
		if(class(GR.SNP) == "GRanges"){
			gr_bg <- GR.SNP
			
		}else{
			
			GR.SNP <- GR.SNP[1]
			
			if(verbose){
				message(sprintf("Load to-be-sampled SNPs from '%s' (%s) ...", GR.SNP, as.character(Sys.time())), appendLF=T)
			}
		
			gr_bg <- xRDataLoader(GR.SNP, RData.location=RData.location, verbose=F)
			if(is.null(gr_bg)){
				GR.SNP <- "dbSNP_GWAS"
				if(verbose){
					message(sprintf("\tinstead, %s will be used", GR.SNP), appendLF=T)
				}
				gr_bg <- xRDataLoader(GR.SNP, RData.location=RData.location, verbose=F)
			}
		}
		
		#####################################
		ind <- which(!is.na(gr_bg$maf))
		df_bg <- GenomicRanges::as.data.frame(gr_bg[ind], row.names=NULL)
		df_best <- GenomicRanges::as.data.frame(gr_best)
		
		if(verbose){
			message(sprintf("%d blocks sampled from %d SNPs over %d times (%s)...", nrow(df_best), nrow(df_bg), num.samples, as.character(Sys.time())), appendLF=T)
			
			tmp_respect <- respect
			if(respect=='both'){
				tmp_respect <- 'maf and distance'
			}
			message(sprintf("\trespect '%s'", tmp_respect), appendLF=T)
			
			if(restrict.chr){
				message(sprintf("\trestrict to the same chromosome"), appendLF=T)
			}
			
			message(sprintf("\tpreserve '%s'", preserve), appendLF=T)
		}
		#####################################
		
		system.time({
		set.seed(seed)
		ls_df <- lapply(1:nrow(df_best), function(i){
			
			if(restrict.chr){
				ind <- which(df_bg$seqnames %in% df_best$seqnames[i])
				df_sample <- df_bg[ind,]
			}else{
				df_sample <- df_bg
			}
			
			if(any(respect %in% c("maf","both"))){
				## maf: +/- 0.0025 (0.5/0.005)
				x_maf <- df_best$maf[i]
				ind_maf <- which(df_sample$maf >= (x_maf-0.0025) & df_sample$maf <= (x_maf+0.0025))
			}
			
			if(any(respect %in% c("distance","both"))){
				## distance: +/- 12500 (2500000/25000)
				x_distance <- df_best$distance[i]
				ind_distance <- which(df_sample$distance >= (x_distance-12500) & df_sample$distance <= (x_distance+12500))
			}
			
			if(respect=="maf"){
				ind <- ind_maf
			}else if(respect=="distance"){
				ind <- ind_distance
			}else{
				ind <- intersect(ind_maf, ind_distance)
			}
			
			## ind_sample
			ind_sample <- base::sample(ind, num.samples, replace=T)
			
			if(preserve=="boundary"){
				y <- df_sample$start[ind_sample]
				res <- data.frame(chr=df_sample$seqnames[ind_sample], start=y+df_best$upstream[i], end=y+df_best$downstream[i], best=rownames(df_best)[i], B=1:num.samples, stringsAsFactors=F)
			}else if(preserve=="exact"){		
				distance_to_best <- grl_block[[i]]$distance_to_best
				ind_sample_rep <- rep(ind_sample, each=length(distance_to_best))
				y <- df_sample$start[ind_sample_rep]+rep(distance_to_best,num.samples)
				res <- data.frame(chr=df_sample$seqnames[ind_sample_rep], start=y, end=y, best=rownames(df_best)[i], B=rep(1:num.samples,each=length(distance_to_best)), stringsAsFactors=F)
			}
		})
		df <- do.call(rbind, ls_df)
		})
		
		system.time({
		#gr <- xGR(df, format='data.frame', add.name=F)
		## construct data GR
		gr <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(df[,1]),
			ranges = IRanges::IRanges(start=df[,2], end=df[,3]),
			strand = S4Vectors::Rle(rep('*',nrow(df)))
		)
		gr$best <- df$best
		gr$B <- df$B
		grl <- GenomicRanges::split(gr, gr$B)
		})
		
	}
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(verbose){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }

	invisible(grl)
}
