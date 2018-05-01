#' Function to conduct region-based enrichment analysis using genomic annotations via sampling
#'
#' \code{xGRviaGenomicAnnoAdv} is supposed to conduct region-based enrichment analysis for the input genomic region data (genome build h19), using genomic annotations (eg active chromatin, transcription factor binding sites/motifs, conserved sites). Enrichment analysis is achieved by comparing the observed overlaps against the expected overlaps which are estimated from the null distribution. The null distribution is generated via sampling, that is, randomly generating samples for data genomic regions from background genomic regions. Background genomic regions can be provided by the user; by default, the annotatable genomic regions will be used. 
#'
#' @param data.file an input data file, containing a list of genomic regions to test. If the input file is formatted as a 'data.frame' (specified by the parameter 'format.file' below), the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. If the format is indicated as "chr:start-end", instead of using the first 3 columns, only the first column will be used and processed. If the file also contains other columns, these additional columns will be ignored. Alternatively, the input file can be the content itself assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns.
#' @param annotation.file an input annotation file containing genomic annotations for genomic regions. If the input file is formatted as a 'data.frame', the first four columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), the ending chromosome position (3rd column), and the genomic annotations (eg transcription factors and histones; 4th column). If the format is indicated as 'bed', the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the format is indicated as "chr:start-end", the first two columns correspond to the chromosome:start-end (1st column) and the genomic annotations (eg transcription factors and histones; 2nd column). If the file also contains other columns, these additional columns will be ignored. Alternatively, the input file can be the content itself assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns.
#' @param background.file an input background file containing a list of genomic regions as the test background. The file format is the same as 'data.file'. By default, it is NULL meaning all annotatable bases (ig non-redundant bases covered by 'annotation.file') are used as background. However, if only one annotation (eg only a transcription factor) is provided in 'annotation.file', the background must be provided.
#' @param format.file the format for input files. It can be one of "data.frame", "chr:start-end", "bed" and "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so).
#' @param background.annotatable.only logical to indicate whether the background is further restricted to annotatable bases (covered by 'annotation.file'). In other words, if the background is provided, the background bases are those after being overlapped with annotatable bases. Notably, if only one annotation (eg only a transcription factor) is provided in 'annotation.file', it should be false.
#' @param num.samples the number of samples randomly generated
#' @param gap.max the maximum distance of background islands to be considered away from data regions. Only background islands no far way from this distance will be considered. For example, if it is 0, meaning that only background islands that overlapp with genomic regions will be considered. By default, it is 50000
#' @param max.distance the maximum distance away from data regions that is allowed when generating random samples. By default, it is NULl meaning no such restriction
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param GR.annotation the genomic regions of annotation data. By default, it is 'NA' to disable this option. Pre-built genomic annotation data are detailed in \code{\link{xDefineGenomicAnno}}. Alternatively, the user can also directly provide a customised GR object (or a list of GR objects)
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. It will depend on whether these two packages "foreach" and "doParallel" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doParallel"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a data frame with 8 columns:
#' \itemize{
#'  \item{\code{name}: the annotation name}
#'  \item{\code{nAnno}: the number of bases covered by that annotation. If the background is provided, they are also restricted by this}
#'  \item{\code{nOverlap}: the number of bases overlapped between input regions and annotation regions. If the background is provided, they are also restricted by this}
#'  \item{\code{fc}: fold change}
#'  \item{\code{zscore}: z-score}
#'  \item{\code{pvalue}: p-value}
#'  \item{\code{adjp}: adjusted p-value. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{nData}: the number of bases covered by input regions}
#'  \item{\code{nBG}: the number of bases covered by background regions}
#' }
#' @note Pre-built genomic annotation data are detailed in \code{\link{xDefineGenomicAnno}}.
#' @export
#' @seealso \code{\link{xDefineGenomicAnno}}
#' @include xGRviaGenomicAnnoAdv.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' 
#' # Enrichment analysis for GWAS SNPs from ImmunoBase
#' ## a) provide input data
#' data.file <- "http://galahad.well.ox.ac.uk/bigdata/ImmunoBase_GWAS.bed"
#' 
#' ## b) perform enrichment analysis using FANTOM expressed enhancers
#' eTerm <- xGRviaGenomicAnnoAdv(data.file=data.file, format.file="bed", GR.annotation="FANTOM5_Enhancer_Cell", num.samples=1000, gap.max=50000, RData.location=RData.location)
#'
#' ## c) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' ## d) barplot of enriched terms
#' bp <- xEnrichBarplot(eTerm, top_num='auto', displayBy="fdr")
#' bp
#'
#' ## e) save enrichment results to the file called 'Regions_enrichments.txt'
#' output <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' utils::write.table(output, file="Regions_enrichments.txt", sep="\t", row.names=FALSE)
#' }

xGRviaGenomicAnnoAdv <- function(data.file, annotation.file=NULL, background.file=NULL, format.file=c("data.frame", "bed", "chr:start-end", "GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), background.annotatable.only=F, num.samples=1000, gap.max=50000, max.distance=NULL, p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), GR.annotation=NA, parallel=TRUE, multicores=NULL, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format.file <- match.arg(format.file)
    build.conversion <- match.arg(build.conversion)
    p.adjust.method <- match.arg(p.adjust.method)
    
    ###################
	if(verbose){
		now <- Sys.time()
		message(sprintf("First, import the files formatted as '%s' (%s) ...", format.file, as.character(now)), appendLF=T)
	}
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\timport the data file (%s) ...", as.character(now)), appendLF=T)
	}
    ## import data file
    if(is.matrix(data.file) | is.data.frame(data.file) | class(data.file)=="GRanges"){
        data <- data.file
    }else if(!is.null(data.file) & any(!is.na(data.file))){
    	if(length(data.file)==1){
			data <- utils::read.delim(file=data.file, header=F, row.names=NULL, stringsAsFactors=F)
			#data <- unique(data[,1])
		}else{
			data <- data.file
		}
    }else{
    	stop("The file 'data.file' must be provided!\n")
    }
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\timport the annotation file (%s) ...", as.character(now)), appendLF=T)
	}
    ## import annotation file
    if(is.matrix(annotation.file) | is.data.frame(annotation.file) | class(annotation.file)=="list"){
        annotation <- annotation.file
    }else if(!is.null(annotation.file)){
		annotation <- utils::read.delim(file=annotation.file, header=F, row.names=NULL, stringsAsFactors=F)
    }else{
    	message("\t\tThe file 'annotation.file' is not provided, so built-in RData will be used instead!")
    	annotation <- NULL
    }
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\timport the background file (%s) ...", as.character(now)), appendLF=T)
	}
	## import background file
    if(is.matrix(background.file) | is.data.frame(background.file) | class(background.file)=="GRanges"){
        background <- background.file
    }else if(!is.null(background.file)){
    	if(length(background.file)==1){
			background <- utils::read.delim(file=background.file, header=F, row.names=NULL, stringsAsFactors=F)
			background <- unique(background[,1])
		}else{
			background <- background.file
		}
    }else{
    	background <- NULL
    }
    
    ###################
	if(verbose){
		now <- Sys.time()
		message(sprintf("Second, construct GenomicRanges object (%s) ...", as.character(now)), appendLF=T)
	}
    
	if(format.file=="data.frame"){
		## construct data GR
		if(ncol(data)>=3){
			data <- data
		}else if(ncol(data)==2){
			data <- cbind(data, data[,2])
		}else{
			stop("Your input 'data.file' is not as expected!\n")
		}
		## make sure positions are numeric
		ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
		data <- data[ind,]
		dGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(data[,1]),
			ranges = IRanges::IRanges(start=as.numeric(data[,2]), end=as.numeric(data[,3])),
			strand = S4Vectors::Rle(rep('*',nrow(data)))
		)
		
		if(!is.null(annotation)){
			## construct annotation GR
			if(ncol(annotation)>=4){
					annotation <- annotation
			}else{
				stop("Your input 'annotation.file' is not as expected!\n")
			}
			anno_ls <- split(x=annotation[,-4], f=annotation[,4])
			aGR <- base::lapply(anno_ls, function(x){
				## make sure positions are numeric
				ind <- suppressWarnings(which(!is.na(as.numeric(x[,2])) & !is.na(as.numeric(x[,3]))))
				x <- x[ind,]
				gr <- GenomicRanges::GRanges(
					seqnames=S4Vectors::Rle(x[,1]),
					ranges = IRanges::IRanges(start=as.numeric(x[,2]), end=as.numeric(x[,3])),
					strand = S4Vectors::Rle(rep('*',nrow(x)))
				)
			})
		}else{
			aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location)
			aGR <- lapply(aGRL, function(x) x)
		}
		
		if(!is.null(background)){
			## construct background GR
			if(ncol(background)>=3){
				background <- background
			}else if(ncol(background)==2){
				background <- cbind(background, background[,2])
			}else{
				stop("Your input 'background.file' is not as expected!\n")
			}
			## make sure positions are numeric
			ind <- suppressWarnings(which(!is.na(as.numeric(background[,2])) & !is.na(as.numeric(background[,3]))))
			background <- background[ind,]
			bGR <- GenomicRanges::GRanges(
				seqnames=S4Vectors::Rle(background[,1]),
				ranges = IRanges::IRanges(start=as.numeric(background[,2]), end=as.numeric(background[,3])),
				strand = S4Vectors::Rle(rep('*',nrow(background)))
			)
		}else{
			bGR <- NULL
		}
		
	}else if(format.file=="chr:start-end"){
		
		## construct data GR
		input <- do.call(rbind, strsplit(data[,1], ":|-"))
		if(ncol(input)>=3){
			data <- input
		}else if(ncol(input)==2){
			data <- cbind(input, input[,2])
		}else{
			stop("Your input 'data.file' does not meet the format 'chr:start-end'!\n")
		}
		## make sure positions are numeric
		ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
		data <- data[ind,]
		dGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(data[,1]),
			ranges = IRanges::IRanges(start=as.numeric(data[,2]), end=as.numeric(data[,3])),
			strand = S4Vectors::Rle(rep('*',nrow(data)))
		)
		
		if(!is.null(annotation)){
			## construct annotation GR
			input <- do.call(rbind, strsplit(annotation[,1], ":|-"))
			if(ncol(input)>=3){
				annotation <- cbind(input[,1:3], annotation[,2])
			}else if(ncol(input)==2){
				annotation <- cbind(input[,c(1,2,2)], annotation[,2])
			}else{
				stop("Your input 'annotation.file' does not meet the format 'chr:start-end'!\n")
			}
			anno_ls <- split(x=annotation[,-4], f=annotation[,4])
			aGR <- base::lapply(anno_ls, function(x){
				## make sure positions are numeric
				ind <- suppressWarnings(which(!is.na(as.numeric(x[,2])) & !is.na(as.numeric(x[,3]))))
				x <- x[ind,]
				gr <- GenomicRanges::GRanges(
					seqnames=S4Vectors::Rle(x[,1]),
					ranges = IRanges::IRanges(start=as.numeric(x[,2]), end=as.numeric(x[,3])),
					strand = S4Vectors::Rle(rep('*',nrow(x)))
				)
			})
		}else{
			aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location)
			aGR <- lapply(aGRL, function(x) x)
		}
		
		if(!is.null(background)){
			## construct background GR
			input <- do.call(rbind, strsplit(background[,1], ":|-"))
			if(ncol(input)>=3){
				background <- input
			}else if(ncol(input)==2){
				background <- cbind(input, input[,2])
			}else{
				stop("Your input 'background.file' does not meet the format 'chr:start-end'!\n")
			}
			## make sure positions are numeric
			ind <- suppressWarnings(which(!is.na(as.numeric(background[,2])) & !is.na(as.numeric(background[,3]))))
			background <- background[ind,]
			bGR <- GenomicRanges::GRanges(
				seqnames=S4Vectors::Rle(background[,1]),
				ranges = IRanges::IRanges(start=as.numeric(background[,2]+1), end=as.numeric(background[,3])),
				strand = S4Vectors::Rle(rep('*',nrow(data)))
			)
		}else{
			bGR <- NULL
		}
		
	}else if(format.file=="bed"){
		## construct data GR
		## make sure positions are numeric
		ind <- suppressWarnings(which(!is.na(as.numeric(data[,2])) & !is.na(as.numeric(data[,3]))))
		data <- data[ind,]
		dGR <- GenomicRanges::GRanges(
			seqnames=S4Vectors::Rle(data[,1]),
			ranges = IRanges::IRanges(start=as.numeric(data[,2])+1, end=as.numeric(data[,3])),
			strand = S4Vectors::Rle(rep('*',nrow(data)))
		)
		
		if(!is.null(annotation)){
			## construct annotation GR
			anno_ls <- split(x=annotation[,-4], f=annotation[,4])
			aGR <- base::lapply(anno_ls, function(x){
				## make sure positions are numeric
				ind <- suppressWarnings(which(!is.na(as.numeric(x[,2])) & !is.na(as.numeric(x[,3]))))
				x <- x[ind,]
				gr <- GenomicRanges::GRanges(
					seqnames=S4Vectors::Rle(x[,1]),
					ranges = IRanges::IRanges(start=as.numeric(x[,2])+1, end=as.numeric(x[,3])),
					strand = S4Vectors::Rle(rep('*',nrow(x)))
				)
			})
		}else{
			aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location)
			aGR <- lapply(aGRL, function(x) x)
		}
		
		if(!is.null(background)){
			## construct background GR
			## make sure positions are numeric
			ind <- suppressWarnings(which(!is.na(as.numeric(background[,2])) & !is.na(as.numeric(background[,3]))))
			background <- background[ind,]
			bGR <- GenomicRanges::GRanges(
				seqnames=S4Vectors::Rle(background[,1]),
				ranges = IRanges::IRanges(start=as.numeric(background[,2])+1, end=as.numeric(background[,3])),
				strand = S4Vectors::Rle(rep('*',nrow(data)))
			)
		}else{
			bGR <- NULL
		}
		
	}else if(format.file=="GRanges"){
		## construct data GR
		dGR <- data
		
		if(!is.null(annotation)){
			## construct annotation GR
			aGR <- annotation
		}else{
			aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location)
			aGR <- lapply(aGRL, function(x) x)
		}
		
		if(!is.null(background)){
			## construct background GR
			bGR <- background
		}else{
			bGR <- NULL
		}
		
	}
	
	#####################################
	## A function to return an GR object storing overlapped regions (ie only overlapped regions!)
	mergeOverlaps <- function(qGR, sGR, out.format=c("GR","counts"), maxgap=-1L, minoverlap=0L){
		out.format <- match.arg(out.format)
		
		hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query=qGR, subject=sGR, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
		qhits <- qGR[hits[,1]]
		shits <- sGR[hits[,2]]
		gr <- IRanges::pintersect(qhits, shits, ignore.strand=T)
		if(out.format=='GR'){
			gr
		}else{
			sum(as.numeric(IRanges::width(gr)))
		}
	}
	
	## A function to return a list of GR objects storing overlapped regions (ie only overlapped regions!)
	mergeOverlaps_GRL <- function(qGR, sGRL, maxgap=-1L, minoverlap=0L, out.format=c("GR", "counts")){
		out.format <- match.arg(out.format)
		
		hits_GRL <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query=qGR, subject=sGRL, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
		qhits_GRL <- hits_GRL[,1]
		shits_GRL <- hits_GRL[,2]
		
		res_ls <- split(x=qhits_GRL, f=shits_GRL)
		names_ls <- as.numeric(names(res_ls))
		if(1){
			res <- lapply(1:length(res_ls), function(i){
				s_ind <- names_ls[i]
				sGR <- sGRL[[s_ind]]
				#q_ind <- res_ls[[i]]
				#qGR_sub <- qGR[q_ind]
			
				mergeOverlaps(qGR=qGR, sGR=sGR, out.format=out.format, maxgap=maxgap, minoverlap=minoverlap)
			})
		}else{
			sGRL_o <- sGRL[names_ls]
			res <- lapply(sGRL_o, function(x){
				mergeOverlaps(qGR=qGR, sGR=x, out.format=out.format, maxgap=maxgap, minoverlap=minoverlap)
			})
		}
		
		if(out.format=='GR'){
			out <- vector('list', length(sGRL))
			names(out) <- names(sGRL)
			out[names_ls] <- res
			out <- lapply(out, function(x){
				if(is.null(x)){
					x <- GenomicRanges::GRanges()
				}else{
					x
				}
			})
		}else{
			out <- as.list(rep(0, length(sGRL)))
			names(out) <- names(sGRL)
			out[names_ls] <- res
		}
		out
	}
	
    # A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d (%s)", i, B, as.character(Sys.time())), appendLF=T)
            }
        }
    }
	#####################################
    
    ############
    # lift over
    ############
    if(!is.na(build.conversion)){
    	## dGR
		if(verbose){
			message(sprintf("\tdata genomic regions: lifted over via genome build conversion `%s`", build.conversion), appendLF=T)
		}
		dGR <- xLiftOver(data.file=dGR, format.file="GRanges", build.conversion=build.conversion, merged=F, verbose=verbose, RData.location=RData.location)
    	## aGR
    	if(!is.null(annotation.file)){
			if(verbose){
				message(sprintf("\tannotation genomic regions: lifted over via genome build conversion `%s`", build.conversion), appendLF=T)
			}
			aGR <- lapply(aGR, function(gr){
				xLiftOver(data.file=gr, format.file="GRanges", build.conversion=build.conversion, merged=F, verbose=verbose, RData.location=RData.location)
			})
		}
    	## bGR
		if(!is.null(bGR)){
			if(verbose){
				message(sprintf("\tbackground genomic regions: lifted over via genome build conversion `%s`", build.conversion), appendLF=T)
			}
			bGR <- xLiftOver(data.file=bGR, format.file="GRanges", build.conversion=build.conversion, merged=F, verbose=verbose, RData.location=RData.location)
		}
	}
    ############
    ############
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Third, define the background (%s) ...", as.character(now)), appendLF=T)
	}
    
    ## get reduced ranges (ie non-overlapping regions)
    ### data GR
    dGR_reduced <- IRanges::reduce(dGR)
    ### annotation GR
	aGR_reduced <- base::lapply(aGR, function(x){
		IRanges::reduce(x)
	})
	### define background GR
	if(is.null(bGR)){
		if(verbose){
			now <- Sys.time()
			message(sprintf("\tall annotatable regions (by default) are used as the background (%s) ...", as.character(now)), appendLF=T)
		}
	
		aGRL <- GenomicRanges::GRangesList(aGR_reduced)
		bGR_reduced <- IRanges::reduce(BiocGenerics::unlist(aGRL))
	}else{
		bGR_reduced <- IRanges::reduce(bGR)
	
		## update annotation GR after considering background
		aGR_reduced <- base::lapply(aGR_reduced, function(gr){
			mergeOverlaps(qGR=gr, sGR=bGR_reduced, maxgap=-1L, minoverlap=0L)
		})
	
		## restrict to the annotatable only?
		if(background.annotatable.only){
			if(verbose){
				now <- Sys.time()
				message(sprintf("\tthe given background regions but restricted to the annotatable are used as the background (%s) ...", as.character(now)), appendLF=T)
			}
		
			## update background GR
			aGRL <- GenomicRanges::GRangesList(aGR_reduced)
			bGR_reduced <- IRanges::reduce(BiocGenerics::unlist(aGRL))
		}else{
			if(verbose){
				now <- Sys.time()
				message(sprintf("\tthe given background regions are used as the background (%s) ...", as.character(now)), appendLF=T)
			}
		}
	}
	
	## update data GR after considering background
	dGR_reduced <- mergeOverlaps(qGR=dGR_reduced, sGR=bGR_reduced, maxgap=-1L, minoverlap=0L)
	
	## find overlap GR between annotation GR and data GR
	aGRL <- GenomicRanges::GRangesList(aGR_reduced)
	oGR_reduced <- mergeOverlaps_GRL(qGR=dGR_reduced, sGRL=aGRL, maxgap=-1L, minoverlap=0L)
	#######################################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("Forth, generate null distribution via doing %d sampling (%s) ...", num.samples, as.character(now)), appendLF=T)
	}
	
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xGRsampling' is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
	sGR_list <- xGRsampling(GR.data=dGR_reduced, GR.background=bGR_reduced, num.samples=num.samples, gap.max=gap.max, max.distance=max.distance, verbose=verbose, RData.location=RData.location)
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xGRsampling' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
	
	#######################################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("Fifth, perform enrichment analysis (%s) ...", as.character(now)), appendLF=T)
	}
	## prepare enrichment analysis
	data_nBases <- sum(as.numeric(IRanges::width(dGR_reduced)))
	annotation_nBases <- base::sapply(aGR_reduced, function(gr){
		sum(as.numeric(IRanges::width(gr)))
	})
	background_nBases <- sum(as.numeric(IRanges::width(bGR_reduced)))
	overlap_nBases <- base::sapply(oGR_reduced, function(gr){
		sum(as.numeric(IRanges::width(gr)))
	})

	if(verbose){
		now <- Sys.time()
		message(sprintf("\tthe number of nucleotides/bases: data (%d), background (%d)", data_nBases, background_nBases), appendLF=T)
		message(sprintf("\tthe number of annotations: %d", length(annotation_nBases)), appendLF=T)
	}

    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
        flag_parallel <- xCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            i <- 1
            b2f <- foreach::`%dopar%` (foreach::foreach(i=1:length(sGR_list), .inorder=T, .combine=rbind), {
				progress_indicate(i, length(sGR_list), 10, flag=T)
				res <- mergeOverlaps_GRL(qGR=sGR_list[[i]], sGRL=aGRL, out.format="counts")
				unlist(res)
            })
        }
    }
    ###### non-parallel computing
    if(flag_parallel==F){
		## find overlap GR between annotation GR and sample GR
		b_overlap_nBases <- base::lapply(1:length(sGR_list), function(i){
			progress_indicate(i, length(sGR_list), 10, flag=T)
			res <- mergeOverlaps_GRL(qGR=sGR_list[[i]], sGRL=aGRL, out.format="counts")
			unlist(res)
		})
		b2f <- do.call(rbind, b_overlap_nBases)
    }
	
	############################################################
    obs <- overlap_nBases	
	exp_mean <- apply(b2f, 2, mean)
	exp_std <- apply(b2f, 2, stats::sd)
	
	## ratio
	ratio <- obs/exp_mean
	
    ## for zscore
    zscore <- (obs-exp_mean)/exp_std

    ## for pvalue
    obs_matrix <- matrix(rep(obs,each=num.samples), nrow=num.samples)
    pvalue <- apply((obs_matrix - b2f)<=0, 2, sum) / num.samples
	####################
	
    zscore[is.na(zscore)] <- 0
    zscore[is.infinite(zscore)] <- max(zscore[!is.infinite(zscore)])
    pvalue[is.na(ratio)] <- 1
    ratio[is.na(ratio)] <- 1
 
	enrichment_df <- data.frame(names(overlap_nBases), annotation_nBases, data_nBases, background_nBases, obs, exp_mean, ratio, zscore, pvalue, row.names=NULL, stringsAsFactors=F)
	colnames(enrichment_df) <- c("name", "nAnno", "nData", "nBG", "nOverlap", "nExpect", "fc", "zscore", "pvalue")

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
    
	res_df <- enrichment_df[, c("name", "nAnno", "nOverlap", "fc", "zscore", "pvalue", "adjp", "nData", "nBG")]
	
	invisible(res_df)
}
