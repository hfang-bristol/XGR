#' Function to conduct region-based enrichment analysis using genomic annotations
#'
#' \code{xEnricherRegions} is supposed to conduct region-based enrichment analysis for the input genomic region dataan, using genomic annotations (eg active chromatin, transcription factor binding sites/motifs, conserved sites). Enrichment analysis is based on either Fisher's exact test or Hypergeometric test for estimating the significance of overlaps at the single base-pair resolution. Test background can be provided; by default, the annotatable bases will be used. 
#'
#' @param data.file an input data file, containing a list of genomic regions to test. If the input file is formatted as a 'data.frame' (specified by the parameter 'format.file' below), the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. If the format is indicated as "chr:start-end", instead of using the first 3 columns, only the first column will be used and processed. If the file also contains other columns, these additional columns will be ignored. Alternatively, the input file can be the content itself assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param annotation.file an input annotation file containing genomic annotations for genomic regions. If the input file is formatted as a 'data.frame', the first four columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), the ending chromosome position (3rd column), and the genomic annotations (eg transcription factors and histones; 4th column). If the format is indicated as 'bed', the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the format is indicated as "chr:start-end", the first two columns correspond to the chromosome:start-end (1st column) and the genomic annotations (eg transcription factors and histones; 2nd column). If the file also contains other columns, these additional columns will be ignored. Alternatively, the input file can be the content itself assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param background.file an input background file containing a list of genomic regions as the test background. The file format is the same as 'data.file'. By default, it is NULL meaning all annotatable bases (ig non-redundant bases covered by 'annotation.file') are used as background. However, if only one annotation (eg only a transcription factor) is provided in 'annotation.file', the background must be provided. 
#' @param format.file the format for input files. It can be one of "data.frame", "chr:start-end", "bed"
#' @param background.annotatable.only logical to indicate whether the background is further restricted to annotatable bases (covered by 'annotation.file'). In other words, if the background is provided, the background bases are those after being overlapped with annotatable bases. Notably, if only one annotation (eg only a transcription factor) is provided in 'annotation.file', it should be false. 
#' @param test the statistic test used. It can be "fisher" for using fisher's exact test, "hypergeo" for using hypergeometric test, or "binomial" for using binomial test. Fisher's exact test is to test the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement). Hypergeometric test is to sample at random (without replacement) from the background containing annotated and non-annotated genes, and thus compare sampling to background. Unlike hypergeometric test, binomial test is to sample at random (with replacement) from the background with the constant probability. In terms of the ease of finding the significance, they are in order: hypergeometric test > binomial test > fisher's exact test. In other words, in terms of the calculated p-value, hypergeometric test < binomial test < fisher's exact test
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param GR.annotation the genomic regions of annotation data. By default, it is 'NA' to disable this option. Predefined annotation data include: ENCODE Transcription Factor ChIP-seq Uniform Peaks data ("Uniform_TFBS"), ENCODE DNaseI Hypersensitivity Uniform Peaks data ("Uniform_DNaseI_HS"), ENCODE Histone Modification ChIP-seq data from different sources ("Broad_Histone","SYDH_Histone","UW_Histone"), FANTOM5 expressed enhancer altas ("FANTOM5_Enhancer_Cell","FANTOM5_Enhancer_Tissue","FANTOM5_Enhancer_Extensive","FANTOM5_Enhancer"), ENCODE combined Genome Segmentation data (ChromHMM and Segway) in 6 cell lines ("Segment_Combined_Gm12878","Segment_Combined_H1hesc","Segment_Combined_Helas3","Segment_Combined_Hepg2","Segment_Combined_Huvec","Segment_Combined_K562"), the human/mouse/rat conserved TFBS ("TFBS_Conserved"), TargetScan miRNA regulatory sites ("TS_miRNA"), TCGA exome mutation data across 12 tumor types ("TCGA"). Beyond these predefined data, the user can specify the customised input. To do so, first save your RData file (a list of GR objects, each is an GR object correponding to an annotation) into your local computer. Then, tell "GR.annotation" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a data frame with 8 columns:
#' \itemize{
#'  \item{\code{name}: the annotation name}
#'  \item{\code{nAnno}: the number of nucleotides/bases covered by that annotation. If the background is provided, they are also restricted by this}
#'  \item{\code{nOverlap}: the number of nucleotides/bases overlapped between input regions and annotation regions. If the background is provided, they are also restricted by this}
#'  \item{\code{zscore}: z-score}
#'  \item{\code{pvalue}: p-value}
#'  \item{\code{adjp}: adjusted p-value. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{expProb}: the probability of expecting nucleotides/bases overlapped between background regions and annotation regions}
#'  \item{\code{obsProb}: the probability of observing nucleotides/bases overlapped between input regions and annotation regions}
#' }
#' @note None
#' @export
#' @seealso \code{\link{xEnrichViewer}}
#' @include xEnricherRegions.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' 
#' # Enrichment analysis for GWAS SNPs from ImmunoBase
#' # a) provide input data
#' data.file <- "http://galahad.well.ox.ac.uk/bigdata/ImmunoBase_GWAS.bed"
#' 
#' # b) perform enrichment analysis using FANTOM expressed enhancers
#' eTerm <- xEnricherRegions(data.file=data.file, format.file="bed", GR.annotation="FANTOM5_Enhancer_Cell")
#'
#' # c) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' # d) save enrichment results to the file called 'Regions_enrichments.txt'
#' output <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' utils::write.table(output, file="Regions_enrichments.txt", sep="\t", row.names=FALSE)
#' }

xEnricherRegions <- function(data.file, annotation.file=NULL, background.file=NULL, format.file=c("data.frame", "bed", "chr:start-end"), background.annotatable.only=T, test=c("hypergeo","fisher","binomial"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), GR.annotation=c(NA,"Uniform_TFBS","Uniform_DNaseI_HS", "Broad_Histone","SYDH_Histone","UW_Histone","FANTOM5_Enhancer_Cell","FANTOM5_Enhancer_Tissue","FANTOM5_Enhancer_Extensive","FANTOM5_Enhancer","Segment_Combined_Gm12878","Segment_Combined_H1hesc","Segment_Combined_Helas3","Segment_Combined_Hepg2","Segment_Combined_Huvec","Segment_Combined_K562","TFBS_Conserved","TS_miRNA","TCGA"), verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/Portal")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format.file <- match.arg(format.file)
    test <- match.arg(test)
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
    if(is.matrix(data.file) | is.data.frame(data.file)){
        data <- data.file
    }else if(!is.null(data.file) & !is.na(data.file)){
		data <- utils::read.delim(file=data.file, header=F, row.names=NULL, stringsAsFactors=F)
    }else{
    	stop("The file 'data.file' must be provided!\n")
    }
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\timport the annotation file (%s) ...", as.character(now)), appendLF=T)
	}
    ## import annotation file
    if(is.matrix(annotation.file) | is.data.frame(annotation.file)){
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
    if(is.matrix(background.file) | is.data.frame(background.file)){
        background <- background.file
    }else if(!is.null(background.file)){
		background <- utils::read.delim(file=background.file, header=F, row.names=NULL, stringsAsFactors=F)
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
			aGR <- lapply(anno_ls, function(x){
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
			if(is.na(GR.annotation)){
				stop("Please specify annotation RData!\n")
			}else{
				if(length(GR.annotation)>1){
					message("\tONLY the first specified annotation RData will be used!\n")
					GR.annotation <- GR.annotation[1]
				}
				aGR <- xRDataLoader(RData.customised=GR.annotation, verbose=verbose, RData.location=RData.location)
				if(is.null(aGR)){
					stop("Your specified annotation RData does not exist!\n")
				}
			}
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
			ranges = IRanges::IRanges(start=as.numeric(data[,2]+1), end=as.numeric(data[,3])),
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
			aGR <- lapply(anno_ls, function(x){
				## make sure positions are numeric
				ind <- suppressWarnings(which(!is.na(as.numeric(x[,2])) & !is.na(as.numeric(x[,3]))))
				x <- x[ind,]
				gr <- GenomicRanges::GRanges(
					seqnames=S4Vectors::Rle(x[,1]),
					ranges = IRanges::IRanges(start=as.numeric(x[,2]+1), end=as.numeric(x[,3])),
					strand = S4Vectors::Rle(rep('*',nrow(x)))
				)
			})
		}else{
			if(is.na(GR.annotation)){
				stop("Please specify annotation RData!\n")
			}else{
				if(length(GR.annotation)>1){
					message("\tONLY the first specified annotation RData will be used!\n")
					GR.annotation <- GR.annotation[1]
				}
				aGR <- xRDataLoader(RData.customised=GR.annotation, verbose=verbose, RData.location=RData.location)
				if(is.null(aGR)){
					stop("Your specified annotation RData does not exist!\n")
				}
			}
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
			aGR <- lapply(anno_ls, function(x){
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
			if(is.na(GR.annotation)){
				stop("Please specify annotation RData!\n")
			}else{
				if(length(GR.annotation)>1){
					message("\tONLY the first specified annotation RData will be used!\n")
					GR.annotation <- GR.annotation[1]
				}
				aGR <- xRDataLoader(RData.customised=GR.annotation, verbose=verbose, RData.location=RData.location)
				if(is.null(aGR)){
					stop("Your specified annotation RData does not exist!\n")
				}
			}
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
		
	}
	
	#####################################
	## A function to return an GR object storing overlapped regions (ie only overlapped regions!)
	mergeOverlaps <- function(qGR, sGR, maxgap=0L, minoverlap=1L){
		hits <- GenomicRanges::findOverlaps(query=qGR, subject=sGR, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)
		qhits <- qGR[GenomicRanges::queryHits(hits)]
		shits <- sGR[GenomicRanges::subjectHits(hits)]

		oGR <- GenomicRanges::pintersect(qhits, shits)
		GenomicRanges::reduce(oGR)
	}
	
    ## Binomial test: sampling at random from the background with the constant probability of having annotated genes (with replacement)
    doBinomialTest <- function(X, K, M, N){
        # X: num of success in sampling
        # K: num of sampling
        # M: num of success in background
        # N: num in background
        p.value <- ifelse(K==0 || M==0 || N==0, 1, stats::pbinom(X,K,M/N, lower.tail=F, log.p=F))
        return(p.value)
    }
    
    ## Fisher's exact test: testing the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not); thus compare sampling to the left part of background (after sampling without replacement)
    doFisherTest <- function(X, K, M, N){
        # X: num of success in sampling
        # K: num of sampling
        # M: num of success in background
        # N: num in background
        ## Prepare a two-dimensional contingency table: #success in sampling, #success in background, #failure in sampling, and #failure in left part
        cTab <- matrix(c(X, K-X, M-X, N-M-K+X), nrow=2, dimnames=list(c("anno", "notAnno"), c("group", "notGroup")))
        p.value <- ifelse(all(cTab==0), 1, stats::fisher.test(cTab, alternative="greater")$p.value)
        return(p.value)
    }
    
    ## Hypergeometric test: sampling at random from the background containing annotated and non-annotated genes (without replacement); thus compare sampling to background
    doHypergeoTest <- function(X, K, M, N){
        # X: num of success in sampling
        # K: num of sampling
        # M: num of success in background
        # N: num in background
        
        x <- X
        m <- M
        n <- N-M # num of failure in background
        k <- K
        p.value <- ifelse(m==0 || k==0, 1, stats::phyper(x,m,n,k, lower.tail=F, log.p=F))
        return(p.value)
    }
	#####################################
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Third, define the background (%s) ...", as.character(now)), appendLF=T)
	}
    
    ## get reduced ranges (ie non-overlapping regions)
    ### data GR
    dGR_reduced <- GenomicRanges::reduce(dGR)
    ### annotation GR
	aGR_reduced <- lapply(aGR, function(x){
		GenomicRanges::reduce(x)
	})
	### define background GR
	if(is.null(bGR)){
		if(verbose){
			now <- Sys.time()
			message(sprintf("\tall annotatable regions (by default) are used as the background (%s) ...", as.character(now)), appendLF=T)
		}
	
		aGRL <- GenomicRanges::GRangesList(aGR_reduced)
		bGR_reduced <- GenomicRanges::reduce(unlist(aGRL))
	}else{
		bGR_reduced <- GenomicRanges::reduce(bGR)
	
		## update annotation GR after considering background
		aGR_reduced <- lapply(aGR_reduced, function(gr){
			mergeOverlaps(qGR=gr, sGR=bGR_reduced, maxgap=0L, minoverlap=1L)
		})
	
		## restrict to the annotatable only?
		if(background.annotatable.only){
			if(verbose){
				now <- Sys.time()
				message(sprintf("\tthe given background regions but restricted to the annotatable are used as the background (%s) ...", as.character(now)), appendLF=T)
			}
		
			## update background GR
			aGRL <- GenomicRanges::GRangesList(aGR_reduced)
			bGR_reduced <- GenomicRanges::reduce(unlist(aGRL))
		}else{
			if(verbose){
				now <- Sys.time()
				message(sprintf("\tthe given background regions are used as the background (%s) ...", as.character(now)), appendLF=T)
			}
		}
	}
	
	## update data GR after considering background
	dGR_reduced <- mergeOverlaps(qGR=dGR_reduced, sGR=bGR_reduced, maxgap=0L, minoverlap=1L)

	## find overlap GR between annotation GR and data GR
	oGR_reduced <- lapply(aGR_reduced, function(gr){
		mergeOverlaps(qGR=gr, sGR=dGR_reduced, maxgap=0L, minoverlap=1L)
	})
	
	#######################################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("Forth, perform region-based enrichment analysis at the single base resolution (%s) ...", as.character(now)), appendLF=T)
	}
	
	## prepare enrichment analysis
	data_nBases <- sum(GenomicRanges::width(dGR_reduced))
	annotation_nBases <- sapply(aGR_reduced, function(gr){
		sum(GenomicRanges::width(gr))
	})
	background_nBases <- sum(GenomicRanges::width(bGR_reduced))
	overlap_nBases <- sapply(oGR_reduced, function(gr){
		sum(GenomicRanges::width(gr))
	})

	if(verbose){
		now <- Sys.time()
		message(sprintf("\tthe number of nucleotides/bases: data (%d), background (%d)", data_nBases, background_nBases), appendLF=T)
		message(sprintf("\tthe number of annotations: %d", length(annotation_nBases)), appendLF=T)
		message(sprintf("\tusing '%s' test", test), appendLF=T)
	}

	## perform enrichment analysis based on the binomial distribution
	res_ls <- lapply(1:length(overlap_nBases), function(i){
		X <- as.numeric(overlap_nBases[i])
		K <- data_nBases
		M <- as.numeric(annotation_nBases[i])
		N <- background_nBases

        p.value <- switch(test,
        	fisher =  doFisherTest(X, K, M, N),
            hypergeo = doHypergeoTest(X, K, M, N),
            binomial = doBinomialTest(X, K, M, N)
        )
 
        ## Z-score based on theoretical calculation
        x.exp <- K*M/N
        var.exp <- K*M/N*(N-M)/N*(N-K)/(N-1)
        if(is.na(var.exp)){
        	z.score <- 0
        }else{
			if(var.exp!=0){
				suppressWarnings(z.score <- (X-x.exp)/sqrt(var.exp))
			}else{
				z.score <- 0
			}
		}
		
		## output
		c(X, K, M, N, X/K, M/N, (X/K)/(M/N), z.score, p.value)
	})
	res_df <- do.call(rbind, res_ls)
	enrichment_df <- data.frame(names(overlap_nBases), res_df, stringsAsFactors=F)
	colnames(enrichment_df) <- c("name", "nOverlap", "nData", "nAnno", "nBG", "obsProb", "expProb", "fc", "zscore", "pvalue")

	## Adjust P-values for multiple comparisons
	p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel")[1]
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
	pvals  <- sapply(pvals, function(x){
		if(x < 0.1 & x!=0){
			as.numeric(format(x,scientific=T))
		}else{
			x
		}
	})
	
	adjpvals <- sapply(adjpvals, function(x){
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
    
	res_df <- enrichment_df[, c("name", "nAnno", "nOverlap", "zscore", "pvalue", "adjp", "expProb", "obsProb")]
	
	invisible(res_df)
}
