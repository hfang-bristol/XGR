#' Function to conduct region-based enrichment analysis using genomic annotations via binomial test
#'
#' \code{xGRviaGenomicAnno} is supposed to conduct region-based enrichment analysis for the input genomic region data (genome build h19), using genomic annotations (eg active chromatin, transcription factor binding sites/motifs, conserved sites). Enrichment analysis is based on binomial test for estimating the significance of overlaps either at the base resolution, at the region resolution or at the hybrid resolution. Test background can be provided; by default, the annotatable will be used. 
#'
#' @param data.file an input data file, containing a list of genomic regions to test. If the input file is formatted as a 'data.frame' (specified by the parameter 'format.file' below), the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. If the format is indicated as "chr:start-end", instead of using the first 3 columns, only the first column will be used and processed. If the file also contains other columns, these additional columns will be ignored. Alternatively, the input file can be the content itself assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns.
#' @param annotation.file an input annotation file containing genomic annotations for genomic regions. If the input file is formatted as a 'data.frame', the first four columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), the ending chromosome position (3rd column), and the genomic annotations (eg transcription factors and histones; 4th column). If the format is indicated as 'bed', the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the format is indicated as "chr:start-end", the first two columns correspond to the chromosome:start-end (1st column) and the genomic annotations (eg transcription factors and histones; 2nd column). If the file also contains other columns, these additional columns will be ignored. Alternatively, the input file can be the content itself assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns.
#' @param background.file an input background file containing a list of genomic regions as the test background. The file format is the same as 'data.file'. By default, it is NULL meaning all annotatable bases (ig non-redundant bases covered by 'annotation.file') are used as background. However, if only one annotation (eg only a transcription factor) is provided in 'annotation.file', the background must be provided.
#' @param format.file the format for input files. It can be one of "data.frame", "chr:start-end", "bed" and "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param resolution the resolution of overlaps being tested. It can be one of "bases" at the base resolution (by default), "regions" at the region resolution, and "hybrid" at the base-region hybrid resolution (that is, data at the region resolution but annotation/background at the base resolution). If regions being analysed are SNPs themselves, then the results are the same even when choosing this parameter as either 'bases' or 'hybrid' or 'regions'
#' @param background.annotatable.only logical to indicate whether the background is further restricted to annotatable bases (covered by 'annotation.file'). In other words, if the background is provided, the background bases are those after being overlapped with annotatable bases. Notably, if only one annotation (eg only a transcription factor) is provided in 'annotation.file', it should be false
#' @param p.tail the tail used to calculate p-values. It can be either "two-tails" for the significance based on two-tails (ie both over- and under-overrepresentation)  or "one-tail" (by default) for the significance based on one tail (ie only over-representation)
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param GR.annotation the genomic regions of annotation data. By default, it is 'NA' to disable this option. Pre-built genomic annotation data are detailed in the section 'Note'. Alternatively, the user can also directly provide a customised GR object (or a list of GR objects)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return 
#' a data frame with following columns (below explanations are based on results at the 'hybrid' resolution):
#' \itemize{
#'  \item{\code{name}: the annotation name}
#'  \item{\code{nAnno}: the number of bases covered by that annotation. If the background is provided, they are also restricted by this}
#'  \item{\code{nOverlap}: the number of regions overlapped between input regions and annotation regions. If the background is provided, they are also restricted by this}
#'  \item{\code{fc}: fold change}
#'  \item{\code{zscore}: z-score}
#'  \item{\code{pvalue}: p-value}
#'  \item{\code{adjp}: adjusted p-value. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{or}: a vector containing odds ratio}
#'  \item{\code{CIl}: a vector containing lower bound confidence interval for the odds ratio}
#'  \item{\code{CIu}: a vector containing upper bound confidence interval for the odds ratio}
#'  \item{\code{expProb}: the probability of expecting bases overlapped between background regions and annotation regions}
#'  \item{\code{obsProb}: the probability of observing regions overlapped between input regions and annotation regions}
#' }
#' @note Pre-built genomic annotation data are detailed in \code{\link{xDefineGenomicAnno}}.
#' @export
#' @seealso \code{\link{xDefineGenomicAnno}}
#' @include xGRviaGenomicAnno.r
#' @examples
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' 
#' \dontrun{
#' # Enrichment analysis for GWAS SNPs from ImmunoBase
#' ## a) provide input data
#' data.file <- "http://galahad.well.ox.ac.uk/bigdata/ImmunoBase_GWAS.bed"
#' 
#' ## b) perform enrichment analysis using FANTOM expressed enhancers
#' ### one-tail p-value calculation (by default)
#' eTerm <- xGRviaGenomicAnno(data.file, format.file="bed", GR.annotation="FANTOM5_Enhancer_Cell", RData.location=RData.location)
#' ### alternatively: two-tails p-value calculation (useful to identify depletions)
#' eTerm_2 <- xGRviaGenomicAnno(data.file, format.file="bed", GR.annotation="FANTOM5_Enhancer_Cell", p.tail="two-tails", RData.location=RData.location)
#'
#' ## c) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' ## d) barplot of enriched terms
#' bp <- xEnrichBarplot(eTerm, top_num='auto', displayBy="fc")
#' bp
#'
#' ## e) forest plot of enriched terms
#' gp <- xEnrichForest(eTerm)
#' gp
#'
#' ## f) save enrichment results to the file called 'Regions_enrichments.txt'
#' output <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' utils::write.table(output, file="Regions_enrichments.txt", sep="\t", row.names=FALSE)
#'
#' ##########################################
#' ### Advanced use: customised GR.annotation
#' ##########################################
#' FANTOM5_CAT_Cell <- xRDataLoader('FANTOM5_CAT_Cell', RData.location=RData.location)
#' ls_gr_lncRNA <- lapply(FANTOM5_CAT_Cell, function(x) x[grep('lncRNA',x$Category)])
#' ls_gr_mRNA <- lapply(FANTOM5_CAT_Cell, function(x) x[grep('coding_mRNA',x$Category)])
#' GR.annotations <- c("ls_gr_lncRNA","ls_gr_mRNA","FANTOM5_CAT_Cell")
#' ls_df <- lapply(1:length(GR.annotations), function(i){
#' 	GR.annotation <- get(GR.annotations[i])
#' 	df <- xGRviaGenomicAnno(data.file=data.file, format.file="bed", GR.annotation=GR.annotation, RData.location=RData.location)
#' 	df$group  <- GR.annotations[i]
#' 	return(df)
#' })
#' df <- do.call(rbind, ls_df)
#' gp <- xEnrichHeatmap(df, fdr.cutoff=0.05, displayBy="zscore")
#' 
#' ##########################################
#' ### Advanced use: customised EpigenomeAtlas_15Segments
#' ##########################################
#' info <- xRDataLoader('EpigenomeAtlas_15Segments_info', RData.location=RData.location)
#' GR.annotations <- paste0('EpigenomeAtlas_15Segments_',names(info))
#' names(GR.annotations) <- info
#' ls_df <- lapply(1:length(GR.annotations), function(i){
#' 	GR.annotation <- GR.annotations[i]
#' 	message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation), as.character(Sys.time())), appendLF=T)
#' 	df <- xGRviaGenomicAnno(data.file=data.file, format.file="bed", GR.annotation=GR.annotation, RData.location=RData.location, verbose=F)
#' 	df$group <- names(GR.annotation)
#' 	return(df)
#' })
#' df <- do.call(rbind, ls_df)
#' gp <- xEnrichHeatmap(df, fdr.cutoff=0.05, displayBy="fdr", reorder="both")
#'
#' }

xGRviaGenomicAnno <- function(data.file, annotation.file=NULL, background.file=NULL, format.file=c("data.frame", "bed", "chr:start-end", "GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), resolution=c("bases","regions","hybrid"), background.annotatable.only=T, p.tail=c("one-tail","two-tails"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), GR.annotation=NA, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format.file <- match.arg(format.file)
    build.conversion <- match.arg(build.conversion)
    resolution <- match.arg(resolution)
    p.adjust.method <- match.arg(p.adjust.method)
    p.tail <- match.arg(p.tail)
    
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
    if(is.matrix(data.file) | is.data.frame(data.file) | is(data.file,"GRanges")){
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
    if(is.matrix(annotation.file) | is.data.frame(annotation.file) | is(annotation.file,"list")){
        annotation <- annotation.file
    }else if(!is.null(annotation.file)){
		annotation <- utils::read.delim(file=annotation.file, header=F, row.names=NULL, stringsAsFactors=F)
    }else{
    	if(verbose){
    		message("\t\tThe file 'annotation.file' is not provided, so built-in RData will be used instead!")
    	}
    	annotation <- NULL
    }
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\timport the background file (%s) ...", as.character(now)), appendLF=T)
	}
	## import background file
    if(is.matrix(background.file) | is.data.frame(background.file) | is(background.file,"GRanges")){
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
			seqnames=data[,1],
			ranges = IRanges::IRanges(start=as.numeric(data[,2]), end=as.numeric(data[,3])),
			strand = rep('*',nrow(data))
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
					seqnames=x[,1],
					ranges = IRanges::IRanges(start=as.numeric(x[,2]), end=as.numeric(x[,3])),
					strand = rep('*',nrow(x))
				)
			})
		}else{
			aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location, guid=guid)
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
				seqnames=background[,1],
				ranges = IRanges::IRanges(start=as.numeric(background[,2]), end=as.numeric(background[,3])),
				strand = rep('*',nrow(background))
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
			seqnames=data[,1],
			ranges = IRanges::IRanges(start=as.numeric(data[,2]), end=as.numeric(data[,3])),
			strand = rep('*',nrow(data))
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
					seqnames=x[,1],
					ranges = IRanges::IRanges(start=as.numeric(x[,2]), end=as.numeric(x[,3])),
					strand = rep('*',nrow(x))
				)
			})
		}else{
			aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location, guid=guid)
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
				seqnames=background[,1],
				ranges = IRanges::IRanges(start=as.numeric(background[,2]), end=as.numeric(background[,3])),
				strand = rep('*',nrow(data))
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
			seqnames=data[,1],
			ranges = IRanges::IRanges(start=as.numeric(data[,2])+1, end=as.numeric(data[,3])),
			strand = rep('*',nrow(data))
		)
		
		if(!is.null(annotation)){
			## construct annotation GR
			anno_ls <- split(x=annotation[,-4], f=annotation[,4])
			aGR <- base::lapply(anno_ls, function(x){
				## make sure positions are numeric
				ind <- suppressWarnings(which(!is.na(as.numeric(x[,2])) & !is.na(as.numeric(x[,3]))))
				x <- x[ind,]
				gr <- GenomicRanges::GRanges(
					seqnames=x[,1],
					ranges = IRanges::IRanges(start=as.numeric(x[,2])+1, end=as.numeric(x[,3])),
					strand = rep('*',nrow(x))
				)
			})
		}else{
			aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location, guid=guid)
			aGR <- lapply(aGRL, function(x) x)
		}
		
		if(!is.null(background)){
			## construct background GR
			## make sure positions are numeric
			ind <- suppressWarnings(which(!is.na(as.numeric(background[,2])) & !is.na(as.numeric(background[,3]))))
			background <- background[ind,]
			bGR <- GenomicRanges::GRanges(
				seqnames=background[,1],
				ranges = IRanges::IRanges(start=as.numeric(background[,2])+1, end=as.numeric(background[,3])),
				strand = rep('*',nrow(data))
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
			aGRL <- xDefineGenomicAnno(GR.annotation, verbose=verbose, RData.location=RData.location, guid=guid)
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
	mergeOverlaps <- function(qGR, sGR, maxgap=-1L, minoverlap=0L){
		hits <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query=qGR, subject=sGR, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
		qhits <- qGR[hits[,1]]
		shits <- sGR[hits[,2]]

		oGR <- IRanges::pintersect(qhits, shits, ignore.strand=T)
		IRanges::reduce(oGR)
	}
	
    ## Binomial test: sampling at random from the background with the constant probability of having annotated genes (with replacement)
    doBinomialTest <- function(X, K, M, N, p.tail){
        # X: num of success in sampling
        # K: num of sampling
        # M: num of success in background
        # N: num in background
        
        N <- max(N, M)
    	#########################
    	if(K==0 || M==0 || N==0){
    		p.value <- 1
    	}else{
    		if(p.tail=='one-tail'){
    			p.value <- stats::pbinom(X,K,M/N, lower.tail=F, log.p=F)
    		}else{
    			if(X>=K*M/N){
    				p.value <- stats::pbinom(X,K,M/N, lower.tail=F, log.p=F)
    			}else{
    				p.value <- stats::pbinom(X,K,M/N, lower.tail=T, log.p=F)
    			}
    		}
    	}
    	#########################
        #p.value <- ifelse(K==0 || M==0 || N==0, 1, stats::pbinom(X,K,M/N, lower.tail=F, log.p=F))

        return(p.value)
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
		dGR <- xLiftOver(data.file=dGR, format.file="GRanges", build.conversion=build.conversion, merged=F, verbose=verbose, RData.location=RData.location, guid=guid)
    	## aGR
    	if(!is.null(annotation.file)){
			if(verbose){
				message(sprintf("\tannotation genomic regions: lifted over via genome build conversion `%s`", build.conversion), appendLF=T)
			}
			aGR <- lapply(aGR, function(gr){
				xLiftOver(data.file=gr, format.file="GRanges", build.conversion=build.conversion, merged=F, verbose=verbose, RData.location=RData.location, guid=guid)
			})
		}
    	## bGR
		if(!is.null(bGR)){
			if(verbose){
				message(sprintf("\tbackground genomic regions: lifted over via genome build conversion `%s`", build.conversion), appendLF=T)
			}
			bGR <- xLiftOver(data.file=bGR, format.file="GRanges", build.conversion=build.conversion, merged=F, verbose=verbose, RData.location=RData.location, guid=guid)
		}
	}
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
	oGR_reduced <- base::lapply(aGR_reduced, function(gr){
		mergeOverlaps(qGR=gr, sGR=dGR_reduced, maxgap=-1L, minoverlap=0L)
	})
	
	#######################################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("Forth, perform enrichment analysis at '%s' resolution with '%s' p-values (%s) ...", resolution, p.tail, as.character(now)), appendLF=T)
	}
	
	
	if(resolution=="bases"){

		## prepare enrichment analysis
		data_nBases <- sum(as.numeric(IRanges::width(dGR_reduced)))
		overlap_nBases <- base::sapply(oGR_reduced, function(gr){
			sum(as.numeric(IRanges::width(gr)))
		})
		annotation_nBases <- base::sapply(aGR_reduced, function(gr){
			sum(as.numeric(IRanges::width(gr)))
		})
		background_nBases <- sum(as.numeric(IRanges::width(bGR_reduced)))

		if(verbose){
			now <- Sys.time()
			message(sprintf("\tthe number of bases: data (%d)", data_nBases), appendLF=T)
			message(sprintf("\tthe number of annotations: %d", length(annotation_nBases)), appendLF=T)
		}
		
	}else if(resolution=="regions"){
	
		## prepare enrichment analysis
		data_nBases <- length(dGR_reduced)
		overlap_nBases <- base::sapply(oGR_reduced, length)
		annotation_nBases <- base::sapply(aGR_reduced, length)
		background_nBases <- max(length(bGR_reduced), max(annotation_nBases))
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("\tthe number of regions: data (%d)", data_nBases), appendLF=T)
			message(sprintf("\tthe number of annotations: %d", length(annotation_nBases)), appendLF=T)
		}
		
	}else if(resolution=="hybrid"){
	
		## prepare enrichment analysis
		data_nBases <- length(dGR_reduced)
		overlap_nBases <- base::sapply(oGR_reduced, length)
		annotation_nBases <- base::sapply(aGR_reduced, function(gr){
			sum(as.numeric(IRanges::width(gr)))
		})
		background_nBases <- sum(as.numeric(IRanges::width(bGR_reduced)))
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("\tthe number of regions: data (%d)", data_nBases), appendLF=T)
			message(sprintf("\tthe number of annotations: %d", length(annotation_nBases)), appendLF=T)
		}
	
	}
	
	## perform enrichment analysis based on the binomial distribution
	res_ls <- base::lapply(1:length(overlap_nBases), function(i){
		X <- as.numeric(overlap_nBases[i])
		K <- data_nBases
		M <- as.numeric(annotation_nBases[i])
		N <- background_nBases
 
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
		
		if(is.na(z.score)){
			z.score <- 0
		}
		
		p.value <- doBinomialTest(X, K, M, N, p.tail)
		
		#### unable to calculate OR for 'hybrid' resolution
		if(all(resolution %in% c("bases","regions"))){
			## odds ratio calculated from Fisher's exact test
			## Prepare a two-dimensional contingency table: #success in sampling, #success in background, #failure in sampling, and #failure in left part
			cTab <- matrix(c(X, K-X, M-X, N-M-K+X), nrow=2, dimnames=list(c("anno", "notAnno"), c("group", "notGroup")))
			
			if(is(suppressWarnings(try(res <- stats::fisher.test(cTab), TRUE)),"try-error")){
				or <- CIl <- CIu <- NA
			}else{			
				res <- stats::fisher.test(cTab)
				or <- as.vector(res$estimate)
				CIl <- as.vector(res$conf.int)[1]
				CIu <- as.vector(res$conf.int)[2]		
			}	
		
		}else{
			or <- CIl <- CIu <- NA
		}
		
		
		## output
		c(X, K, M, N, X/K, M/N, (X/K)/(M/N), z.score, p.value, or, CIl, CIu)
	})
	res_df <- do.call(rbind, res_ls)
	enrichment_df <- data.frame(names(overlap_nBases), res_df, stringsAsFactors=F)
	colnames(enrichment_df) <- c("name", "nOverlap", "nData", "nAnno", "nBG", "obsProb", "expProb", "fc", "zscore", "pvalue", "or", "CIl", "CIu")

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
    
	res_df <- enrichment_df[, c("name", "nAnno", "nOverlap", "fc", "zscore", "pvalue", "adjp", "or", "CIl", "CIu", "expProb", "obsProb")]
	eTerm <- xEnrichViewer(res_df, top_num=nrow(res_df), sortBy='zscore')
	
	invisible(eTerm)
}
