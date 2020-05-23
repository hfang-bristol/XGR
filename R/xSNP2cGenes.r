#' Function to define HiC genes given a list of SNPs
#'
#' \code{xSNP2cGenes} is supposed to define HiC genes given a list of SNPs. The HiC weight is calcualted as Cumulative Distribution Function of HiC interaction scores. 
#'
#' @param data an input vector containing SNPs. SNPs should be provided as dbSNP ID (ie starting with rs) or in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for example, 'chr16:28525386'. Alternatively, it can be other formats/entities (see the next parameter 'entity')
#' @param entity the data entity. By default, it is "SNP". For general use, it can also be one of "chr:start-end", "data.frame", "bed" or "GRanges"
#' @param include.HiC genes linked to input SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, those genes linked to SNPs will be included according to Promoter Capture HiC (PCHiC) datasets. Pre-built HiC datasets are detailed in \code{\link{xDefineHIC}}
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param cdf.function a character specifying a Cumulative Distribution Function (cdf). It can be one of 'exponential' based on exponential cdf, 'empirical' for empirical cdf
#' @param plot logical to indicate whether the histogram plot (plus density or CDF plot) should be drawn. By default, it sets to false for no plotting
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return
#' a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: SNP-interacting genes caputured by HiC}
#'  \item{\code{SNP}: SNPs}
#'  \item{\code{Sig}: the interaction score (the higher stronger)}
#'  \item{\code{Weight}: the HiC weight}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xDefineHIC}}, \code{\link{xSymbol2GeneID}}
#' @include xSNP2cGenes.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' data(ImmunoBase)
#' data <- names(ImmunoBase$AS$variants)
#'
#' # b) define HiC genes
#' df_cGenes <- xSNP2cGenes(data, include.HiC="Monocytes", RData.location=RData.location)
#' }

xSNP2cGenes <- function(data, entity=c("SNP","chr:start-end","data.frame","bed","GRanges"), include.HiC=NA, GR.SNP=c("dbSNP_GWAS","dbSNP_Common"), cdf.function=c("empirical","exponential"), plot=FALSE, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
	entity <- match.arg(entity)
    cdf.function <- match.arg(cdf.function)
    
    ######################################################
    # Link to targets based on HiC
    ######################################################
	
	## all
	df_FTS <- xDefineHIC(data=NULL, include.HiC=include.HiC, verbose=verbose, RData.location=RData.location, guid=guid)
	
	##################
	if(is.null(data)){
		return(NULL)
	}
	##################
		
	## only data
	PCHiC <- xDefineHIC(data=data, entity=entity, include.HiC=include.HiC, GR.SNP=GR.SNP, verbose=verbose, RData.location=RData.location, guid=guid)
	df_data <- PCHiC$df
	
	if(!is.null(df_FTS)){
		
		## all
		uid <- paste(df_FTS[,1], df_FTS[,2], sep='_')
		df <- cbind(uid, df_FTS)
		res_list <- split(x=df$score, f=df$uid)
		raw_score <- unlist(lapply(res_list, max))
		
		## only data
		uid_data <- paste(df_data[,1], df_data[,2], sep='_')
		
		if(cdf.function == "exponential"){
			##  fit raw_score to the cumulative distribution function (CDF; depending on exponential empirical distributions)
			lambda <- MASS::fitdistr(raw_score, "exponential")$estimate
			
			## HiC weight for input SNPs
			## weights according to HiC
			wE <- stats::pexp(df_data$score, rate=lambda)
			
			#########
			if(nrow(df_data)==0){
				df_cGenes <- NULL
			}else{
				Gene <- sapply(1:nrow(df_data), function(i){
					if(df_data$SNP_end[i]=='bait/from'){
						df_data$to_genes[i]
					}else{
						df_data$from_genes[i]
					}
				})
				
				df_cGenes_ori <- data.frame(Gene=Gene, SNP=df_data$SNP, Sig=df_data$score, Weight=wE, row.names=NULL, stringsAsFactors=FALSE)
				
				ls_tmp <- strsplit(df_cGenes_ori$Gene, ';')
				res_ls <- lapply(1:length(ls_tmp), function(i){
					x <- ls_tmp[[i]]
					x <- x[x!='.']
					data.frame(Gene=x, df_cGenes_ori[rep(i,length(x)),-1], row.names=NULL, stringsAsFactors=FALSE)
				})
				df_cGenes <- do.call(rbind, res_ls)
				
			}
			#########
			
			if(plot){
				hist(raw_score, breaks=1000, freq=FALSE, col="grey", xlab="Score", main="")
				curve(stats::dexp(x=raw_score,rate=lambda), 5:max(raw_score), col=2, add=TRUE)
			}
			
			if(verbose){
				now <- Sys.time()
				message(sprintf("HiC weights are CDF of exponential empirical distributions (parameter lambda=%f)", lambda), appendLF=TRUE)
			}
			
		}else if(cdf.function == "empirical"){
			## Compute an empirical cumulative distribution function
			my.CDF <- stats::ecdf(raw_score)
			
			## HiC weight for input SNPs
			## weights according to HiC
			wE <- my.CDF(df_data$score)
			
			#########
			if(nrow(df)==0){
				df_cGenes <- NULL
			}else{
				Gene <- sapply(1:nrow(df_data), function(i){
					if(df_data$SNP_end[i]=='bait/from'){
						df_data$to_genes[i]
					}else{
						df_data$from_genes[i]
					}
				})
				
				df_cGenes_ori <- data.frame(Gene=Gene, SNP=df_data$SNP, Sig=df_data$score, Weight=wE, row.names=NULL, stringsAsFactors=FALSE)
				
				ls_tmp <- strsplit(df_cGenes_ori$Gene, ';')
				res_ls <- lapply(1:length(ls_tmp), function(i){
					x <- ls_tmp[[i]]
					x <- x[x!='.']
					data.frame(Gene=x, df_cGenes_ori[rep(i,length(x)),-1], row.names=NULL, stringsAsFactors=FALSE)
				})
				df_cGenes <- do.call(rbind, res_ls)
				
				df_cGenes <- df_cGenes[order(df_cGenes$Gene,-df_cGenes$Sig,df_cGenes$SNP,decreasing=FALSE),]
			}
			#########
			
			if(plot){
				#plot.ecdf
				plot(my.CDF, xlab="Score", ylab="Empirical CDF (HiC weights)", main="", xlim=c(0,50))
			}
			
			if(verbose){
				now <- Sys.time()
				message(sprintf("HiC weights are CDF of empirical distributions"), appendLF=TRUE)
			}
			
		}
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("%d cGenes are defined involving %d SNP", length(unique(df_cGenes$Gene)), length(unique(df_cGenes$SNP))), appendLF=TRUE)
		}
	
	}else{
		df_cGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No HiC genes are defined"), appendLF=TRUE)
		}
	}
	
	####################################
	# only keep those genes with GeneID
	####################################
	if(!is.null(df_cGenes)){
		ind <- xSymbol2GeneID(df_cGenes$Gene, details=FALSE, verbose=verbose, RData.location=RData.location, guid=guid)
		df_cGenes <- df_cGenes[!is.na(ind), ]
		
		if(nrow(df_cGenes)==0){
			df_cGenes <- NULL
		}else{
			if(1){
				#################################
				# remove HLA genes and histone genes
				#ind <- which(!grepl('^HLA-|^HIST', df_cGenes$Gene))
				# remove histone genes
				ind <- which(!grepl('^HIST', df_cGenes$Gene))
				df_cGenes <- df_cGenes[ind,]
				#################################
				if(nrow(df_cGenes)==0){
					df_cGenes <- NULL
				}
			}
		}
	}
	####################################
	
    invisible(df_cGenes)
}
