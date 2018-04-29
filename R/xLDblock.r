#' Function to obtain LD blocks
#'
#' \code{xLDblock} is supposed to obtain LD blocks for a list of Lead SNPs together with the significance level.
#'
#' @param data a named input vector containing the significance level for nodes (dbSNP). For this named vector, the element names are dbSNP (starting with rs or in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for example, 'chr16:28525386'), the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for dbSNP, 2nd column for the significance level. 
#' @param include.LD additional SNPs in LD with Lead SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, LD SNPs will be included based on one or more of 26 populations and 5 super populations from 1000 Genomics Project data (phase 3). The population can be one of 5 super populations ("AFR", "AMR", "EAS", "EUR", "SAS"). Explanations for population code can be found at \url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
#' @param LD.customised a user-input matrix or data frame with 3 compulsory columns: 1st column for Lead SNPs, 2nd column for LD SNPs, and 3rd for LD r2 value. The recommended columns are 'maf', 'distance' (to the nearest gene) and 'cadd'. It is designed to allow the user analysing their precalcuated LD info. This customisation (if provided) has the high priority over built-in LD SNPs
#' @param LD.r2 the LD r2 value. By default, it is 0.8, meaning that SNPs in LD (r2>=0.8) with input SNPs will be considered as LD SNPs. It can be any value from 0.1 to 1
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'LDblock_GR', that is, SNPs from dbSNP (version 150) restricted to GWAS SNPs and their LD SNPs (hg19). Note: you can also load your customised GR object directly
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "bLD", a list with following components:
#' \itemize{
#'  \item{\code{best}: a GR object. It has optional meta-columns 'maf', 'distance' (to the nearest gene) and 'cadd', and compulsory meta-columns 'pval', 'score' (-log10(pval)),  'upstream' (the lower boundary away from the best SNP, non-positive value), 'downstream' (the upper boundary away from the best SNP, non-negative value) and 'num' (the number of SNPs in the block)}
#'  \item{\code{block}: a GRL object, each element corresponding to a block for the best SNP with optional meta-columns 'maf', 'distance' (to the nearest gene) and 'cadd', and compulsory meta-columns 'pval', 'score' (-log10(pval)*R2, based on pval for its lead SNP), 'best' (the best SNP) and 'distance_to_best' (to the best SNP)}
#' }
#' @note None
#' @export
#' @seealso \code{\link{xLDblock}}
#' @include xLDblock.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#'
#' \dontrun{
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' data(ImmunoBase)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' data <- GenomicRanges::mcols(gr)[,c(1,3)]
#'
#' # b) get LD block (EUR population)
#' bLD <- xLDblock(data, include.LD="EUR", LD.r2=0.8, RData.location=RData.location)
#' 
#' # c1) manhattan plot of the best
#' best <- bLD$best
#' best$value <- best$score
#' gp <- xGRmanhattan(best, top=length(best))
#' gp
#' # c2) manhattan plot of all LD block
#' grl_block <- bLD$block
#' gr_block <- BiocGenerics::unlist(grl_block,use.names=F)
#' gr_block$value <- gr_block$score
#' top.label.query <- names(gr_block)[!is.na(gr_block$pval)]
#' #gr_block <- gr_block[as.character(GenomicRanges::seqnames(gr_block)) %in% c('chr1','chr2')]
#' gp <- xGRmanhattan(gr_block, top=length(gr_block), top.label.query=top.label.query)
#' 
#' # d) track plot of 1st LD block
#' gr_block <- bLD$block[[1]]
#' cnames <- c('score','maf','cadd')
#' ls_gr <- lapply(cnames, function(x) gr_block[,x])
#' names(ls_gr) <- cnames
#' ls_gr$score$Label <- names(gr_block)
#' ls_gr$score$Label[is.na(gr_block$pval)] <-''
#' GR.score.customised <- ls_gr
#' ## cse.query
#' df_block <- as.data.frame(gr_block)
#' chr <- unique(df_block$seqnames)
#' xlim <- range(df_block$start)
#' cse.query <- paste0(chr,':',xlim[1],'-',xlim[2])
#' #cse.query <- paste0(chr,':',xlim[1]-1e4,'-',xlim[2]+1e4)
#' ## xGRtrack
#' tks <- xGRtrack(cse.query=cse.query, GR.score="RecombinationRate", GR.score.customised=GR.score.customised, RData.location=RData.location)
#' tks
#' 
#' ###############
#' # Advanced use: get LD block (based on customised LD and SNP data)
#' ###############
#' LD.customised <- xRDataLoader('LDblock_EUR', RData.location=RData.location)
#' GR.SNP <- xRDataLoader('LDblock_GR', RData.location=RData.location)
#' bLD <- xLDblock(data, LD.customised=LD.customised, LD.r2=0.8, GR.SNP=GR.SNP, RData.location=RData.location)
#' }

xLDblock <- function(data, include.LD=c("AFR","AMR","EAS","EUR","SAS"), LD.customised=NULL, LD.r2=0.8, GR.SNP="LDblock_GR", verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    startT <- Sys.time()
    if(verbose){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }
	####################################################################################
	
    if(is.null(data)){
        stop("The input data must be not NULL.\n")
    }else{
    	
    	if(class(data)=='DataFrame'){
    		data <- as.data.frame(data)
    	}
    
		if (is.vector(data)){
			if(length(data)>1){
				# assume a vector
				if(is.null(names(data))){
					stop("The input data must have names with attached dbSNP ID.\n")
				}
			}else{
				# assume a file
				data <- utils::read.delim(file=data, header=F, row.names=NULL, stringsAsFactors=F)
			}
		}
		
		if (is.vector(data)){
			pval <- data
		}else if(is.matrix(data) | is.data.frame(data)){
			data <- as.data.frame(data, stringsAsFactors=F)
			colnames(data) <- c('snp','pval')
			data$pval <- as.numeric(data$pval)
			data <- data[!is.na(data$pval),]

			snp <- pval <- NULL
			data <- as.data.frame(data %>% dplyr::group_by(snp) %>% dplyr::summarise(min(pval)))
			pval <- data[,2]
			names(pval) <- data[,1]
		}
		
		# force those zeros to be miminum of non-zeros
		#tmp <- as.numeric(format(.Machine)['double.xmin'])
		tmp <- min(pval[pval!=0])
		pval[pval < tmp] <- tmp
	}
	
	## replace '_' with ':'
	tmp <- names(pval)
	tmp <- gsub("_", ":", tmp, perl=T)
	## replace 'imm:' with 'chr'
	names(pval) <- gsub("imm:", "chr", tmp, perl=T)
	
	pval <- sort(pval)
	Lead_Sig <- data.frame(SNP=names(pval), Sig=pval, row.names=NULL, stringsAsFactors=F)
	leads <- Lead_Sig[,1]
	sigs <- Lead_Sig[,2]

	if(verbose){
		message(sprintf("A total of %d Lead SNPs are input (%s)", length(leads), as.character(Sys.time())), appendLF=T)
	}

	###########################
	## include additional SNPs that are in LD with input SNPs
	if(LD.r2>=0.1 & LD.r2<=1){
		default.include.LD <- c("AFR","AMR","EAS","EUR","SAS")
		ind <- match(default.include.LD, include.LD)
		include.LD <- default.include.LD[!is.na(ind)]
	}else{
		include.LD <- NULL
	}
	
	LLR <- NULL
	if(length(include.LD) > 0 & is.null(LD.customised)){
		
		#############################
		## LLR contains all (maximum) LDs for leads in any populations
		#############################
		
		if(verbose){
			message(sprintf("LD SNPs are defined based on population '%s' (%s):", paste(include.LD, collapse=','), as.character(Sys.time())), appendLF=T)
		}

		res_list <- lapply(include.LD, function(x){
		
			if(verbose){
				message(sprintf("\tpopulation '%s' (%s) ...", x, as.character(Sys.time())), appendLF=T)
			}
		
			data_ld <- xRDataLoader(paste0("LDblock_", x), RData.location=RData.location, verbose=F)
			ind <- match(data_ld$Lead, leads)
			ind_lead <- which(!is.na(ind))
			
			if(length(ind_lead) >= 1){
				res <- as.data.frame(data_ld %>% dplyr::slice(ind_lead))
			}else{
				NULL
			}
		})
		## get data frame (Lead LD R2 distance)
		LLR <- do.call(rbind, res_list)
		
		Lead <- LD <- R2 <- NULL
		LLR <- as.data.frame(LLR %>% dplyr::group_by(Lead, LD) %>% dplyr::summarise(R2=max(R2)))
		
	}else if(!is.null(LD.customised)){
	
		if (is.vector(LD.customised)){
			# assume a file
			LLR <- utils::read.delim(file=LD.customised, header=F, row.names=NULL, stringsAsFactors=F)
		}else if(is.matrix(LD.customised) | is.data.frame(LD.customised)){
			LLR <- LD.customised
		}
		
		if(!is.null(LLR)){
			flag <- LLR[,3]>=LD.r2
			if(sum(flag)>0){
				colnames(LLR)[1:3] <- c("Lead", "LD", "R2")
				
				ind <- match(LLR$Lead, leads)
				ind_lead <- which(!is.na(ind))
			
				if(length(ind_lead) >= 1){
					LLR <- as.data.frame(LLR %>% dplyr::slice(ind_lead))
					
					if(verbose){
						message(sprintf("LD SNPs are defined based on customised data (%s):", as.character(Sys.time())), appendLF=T)
					}
					
				}else{
					LLR <- NULL
				}

			}else{
				LLR <- NULL
			}
		}
		
	}
	
	R2 <- Lead <- LD <- pval <- maf <- score <- index <- distance_to_best <- NULL
	
	gr_best <- NULL
	grl_block <- NULL
	if(!is.null(LLR)){
	
		if(verbose){
			message(sprintf("Define LD blocks (%s) ...", as.character(Sys.time())), appendLF=T)
		}

		if(verbose){
			message(sprintf("\tload positional information for SNPs (%s) ...", as.character(Sys.time())), appendLF=T)
		}
		if(class(GR.SNP) == "GRanges"){
			LDblock_GR <- GR.SNP
		}else{
			LDblock_GR <- xRDataLoader(GR.SNP[1], RData.location=RData.location, verbose=F)
			if(is.null(LDblock_GR)){
				GR.SNP <- "LDblock_GR"
				if(verbose){
					message(sprintf("\tinstead, %s will be used", GR.SNP), appendLF=T)
				}
				LDblock_GR <- xRDataLoader(GR.SNP, RData.location=RData.location, verbose=F)
			}
		}

		## LDblock_GR_lead
		ind <- match(names(LDblock_GR), unique(LLR$Lead))
		LDblock_GR_lead <- LDblock_GR[!is.na(ind)]
		## LDblock_GR_ld
		ind <- match(names(LDblock_GR), unique(LLR$LD))
		LDblock_GR_ld <- LDblock_GR[!is.na(ind)]
	
		## LDblock_GR_both
		ind <- match(names(LDblock_GR), unique(c(LLR$Lead,LLR$LD)))
		LDblock_GR_both <- LDblock_GR[!is.na(ind)]
	
		###########################	
		## construct igraph for leads
		df_nodes <- data.frame(name=unique(LLR$Lead), stringsAsFactors=F)
		ind <- match(LLR$LD, leads)
		ind <- which(!is.na(ind))
		df_edge <- as.data.frame(LLR %>% dplyr::slice(ind) %>% dplyr::filter(R2>=LD.r2))
		g <- igraph::graph.data.frame(d=df_edge, directed=TRUE, vertices=df_nodes)
		ls_block <- igraph::groups(igraph::components(g, mode=c("weak","strong")))
	
		if(verbose){
			message(sprintf("\tLD blocks: %d (%s)", length(ls_block), as.character(Sys.time())), appendLF=T)
		}
	
		ls_best_block <- lapply(ls_block, function(x){
			ind <- match(Lead_Sig$SNP, x)
			y <- Lead_Sig[!is.na(ind),]
			
			## gr for all input SNPs
			ind <- match(names(LDblock_GR_lead), y$SNP)
			gr <- LDblock_GR_lead[!is.na(ind)]
			gr$pval <- y$Sig[ind[!is.na(ind)]]
			
			## best_gr
			best_ind <- which(gr$pval==min(gr$pval))[1]
			best_gr <- gr[best_ind]
			best_gr$score <- -log10(best_gr$pval)
			
			## df_block
			ind <- match(LLR$Lead, x)
			ind <- which(!is.na(ind))
			df_block <- as.data.frame(LLR %>% dplyr::slice(ind) %>% dplyr::filter(R2>=LD.r2) %>% dplyr::arrange(Lead,-R2))[,c('Lead','LD','R2')]
			df_tmp <- data.frame(Lead=x, LD=x, R2=rep(1,length(x)), stringsAsFactors=1)
			if(nrow(df_block) == 0){
				df_block <- df_tmp
			}else{
				df_block <- rbind(df_tmp, df_block)
			}
			## append 'pval' and 'score'
			ind <- match(df_block$Lead, y$SNP)
			df_block$pval <- y$Sig[ind]
			df_block$score <- -log10(df_block$pval) * df_block$R2
			
			if(nrow(df_block) > 0){
				
				## maximum score
				df_block_LD <- as.data.frame(df_block %>% dplyr::group_by(LD) %>% dplyr::summarise(score=max(score)))
				
				## ld_gr
				ind <- match(names(LDblock_GR_both), df_block_LD$LD)
				ld_gr <- LDblock_GR_both[!is.na(ind)]
				
				####
				# add 'pval'
				ind <- match(names(ld_gr), names(gr))
				ld_gr$pval <- gr$pval[ind]
				####

				####
				# add 'score'
				ind <- match(names(ld_gr), df_block_LD$LD)
				ld_gr$score <- df_block_LD$score[ind[!is.na(ind)]]
				####
				
				####
				# add 'best'
				ld_gr$best <- names(best_gr)
				####
				
				####
				# add 'distance_to_best'
				ld_gr$distance_to_best <- GenomicRanges::start(best_gr) - GenomicRanges::start(ld_gr)
				####
				
				## calculate upstream and downstream
				updown_dist <- range(ld_gr$distance_to_best <- GenomicRanges::start(best_gr) - GenomicRanges::start(ld_gr))
				best_gr$upstream <- min(0, updown_dist[1])
				best_gr$downstream <- max(0, updown_dist[2])
			
			}else{
				ld_gr <- best_gr
				ld_gr$best <- names(best_gr)
				ld_gr$distance_to_best <- 0
				
				best_gr$upstream <- 0
				best_gr$downstream <- 0
			}
			
			## sort: pval, score, distance_to_best
			df_sort <- as.data.frame(ld_gr)
			df_sort$index <- 1:nrow(df_sort)
			df_sort <- as.data.frame(df_sort %>% dplyr::arrange(pval,-score,distance_to_best))
			ld_gr <- ld_gr[df_sort$index]
			
			res <- list(best=best_gr, block=ld_gr)
			
			return(res)
		})
		
		# gr_best
		ls_best_gr <- lapply(ls_best_block, function(x) x$best)
		grl <- GenomicRanges::GRangesList(ls_best_gr)
		gr_best <- BiocGenerics::unlist(grl,use.names=F)
		
		# lgr_block
		lgr_block <- lapply(ls_best_block, function(x) x$block)
		names(lgr_block) <- names(gr_best)
		
		## add 'num'
		gr_best$num <- sapply(lgr_block, length)
		
		########################
		# order: by pval
		ind <- order(gr_best$pval)
		gr_best <- gr_best[ind]
		lgr_block <- lgr_block[ind]
		grl_block <- GenomicRanges::GRangesList(lgr_block)
		########################

	}
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(verbose){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
	
    bLD <- list(best = gr_best,
    			block = grl_block
                 )
    class(bLD) <- "bLD"
	
    return(bLD)
}
