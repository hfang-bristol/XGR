#' Function to priorise genes given a list of seed SNPs together with the significance level (e.g. GWAS reported p-values)
#'
#' \code{xPrioritiserSNPs} is supposed to priorise genes given a list of seed SNPs together with the significance level. To priorise genes, it first defines seed genes and their weights that take into account the distance to and the significance of seed SNPs. With seed genes and weights, it then uses Random Walk with Restart (RWR) to calculate the affinity score of all nodes in the input graph to the seed genes. The priority score is the affinity score. Parallel computing is also supported for Linux or Mac operating systems. It returns an object of class "pNode".
#'
#' @param data a named input vector containing the sinificance level for nodes (dbSNP). For this named vector, the element names are dbSNP, the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for dbSNP, 2nd column for the significance level
#' @param include.LD additional SNPs in LD with Lead SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, LD SNPs will be included based on one or more of 26 populations and 5 super populations from 1000 Genomics Project data (phase 3). The population can be one of 5 super populations ("AFR", "AMR", "EAS", "EUR", "SAS"), or one of 26 populations ("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"). Explanations for population code can be found at \url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
#' @param LD.r2 the LD r2 value. By default, it is 0.8, meaning that SNPs in LD (r2>=0.8) with input SNPs will be considered as LD SNPs. It can be any value from 0.8 to 1
#' @param include.eQTL genes modulated by eQTL (also Lead SNPs or in LD with Lead SNPs) are also included. By default, it is 'NA' to disable this option. Otherwise, those genes modulated either by cis-eQTLs ('JKscience_TS2B') or trans-eQTLs ('JKscience_TS3A') will be inlcuded according to this work by Fairfax et al. Science 2014, 343(6175):1246949
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathways Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), and "STRING_medium" for interactions with medium confidence (confidence scores>=400). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addtion to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of SNPs into the significance-component weights. If given, those SNPs below this are considered significant and thus weighted positively. Instead, those above this are considered insigificant and thus receive no weight
#' @param distance.max the maximum distance between genes and SNPs. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby SNPs per gene
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for Random Walk with Restart (RWR). The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "pNode", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a matrix of nNode X 4 containing node priority information, where nNode is the number of nodes in the input graph, and the 4 columns are "name" (node names), "seed" (1 for seeds, 0 for non-seeds), "weight" (weight values for seeds),  "priority" (the priority scores that are rescaled to the range [0,1]), "rank" (ranks of the priority scores)}
#'  \item{\code{g}: an input "igraph" object}
#'  \item{\code{SNP}: a data frame of nSNP X 3 containing input SNPs and/or LD SNPs info, where nSNP is the number of input SNPs and/or LD SNPs, and the 3 columns are "SNP" (dbSNP), "Pval" (the SNP p-value), "Score" (the SNP score)}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The search procedure is heuristic to find the subgraph with the maximum score:
#' \itemize{
#' \item{i) transform the significance level of SNPs into the significance-component weights (noded as 'wS'). If the intolerable significance threshold is given, those SNPs below this are considered significant and thus weighted positively. Instead, those above this are considered insigificant and thus receive no weight.}
#' \item{ii) find genes located away from seed SNPs within the certain range (by default 500kb) and, for nearby SNPs per gene, calculate the distance-component weights (noded as 'wD').}
#' \item{iii) define seed genes as those found in ii) and their weights as the maximum of 'wS * wD'.}
#' \item{iv) \code{\link{xPrioritiserGenes}} used to prioritise genes using an input graph and a list of seed genes weighted from iii). The priority score is the affinity score estimated by Random Walk with Restart (RWR), measured as the affinity of all nodes in the graph to the seeds.}
#' }
#' @export
#' @import MASS
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xPrioritiser}}, \code{\link{xPrioritiserGenes}}, \code{\link{xPrioritiserPathways}}
#' @include xPrioritiserSNPs.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' library(dnet)
#' library(ggbio)
#'
#' # a) provide the SNPs with the significance info
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' AS <- read.delim(file.path(path.package("XGR"),"AS.txt"), stringsAsFactors=FALSE)
#'
#' # b) perform priority analysis
#' pNode <- xPrioritiserSNPs(data=AS, network="PCommonsUN_medium",restart=0.7)
#'
#' # c) save to the file called 'SNPs_priority.txt'
#' write.table(pNode$priority, file="SNPs_priority.txt", sep="\t", row.names=FALSE)
#' 
#' # d) manhattan plot
#' mp <- xPrioritiserManhattan(pNode, highlight.top=10)
#' #pdf(file="Gene_manhattan.pdf", height=6, width=12, compress=TRUE)
#' print(mp)
#' #dev.off()
#' }

xPrioritiserSNPs <- function(data, include.LD=NA, LD.r2=0.8, include.eQTL=NA, network=c("STRING_highest","STRING_high","STRING_medium","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD"), network.customised=NULL, significance.threshold=5e-5, distance.max=200000, normalise=c("laplacian","row","column","none"), restart=0.75, normalise.affinity.matrix=c("none","quantile"), parallel=TRUE, multicores=NULL, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    if(is.null(data)){
        stop("The input data must be not NULL.\n")
    }else{
		if (is.vector(data)){
			if(length(data)>1){
				# assume a vector
				if(is.null(names(data))){
					stop("The input data must have names with attached gene symbols.\n")
				}
			}else{
				# assume a file
				data <- utils::read.delim(file=data, header=F, row.names=NULL, stringsAsFactors=F)
			}
		}
		
		if (is.vector(data)){
			pval <- data
		}else if(is.matrix(data) | is.data.frame(data)){
			data <- as.matrix(data)
			data_list <- split(x=data[,2], f=as.character(data[,1]))
			res_list <- lapply(data_list, function(x){
				x <- as.numeric(x)
				x <- x[!is.na(x)]
				if(length(x)>0){
					min(x)
				}else{
					NULL
				}
			})
			pval <- unlist(res_list)
		}
		
		# force those zeros to be miminum of non-zeros
		#tmp <- as.numeric(format(.Machine)['double.xmin'])
		tmp <- min(pval[pval!=0])
		pval[pval < tmp] <- tmp
	}
	
	Lead_Sig <- data.frame(SNP=names(pval), Sig=pval, row.names=NULL, stringsAsFactors=F)
	leads <- Lead_Sig[,1]
	sigs <- Lead_Sig[,2]

	###########################
	## include additional SNPs that are in LD with input SNPs
	if(LD.r2>=0.8 & LD.r2<=1){
		default.include.LD <- c("ACB","AFR","AMR","ASW","BEB","CDX","CEU","CHB","CHS","CLM","EAS","ESN","EUR","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","SAS","STU","TSI","YRI")
		ind <- match(default.include.LD, include.LD)
		include.LD <- default.include.LD[!is.na(ind)]
	}
	
	if(length(include.LD) > 0){
		GWAS_LD <- xRDataLoader(RData.customised='GWAS_LD', RData.location=RData.location, verbose=verbose)
		res_list <- lapply(include.LD, function(x){
			data_ld <- ''
			eval(parse(text=paste("data_ld <- GWAS_LD$", x, sep="")))
			ind <- match(rownames(data_ld), leads)
			ind_lead <- which(!is.na(ind))
			ind_ld <- which(Matrix::colSums(data_ld[ind_lead,]>=LD.r2)>0)
		
			sLL <- data_ld[ind_lead, ind_ld]
			summ <- summary(sLL)
			res <- data.frame(Lead=rownames(sLL)[summ$i], LD=colnames(sLL)[summ$j], R2=summ$x, stringsAsFactors=F)
		})
		## get data frame (Lead LD R2)
		LLR <- do.call(rbind, res_list)
		
		###########################
		## also based on ImmunoBase
		if(1){
			ImmunoBase_LD <- xRDataLoader(RData.customised='ImmunoBase_LD', RData.location=RData.location, verbose=verbose)
			res_list <- lapply(include.LD, function(x){
				data_ld <- ''
				eval(parse(text=paste("data_ld <- ImmunoBase_LD$", x, sep="")))
				ind <- match(rownames(data_ld), leads)
				ind_lead <- which(!is.na(ind))
				ind_ld <- which(Matrix::colSums(data_ld[ind_lead,]>=LD.r2)>0)
		
				sLL <- data_ld[ind_lead, ind_ld]
				summ <- summary(sLL)
				res <- data.frame(Lead=rownames(sLL)[summ$i], LD=colnames(sLL)[summ$j], R2=summ$x, stringsAsFactors=F)
			})
			## get data frame (Lead LD R2)
			LLR_tmp <- do.call(rbind, res_list)
			LLR <- rbind(LLR, LLR_tmp)
		}
		###########################
				
		## get data frame (LD Sig)
		ld_list <- split(x=LLR[,-2], f=LLR[,2])
		res_list <- lapply(ld_list, function(x){
			ind <- match(x$Lead, leads)
			## power transformation of p-values X by R2, then keep the min
			min(sigs[ind] ^ x$R2)
		})
		vec <- unlist(res_list)
		LD_Sig <- data.frame(SNP=names(vec), Sig=vec, row.names=NULL, stringsAsFactors=F)
	
		## merge Lead and LD
		df <- rbind(Lead_Sig, as.matrix(LD_Sig))
		res_list <- split(x=df$Sig, f=df$SNP)
		res <- lapply(res_list, function(x){
			min(x)
		})
		vec <- unlist(res)
		SNP_Sig <- data.frame(SNP=names(vec), FDR=vec, row.names=NULL, stringsAsFactors=F)
	
	}else{
		SNP_Sig <- Lead_Sig
	}
	###########################
	
	pval <- as.numeric(SNP_Sig[,2])
	names(pval) <- SNP_Sig[,1]
	
	# transformed into scores according to log-likelihood ratio between the true positives and the false positivies
	scores <- log10(2) * dFDRscore(pval, fdr.threshold=significance.threshold, scatter=F)
	scores[scores<0] <- 0
	seeds.snps <- scores
    
    #########
    # get a data frame: SNPs pvalue score
    df_snp <- data.frame(SNP=names(pval), Pval=pval, Score=seeds.snps, row.names=NULL, stringsAsFactors=F)
    #########
    
    ######################################################
    # Link to targets based on genomic distance
    ######################################################
    
  	## load positional information
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for SNPs (%s) ...", as.character(now)), appendLF=T)
	}
  	pos_SNP <- xRDataLoader(RData.customised="RegulomeDB_SNPs", RData.location=RData.location, verbose=verbose)
  	ind <- match(names(seeds.snps), names(pos_SNP))
  	ind <- ind[!is.na(ind)]
  	if(length(ind)){
  		gr_SNP <- pos_SNP[ind,]
  		
  		## append p-value weight
  		wS <- seeds.snps[names(gr_SNP)]
  		GenomicRanges::mcols(gr_SNP) <- data.frame(mcols(gr_SNP), wS=wS)
  		
		if(verbose){
			now <- Sys.time()
			message(sprintf("\t %d SNPs are used as seeds", length(gr_SNP)), appendLF=T)
		}
  	}
  	
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=T)
	}
  	gr_Gene <- xRDataLoader(RData.customised="UCSC_genes", RData.location=RData.location, verbose=verbose)
    
	# genes: get all UCSC genes within 500k away from variants
	maxgap <- distance.max
	minoverlap <- 1L # 1b overlaps
	subject <- gr_Gene
	query <- gr_SNP
	q2r <- as.matrix(suppressWarnings(GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
	
	list_gene <- split(x=q2r[,1], f=q2r[,2])
	ind_gene <- as.numeric(names(list_gene))
	res_list <- lapply(1:length(ind_gene), function(i){
		x <- subject[ind_gene[i],]
		y <- query[list_gene[[i]],]
		dists <- GenomicRanges::distance(x, y, select="all", ignore.strand=T)
		
		## weights according to distance away from lead SNPs
		#wD <- 1- dists/maxgap
		wD <- 10^(-1*dists/maxgap)
		## weights according to P-values
		wS <- mcols(y)$wS
		
		## seeds weights according to wD and WP
		res <- max(wD * wS)
		names(res) <- mcols(x)$Symbol
		res
	})
	dist.seeds.genes <- unlist(res_list)
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\t %d Genes are defined as seeds according to genomic distance", length(dist.seeds.genes)), appendLF=T)
	}
    
    ######################################################
    # Link to targets based on eQTL
    ######################################################
    
    default.include.eQTL <- c("JKscience_TS2B","JKscience_TS3A")
	ind <- match(default.include.eQTL, include.eQTL)
	include.eQTL <- default.include.eQTL[!is.na(ind)]
    
    if(length(include.eQTL) > 0){
    
		res_list <- lapply(include.eQTL, function(x){
		
			if(x=='JKscience_TS2B'){
				# cis-eQTL
				cis <- xRDataLoader(RData.customised='JKscience_TS2B', RData.location=RData.location, verbose=verbose)
				minFDR <- apply(cis[,c(9:12)], 1, min, na.rm=T)
				df <- data.frame(SNP=cis[,1], Gene=cis[,4], FDR=minFDR, stringsAsFactors=F)
			}else if(x=='JKscience_TS3A'){
				# trans-eQTL
				trans <- xRDataLoader(RData.customised='JKscience_TS3A', RData.location=RData.location, verbose=verbose)
				minFDR <- apply(trans[,c(9:12)], 1, min, na.rm=T)
				df <- data.frame(SNP=trans[,1], Gene=trans[,4], FDR=minFDR, stringsAsFactors=F)
			}
			df
		})
		## get data frame (SNP Gene FDR)
		e_df <- do.call(rbind, res_list)
		uid <- paste(e_df[,1], e_df[,2], sep='_')
		df <- cbind(uid, e_df)
		res_list <- split(x=df$FDR, f=df$uid)
		res <- lapply(res_list, function(x){
			min(x)
		})
		vec <- unlist(res)
		raw_score <- -1*log10(vec)
		##  fit all FDR to the cumulative distribution function (CDF; depending on exponential empirical distributions)
		lamda <- MASS::fitdistr(raw_score, "exponential")$estimate
	
		## transformed score
		ind <- match(e_df[,1], names(seeds.snps))
		df <- data.frame(e_df[!is.na(ind),], wS=seeds.snps[ind[!is.na(ind)]])
		## weights according to eQTL (either lead SNPs or ld SNPs)
		wE <- stats::pexp(-log10(df$FDR), rate=lamda)
		## seeds weights according to wE and WP
		gw <- cbind(Gene=df$Gene, Weight=wE * df$wS)
		res_list <- split(x=gw[,2], f=gw[,1])
		res <- lapply(res_list, function(x){
			max(as.numeric(x))
		})
		eqtl.seeds.genes <- unlist(res)
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("\t %d Genes are defined as seeds according to eQTL mapping", length(eqtl.seeds.genes)), appendLF=T)
		}
		
		# merge both seed genes
		s_both <- c(dist.seeds.genes, eqtl.seeds.genes)
		df <- data.frame(names(s_both), s_both)
		res_list <- split(x=df[,2], f=df[,1])
		res <- lapply(res_list, function(x){
			sum(as.numeric(x))
		})
		seeds.genes <- unlist(res)
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("\t %d Genes are defined as seeds according to genomic distance and eQTL mapping", length(seeds.genes)), appendLF=T)
		}
	
	}else{
		seeds.genes <- dist.seeds.genes
	}
    
    
    ######################################################################################
    
    pNode <- suppressMessages(xPrioritiserGenes(data=seeds.genes, network=network, network.customised=network.customised, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose, RData.location=RData.location))
    
    pNode[['SNP']] <- df_snp
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(pNode)
}
