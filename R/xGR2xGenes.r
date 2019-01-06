#' Function to define genes from an input list of genomic regions given the crosslink info
#'
#' \code{xGR2xGenes} is supposed to define genes crosslinking to an input list of genomic regions (GR). Also required is the crosslink info with a score quantifying the link of a GR to a gene. Currently supported built-in crosslink info is enhancer genes, eQTL genes, conformation genes and nearby genes (purely), though the user can customise it via 'crosslink.customised'; if so, it has priority over the built-in data.
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. It can be one of 'genehancer' (enhancer genes; PMID:28605766), 'nearby' (nearby genes; if so, please also specify the relevant parameters 'nearby.distance.max', 'nearby.decay.kernel' and 'nearby.decay.exponent' below), 'PCHiC_combined' (conformation genes; PMID:27863249), 'GTEx_V6p_combined' (eQTL genes; PMID:29022597), 'eQTL_scRNAseq_combined' (eQTL genes; PMID:29610479), 'eQTL_jpRNAseq_combined' (eQTL genes; PMID:28553958), 'eQTL_ImmuneCells_combined' (eQTL genes; PMID:24604202,22446964,26151758,28248954,24013639)
#' @param crosslink.customised the crosslink info with a score quantifying the link of a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if not provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param cdf.function a character specifying how to transform the input crosslink score. It can be one of 'original' (no such transformation), and 'empirical' for looking at empirical Cumulative Distribution Function (cdf; as such it is converted into pvalue-like values [0,1])
#' @param scoring logical to indicate whether gene-level scoring will be further calculated. By default, it sets to false
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param scoring.rescale logical to indicate whether gene scores will be further rescaled into the [0,1] range. By default, it sets to false
#' @param nearby.distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param nearby.decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param nearby.decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' If scoring sets to false, a data frame with following columns:
#' \itemize{
#'  \item{\code{GR}: genomic regions}
#'  \item{\code{Gene}: crosslinked genes}
#'  \item{\code{Score}: the original score between the gene and the GR (if cdf.function is 'original'); otherwise cdf (based on the whole crosslink inputs)}
#'  \item{\code{Context}: the context}
#' }
#' If scoring sets to true, a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: crosslinked genes}
#'  \item{\code{Score}: gene score summarised over its list of crosslinked GR}
#'  \item{\code{Pval}: p-value-like significance level transformed from gene scores}
#'  \item{\code{Context}: the context}
#' }
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xGR}}
#' @include xGR2xGenes.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#'
#' # 1) provide the genomic regions
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' names(gr) <- NULL
#' dGR <- xGR(gr, format="GRanges")
#'
#' # 2) using built-in crosslink info
#' ## enhancer genes
#' df_xGenes <- xGR2xGenes(dGR, format="GRanges", crosslink="genehancer", RData.location=RData.location)
#' ## conformation genes
#' df_xGenes <- xGR2xGenes(dGR, format="GRanges", crosslink="PCHiC_combined", RData.location=RData.location)
#' ## eQTL genes
#' df_xGenes <- xGR2xGenes(dGR, format="GRanges", crosslink="GTEx_V6p_combined", RData.location=RData.location)
#' ## nearby genes (50kb, decaying rapidly)
#' df_xGenes <- xGR2xGenes(dGR, format="GRanges", crosslink="nearby", nearby.distance.max=50000, nearby.decay.kernel="rapid", RData.location=RData.location)
#'
#' # 3) advanced use
#' # 3a) provide crosslink.customised
#' ## illustration purpose only (see the content of 'crosslink.customised')
#' df <- xGR2nGenes(dGR, format="GRanges", RData.location=RData.location)
#' crosslink.customised <- data.frame(GR=df$GR, Gene=df$Gene, Score=df$Weight, Context=rep('C',nrow(df)), stringsAsFactors=F)
#' #crosslink.customised <- data.frame(GR=df$GR, Gene=df$Gene, Score=df$Weight, stringsAsFactors=F)
#' # 3b) define crosslinking genes
#' # without gene scoring
#' df_xGenes <- xGR2xGenes(dGR, format="GRanges", crosslink.customised=crosslink.customised, RData.location=RData.location)
#' # with gene scoring
#' df_xGenes <- xGR2xGenes(dGR, format="GRanges", crosslink.customised=crosslink.customised, scoring=T, scoring.scheme="max", RData.location=RData.location)
#' }

xGR2xGenes <- function(data, format=c("chr:start-end","data.frame","bed","GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("genehancer","PCHiC_combined","GTEx_V6p_combined","nearby"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring=F, scoring.scheme=c("max","sum","sequential"), scoring.rescale=F, nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, verbose=T, silent=F, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    startT <- Sys.time()
    if(!silent){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }else{
    	verbose <- FALSE
    }
    ####################################################################################
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    #crosslink <- match.arg(crosslink)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
	
	if(class(data)=='GRanges'){
		names(data) <- NULL
	}
	
	dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
	
	###################
	if(is.null(dGR)){
		return(NULL)
	}
	###################
	####################################################################################
	df_SGS_customised <- NULL
    if(!is.null(crosslink.customised)){
    
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the customised crosslink (%s) ...", as.character(now)), appendLF=T)
		}

		###########################	
		# customised df_SGS
		###########################
		df <- NULL
		if(is.vector(crosslink.customised)){
			# assume a file
			df <- utils::read.delim(file=crosslink.customised, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
		}else if(is.matrix(crosslink.customised) | is.data.frame(crosslink.customised)){
			df <- crosslink.customised
		}
		
		if(!is.null(df) && (ncol(df)==4 | ncol(df)==3)){
		
			if(ncol(df)==4){
				SGS_customised <- df
			}else{
				SGS_customised <- df
				SGS_customised$Context <- 'C'
			}
			colnames(SGS_customised) <- c("GR", "Gene", "Score", "Context")

			############################
			# remove Gene if NA
			# remove GR if NA
			# remove Score if NA
			df_SGS_customised <- SGS_customised[!is.na(SGS_customised[,1]) & !is.na(SGS_customised[,2]) & !is.na(SGS_customised[,3]),]
			############################
			
			if(verbose){
				message(sprintf("Genes (%d) and genomic regions (%d) are considered based on customised contexts (%d) (%s) ...", length(unique(df_SGS_customised[,2])), length(unique(df_SGS_customised[,1])), length(unique(df_SGS_customised[,4])), as.character(Sys.time())), appendLF=TRUE)
			}
		}
	
		#########################################
		if(!is.null(df_SGS_customised)){
			############################
			# remove Gene if ''
			# remove GR if ''
			df_SGS_customised <- df_SGS_customised[df_SGS_customised[,1]!='' & df_SGS_customised[,2]!='',]
			############################
		}
		
	}

	if(is.null(df_SGS_customised)){
		
		default.crosslink <- c("genehancer","PCHiC_combined","PCHiC_combined_PE","GTEx_V6p_combined","eQTL_ImmuneCells_combined","eQTL_eQTLGen","eQTL_scRNAseq_combined","eQTL_jpRNAseq_combined","FANTOM5_Cell","FANTOM5_Tissue", "REG_lncRNA","REG_enhancer", "nearby", "PCHiC_Monocytes","PCHiC_Macrophages_M0","PCHiC_Macrophages_M1","PCHiC_Macrophages_M2","PCHiC_Neutrophils","PCHiC_Megakaryocytes","PCHiC_Endothelial_precursors","PCHiC_Erythroblasts","PCHiC_Fetal_thymus","PCHiC_Naive_CD4_T_cells","PCHiC_Total_CD4_T_cells","PCHiC_Activated_total_CD4_T_cells","PCHiC_Nonactivated_total_CD4_T_cells","PCHiC_Naive_CD8_T_cells","PCHiC_Total_CD8_T_cells","PCHiC_Naive_B_cells","PCHiC_Total_B_cells", "PCHiC_PE_Monocytes","PCHiC_PE_Macrophages_M0","PCHiC_PE_Macrophages_M1","PCHiC_PE_Macrophages_M2","PCHiC_PE_Neutrophils","PCHiC_PE_Megakaryocytes","PCHiC_PE_Erythroblasts","PCHiC_PE_Naive_CD4_T_cells","PCHiC_PE_Naive_CD8_T_cells", "GTEx_V6p_Adipose_Subcutaneous","GTEx_V6p_Adipose_Visceral_Omentum","GTEx_V6p_Adrenal_Gland","GTEx_V6p_Artery_Aorta","GTEx_V6p_Artery_Coronary","GTEx_V6p_Artery_Tibial","GTEx_V6p_Brain_Anterior_cingulate_cortex_BA24","GTEx_V6p_Brain_Caudate_basal_ganglia","GTEx_V6p_Brain_Cerebellar_Hemisphere","GTEx_V6p_Brain_Cerebellum","GTEx_V6p_Brain_Cortex","GTEx_V6p_Brain_Frontal_Cortex_BA9","GTEx_V6p_Brain_Hippocampus","GTEx_V6p_Brain_Hypothalamus","GTEx_V6p_Brain_Nucleus_accumbens_basal_ganglia","GTEx_V6p_Brain_Putamen_basal_ganglia","GTEx_V6p_Breast_Mammary_Tissue","GTEx_V6p_Cells_EBVtransformed_lymphocytes","GTEx_V6p_Cells_Transformed_fibroblasts","GTEx_V6p_Colon_Sigmoid","GTEx_V6p_Colon_Transverse","GTEx_V6p_Esophagus_Gastroesophageal_Junction","GTEx_V6p_Esophagus_Mucosa","GTEx_V6p_Esophagus_Muscularis","GTEx_V6p_Heart_Atrial_Appendage","GTEx_V6p_Heart_Left_Ventricle","GTEx_V6p_Liver","GTEx_V6p_Lung","GTEx_V6p_Muscle_Skeletal","GTEx_V6p_Nerve_Tibial","GTEx_V6p_Ovary","GTEx_V6p_Pancreas","GTEx_V6p_Pituitary","GTEx_V6p_Prostate","GTEx_V6p_Skin_Not_Sun_Exposed_Suprapubic","GTEx_V6p_Skin_Sun_Exposed_Lower_leg","GTEx_V6p_Small_Intestine_Terminal_Ileum","GTEx_V6p_Spleen","GTEx_V6p_Stomach","GTEx_V6p_Testis","GTEx_V6p_Thyroid","GTEx_V6p_Uterus","GTEx_V6p_Vagina","GTEx_V6p_Whole_Blood", "eQTL_ImmuneCells_bcell","eQTL_ImmuneCells_Blood","eQTL_ImmuneCells_CD4","eQTL_ImmuneCells_CD8","eQTL_ImmuneCells_JKscience_CD14","eQTL_ImmuneCells_JKscience_IFN","eQTL_ImmuneCells_JKscience_LPS2","eQTL_ImmuneCells_JKscience_LPS24","eQTL_ImmuneCells_mono","eQTL_ImmuneCells_Neutrophils","eQTL_ImmuneCells_NK", "eQTL_scRNAseq_Bcell","eQTL_scRNAseq_PBMC","eQTL_scRNAseq_NK","eQTL_scRNAseq_Mono","eQTL_scRNAseq_DC","eQTL_scRNAseq_CD8","eQTL_scRNAseq_CD4", "eQTL_jpRNAseq_combined","eQTL_jpRNAseq_Bcell","eQTL_jpRNAseq_CD4","eQTL_jpRNAseq_CD8","eQTL_jpRNAseq_Mono","eQTL_jpRNAseq_NK","eQTL_jpRNAseq_PBMC", "PCHiC_PMID25938943_GM12878","PCHiC_PMID25938943_CD34")
		ind <- match(default.crosslink, crosslink)
		crosslink <- default.crosslink[!is.na(ind)]
		if(length(crosslink)==0){
			return(NULL)
		}else{
			## only keep the first one
			crosslink <- crosslink[1]
		}
		
		if(crosslink!='nearby'){
			if(verbose){
				now <- Sys.time()
				message(sprintf("Load the built-in crosslink '%s' (%s) ...", crosslink, as.character(now)), appendLF=T)
			}
		
			if(crosslink=="genehancer"){
				if(0){
					ig <- xRDataLoader('ig.genehancer', verbose=F, RData.location=RData.location)
					V(ig)$name <- V(ig)$id
					df_edges <- get.data.frame(ig, what="edges")
					df_nodes <- get.data.frame(ig, what="vertices")
					df_nodes <- subset(df_nodes, df_nodes$type=='enhancer')
					ind <- match(df_edges$from, df_nodes$id)
					df_edges$GR_score <- as.numeric(df_nodes$score[ind])
					## final score: Snode * Slink
					crosslink.customised <- data.frame(GR=df_edges$from, Gene=df_edges$to, Score=df_edges$score * df_edges$GR_score, Context=rep('genehancer',nrow(df_edges)), stringsAsFactors=F)
					df_SGS_customised <- crosslink.customised
				}else{
					df_SGS_customised <- xRDataLoader('crosslink.customised.genehancer', verbose=verbose, RData.location=RData.location)
				}
			
			}else if(sum(grep("PCHiC_",crosslink,perl=TRUE)) > 0){
				rdata <- paste0('crosslink.customised.', crosslink)
				df_SGS_customised <- xRDataLoader(rdata, verbose=verbose, RData.location=RData.location)
				
			}else if(sum(grep("GTEx_V6p_",crosslink,perl=TRUE)) > 0){
				rdata <- paste0('crosslink.customised.', crosslink)
				df_SGS_customised <- xRDataLoader(rdata, verbose=verbose, RData.location=RData.location)
					
			}else if(sum(grep("FANTOM5_",crosslink,perl=TRUE)) > 0){
				rdata <- paste0('crosslink.customised.', crosslink)
				df_SGS_customised <- xRDataLoader(rdata, verbose=verbose, RData.location=RData.location)
					
			}else{
				## general use
				rdata <- paste0('crosslink.customised.', crosslink)
				df_SGS_customised <- xRDataLoader(rdata, verbose=verbose, RData.location=RData.location)
			}

			############################
			# remove Gene if NA
			# remove GR if NA
			# remove Score if NA
			df_SGS_customised <- df_SGS_customised[!is.na(df_SGS_customised[,1]) & !is.na(df_SGS_customised[,2]) & !is.na(df_SGS_customised[,3]),]
			############################

			if(!is.null(df_SGS_customised)){

				if(verbose){
					message(sprintf("\tGenes (%d) and genomic regions (%d) are considered based on built-in '%s' (%s) ...", length(unique(df_SGS_customised[,2])), length(unique(df_SGS_customised[,1])), unique(df_SGS_customised[,4]), as.character(Sys.time())), appendLF=TRUE)
				}
			}
		}		
	}


	####################################################################################
	####################################################################################

	if(!is.null(df_SGS_customised)){
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("Define crosslinked genes from an input list of %d genomic regions (%s) ...", length(dGR), as.character(now)), appendLF=T)
		}
		
		Gene <- Weight <- Score <- NULL
		dgr <- NULL
		
		ls_df_SGS <- split(x=df_SGS_customised, f=df_SGS_customised$Context)
		ls_res_df <- lapply(1:length(ls_df_SGS), function(j){
			df_SGS <- ls_df_SGS[[j]]
			
			if(cdf.function=="empirical"){
				## Compute an empirical cumulative distribution function
				my.CDF <- stats::ecdf(df_SGS$Score)
				df_SGS$Weight <- my.CDF(df_SGS$Score)
				
			}else{
				df_SGS$Weight <- df_SGS$Score
			}
			
			gr <- xGR(df_SGS$GR, format="chr:start-end", verbose=verbose, RData.location=RData.location)
			
			q2r <- as.data.frame(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=0L, type="any", select="all", ignore.strand=TRUE))
			q2r$gr <- names(gr[q2r[,2]])
			q2r$dgr <- names(dGR[q2r[,1]])
			
			#################################
			## to reduce runtime significantly
			ind <- match(df_SGS$GR, q2r$gr)
			df_found <- df_SGS[!is.na(ind), ]
			#################################
			
			if(1){
				system.time({
				
				## very fast
				ls_df_found <- split(x=df_found[,c('GR','Gene','Weight')], f=df_found$GR)
				
				### df_found_reorder
				ind <- match(q2r$gr, names(ls_df_found))
				ls_df_found_reorder <- ls_df_found[ind]
				df_found_reorder <- do.call(rbind, ls_df_found_reorder)
				
				### df_q2r
				vec_nrow <- sapply(ls_df_found_reorder, nrow)
				ind_times <- rep(1:nrow(q2r),times=vec_nrow)
				df_q2r <- q2r[ind_times,]
				
				### df
				df <- cbind(df_q2r, df_found_reorder)
				
				#################################
				## keep maximum weight per gene and dgr if there are many overlaps
				#################################
				df_xGenes <- as.data.frame(df %>% dplyr::group_by(dgr,Gene) %>% dplyr::summarize(Score=max(Weight)))
				colnames(df_xGenes) <- c('GR','Gene','Score')
				
				})
				
			}else{

				system.time({
			
				## very slow
				ls_dgr <- split(x=q2r$gr, f=q2r$dgr)
				ls_df <- lapply(1:length(ls_dgr), function(i){
				
					#################
					## very important
					#################
					if(0){
					ind <- match(df_SGS$GR, ls_dgr[[i]])
					df <- df_SGS[!is.na(ind), ]
					}else{
					ind <- match(df_found$GR, ls_dgr[[i]])
					df <- df_found[!is.na(ind), ]
					}
				
					#################################
					## keep maximum weight per gene if there are many overlaps
					#################################
					df <- as.data.frame(df %>% dplyr::group_by(Gene) %>% dplyr::summarize(Score=max(Weight)))
					data.frame(GR=rep(names(ls_dgr)[i],nrow(df)), df, stringsAsFactors=F)
				})
				df_xGenes <- do.call(rbind, ls_df)

				})
			}
			
			########################################
			# check gene (make sure official symbol)
			ind <- !is.na(XGR::xSymbol2GeneID(df_xGenes$Gene, details=TRUE, verbose=FALSE, RData.location=RData.location)$Symbol)
			df_xGenes <- df_xGenes[ind,]
			########################################
	
			if(verbose){
				message(sprintf("\t%d xGenes (%d genomic regions) are defined based on built-in '%s' (%s)", length(unique(df_xGenes$Gene)), length(unique(df_xGenes$GR)), names(ls_df_SGS)[j], as.character(Sys.time())), appendLF=T)
			}
		
			############################################
			## whether gene scoring
			if(scoring){
				
				## summaryFun is applied over GR for each seed gene
				
				## calculate genetic influence score under a set of GR for each seed gene
				if(scoring.scheme=="max"){
					summaryFun <- max
				}else if(scoring.scheme=="sum"){
					summaryFun <- sum
				}else if(scoring.scheme=="sequential"){
					summaryFun <- function(x){
						base::sum(x / base::rank(-x,ties.method="min"))
					}
				}
				
				df_xGenes <- as.data.frame(df_xGenes %>% dplyr::group_by(Gene) %>% dplyr::summarise(Score=summaryFun(Score)))

				if(verbose){
					now <- Sys.time()
					message(sprintf("\t%d xGenes are scored using '%s' scoring scheme (%s)", length(unique(df_xGenes$Gene)), scoring.scheme, as.character(now)), appendLF=T)
				}

				if(scoring.rescale){
					if(verbose){
						now <- Sys.time()
						message(sprintf("\talso rescale score into the [0,1] range (%s)", as.character(now)), appendLF=T)
					}
					# rescale to [0 1]
					rescaleFun <- function(x){
						(x - min(x))/(max(x) - min(x))
					}
					
					df_xGenes$Score <- rescaleFun(df_xGenes$Score)
				}
				
				##############
				## df_xGenes$Pval
				##############
				seeds.genes <- df_xGenes$Score
				names(seeds.genes) <- df_xGenes$Gene
				# rescale to [0.100001 1]
				rescaleFun <- function(x){
					0.100001 + 0.9*0.99999888888*(x - min(x))/(max(x) - min(x))
				}
	
				x <- rescaleFun(seeds.genes)
				# convert into pvalue by 10^(-x*10)
				# [1e-10, 0.0999977]
				pval <- 10^(-x*10)
				
				df_xGenes$Pval <- pval
				##############
		
			}
	
			data.frame(df_xGenes, Context=rep(names(ls_df_SGS)[j],nrow(df_xGenes)), stringsAsFactors=F)
		})
		res_df <- do.call(rbind, ls_res_df)
	
	}else{
		
		## only for the option 'nearby'
		if(crosslink=='nearby'){
			df <- xGR2nGenes(data=dGR, format="GRanges", distance.max=nearby.distance.max, decay.kernel=nearby.decay.kernel, decay.exponent=nearby.decay.exponent, GR.Gene="UCSC_knownGene", scoring=scoring, scoring.scheme=scoring.scheme, scoring.rescale=scoring.rescale, verbose=F, RData.location=RData.location)
			
			context <- paste0('nearby_',nearby.distance.max,'_',nearby.decay.kernel)
			if(scoring){
				res_df <- data.frame(df, Context=rep(context,nrow(df)), stringsAsFactors=F)
			}else{
				res_df <- data.frame(GR=df$GR, Gene=df$Gene, Score=df$Weight, Context=rep(context,nrow(df)), stringsAsFactors=F)
			}
			
			if(verbose){
				if(scoring){
					message(sprintf("\t%d xGenes (out of %d genomic regions) are defined as nearby genes within %d(bp) genomic distance window using '%s' decay kernel (%s)", length(unique(res_df$Gene)), length(dGR), nearby.distance.max, nearby.decay.kernel, as.character(Sys.time())), appendLF=T)
				}else{
					message(sprintf("\t%d xGenes from an input list of %d (out of %d) genomic regions are defined as nearby genes within %d(bp) genomic distance window using '%s' decay kernel (%s)", length(unique(res_df$Gene)), length(unique(res_df$GR)), length(dGR), nearby.distance.max, nearby.decay.kernel, as.character(Sys.time())), appendLF=T)
				}
			}
			
		}
	
	}
	
	df_xGenes <- res_df
	
	##### order
	Context <- GR <- Score <- NULL
	if(scoring){
		df_xGenes <- df_xGenes %>% dplyr::arrange(Context, -Score)
	}else{
		df_xGenes <- df_xGenes %>% dplyr::arrange(-Score)
		#### sort GR by chromosome, start and end
		ind <- xGRsort(df_xGenes$GR)
		if(!is.null(ind)){
			df_xGenes <- df_xGenes[ind,]
			####
			df_xGenes <- df_xGenes %>% dplyr::arrange(Context)
		}
	}
	
	####################################
	## also output igraph (genes with genomic location)
	if(0){
		GR.Gene <- "UCSC_knownGene"
		gr_Gene <- xRDataLoader(RData.customised=GR.Gene, verbose=FALSE, RData.location=RData.location)
		
		tmp_df <- df_xGenes
		ind <- match(tmp_df$Gene, names(gr_Gene))
		tmp_df <- tmp_df[!is.na(ind),]
		tmp <- as.data.frame(gr_Gene)[ind[!is.na(ind)],]
		tmp_df$Gene_gr <- paste0(tmp$seqnames,':',tmp$start,'-',tmp$end)
		
		####################################
		## also output igraph
		context_ls <- split(x=tmp_df, f=tmp_df$Context)
		ls_ig <- lapply(1:length(context_ls), function(i){
			df <- context_ls[[i]]
			
			## edges
			relations <- df[,1:3]
			
			## nodes
			df_gr <- data.frame(df[,c("GR","GR")], type=rep('GR',nrow(df)), stringsAsFactors=F)
			df_gene <- data.frame(df[,c("Gene","Gene_gr")], type=rep('Gene',nrow(df)), stringsAsFactors=F)
			nodes <- base::as.data.frame(rbind(as.matrix(df_gr), as.matrix(df_gene)), stringsAsFactors=FALSE)
			nodes <- nodes[!duplicated(nodes),]
			colnames(nodes) <- c("name","id","type")
			ig <- graph.data.frame(d=relations, directed=TRUE, vertices=nodes)
			
			## Circos plot
			if(0){
				GR.SNP <- xGR(data=V(ig)$id[V(ig)$type=='GR'], format="chr:start-end", verbose=FALSE, RData.location=RData.location)
				names(GR.SNP) <- V(ig)$name[V(ig)$type=='GR']
				GR.Gene <- xGR(data=V(ig)$id[V(ig)$type=='Gene'], format="chr:start-end", verbose=FALSE, RData.location=RData.location)
				names(GR.Gene) <- V(ig)$name[V(ig)$type=='Gene']
				#library(RCircos)
				xCircos(g=ig, entity="Both", top_num=10, GR.SNP=GR.SNP, GR.Gene=GR.Gene, colormap=c("orange-darkred"), ideogram=F, chr.exclude="auto", entity.label.cex=0.6, entity.label.side="out", RData.location=RData.location)
			}
	
			return(ig)
			
		})
		names(ls_ig) <- names(context_ls)
		####################################
	}
	
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xGR2xGenes): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
  	
	invisible(df_xGenes)
}
