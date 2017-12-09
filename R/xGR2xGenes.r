#' Function to define genes from an input list of genomic regions given the crosslink info
#'
#' \code{xGR2xGenes} is supposed to define genes crosslinking to an input list of genomic regions (GR). Also required is the crosslink info with a score quantifying the link of a GR to a gene. Currently supported built-in crosslink info is enhancer genes and nearby genes (purely), though the user can customise it via 'crosslink.customised'; if so, it has priority over the built-in data.
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. It can be one of 'genehancer' (enhancer genes; PMID:28605766) or 'nearby' (nearby genes; if so, please also specify the relevant parameters 'nearby.distance.max', 'nearby.decay.kernel' and 'nearby.decay.exponent' below)
#' @param crosslink.customised the crosslink info with a score quantifying the link of a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if nor provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param cdf.function a character specifying how to transform the input crosslink score. It can be one of 'original' (no such transformation), and 'empirical'  for looking at empirical Cumulative Distribution Function (cdf; as such it is converted into pvalue-like values [0,1])
#' @param scoring logical to indicate whether gene-level scoring will be further calculated. By default, it sets to false
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param scoring.rescale logical to indicate whether gene scores will be further rescaled into the [0,1] range. By default, it sets to false
#' @param nearby.distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param nearby.decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param nearby.decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' If scoring sets to false, a data frame with following columns:
#' \itemize{
#'  \item{\code{GR}: genomic regions}
#'  \item{\code{Gene}: crosslinked genes}
#'  \item{\code{Score}: the score between the gene and the GR}
#'  \item{\code{Context}: the context}
#' }
#' If scoring sets to true, a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: crosslinked genes}
#'  \item{\code{Score}: gene score summarised over its list of crosslinked GR}
#'  \item{\code{Context}: the context}
#' }
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xGR}}
#' @include xGR2xGenes.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
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
#' # with their scores
#' df_xGenes <- xGR2xGenes(dGR, format="GRanges", crosslink.customised=crosslink.customised, scoring=T, scoring.scheme="max", RData.location=RData.location)
#' }

xGR2xGenes <- function(data, format=c("chr:start-end","data.frame","bed","GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("genehancer","nearby"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring=F, scoring.scheme=c("max","sum","sequential"), scoring.rescale=F, nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    crosslink <- match.arg(crosslink)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
	
	dGR <- xGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, RData.location=RData.location)
	
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
		if(is.vector(crosslink.customised)){
			# assume a file
			df <- utils::read.delim(file=crosslink.customised, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
		}else if(is.matrix(crosslink.customised) | is.data.frame(crosslink.customised)){
			df <- crosslink.customised
		}
		
		if(!is.null(df) & (ncol(df)==4 | ncol(df)==3)){
		
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
			df_SGS_customised <- SGS_customised[!is.na(SGS_customised[,1]) & !is.na(SGS_customised[,2]),]
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
		
		if(crosslink!='nearby'){
			if(verbose){
				now <- Sys.time()
				message(sprintf("Load the built-in crosslink '%s' (%s) ...", crosslink, as.character(now)), appendLF=T)
			}
		
			if(crosslink=="genehancer"){
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
			
			}

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
				message(sprintf("\t%d xGenes from an input list of %d (out of %d) genomic regions are defined as nearby genes within %d(bp) genomic distance window using '%s' decay kernel (%s)", length(unique(res_df$Gene)), length(unique(res_df$GR)), length(dGR), nearby.distance.max, nearby.decay.kernel, as.character(Sys.time())), appendLF=T)
			}
			
		}
	
	}
	
	df_xGenes <- res_df
	
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
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
  	
	invisible(df_xGenes)
}
