#' Function to visualise genes within a genomic region using track plot
#'
#' \code{xGRtrack} is supposed to visualise genes within a genomic region using track plot. Genes in query within a genomic region are displayed on the gene model track along with nearby genes of desired window or number. If scores for genomic region are also provided, the genomic score track will be also displayed at the top.
#'
#' @param cse.query a genomic region in query. By default it is NULL; otherwise provided as 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If provided, it will overwrite the parameter 'gene.query' below
#' @param gene.query which gene in query will be visualised. By default it is NULL
#' @param window the maximum distance defining nearby genes around the gene in query. By default it is 1e5
#' @param nearby the maximum number defining nearby genes around the gene in query. By default it is NULL. If not NULL, it will overwrite the parameter 'window' above
#' @param name.scoretrack the name for the score track. By default, it is "Genomic scores"
#' @param gene.model the genomic regions of the gene model. By default, it is 'UCSC_knownGene_model', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical_model', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
#' @param GR.score the genomic regions together with score data. By default, it is 'NA' to disable this option. Pre-built genomic score data: 'RecombinationRate' (recombintion rate, \url{http://www.ncbi.nlm.nih.gov/pubmed/17943122})), 'phastCons100way', 'phyloP100way', 'GERP'.
#' @param GR.score.customised the customised genomic score data. By default, it is NA to disable this option; otherwise load your customised GR object directly (with the first meta column for scores; if not provided, it will be valued at 1). If provided, it will be appended to 'GR.score' above. Also supported is the labelling by providing a meta-column called 'Label'. It can be also a list of GR objects
#' @param name.customised the name for customised genomic score data. By default, it is "Customised"
#' @param type.customised how to plot customised genomic score data. It can be "point" or "line"
#' @param label.size label size
#' @param label.col label color
#' @param label.force the repelling force between overlapping labels
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return a Tracks object.
#' @note none
#' @export
#' @seealso \code{\link{xGRoverlap}}
#' @include xGRtrack.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' ## given a query gene
#' tks <- xGRtrack(gene.query='TNF', nearby=10, gene.model="UCSC_knownGene_model", GR.score=c("RecombinationRate","phastCons100way"), RData.location=RData.location)
#' tks
#' ## given a query genomic region
#' tks <- xGRtrack(cse.query='chr6:31497996-31584798', gene.model="UCSC_knownGene_model", GR.score=c("RecombinationRate","phastCons100way"), RData.location=RData.location)
#' 
#' ### GWAS catalog
#' GWAScatalog <- xRDataLoader('GWAScatalog', RData.location=RData.location)
#' gwas <- xGR(GWAScatalog$cse_hg19, format="chr:start-end")
#' ind <- match(names(gwas), GWAScatalog$cse_hg19)
#' gwas$pvalue <- -log10(GWAScatalog$pvalue[ind])
#' tks <- xGRtrack(gene.query='TNF', nearby=10, gene.model="UCSC_knownGene_model", GR.score="RecombinationRate", GR.score.customised=gwas, RData.location=RData.location)
#' tks
#' 
#' ##########################
#' ## Advanced use: customised GR.score
#' ##########################
#' gene.model <- xRDataLoader("UCSC_knownGene_model", RData.location=RData.location)
#' 
#' ### LDblock_GR
#' LDblock_GR <- xRDataLoader("LDblock_GR", RData.location=RData.location)
#' gr <- LDblock_GR
#' GenomicRanges::mcols(gr) <- NULL
#' maf <- distance <- cadd <- gr
#' maf$maf <- LDblock_GR$maf
#' distance$distance <- LDblock_GR$distance
#' cadd$cadd <- LDblock_GR$cadd
#' ### GR.score.customised as a list of GR objects
#' GR.score.customised <- list(maf=maf, distance=distance, cadd=cadd)
#' tks <- xGRtrack(gene.query='TNF', window=1e0, gene.model=gene.model, GR.score=NA, GR.score.customised=GR.score.customised, type.customised='point', RData.location=RData.location)
#' tks
#'
#' ### the built-in provided as the customised
#' customised <- c("RecombinationRate","phastCons100way","phyloP100way","GERP","dbSNP_GWAS")
#' GR.score.customised <- lapply(customised, function(x) xRDataLoader(x, RData.location=RData.location))
#' tks <- xGRtrack(gene.query='TNF', nearby=10, gene.model=gene.model, GR.score=NA, GR.score.customised=GR.score.customised, type.customised='line', RData.location=RData.location)
#' tks
#' }

xGRtrack <- function(cse.query=NULL, gene.query=NULL, window=1e5, nearby=NULL, name.scoretrack="Genomic scores", gene.model=c("UCSC_knownGene_model","UCSC_knownCanonical_model"), GR.score=c(NA, "RecombinationRate","phastCons100way","phyloP100way","GERP"), GR.score.customised=NULL, name.customised="Customised", type.customised=c('point','line'), label.size=2, label.col="black", label.force=0.05, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    startT <- Sys.time()
    if(verbose){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }
	####################################################################################

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    type.customised <- match.arg(type.customised)
	
	gene.model <- gene.model
	if(class(gene.model) == "GRanges"){
		if(verbose){
			message(sprintf("Load customised gene model (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		gr_gene <- gene.model
		
	}else{
		if(verbose){
			message(sprintf("Load gene model '%s' (%s) ...", gene.model[1], as.character(Sys.time())), appendLF=TRUE)
		}
		gr_gene <- xRDataLoader(RData.customised=gene.model[1], verbose=F, RData.location=RData.location)
		if(is.null(gr_gene)){
			gene.model <- "UCSC_knownGene_model"
			if(verbose){
				message(sprintf("Instead, %s will be used", gene.model), appendLF=TRUE)
			}
			gr_gene <- xRDataLoader(RData.customised=gene.model, verbose=verbose, RData.location=RData.location)
		}
    }
    
    grl_gene <- GenomicRanges::split(gr_gene, gr_gene$Symbol)
    
	########################################
	if(is.null(cse.query)){
		if(!is.null(gene.query)){
			ind <- match(gene.query, names(grl_gene))
			
			if(length(ind) == 0 | is.na(ind)){
				warning(sprintf("\tNo found for queried %s", gene.query), appendLF=TRUE)
				return(NULL)
			}else{
				gr_tmp <- BiocGenerics::unlist(grl_gene[ind[!is.na(ind)]])
				
				if(is.null(nearby)){
				
					if(verbose){
						message(sprintf("View for the window %d centered on the gene %s (%s) ...", window, gene.query, as.character(Sys.time())), appendLF=TRUE)
					}
				
					## a window (eg 1e6) of upstream and downstream from the gene in query
					q2r <- as.matrix(as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=gr_gene, subject=gr_tmp, maxgap=window-1, minoverlap=0L, type="any", select="all", ignore.strand=TRUE))))
					genes_tmp <- unique(gr_gene$Symbol[q2r[,1]])
					
				}else{
					
					if(verbose){
						message(sprintf("View for the gene %s and its %d nearby genes (%s) ...", gene.query, nearby, as.character(Sys.time())), appendLF=TRUE)
					}
					
					Symbol <- Dist <- NULL
					
					dists <- as.matrix(as.data.frame(suppressWarnings(GenomicRanges::distance(x=gr_gene, y=gr_tmp, select="all", ignore.strand=TRUE))))
					gr_gene_tmp <- gr_gene[!is.na(dists[,1])]
					gr_gene_tmp$dist <- dists[!is.na(dists[,1])]
					df <- data.frame(Symbol=gr_gene_tmp$Symbol, Dist=gr_gene_tmp$dist, stringsAsFactors=F)
					df <- df %>% dplyr::group_by(Symbol) %>% dplyr::summarize(Dist=min(Dist)) %>% dplyr::arrange(Dist)
					genes_tmp <- df$Symbol[1:min(nearby+1,nrow(df))]
				}

				grl_subset <- grl_gene[genes_tmp]
				gr_subset <- BiocGenerics::unlist(grl_subset)
				
				chr <- unique(as.character(GenomeInfoDb::seqnames(gr_subset)))
				xlim <- c(min(BiocGenerics::start(gr_subset)), max(BiocGenerics::end(gr_subset)))
				cse.query <- paste0(chr,':',xlim[1],'-',xlim[2])
				gr_cse <- xGR(cse.query, format='chr:start-end')

				#seqlevels(grl_subset) <- chr
				#seqlengths(grl_subset) <- width(gr_cse)

			}
		}else{
			warning(sprintf("\tNone is queried %s", gene.query), appendLF=TRUE)
			return(NULL)
		}
	
	}else{
	
		gr_cse <- xGR(cse.query, format='chr:start-end')
		
		if(!is.null(gr_cse)){
			
			if(verbose){
				message(sprintf("View for %s (%s) ...", cse.query, as.character(Sys.time())), appendLF=TRUE)
			}
		
			q2r <- as.matrix(as.data.frame(GenomicRanges::findOverlaps(query=gr_gene, subject=gr_cse, maxgap=-1L, minoverlap=0L, type="any", select="all", ignore.strand=T)))
			
			if(length(q2r[,1])>0){
				genes_tmp <- unique(gr_gene$Symbol[q2r[,1]])
				grl_subset <- grl_gene[genes_tmp]
				gr_subset <- BiocGenerics::unlist(grl_subset)
			
				# update gr_cse
				chr <- unique(as.character(GenomeInfoDb::seqnames(gr_subset)))
				xlim <- c(min(BiocGenerics::start(gr_subset)), max(BiocGenerics::end(gr_subset)))
				cse.query <- paste0(chr,':',min(BiocGenerics::start(gr_cse),xlim[1]),'-',max(BiocGenerics::end(gr_cse),xlim[2]))
				gr_cse <- xGR(cse.query, format='chr:start-end')
				
			}else{
				grl_subset <- gr_subset <- NULL
			}
		}else{
			warning(sprintf("\tNone is queried %s", gr_cse), appendLF=TRUE)
			return(NULL)
		}
	}

	########################################
	## gp_model
	gp_model <- NULL
	
	if(length(gr_subset)>0){
	
		if(verbose){
			message(sprintf("Plot gene model '%s' (%s) ...", gene.model[1], as.character(Sys.time())), appendLF=TRUE)
		}
	
		Symbol <- NULL
		## gp_model
		if(0){
			gp <- ggplot(grl_subset) + ggbio::geom_alignment(facets=Symbol~., gap.geom="arrow", label=T, label.color="blue", length=unit(0.1,"cm"), color='lightblue', fill='lightblue')
		}else{
			gp <- ggbio::autoplot(grl_subset, aes(fill=Symbol,color=Symbol), gap.geom="arrow", label=T, label.color="grey20", length=unit(0.1,"cm"))
		}
		gp_model <- suppressWarnings(suppressMessages(gp + ggbio::theme_alignment(base_size=8)))
	
	}
	
	########################################
	
	## gp_score
	gp_score <- NULL
	
	default.GR.score <- c("RecombinationRate", "phastCons100way", "phyloP100way")
	names(default.GR.score) <- c('RecombRate', 'PhastCons', 'PhyloP')
	ind <- match(default.GR.score, GR.score)
	GR.score <- default.GR.score[!is.na(ind)]
	ls_gr_anno <- NULL
	if(length(GR.score) > 0){
	
		if(verbose){
			message(sprintf("Plot genomic score '%s' (%s) ...", paste(GR.score,collapse=','), as.character(Sys.time())), appendLF=TRUE)
		}

		ls_gr_anno <- lapply(1:length(GR.score), function(i){
			x <- GR.score[i]
			
			if(verbose){
				message(sprintf("\tcalculating '%s' genomic score (%s) ...", x, as.character(Sys.time())), appendLF=TRUE)
			}
			
			gr <- xGRoverlap(data=gr_cse, format="GRanges", GR.score=x, verbose=F, RData.location=RData.location)
			gr$Anno <- names(x)
			GenomicRanges::mcols(gr) <- GenomicRanges::mcols(gr)[,c('GScore','Anno')]
			gr$Label <- ''
			gr
		})
		names(ls_gr_anno) <- names(GR.score)
	}
	
	if(class(GR.score.customised)=='GRanges'){
		GR.score.customised <- list(GR.score.customised)
		names(GR.score.customised) <- name.customised
	}
	
	if(class(GR.score.customised)=='list'){
		
		if(is.null(names(GR.score.customised))){
			names(GR.score.customised) <- paste0("C", 1:length(GR.score.customised))
		}
		
		if(verbose){
			message(sprintf("Plot %d customised genomic score '%s' (%s) ...", length(GR.score.customised), paste(names(GR.score.customised),collapse=','), as.character(Sys.time())), appendLF=TRUE)
		}
		
		for(i in 1:length(GR.score.customised)){
			x <- GR.score.customised[[i]]
			if(class(x)=='GRanges'){
				gr <- xGRoverlap(data=gr_cse, format="GRanges", GR.score=x, verbose=F, RData.location=RData.location)
				if(length(gr)>0){
					## append 'Anno'
					gr$Anno <- names(GR.score.customised)[i]
					## append 'Label'
					if(is.null(gr$Label)){
						gr$Label <- ''
					}
					## extract 'GScore','Anno','Label'
					GenomicRanges::mcols(gr) <- GenomicRanges::mcols(gr)[,c('GScore','Anno','Label')]
			
					## just in case GR.score is NA
					if(length(ls_gr_anno)==0){
						ls_gr_anno <- list(gr)
						names(ls_gr_anno) <- names(GR.score.customised)[i]
					}else{
						ls_gr_anno[[length(ls_gr_anno)+1]] <- gr
						names(ls_gr_anno)[length(ls_gr_anno)] <- names(GR.score.customised)[i]
					}
			
				}else{
					if(verbose){
						message(sprintf("\tnone found in the region"), appendLF=TRUE)
					}
				}	
			}	
		}
	}
	
	if(!is.null(ls_gr_anno)){
		gr_anno <- BiocGenerics::unlist(GenomicRanges::GRangesList(ls_gr_anno))
		
		x <- y <- anno <- score <- label <- NULL

		if(1){
			df <- as.data.frame(gr_anno, row.names=NULL)
			df_start <- df[,c('start','GScore','Anno','Label')]
			colnames(df_start) <- c('x','y','anno','label')
			df_end <- df[,c('end','GScore','Anno','Label')]
			colnames(df_end) <- c('x','y','anno','label')
			data <- rbind(df_start, df_end)
			data <- unique(data)
			
			#######
			## order by Anno
			data$anno <- factor(data$anno, levels=names(ls_gr_anno))
			#######
						
			## customised
			gp_score <- ggplot(data=subset(data,anno %in% names(GR.score.customised)), aes(x=x,y=y,color=anno)) + geom_point(size=0.5)
			if(type.customised=='line'){
				gp_score <- gp_score + geom_line(size=0.8, alpha=0.5)
			}
			
			## build in
			gp_score <- gp_score + geom_line(data=subset(data,!(data$anno %in% names(GR.score.customised))), aes(x=x,y=y,color=anno), size=0.8, alpha=0.5)
			
			## facet_grid
			gp_score <- gp_score + xlab('') + ylab(name.scoretrack) + facet_grid(anno~., scales="free_y")
			
			## label for customised
			if(sum(data$label!="")){
				gp_score <- gp_score + ggrepel::geom_text_repel(data=subset(data,label!=''), aes(x=x,y=y,label=label), size=label.size, color=label.col, force=label.force, fontface='bold.italic', point.padding=unit(0.1,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')))	
			}

		}else{
			gp_score <- ggbio::autoplot(gr_anno, geom="bar", aes(fill=score)) + ggbio::scale_fill_fold_change() + ggbio::theme_alignment() + theme(legend.position="none") + xlim(xlim)
		}
	}
		
	####################
	alpha_my <- function(colour, alpha=NA){
		col <- grDevices::col2rgb(colour, TRUE)/255
		if (length(colour) != length(alpha)) {
			if (length(colour) > 1 && length(alpha) > 1) {
				stop("Only one of colour and alpha can be vectorised")
			}
			if (length(colour) > 1) {
				alpha <- rep(alpha, length.out = length(colour))
			}
			else if (length(alpha) > 1) {
				col <- col[, rep(1, length(alpha)), drop = FALSE]
			}
		}
		alpha[is.na(alpha)] <- col[4, ][is.na(alpha)]
		new_col <- grDevices::rgb(col[1, ], col[2, ], col[3, ], alpha)
		new_col[is.na(colour)] <- NA
		new_col
	}
	
	theme_tracks_my <- function(bg="#fffedb", alpha=1, ...){
			res <- ggbio::theme_clear(grid.x.major=FALSE, ...)
			attr(res, "track.plot.color") <- sapply(bg, alpha_my, alpha) 
			attr(res, "track.bg.color") <- bg
			attr(res, "label.bg.fill") <- "transparent"
			attr(res, "label.text.color") <- "black"
			attr(res, "label.text.cex") <- 0.8
			res
		}
	
	####################
	xlim <- c(min(GenomicRanges::start(gr_cse)), max(GenomicRanges::end(gr_cse)))
	
	tks <- NULL
	if(!is.null(gp_model) | !is.null(gp_score)){
		if(is.null(gp_model)){
			tks <- suppressWarnings(suppressMessages(ggbio::tracks(gp_score, xlab=cse.query, xlim=xlim)))
			
		}else{
			if(is.null(gp_score)){
				tks <- suppressWarnings(suppressMessages(ggbio::tracks(gp_model, xlab=cse.query, xlim=xlim)))
			}else{
				tks <- suppressWarnings(suppressMessages(ggbio::tracks(gp_model, gp_score, xlab=cse.query, xlim=xlim)))
			}
		}
		tks <- tks + theme_tracks_my(bg='transparent') + theme(legend.position="none",axis.title.y=element_text(size=8,color="black"), axis.text.y=element_text(size=6,color="black"),strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=7,face="bold.italic"))
	}

    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(verbose){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }

    invisible(tks)
}


