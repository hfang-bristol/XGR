#' Function to visualise genomic regions using karyogram plot
#'
#' \code{xGRkaryogram} is supposed to visualise genomic regions using manhattan plot. It returns an object of class "ggplot".
#'
#' @param gr a GenomicRange object. If the meta-column 'label' is not provided, it will the name of this object
#' @param cytoband logical to indicate whether cytoband will be displayed. By default, it sets to false
#' @param color the rect color. By default it is 'royalblue'
#' @param size the rect size
#' @param label logical to indicate whether to label the rect. By default, it sets to false
#' @param label.size the label size
#' @param label.col the label color ('magenta' by default)
#' @param label.force the repelling force between overlapping labels
#' @param label.query only query will be labelled. By default, it sets to NULL meaning all will be displayed. If labels in query can not be found, then all will be displayed
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return a ggplot object.
#' @note none
#' @export
#' @seealso \code{\link{xGRkaryogram}}
#' @include xGRkaryogram.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' ### GWAS catalog
#' GWAScatalog <- xRDataLoader('GWAScatalog', RData.location=RData.location)
#' gwas <- xGR(GWAScatalog$cse_hg19, format="chr:start-end")
#' ind <- match(names(gwas), GWAScatalog$cse_hg19)
#' names(gwas) <- GWAScatalog$snp_id_current[ind]
#' gwas$label <- names(gwas)
#' gp <- xGRkaryogram(gwas)
#' gp
#' }

xGRkaryogram <- function(gr, cytoband=F, color="royalblue", size=0.5, label=F, label.size=2, label.col="magenta", label.force=0.05, label.query=NULL, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

	if(verbose){
		message(sprintf("First, load hg19 ideogram (%s) ...", as.character(Sys.time())), appendLF=T)
	}
    hg19_ideogram <- xRDataLoader(RData.customised="hg19_ideogram", verbose=verbose, RData.location=RData.location)
	gr_cytoband <- hg19_ideogram$cytoband
	gr_ideogram <- hg19_ideogram$ideogram
	
	if(verbose){
		message(sprintf("Second, add seqlengths and sort seqlevels (%s) ...", as.character(Sys.time())), appendLF=T)
	}
	## add 'seqlengths' and use the sorted seqlevels
	GenomeInfoDb::seqlengths(gr) <- GenomeInfoDb::seqlengths(gr_ideogram)[names(GenomeInfoDb::seqlengths(gr))]
	seqnames <- paste0("chr", c(1:22,"X","Y"))
	ind <- match(seqnames, names(GenomeInfoDb::seqlengths(gr)))
	seqnames <- seqnames[!is.na(ind)]
	gr <- GenomeInfoDb::keepSeqlevels(gr, seqnames)

	if(verbose){
		if(cytoband){
			message(sprintf("Third, visualise karyogram with cytoband (%s) ...", as.character(Sys.time())), appendLF=T)
		}else{
			message(sprintf("Third, visualise karyogram without cytoband (%s) ...", as.character(Sys.time())), appendLF=T)
		}
	}
	if(cytoband){
		gp <- ggplot(gr_cytoband) + ggbio::layout_karyogram(cytoband=T)
	}else{
		gp <- ggbio::autoplot(GenomeInfoDb::seqinfo(gr), layout="karyogram")
	}
	
	gp <- gp + theme(strip.background=element_rect(fill="transparent",color="transparent"))
	
	gp <- suppressMessages(gp + ggbio::layout_karyogram(gr, geom=c("rect","point")[1], ylim=c(10,40), color=color, size=size))
	#gp <- suppressMessages(gp + ggbio::layout_karyogram(gr, aes(x=start, y=num), geom=c("rect","point")[1], ylim=c(10,40), color=color, size=size))
	
	gp <- gp + theme(legend.position='none', axis.title.x=element_blank())
	
	if(label){
		df_highlight <- as.data.frame(gr)
		if(is.null(df_highlight$label)){
			df_highlight$label <- rownames(df_highlight)
		}
		#############################		
		## restrict to top in query for labels
		if(!is.null(label.query)){
			ind <- match(df_highlight$label, label.query)
			if(sum(!is.na(ind)) >= 1){
				df_highlight <- df_highlight[!is.na(ind), ]
			}else{
				df_highlight <- NULL
			}
		}
		#############################	
		
		start <- label <- NULL
		if(!is.null(df_highlight)){
			gp <- gp + ggrepel::geom_text_repel(data=df_highlight, aes(x=start,y=40,label=label), size=label.size, color=label.col, force=label.force, fontface='bold.italic', point.padding=unit(0.2,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')))
		}
	}
	
    kp <- gp
    kp$gr <- gr

    invisible(kp)
}


