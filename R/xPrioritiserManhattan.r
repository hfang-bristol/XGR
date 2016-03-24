#' Function to visualise prioritised genes using manhattan plot
#'
#' \code{xPrioritiserManhattan} is supposed to visualise prioritised genes using manhattan plot. 
#'
#' @param pNode an object of class "pNode"
#' @param color a character vector for point colors to alternate
#' @param cex a numeric value for point size
#' @param highlight.top the number of the top targets to be highlighted
#' @param highlight.col the highlight colors
#' @param highlight.label.size the highlight label size
#' @param highlight.label.offset the highlight label offset
#' @param highlight.label.col the highlight label color
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @import ggbio
#' @importFrom ggplot2 theme element_text element_rect Position
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xPrioritiser}}, \code{\link{xPrioritiserSNPs}}, \code{\link{xPrioritiserGenes}}, \code{\link{xPrioritiserPathways}}
#' @include xPrioritiserManhattan.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' library(dnet)
#' library(GenomicRanges)
#' #library(qqman)
#' library(ggbio)
#'
#' RData.location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0"
#' # a) provide the seed nodes/genes with the weight info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get genes within 500kb away from AS GWAS lead SNPs
#' seeds.genes <- ImmunoBase$AS$genes_variants
#' ## seeds weighted according to distance away from lead SNPs
#' data <- 1- seeds.genes/500000
#'
#' # b) perform priority analysis
#' pNode <- xPrioritiserGenes(data=data, network="PCommonsDN_medium",restart=0.7, RData.location=RData.location)
#' 
#' # c) manhattan plot
#' mp <- xPrioritiserManhattan(pNode, highlight.top=10, RData.location=RData.location)
#' #pdf(file="Gene_manhattan.pdf", height=6, width=12, compress=TRUE)
#' print(mp)
#' #dev.off()
#' }

xPrioritiserManhattan <- function(pNode, color=c("darkred","darkgreen"), cex=0.5, highlight.top=10, highlight.col="deepskyblue", highlight.label.size=2, highlight.label.offset=0.02, highlight.label.col="darkblue", verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
{
    if (class(pNode) != "pNode" ){
        stop("The function must apply to a 'pNode' object.\n")
    }
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=T)
	}
  	gr_Gene <- xRDataLoader(RData.customised="UCSC_genes", RData.location=RData.location)
    
    ## ONLY restricted to genes with genomic locations
	ind <- match(rownames(pNode$priority), mcols(gr_Gene)$Symbol)
	p_gr <- gr_Gene[ind[!is.na(ind)],]
	p_matrix <- pNode$priority[!is.na(ind),]
	
	## append genomic locations to GR object
	gr <- p_gr
	GenomicRanges::mcols(gr) <- cbind(mcols(gr), p_matrix[,2:4])
	## for sorting
	chrlabs <- paste('chr', as.character(c(1:22,'X','Y')), sep='')
	eval(parse(text=paste("seqlevels(gr) <- chrlabs",sep="")))
	
	## seed only
	if(0){
		ind <- which(mcols(gr)$seed==1)
		gr <- gr[ind,]
	}
	
	## highlight points
	highlight.top <- as.integer(highlight.top)
    if ( highlight.top > length(gr) ){
        highlight.top <- length(gr)
    }
	df <- data.frame(index=1:length(gr), val=mcols(gr)$priority)
	ind_o <- df[with(df,order(-val))[1:highlight.top],1]
	gro <- gr[ind_o,]
	names(gro) <- mcols(gro)$Symbol
	
	## draw plot
    suppressWarnings(
    mp <- plotGrandLinear(gr, eval(parse(text=paste("aes(y=priority)",sep=""))), color=color, spaceline=T, cex=cex, ylab='Priority', highlight.gr=gro, highlight.col=highlight.col, highlight.label=F, highlight.label.size=highlight.label.size, highlight.label.offset=highlight.label.offset, highlight.label.col=highlight.label.col) + theme(axis.text.x=element_text(angle=45, hjust=1,color="black",size=12), panel.background=element_rect(fill=rgb(0.95,0.95,0.95,1)))
    )
	
	## use qqman
	if(0){
	## prepare data.frame for qqman
	Gene <- mcols(gr)$Symbol
	CHR <- as.character(seqnames(gr))
	CHR <- gsub('X', '23', CHR, perl=T)
	CHR <- gsub('Y', '24', CHR, perl=T)
	CHR <- as.numeric(gsub('^chr', '', CHR, perl=T))
	BP <- start(gr)
	P <- mcols(gr)$priority
	df <- data.frame(Gene=Gene, CHR=CHR, BP=BP, P=P, stringsAsFactors=F)
	chrlabs <- c(1:22,"X","Y")
	#manhattan(df, snp='Gene', chrlabs=NULL, suggestiveline=F, genomewideline=F, logp=F, ylim=c(0,max(df$P)), ylab='Priority', cex=0.5, cex.axis=0.6, col=c("blue4","orange3"))
	}
	
    invisible(mp)
}
