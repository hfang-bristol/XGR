#' Function to visualise genomic regions using manhattan plot
#'
#' \code{xGRmanhattan} is supposed to visualise genomic regions using manhattan plot. It returns an object of class "ggplot".
#'
#' @param gr a GenomicRange object with a meta-column 'value'. If the meta-column 'label' is not provided, it will the name of this object
#' @param chromosome.only logical to indicate whether only those from input data will be displayed. By default, it sets to TRUE
#' @param color a character vector for colors to alternate chromosome colorings. If NULL, ggplot2 default colors will be used. If a single character is provided, it can be "jet" (jet colormap) or "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta)
#' @param y.scale how to transform the y scale. It can be "normal" for no transformation, "sqrt" for square root transformation, and "log" for log-based transformation
#' @param y.lab the y labelling. If NULL (by default), it shows the column of input data
#' @param top the number of the top targets to be labelled/highlighted
#' @param top.label.type how to label the top targets. It can be "box" drawing a box around the labels , and "text" for the text only
#' @param top.label.size the highlight label size
#' @param top.label.col the highlight label color
#' @param top.label.force the repelling force between overlapping labels
#' @param top.label.query which top genes in query will be labelled. By default, it sets to NULL meaning all top genes will be displayed. If labels in query can not be found, then all will be displayed
#' @param label.query.only logical to indicate whether only those in query will be displayed. By default, it sets to FALSE. It only works when labels in query are enabled/found
#' @param top.label.chr logical to indicate whether the top hit per chromosome will be displayed. By default, it sets to TRUE. It only works when the parameter 'top' is null
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @return a ggplot object.
#' @note none
#' @export
#' @seealso \code{\link{xGRmanhattan}}
#' @include xGRmanhattan.r
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
#' gwas$value <- -log10(GWAScatalog$pvalue[ind])
#' names(gwas) <- GWAScatalog$snp_id_current[ind]
#' gwas$label <- names(gwas)
#' gp <- xGRmanhattan(gwas)
#' gp
#' }

xGRmanhattan <- function(gr, chromosome.only=TRUE, color=c("royalblue","sandybrown"), y.scale=c("normal","sqrt","log"), y.lab=NULL, top=NULL, top.label.type=c("text","box"), top.label.size=2, top.label.col="black", top.label.force=0.05, top.label.query=NULL, label.query.only=FALSE, top.label.chr=T, verbose=TRUE)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    y.scale <- match.arg(y.scale)
    top.label.type <- match.arg(top.label.type)

	if(!all(c("value") %in% colnames(GenomicRanges::mcols(gr)))){
		return(NULL)
	}
	
	if(is.null(names(gr))){
		gr <- xGR(gr, format='GRanges')
	}
	
	if(is.null(gr$label)){
		gr$label <- names(gr)
	}

	########################################
	if(label.query.only){
		if(!is.null(top.label.query)){
			top.label.query <- as.vector(t(top.label.query)) # just in case converting data.frame to vector
			ind <- match(names(gr), top.label.query)
			if(sum(!is.na(ind)) >= 1){
				gr <- gr[!is.na(ind)]
			}
		}
	}
	########################################

	## for sorting
	chrlabs <- paste('chr', as.character(c(1:22,'X','Y')), sep='')
	#######
	#if(chromosome.only){
	if(1){
		ind <- chrlabs %in% unique(as.character(gr@seqnames@values))
		chrlabs <- chrlabs[ind]
	}
	#######	
	#eval(parse(text="seqlevels(gr) <- chrlabs"))
	GenomeInfoDb::seqlevels(gr) <- chrlabs

	value <- seqnames <- NULL
	###############################
	## calling ggbio::autoplot
	suppressMessages(ggp <- ggbio::autoplot(object=gr, aes(y=value,color=seqnames), coord=c("genome","default")[1], geom='point', layout=c("linear","circle")[1], space.skip=0.01))
	
	## extract ggplot
	bp <- ggp@ggplot
	df_data <- bp$data
	df_data <- df_data %>% dplyr::arrange(-value)

	## alternative colors
	if(!is.null(color)){
		if(length(color)>=2){
			alternative_colors <- color
			chrs <- levels(df_data[,1])
			N <- length(chrs)
			cols <- rep(alternative_colors, round(N/length(alternative_colors)) + 1)[1:N]
			names(cols) <- chrs
			bp <- bp + scale_color_manual(values=cols) + theme(legend.position="none")
		}else if(length(color)==1){
			chrs <- levels(df_data[,1])
			N <- length(chrs)
			cols <- xColormap(color)(N)
			names(cols) <- chrs
			bp <- bp + scale_color_manual(values=cols) + theme(legend.position="none")
		}
	}else{
		bp <- bp + theme(legend.position="none")
	}

	## vline
  	if(TRUE){
		vline.df <- df_data
		vline.df <- do.call(rbind, by(vline.df, vline.df$seqnames, function(dd){
			data.frame(start=min(dd$start), end=max(dd$end))
		}))
		## compute gap
		gap <- (vline.df$start[-1] + vline.df$end[-nrow(vline.df)])/2
		bp <- bp + geom_vline(xintercept=gap, alpha=0.5, color='gray70') + theme(panel.grid=element_blank())
  	}
	#################################################################################

	## highlight points
	if(!is.null(top)){
		top <- as.integer(top)
		if(top > length(gr)){
			top <- length(gr)
		}
		
	}else{
		## the top per chromosome
		if(top.label.chr){
			top <- length(gr)
			
			seqnames <- value <- NULL
			
			df_chr <- df_data %>% dplyr::arrange(seqnames,-value)
			ind <- which(!duplicated(df_chr$seqnames))
			top.label.query <- df_chr$label[ind]
			
		}
		
	}
  	
	############
	## highlight top label
	############
	if(!is.null(top)){
		df_highlight <- df_data[1:top,]
		
		#############################		
		## restrict to top in query for labels
		if(!is.null(top.label.query)){
			ind <- match(df_highlight$label, top.label.query)
			if(sum(!is.na(ind)) >= 1){
				df_highlight <- df_highlight[!is.na(ind), ]
			}else{
				df_highlight <- NULL
			}
		}
		#############################		
		
		###########
		## potentially controlling only labels those in specific chromosome
		if(FALSE){
			ind <- match(df_highlight$seqnames, "chr1")
			if(sum(!is.na(ind)) >= 1){
				df_highlight <- df_highlight[!is.na(ind), ]
			}else{
				df_highlight <- NULL
			}
		}
		###########
		
		midpoint <- value <- label <- NULL
		if(!is.null(df_highlight)){
			if(top.label.type=="text"){
				bp <- bp + ggrepel::geom_text_repel(data=df_highlight, aes(x=midpoint,y=value,label=label), size=top.label.size, color=top.label.col, force=top.label.force, fontface='bold.italic', point.padding=unit(0.2,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')))
			}else if(top.label.type=="box"){
				bp <- bp + ggrepel::geom_label_repel(data=df_highlight, aes(x=midpoint,y=value,label=label), size=top.label.size, color=top.label.col, force=top.label.force, fontface='bold.italic', box.padding=unit(0.2,"lines"), point.padding=unit(0.2,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')))
			}
		}
	}
	
	#################################################################################
	
	## y scale
    if(y.scale=="sqrt"){
    	x <- NULL
    	bp <- bp + scale_y_continuous(trans=scales::sqrt_trans(), breaks=scales::trans_breaks("log10", function(x) 10^x, n=2))
    }else if(y.scale=="log"){
    	x <- NULL
    	bp <- bp + scale_y_continuous(trans=scales::log_trans(), breaks=scales::trans_breaks("log10", function(x) 10^x, n=2)) + annotation_logticks(sides='l')
    }
	
	if(!is.null(y.lab)){
		bp <- bp + ylab(y.lab)
	}
	
	bp <- bp + theme(axis.title.y=element_text(size=12), axis.text.y=element_text(color="black",size=8), axis.text.x=element_text(angle=45, hjust=1,color="black",size=10), panel.background=element_rect(fill=rgb(0.95,0.95,0.95,1)))
	
	## put arrows on y-axis and x-axis
	bp <- bp + theme(axis.line.y=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")), axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))

    mp <- bp
    mp$gr <- gr
  	
  	#mp + ggforce::facet_zoom(x=(seqnames=="chr2"))
  	
    invisible(mp)
}


