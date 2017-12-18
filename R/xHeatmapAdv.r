#' Function to draw heatmap together with sidebars on rows using ggplot2
#'
#' \code{xHeatmapAdv} is supposed to draw heatmap together with sidebars on rows using ggplot2.
#'
#' @param data.main a data frame/matrix for main heatmap. The coloring can be continuous (numeric matrix) or discrete (factor matrix)
#' @param data.meta a data frame/matrix for metadata visualisation. The per-column coloring can be continuous (numeric) or discrete (factor)
#' @param reorder how to reorder rows and columns. It can be "none" for no reordering, "row" for reordering rows according to number of sharings (by default), "col" for reordering columns, and "both" for reordering rows and columns
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to the range of the finite values of displayed matrix
#' @param barwidth the width of the colorbar. Default value is 'legend.key.width' or 'legend.key.size' in 'theme' or theme
#' @param barheight the height of the colorbar. Default value is 'legend.key.height' or 'legend.key.size' in 'theme' or theme
#' @param nbin the number of bins for drawing colorbar 
#' @param legend.title the title of the colorbar. By default, it is ''
#' @param x.rotate the angle to rotate the x tick labelings. By default, it is 60
#' @param x.text.size the text size of the x tick labelings. By default, it is 6
#' @param y.text.size the text size of the y tick labelings. By default, it is 6
#' @param legend.text.size the text size of the legend tick labelings. By default, it is 5
#' @param legend.title.size the text size of the legend titles. By default, it is 6
#' @param shape the number specifying the shape. By default, it is 19
#' @param size the number specifying the shape size. By default, it is 2
#' @param plot.margin the margin (t, r, b, l) around plot. By default, it is unit(c(5.5,5.5,5.5,5.5),"pt")
#' @param font.family the font family for texts
#' @param na.color the color for NAs. By default, it is 'grey80'
#' @param data.label a data frame/matrix used for the labelling
#' @param label.size the label size
#' @param label.color the label color
#' @param meta.colormap the colormap for metadata
#' @param meta.x.rotate the angle to rotate the x tick labelings for the metadata. By default, it is 90
#' @param meta.shape.continuous the number specifying the shape for continuous metadata. By default, it is 15
#' @param meta.shape.discrete the number specifying the shape for discrete metadata. By default, it is 95
#' @param meta.size the number specifying the shape size for metadata. By default, it is 2
#' @param meta.location the location of metadata. It can be "right" or "left"
#' @param meta.width the width for each column in metadata. By default, it is 0.5 (relative to each column in main data)
#' @param gap.width the width for the gap between panels. By default, it is 0.5 (relative to each column in main data)
#' @param legend.width the width for the legend. By default, it is NULL automatically determined (which can be used as a reference to define later for the better visualisation)
#' @param legend.direction the direction of the legend. It can be "vertical" or "horizontal"
#' @param legend.nrow the row number for legends. By default, it is 3 for the vertical direction of the legend; 6 for the horizontal direction of the legend
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param ... additional graphic parameters for supraHex::visTreeBootstrap
#' @return 
#' a gtable object
#' @note none
#' @export
#' @seealso \code{\link{xHeatmap}}
#' @include xHeatmapAdv.r
#' @examples
#' \dontrun{
#' # Load the XGR package
#' library(XGR)
#' data(mtcars)
#' data.main <- mtcars[,1:6]
#' data.meta <- mtcars[,7:11]
#' gt <- xHeatmapAdv(data.main, data.meta, barwidth=0.3, barheight=2.5, meta.location="right", legend.nrow=3, meta.width=0.4, gap.width=0.2, legend.width=NULL)
#' gt <- xHeatmapAdv(data.main, data.meta, barwidth=0.3, barheight=4, meta.location="right", legend.nrow=6, meta.width=0.4, gap.width=0.2, legend.width=4)
#' dev.new(); grid.draw(gt)
#' }

xHeatmapAdv <- function(data.main, data.meta, reorder=c("none","row","col","both"), colormap="spectral", ncolors=64, zlim=NULL, barwidth=0.3, barheight=4, nbin=64, legend.title="Main", x.rotate=60, x.text.size=6, y.text.size=6, legend.text.size=5, legend.title.size=6, shape=19, size=2, plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt"), font.family="sans", na.color='grey80', data.label=NULL, label.size=1, label.color="black", meta.colormap="spectral",meta.x.rotate=75,meta.shape.continuous=15,meta.shape.discrete=95,meta.size=2,meta.location=c("right","left"), meta.width=0.5, gap.width=0.5, legend.width=NULL, legend.direction=c("vertical","horizontal"), legend.nrow=NULL, verbose=TRUE, ...)
{

    reorder <- match.arg(reorder)
    meta.location <- match.arg(meta.location)
    legend.direction <- match.arg(legend.direction)
    
    if(is.null(legend.nrow)){
    	if(legend.direction == "vertical"){
    		legend.nrow <- 3
    	}else if(legend.direction == "horizontal"){
    		legend.nrow <- 6
    	}
    }
    
    ## main gp
	gp_main <- xHeatmap(data.main, reorder=reorder, colormap=colormap, ncolors=ncolors, zlim=zlim, barwidth=barwidth, barheight=barheight, nbin=nbin, legend.title=legend.title, x.rotate=x.rotate, x.text.size=x.text.size, y.text.size=y.text.size, legend.text.size=legend.text.size, legend.title.size=legend.title.size, shape=shape, size=size, plot.margin=plot.margin, font.family=font.family, na.color=na.color, data.label=data.label, label.size=label.size, label.color=label.color, ...)
    
    if(legend.direction == "horizontal"){
		gp_main <- gp_main + theme(legend.position="right",legend.direction="horizontal") + guides(color=guide_colorbar(title=legend.title,title.position="top",barwidth=barheight,barheight=barwidth))
    }
    
    gp_main_void <- gp_main + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())
    
	#######################################
	if(is.null(rownames(data.meta))){
		rownames(data.meta) <- paste('R', 1:nrow(data.meta), sep=' ')
	}
	if(is.null(colnames(data.meta))){
		colnames(data.meta) <- paste('C', 1:ncol(data.meta), sep=' ')
	}
	
	## rows reordered according to gp_main
	main_y_texts <- rev(levels(gp_main$data$gene))
	main_y_ind <- match(main_y_texts, rownames(data.meta))
	meta <- data.meta %>% dplyr::arrange(main_y_ind)
	rownames(meta) <- rownames(data.meta)[main_y_ind]
	#######################################

    ## list of meta gp
	ls_gp_meta <- lapply(1:ncol(meta), function(i){
		data <- meta %>% dplyr::select(i)
		meta.legend.title <- colnames(meta)[i]
		colormap <- meta.colormap
		meta.shape <- meta.shape.continuous
		#meta.size <- 2.5
		#meta.plot.margin <- unit(c(5.5,0,5.5,0),"pt")
		meta.plot.margin <- unit(c(plot.margin[[1]],0,plot.margin[[3]],0),"pt")
		if(is.factor(unlist(data))){
			#shape_discrete <- c(4,8,95)
			meta.shape <- meta.shape.discrete
			if(meta.shape>25){
				meta.size <- meta.size * 1.5
			}
		}
		gp <- xHeatmap(data, reorder="none", colormap=meta.colormap, ncolors=ncolors, zlim=NULL, barwidth=barwidth, barheight=barheight, nbin=nbin, legend.title=meta.legend.title, x.rotate=meta.x.rotate, x.text.size=x.text.size, legend.text.size=legend.text.size, legend.title.size=legend.title.size, shape=meta.shape, size=meta.size, plot.margin=meta.plot.margin, font.family=font.family, na.color=na.color)
		
		if(legend.direction == "horizontal"){
			gp <- gp + theme(legend.position="right",legend.direction="horizontal") + guides(color=guide_colorbar(title=meta.legend.title,title.position="top",barwidth=barheight,barheight=barwidth))
		}
		
		return(gp)
	})
    
	ls_gp_meta_void <- lapply(ls_gp_meta, function(gp){
		gp_void <- gp + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.border=element_blank())

	})
	
	#######################################
	## define a function to extract the legend
	extract_gp_legend <- function(gp){
		gt <- ggplot_gtable(ggplot_build(gp)) 
		ind <- which(sapply(gt$grobs, function(x) x$name)=="guide-box")
		legend <- gt$grobs[[ind]]
		#grid.draw(legend)
		invisible(legend)
	}
	
	## extract legends
	ls_gp <- c(list(gp_main), ls_gp_meta)
	ls_gp_legend <- lapply(ls_gp, function(gp){
		extract_gp_legend(gp)
	})
	#gt_legend <- gridExtra::grid.arrange(grobs=ls_gp_legend, nrow=legend.nrow, as.table=TRUE)
	gt_legend <- gridExtra::arrangeGrob(grobs=ls_gp_legend, nrow=legend.nrow, as.table=TRUE, top='', padding=unit(c(1,0,0,0),"line"))
	#grid::grid.draw(gt_legend)
	
	#######################################
	## aligning axes (keep the same heights)
	ls_gp_void <- c(list(gp_main_void), ls_gp_meta_void)
	ls_gt_void <- lapply(ls_gp_void, function(gp){
		gt <- ggplot_gtable(ggplot_build(gp))
	})
	ls_heights <- lapply(ls_gt_void, function(gt){
		gt$heights[4:5]
		#gt$widths[2:3]
	})
	maxHeight <- do.call(grid::unit.pmax, ls_heights)
	ls_gt_void_aligned <- lapply(ls_gt_void, function(gt){
		gt$heights[4:5] <- maxHeight
		invisible(gt)
	})

	#######################################
	## piece together
	n1 <- ncol(data.main)
	n2 <- ncol(data.meta)
	gt_gap <- grid::rectGrob(gp=grid::gpar(col="transparent"))
	
	if(is.null(legend.width)){
		legend.width <- 3*meta.width*ceiling((n2+1)/legend.nrow)
	}
	
	if(meta.location=='right'){
		### main gap meta blank gap
		grobs <- c(ls_gt_void_aligned[1], list(gt_gap), ls_gt_void_aligned[-1], list(gt_gap), list(gt_legend))
		widths <- c(n1, gap.width, rep(meta.width,n2), gap.width, legend.width)

	}else if(meta.location=='left'){
		### gap meta main gap legend
		grobs <- c(list(gt_gap), ls_gt_void_aligned[-1], ls_gt_void_aligned[1], list(gt_gap), list(gt_legend))
		widths <- c(gap.width, rep(meta.width,n2), n1, gap.width, legend.width)

	}
	
	#gt_heatmap <- gridExtra::grid.arrange(grobs=grobs, nrow=1, widths=widths, padding=unit(0, "line"))
	gt_heatmap <- gridExtra::arrangeGrob(grobs=grobs, nrow=1, widths=widths)
	#grid.draw(gt_heatmap)
	
	if(verbose){
		message(sprintf("data.main: %d rows X %d columns", nrow(data.main), ncol(data.main)), appendLF=TRUE)
		n <- sum(sapply(data.meta, is.factor))
		m <- ncol(data.meta)-n 
		message(sprintf("data.meta: %d rows X %d columns (%d continuous + %d discrete)", nrow(data.meta), ncol(data.meta), m, n), appendLF=TRUE)
		message(sprintf("widths: %.2f (meta), %.2f (blank), %.2f (legend)", meta.width, gap.width, legend.width), appendLF=TRUE)
	}
	
	invisible(gt_heatmap)
}
