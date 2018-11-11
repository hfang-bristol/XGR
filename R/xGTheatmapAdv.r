#' Function to draw heatmap on tree tips
#'
#' \code{xGTheatmapAdv} is supposed to draw heatmap on tree tips.
#'
#' @param gp a ggplot object resulting from ggtree
#' @param data a data frame/matrix for coloring
#' @param ratio.width the width of the heatmap relative to the tree. Be default this ratio is 1 (that is, the equal width as the tree)
#' @param gap.width the gap between the heatmap and the tree. Be default it is 0.1 multiplied tree width
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to the range of the finite values of displayed matrix
#' @param barwidth the width of the colorbar. Default value is 'legend.key.width' or 'legend.key.size' in 'theme' or theme
#' @param barheight the height of the colorbar. Default value is 'legend.key.height' or 'legend.key.size' in 'theme' or theme
#' @param nbin the number of bins for drawing colorbar 
#' @param legend.title the title of the colorbar. By default, it is ''
#' @param legend.text.size the text size of the legend tick labelings. By default, it is 5
#' @param legend.title.size the text size of the legend titles. By default, it is 6
#' @param y.text.size the text size of the y tick labelings. By default, it is 6
#' @param x.rotate the angle to rotate the x tick labelings. By default, it is 60
#' @param x.text.size the text size of the x tick labelings. By default, it is 6
#' @param x.text.hjust the hjust of the x tick labelings. By default, it is 0.5
#' @param ... additional graphic parameters for ggtree::gheatmap
#' @return 
#' a ggplot object
#' @note none
#' @export
#' @seealso \code{\link{xHeatmap}}
#' @include xHeatmapAdv.r
#' @examples
#' \dontrun{
#' # Load the XGR package
#' library(XGR)
#' set.seed(825)
#' tree <- ape::rtree(50)
#' gp <- xGGtree(tree)
#' # heatmap data
#' x <- matrix(rnorm(length(tree$tip.label)*5), ncol=5)
#' colnames(x) <- paste0('C',1:5)
#' rownames(x) <- tree$tip.label
#' gp_heatmap <- xGTheatmapAdv(gp, x, ratio.width=0.5)
#' gp_heatmap
#' gp_heatmap + coord_polar(theta="y")
#' }

xGTheatmapAdv <- function(gp, data, ratio.width=1, gap.width=NULL, colormap="spectral", ncolors=64, zlim=NULL, barwidth=0.5, barheight=4, nbin=64, legend.title="", legend.text.size=5, legend.title.size=6, y.text.size=2, x.rotate=0, x.text.size=2, x.text.hjust=0.5, ...)
{
    
    if(!all(class(gp) %in% c("ggtree","gg","ggplot"))){
    	warnings("The function must apply to a 'ggplot' object.\n")
        return(NULL)
    }
	
    if(any(class(data) %in% c("matrix","data.frame"))){
    	data <- as.matrix(data)
    }else{
    	return(NULL)
    }
	
	if(is.null(gap.width)){
		gap.width <- max(gp$data$x) * 0.1
	}
	
	isTip <- y <- NULL
	
	df_tips <- gp$data %>% dplyr::filter(isTip) %>% dplyr::arrange(-y)
	ind <- match(df_tips$label, rownames(data))
	mat <- data[ind,]
    
    gp <- gp + ggtree::geom_tiplab(size=y.text.size,align=T)
    
	gp_heatmap <- gp %>% ggtree::gheatmap(mat, width=ratio.width, offset=gap.width, font.size=x.text.size, colnames_angle=x.rotate, hjust=x.text.hjust, ...) %>% ggtree::scale_x_ggtree() + scale_fill_gradientn(colors=xColormap(colormap)(ncolors), limits=zlim, guide=guide_colorbar(title=legend.title,title.position="top",barwidth=barwidth,barheight=barheight,nbin=nbin,draw.ulim=FALSE,draw.llim=FALSE)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.title=element_text(face="bold",color="black",size=legend.title.size),legend.text=element_text(face="bold",color="black",size=legend.text.size),legend.title.align=0.5)
	
	invisible(gp_heatmap)
}
