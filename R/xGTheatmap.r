#' Function to plot heatmap given associated data on tree tips
#'
#' \code{xGTheatmap} is supposed to plot heatmap given associated data on tree tips. It returns an object of class "ggplot".
#'
#' @param gp a ggplot object resulting from ggtree
#' @param data a data frame/matrix for coloring
#' @param reorder how to reorder rows and columns. It can be "none" for no reordering, "col" for reordering columns
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to the range of the finite values of displayed matrix
#' @param barwidth the width of the colorbar. Default value is 'legend.key.width' or 'legend.key.size' in 'theme' or theme
#' @param barheight the height of the colorbar. Default value is 'legend.key.height' or 'legend.key.size' in 'theme' or theme
#' @param nbin the number of bins for drawing colorbar 
#' @param legend.title the title of the colorbar. By default, it is ''
#' @param x.rotate the angle to rotate the x tick labelings. By default, it is 60
#' @param x.text.size the text size of the x tick labelings. By default, it is 6
#' @param x.text.hjust the hjust of the x tick labelings. By default, it is 0.5
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
#' @param ... additional graphic parameters for xHeatmap
#' @return 
#' a ggplot object
#' @note none
#' @export
#' @seealso \code{\link{xGT}}
#' @include xGTheatmap.r
#' @examples
#' \dontrun{
#' # Load the XGR package
#' library(XGR)
#' set.seed(825)
#' tree <- ape::rtree(50)
#' gp <- xGT(tree) + ggtree::geom_tiplab(size=2)
#' # heatmap data
#' x <- matrix(rnorm(length(tree$tip.label)*5), ncol=5)
#' colnames(x) <- paste0('C',1:5)
#' rownames(x) <- tree$tip.label
#' gp_heatmap <- xGTheatmap(gp, x)
#' }

xGTheatmap <- function(gp, data, reorder=c("none","col"), colormap="spectral", ncolors=64, zlim=NULL, barwidth=0.3, barheight=4, nbin=64, legend.title="", x.rotate=60, x.text.size=6, x.text.hjust=0.5, y.text.size=6, legend.text.size=5, legend.title.size=6, shape=19, size=2, plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt"), font.family="sans", na.color='grey80', data.label=NULL, label.size=1, label.color="black", ...)
{

    reorder <- match.arg(reorder)
    
    if(!all(class(gp) %in% c("ggtree","gg","ggplot"))){
    	warnings("The function must apply to a 'ggplot' object.\n")
        return(NULL)
    }
	
    if(any(class(data) %in% c("matrix","data.frame"))){
    	data <- as.matrix(data)
    }else{
    	return(NULL)
    }
	
	isTip <- y <- NULL
	
	df_tips <- gp$data %>% dplyr::filter(isTip) %>% dplyr::arrange(-y)
	ind <- match(df_tips$label, rownames(data))
	mat <- data[ind,]
    
	gp_heatmap <- xHeatmap(mat, reorder=reorder, colormap=colormap, ncolors=ncolors, zlim=zlim, barwidth=barwidth, barheight=barheight, nbin=nbin, legend.title=legend.title, x.rotate=x.rotate, x.text.size=x.text.size, x.text.hjust=x.text.hjust, y.text.size=y.text.size, legend.text.size=legend.text.size, legend.title.size=legend.title.size, shape=shape, size=size, plot.margin=plot.margin, font.family=font.family, na.color=na.color, data.label=data.label, label.size=label.size, label.color=label.color, ...)
	
	invisible(gp_heatmap)
}
