#' Function to plot tree with associated data on tips
#'
#' \code{xGTplot} is supposed to plot tree with associated data on tips. It returns an object of class "ggplot".
#'
#' @param gp a ggplot object resulting from ggtree
#' @param data a data frame/matrix for coloring
#' @param type the plot type. It can be "bar" and "ridge"
#' @param combined logical to indicate whether the tree is combined. By default, it sets to true
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param y.text.size the text size of the y tick labelings. By default, it is 6
#' @param ... additional graphic parameters for ggridges::geom_density_ridges_gradient (such as scale=1)
#' @return 
#' a ggplot object
#' @note none
#' @export
#' @seealso \code{\link{xGT}}
#' @include xGTplot.r
#' @examples
#' \dontrun{
#' # Load the XGR package
#' library(XGR)
#' set.seed(825)
#' tree <- ape::rtree(50)
#' gp <- xGT(tree)
#' 
#' # barplot plot
#' data <- data.frame(id=tree$tip.label, val=rnorm(length(tree$tip.label)))
#' gp_bar <- xGTplot(gp, data, type="bar")
#' gp_bar
#' 
#' # ridge plot
#' id <- rep(tree$tip.label,each=20)
#' val <- sapply(1:length(tree$tip.label), function(i) rnorm(20,mean=i))
#' data <- data.frame(id=id, val=as.vector(val))
#' gp_ridge <- xGTplot(gp, data, type="ridge")
#' gp_ridge
#' 
#' # combine gp_bar and gp_ridge
#' g1 <- gtable::ggplotGrob(gp_bar)
#' g2 <- gtable::ggplotGrob(gp_ridge)
#' g <- cbind(g1, g2, size = "first")
#' g$heights <- grid::unit.pmax(g1$heights, g2$heights)
#' grid::grid.draw(g)
#' }

xGTplot <- function(gp, data, type=c("bar","ridge"), combined=T, colormap="spectral", ncolors=64, y.text.size=2, ...)
{

    type <- match.arg(type)
    
    if(!all(class(gp) %in% c("ggtree","gg","ggplot"))){
    	warnings("The function must apply to a 'ggplot' object.\n")
        return(NULL)
    }
	
    if(any(class(data) %in% c("matrix","data.frame"))){
    	data <- data.frame(id=data[,1], value=data[,2], stringsAsFactors=F)
    }else{
    	return(NULL)
    }
	
	isTip <- y <- value <- x <- label <- panel <- ..density.. <- NULL
	
	## df: bind dat to p$data
	df_tips <- gp$data %>% dplyr::filter(isTip) %>% dplyr::arrange(-y)
	ind <- match(data$id, df_tips$label)
	df <- data.frame(df_tips[ind[!is.na(ind)],], data[!is.na(ind),], stringsAsFactors=F) %>% dplyr::arrange(-y) 
	
	if(type=="bar"){
		if(combined){
			## add panel to df and gp$data
			df$panel <- factor('Bar', levels=c('Tree','Bar'))
			gp$data$panel <- factor('Tree', levels=c('Tree','Bar'))
			gp <- gp + ggtree::geom_tiplab(size=y.text.size,align=T)
			gp2 <- gp + geom_segment(data=df, aes(x=0, y=y, xend=value, yend=y,color=value), size=1) + scale_color_gradientn(colors=xColormap(colormap)(ncolors), guide=F) + facet_grid(.~panel, scales="free_x")
		}else{
			gp2 <- ggplot(df) + geom_segment(aes(x=0, y=y, xend=value, yend=y, color=value), size=1) + scale_color_gradientn(colors=xColormap(colormap)(ncolors), guide=F) + scale_y_continuous(breaks=1:nrow(df_tips),labels=rev(df_tips$label),limits=c(1,nrow(df_tips)),expand=c(0,0.5)) + theme_bw() + theme(legend.position="none",axis.text.y=element_text(size=y.text.size*4),axis.title.y=element_blank(), axis.title.x=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank())
		}
		
	}else if(type=="ridge"){
		if(combined){
			## add panel to df and gp$data
			df$panel <- factor('Ridge', levels=c('Tree','Ridge'))
			gp$data$panel <- factor('Tree', levels=c('Tree','Ridge'))
			gp <- gp + ggtree::geom_tiplab(size=y.text.size,align=T)
			gp2 <- gp + ggridges::geom_density_ridges_gradient(data=df, aes(x=value, group=label, fill=..density..), lwd=0.2, jittered_points=T, position=ggridges::position_points_jitter(width=0.05,height=0), point_shape='|', point_size=1, point_alpha=0.5, ...) + scale_fill_gradientn(colors=xColormap(colormap)(ncolors), guide=F) + facet_grid(.~panel, scales="free_x")
		}else{
			gp2 <- ggplot(df) + ggridges::geom_density_ridges_gradient(aes(x=value, y=y, group=label, fill=..density..), lwd=0.2, jittered_points=T, position=ggridges::position_points_jitter(width=0.05,height=0), point_shape='|', point_size=1, point_alpha=0.5, ...) + scale_fill_gradientn(colors=xColormap(colormap)(ncolors), guide=F) + scale_y_continuous(breaks=1:nrow(df_tips),labels=rev(df_tips$label),expand=c(0,0.5)) + theme_bw() + theme(legend.position="none",axis.text.y=element_text(size=y.text.size*4),axis.title.y=element_blank(), axis.title.x=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank())
		}

	}
	
	invisible(gp2)
}
