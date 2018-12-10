#' Function to draw multiple component planes reorded within a sheet-shape rectangle grid using ggplot2
#'
#' \code{xSHreorder} is supposed to draw multiple component planes reorded within a sheet-shape rectangle grid using ggplot2.
#'
#' @param sMap an object of class "sMap"
#' @param sReorder an object of class "sReorder"
#' @param title.size the title size
#' @param embed.width the width/height of embedded ggplot obejcts
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to the range of the finite values of displayed codebook matrix
#' @param border.color the border color for each hexagon
#' @param barwidth the width of the colorbar. Default value is 'legend.key.width' or 'legend.key.size' in 'theme' or theme
#' @param barheight the height of the colorbar. Default value is 'legend.key.height' or 'legend.key.size' in 'theme' or theme
#' @param nbin the number of bins for drawing colorbar 
#' @param legend.title the title of the colorbar. By default, it is ''
#' @param legend.text.size the text size of the legend tick labelings
#' @param legend.title.size the text size of the legend titles
#' @return 
#' a ggplot2 object
#' @note none
#' @export
#' @seealso \code{\link{xSHreorder}}
#' @include xSHreorder.r
#' @examples
#' \dontrun{
#' library(XGR)
#' 
#' gp <- xSHreorder(sMap, sReorder)
#' gp
#' ls_gp <- gp$ls_gp
#' }

xSHreorder <- function(sMap, sReorder, title.size=8, embed.width=1, colormap="spectral", zlim=NULL, border.color="transparent", barwidth=0.4, barheight=NULL, nbin=64, legend.title='', legend.text.size=6, legend.title.size=8)
{

    if(class(sMap) != "sMap"){
        stop("The funciton must apply to 'sMap' object.\n")
    }
    
    codebook <- sMap$codebook
    cnames <- colnames(codebook)
    if(is.null(cnames)){
        cnames <- seq(1,ncol(codebook))
    }

	if(is.null(zlim)){
		zlim <- c(floor(min(codebook,na.rm=T)*10)/10, ceiling(max(codebook,na.rm=T)*10)/10)
		
		if(zlim[1]==zlim[2]){
			zlim[1] <- zlim[2]/2
		}
	}
	codebook[codebook<=zlim[1]] <- zlim[1]
	codebook[codebook>=zlim[2]] <- zlim[2]

    ## list of gp
    ls_gp <-lapply(1:length(cnames), function(i){
    	xSH(sMap, customised.comp=codebook[,i], colormap=colormap, zlim=zlim, border.color=border.color, barwidth=barwidth, barheight=barheight, nbin=nbin, legend.title=legend.title, legend.text.size=legend.text.size, legend.title.size=legend.title.size) + labs(title=cnames[i]) + theme(plot.title=element_text(hjust=0.5,size=title.size), plot.margin=unit(rep(0,4),rep("lines",4)))
    	#+ theme(legend.position="none")
    })
    names(ls_gp) <- cnames

    if(class(sReorder) != "sReorder"){
        stop("The funciton must apply to 'sReorder' object.\n")
    }

	##################
	# df: geom_subview
	##################
	df <- data.frame(x=sReorder$coord[,1], y=sReorder$coord[,2], stringsAsFactors=F)
	## append ls_gp
	df_max <- df[,1:2] %>% dplyr::arrange(-x,y)
	ind <- which(df[,1]==df_max[1,1] & df[,2]==df_max[1,2])
	df$ls_gp <- lapply(1:length(ls_gp), function(i){
		if(i!=ind){
			ls_gp[[i]] + theme(legend.position="none")
		}else{
			ls_gp[[i]] + theme(legend.position="none")
		}
	})
	
	##################
	## data: initialise the plot with 3 columns: x, y, fill
	##################
	if(1){
		xdim <- sReorder$xdim
		ydim <- sReorder$ydim
	}else{
		xdim <- max(df$x)
		ydim <- max(df$y)
	}
	data <- data.frame(x=c(0.5,xdim+0.5), y=c(0.5,ydim+0.5), stringsAsFactors=F)
	data$fill <- zlim
	
	x <- y <- fill <- NULL
	
	## ggplot
	#gp <- ggplot(data, aes(x,y)) + geom_blank()
	
	gp <- ggplot(data, aes(x,y)) + geom_point(aes(fill=fill), alpha=0) + scale_fill_gradientn(colors=xColormap(colormap)(nbin), limits=zlim, guide=guide_colorbar(title=legend.title,title.position="top",barwidth=barwidth,barheight=barheight,nbin=nbin,draw.ulim=FALSE,draw.llim=FALSE))
	
	gp <- gp + ggimage::geom_subview(data=df, aes(x=x, y=y, subview=ls_gp), width=embed.width, height=embed.width) + coord_fixed(ratio=1)
	
	gp <- gp + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), , panel.background=element_rect(fill=rgb(0.95,0.95,0.95,0.95)))
	gp <- gp + theme(legend.position="right",legend.box="vertical",legend.justification="bottom",axis.title.x=element_blank(), axis.title.y=element_blank())
	# gp + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
	
	## a sheet-shape rectangle grid
	gp <- gp + scale_x_continuous(breaks=1:xdim,labels=1:xdim,limits=c(0.5,xdim+0.5),expand=c(0,0))
	gp <- gp + scale_y_continuous(breaks=1:ydim,labels=1:ydim,limits=c(0.5,ydim+0.5),expand=c(0,0))
	
	gp$ls_gp <- ls_gp
	
	invisible(gp)
}
