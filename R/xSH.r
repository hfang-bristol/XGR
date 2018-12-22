#' Function to draw multiple component planes of a supra-hexagonal grid using ggplot2
#'
#' \code{xSH} is supposed to draw multiple component planes of a supra-hexagonal grid using ggplot2.
#'
#' @param sMap an object of class "sMap"
#' @param which.components an integer vector specifying which compopnets will be visualised. By default, it is NULL meaning all components will be visualised
#' @param customised.comp the customised codebook matrix. It has a high priority over the built one. It can be a vector
#' @param ncolumns an integer specifying the number of columns for a rectangle grid wherein the component planes are placed. By defaul, it is NULL (decided on according to the number of component planes that will be visualised)
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
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
#' @seealso \code{\link{xSH}}
#' @include xSH.r
#' @examples
#' \dontrun{
#' library(XGR)
#' 
#' gp <- xSH(sMap)
#' gp <- xSH(sMap, ncolumns=8)
#' gp <- xSH(sMap, which.components=1)
#' gp + theme(legend.position="none")
#' 
#' ## advanced use
#' gp <- xSH(sMap, customised.comp=sMap$codebook[,1:2])
#' 
#' ## list of gp
#' ls_gp <-lapply(1:5, function(i) xSH(sMap, which.components=i) + theme(legend.position="none"))
#' 
#' ## save into a gif file via magick
#' ls_gp <-lapply(1:ncol(sMap$codebook), function(i) xSH(sMap, which.components=i,zlim=c(0,5))+ggtitle(colnames(sMap$codebook)[i]))
#' img <- magick::image_graph(480, 480, res=96)
#' ls_gp; dev.off()
#' animation <- magick::image_animate(img, fps=2)
#' magick::image_write(animation, "xSH.gif")
#' 
#' ## save into a gif file via magick and tweenr
#' cnames <- colnames(sMap$codebook)
#' data <- lapply(1:length(cnames), function(i) data.frame(x=sMap$codebook[,i]))
#' k <- 8
#' data_tween <- tweenr::tween_states(c(data,data[1]), tweenlength=1, statelength=1, ease='cubic-in-out', nframes=length(cnames) * k)
#' ls_gp <- lapply(1:(length(cnames)*k), function(i) {
#'    p <- xSH(sMap, customised.comp=subset(data_tween,.frame==i)$x, zlim=c(0,5)) + ggtitle(cnames[ceiling(i/k)])
#' })
#' img <- magick::image_graph(480, 480, res=96)
#' ls_gp; dev.off()
#' animation <- magick::image_animate(img, fps=5)
#' magick::image_write(animation, "xSH_advanced.gif")
#' }

xSH <- function(sMap, which.components=NULL, customised.comp=NULL, ncolumns=NULL, colormap="spectral", ncolors=64, zlim=NULL, border.color="transparent", barwidth=0.4, barheight=NULL, nbin=64, legend.title='', legend.text.size=6, legend.title.size=8)
{

    if (class(sMap) != "sMap"){
        stop("The funciton must apply to 'sMap' object.\n")
    }
    codebook <- sMap$codebook
    cnames <- colnames(codebook)
    if(is.null(cnames)){
        cnames <- seq(1,ncol(codebook))
    }
   	
   	# customised comp
   	if(!is.null(customised.comp)){
   		if(is.vector(customised.comp)){
   			codebook <- matrix(customised.comp, ncol=1)
   		}else{
   			codebook <- customised.comp
   		}
   		
   	}else{
		if(all(!is.null(which.components))){
			which.components <- as.integer(which.components)
			if(all(which.components>=1 & which.components<=length(cnames))){
				codebook <- matrix(codebook[,which.components], ncol=length(which.components))
				cnames <- cnames[which.components]
				colnames(codebook) <- cnames
			}
		}
   	
   	}
   
	if(is.null(zlim)){
		zlim <- c(floor(min(codebook,na.rm=T)*10)/10, ceiling(max(codebook,na.rm=T)*10)/10)
			
		if(zlim[1]==zlim[2]){
			zlim[1] <- zlim[2]/2
		}
	}
	codebook[codebook<=zlim[1]] <- zlim[1]
	codebook[codebook>=zlim[2]] <- zlim[2]
	
	# df_polygon
	df_polygon <- sMap$polygon
	ls_df <- lapply(1:ncol(codebook),function(i){
		df <- df_polygon
		df$component <- colnames(codebook)[i]
		df$value <- codebook[df_polygon$index, i]
		df
	})
	df <- do.call(rbind, ls_df)
	
	group <- value <- component <- NULL
	x <- y <- index <- NULL
	
	# ggplot
	gp <- ggplot(data=df, aes(x,y)) + geom_polygon(aes(fill=value,group=index),color=border.color)
	
	gp <- gp + scale_fill_gradientn(colors=xColormap(colormap)(ncolors), limits=zlim, guide=guide_colorbar(title=legend.title,title.position="top",barwidth=barwidth,barheight=barheight,nbin=nbin,draw.ulim=FALSE,draw.llim=FALSE))
	
	gp <- gp + coord_fixed(ratio=1) + theme_void() + theme(legend.position="right") + theme(legend.title=element_text(face="bold",color="black",size=legend.title.size),legend.text=element_text(face="bold",color="black",size=legend.text.size),legend.title.align=0.5) + theme(legend.background=element_rect(fill="transparent"))
	
	# facet_wrap: partitions a plot into a matrix of panels
	if(ncol(codebook)>1){
		if(is.null(ncolumns)){
			ncolumns <- ceiling(sqrt(ncol(codebook)))
		}
		gp <- gp + facet_wrap(~component, ncol=ncolumns)
	}
	
	invisible(gp)
}
