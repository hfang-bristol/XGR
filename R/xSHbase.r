#' Function to visualise bases partitioned from a supra-hexagonal grid using ggplot2
#'
#' \code{xSHbase} is supposed to visualise bases partitioned from a supra-hexagonal grid using ggplot2.
#'
#' @param sMap an object of class "sMap"
#' @param sBase an object of class "sBase". It can be an integer vector specifying clusters
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param border.color the border color for each hexagon
#' @param legend.title the title of the colorbar. By default, it is ''
#' @param legend.text.size the text size of the legend tick labelings. By default, it is 5
#' @param legend.title.size the text size of the legend titles. By default, it is 6
#' @return 
#' a ggplot2 object
#' @note none
#' @export
#' @seealso \code{\link{xSHbase}}
#' @include xSHbase.r
#' @examples
#' \dontrun{
#' library(XGR)
#' 
#' gp <- xSHbase(sMap, sBase)
#' gp + theme(legend.position="none")
#' 
#' ## advanced use
#' # steps
#' df <- sMap$polygon %>% dplyr::group_by(index) %>% dplyr::summarise(val=unique(stepCentroid))
#' sBase <- df$val
#' gp <- xSHbase(sMap, sBase, legend.title='Steps')
#' # further labelled by hits
#' df_coord <- data.frame(sMap$coord, hits=sMap$hits, stringsAsFactors=F)
#' gp + geom_text(data=df_coord, aes(x,y,label=hits))
#' 
#' # angles
#' df <- sMap$polygon %>% dplyr::group_by(index) %>% dplyr::summarise(val=unique(angleCentroid))
#' sBase <- ceiling(180*(df$val/3.14))
#' gp <- xSHbase(sMap, sBase, legend.title='Angles')
#' }

xSHbase <- function(sMap, sBase, colormap="spectral", border.color="white", legend.title="", legend.text.size=6, legend.title.size=8)
{

    if (class(sMap) != "sMap"){
        stop("The funciton must apply to 'sMap' object.\n")
    }
	
	gp <- NULL
	
	vec_labels <- NULL
	
    if(class(sBase) == "sBase"){
    	## vec_bases
        vec_bases <- sBase$bases
        
    	## vec_labels
        vec_labels <- rep("", length(vec_bases))
        vec_labels[sBase$seeds] <- as.character(seq(1,length(sBase$seeds)))
          
    }else if(is.vector(sBase)){
    	vec_bases <- sBase
    	
    }

	if(is.vector(vec_bases)){
		
		# my_colors
        tmp <- sort(unique(vec_bases))
        my_colors <- xColormap(colormap, data=tmp)
        names(my_colors) <- tmp
        
		# df_polygon
		df_polygon <- sMap$polygon
		df_polygon$base <- vec_bases[df_polygon$index]
	
		df_polygon$base <- factor(df_polygon$base, levels=names(my_colors))
	
		group <- base <- NULL
		x <- y <- index <- NULL
		
		# ggplot
		gp <- ggplot(data=df_polygon, aes(x,y)) + geom_polygon(aes(fill=base,group=index),color=border.color) + scale_fill_manual(values=my_colors) + guides(fill=guide_legend(title=legend.title,keywidth=0.8, keyheight=0.8))
	
		gp <- gp + coord_fixed(ratio=1) + theme_void() + theme(legend.position="bottom") + theme(legend.title=element_text(face="bold",color="black",size=legend.title.size),legend.text=element_text(face="bold",color="black",size=legend.text.size),legend.title.align=0.5) + theme(legend.background=element_rect(fill="transparent"))

		if(!is.null(vec_labels)){
			label <- NULL

			df_coord <- data.frame(sMap$coord, label=vec_labels, stringsAsFactors=F)
			
			gp <- gp + geom_text(data=df_coord, aes(x,y,label=label))
			#gp + theme(legend.position="none")
		}

	}

	invisible(gp)
}
