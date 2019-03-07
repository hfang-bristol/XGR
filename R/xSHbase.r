#' Function to visualise bases partitioned from a supra-hexagonal grid using ggplot2
#'
#' \code{xSHbase} is supposed to visualise bases partitioned from a supra-hexagonal grid using ggplot2.
#'
#' @param sMap an object of class "sMap"
#' @param sBase an object of class "sBase". It can be an integer vector specifying clusters or NULL (map shown without coloring/legend)
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param border.color the border color for each hexagon
#' @param legend.title the title of the colorbar. By default, it is ''
#' @param legend.text.size the text size of the legend tick labelings
#' @param legend.title.size the text size of the legend titles
#' @param boundary logical to indicate whether the base boundary is added
#' @param boundary.color the boundary color for each base
#' @param boundary.type the type of boundary. It can be 'line' (polygon; likely problematic due to the ordering issue), 'point' (simply dots)
#' @return
#' a ggplot2 object
#' @note Illustrator: cut a shape of points into each segmented line: 1) select the shape; 2) Direct Select Tool (A); 3) hold Shift and select an anchor point; 4) select Cut path of selected anchor points in the Anchor Point toolbar (The icon with the segmented line and scissors).
#' @export
#' @seealso \code{\link{xSHbase}}
#' @include xSHbase.r
#' @examples
#' \dontrun{
#' library(XGR)
#' 
#' gp <- xSHbase(sMap)
#' 
#' gp <- xSHbase(sMap, sBase)
#' gp + theme(legend.position="none")
#' 
#' ## advanced use
#' # steps to the centroid
#' df <- unique(sMap$polygon[,c("index","stepCentroid")]) %>% dplyr::arrange(index)
#' gp <- xSHbase(sMap, df$stepCentroid, legend.title='Steps')
#' # further labelled by hits
#' df_coord <- data.frame(sMap$coord, hits=sMap$hits, index=1:length(sMap$hits), stringsAsFactors=F)
#' gp + geom_text(data=df_coord, aes(x,y,label=hits))
#' 
#' # angles to the centroid
#' df <- unique(sMap$polygon[,c("index","angleCentroid")]) %>% dplyr::arrange(index)
#' vec <- ceiling(180*(df$angleCentroid/3.14))
#' gp <- xSHbase(sMap, vec, legend.title='Angles')
#' 
#' ##################
#' # labelled by data
#' gp <- xSHbase(sMap)
#' # df_output
#' df_output <- sWriteData(sMap, data) %>% dplyr::group_by(Hexagon_index) %>% dplyr::summarise(IDs=paste(ID,collapse="\n"))
#' # df_coord
#' df_coord <- data.frame(sMap$coord, hits=sMap$hits, index=1:length(sMap$hits), stringsAsFactors=F)
#' # df_output_coord
#' ind <- match(df_output$Hexagon_index, df_coord$index)
#' df <- data.frame(df_output, df_coord[ind,])
#' gp + geom_point(data=df,aes(x,y,size=hits),alpha=0.1) + geom_text(data=df,aes(x,y,label=IDs),size=2)
#' gp + geom_point(data=df,aes(x,y,size=hits),alpha=0.1) + ggrepel::geom_text_repel(data=df,aes(x,y,label=IDs),size=1.5)
#' ##################
#' 
#' # add boundary
#' gp <- xSHbase(sMap, sBase, boundary=T, colormap="rainbow_hcl")
#' gp <- xSHbase(sMap, sBase, boundary=T, colormap="transparent-transparent")
#' gp + theme(legend.position="none")
#' 
#' # interactive (ggiraph)
#' df_polygon <- sMap$polygon
#' df_polygon$onclick <- 'alert(this.getAttribute("data-id"))'
#' df_polygon$onclick <- paste0('window.open("https://en.wikipedia.org/wiki/', df_polygon$stepCentroid,'")')
#' gg <- ggplot(df_polygon, aes(x, y)) + ggiraph::geom_polygon_interactive(aes(fill=index, group=index, tooltip=stepCentroid, data_id=stepCentroid, onclick=onclick)) + coord_fixed(ratio=1) + theme_void()
#' gr <- ggiraph::ggiraph(code=print(gg), width=0.6)
#' gr <- ggiraph::girafe_options(gr, ggiraph::opts_tooltip(use_fill=T), ggiraph::opts_hover(css="fill:orange"), ggiraph::opts_toolbar(position="topright", saveaspng=F))
#' htmlwidgets::saveWidget(gr, file="out.html", background="white", selfcontained=T)
#' }

xSHbase <- function(sMap, sBase=NULL, colormap="rainbow_hcl", border.color="grey", legend.title="", legend.text.size=6, legend.title.size=8, boundary=F, boundary.color="black", boundary.type=c("line","point"))
{

    boundary.type <- match.arg(boundary.type)

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
    	
    }else if(is.null(sBase)){
    	vec_bases <- 1:nrow(sMap$coord)
    	colormap <- NA
    }

	#########################################
	if(is.na(colormap)){
		colormap <- "transparent-transparent"
	}
	#########################################

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
	
		gp <- gp + coord_fixed(ratio=1) + theme_void() + theme(legend.position="right",legend.box="vertical") + theme(legend.title=element_text(face="bold",color="black",size=legend.title.size),legend.text=element_text(face="bold",color="black",size=legend.text.size),legend.title.align=0.5) + theme(legend.background=element_rect(fill="transparent"))
		
		if(colormap=="transparent-transparent"){
			gp <- gp + theme(legend.position="none")
		}
		
		if(!is.null(vec_labels)){
			label <- NULL

			df_coord <- data.frame(sMap$coord, label=vec_labels, stringsAsFactors=F)
			
			gp <- gp + geom_text(data=df_coord, aes(x,y,label=label))
			#gp + theme(legend.position="none")
		}
		
		###################################
		###################################
		if(boundary){
			x <- y <- angle <- base <- NULL
			
			######################
			# add base boundary 
			#gp_main <- xSHbase(sMap, sBase, colormap="lightblue-darkorange")
			######################
		
			## very important! because of precision
			df_polygon$y <- round(df_polygon$y, digits=2)
		
			## df1: insider points (>= 2 bases)
			df <- df_polygon
			df <- unique(df[,c('x','y','base')])
			df1 <- as.data.frame(df %>% dplyr::group_by(x,y) %>% dplyr::group_by(n=n(),add=T) %>% dplyr::filter(n!=1))[,c('x','y','base')]
			## df2: outsider points (<= 2 index)
			df <- df_polygon
			df2 <- as.data.frame(df %>% dplyr::group_by(x,y) %>% dplyr::group_by(n=n(),add=T) %>% dplyr::filter(n<=2))[,c('x','y','base')]
			## df: both
			df <- unique(rbind(df1, df2))
			#ggplot(data=df, aes(x,y)) + geom_point()
			#gp_main + geom_point(data=df, aes(x,y))		
			
			## add angle per base, relative to its centroid point
			ls_df <- split(x=df, f=df$base)
			ls_res <- lapply(1:length(ls_df), function(k){
				z <- ls_df[[k]]
				a <- c(1,0) # according to the angle with the x-axis
				z$angle <- sapply(1:nrow(z), function(i){
					b <- z[i,c('x','y')] - c(mean(z$x),mean(z$y)) # move the coordinate system
					res <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
					## convert into [0, 2*pi]
					if(b[2] < 0){
						res <- 2*pi - res
					}
					res
				})
				z <- z %>% dplyr::arrange(angle)
				#gp_main + geom_point(data=z, aes(x,y,color=angle)) + scale_color_gradientn(colors=xColormap('jet')(64))
				#gp_main + geom_polygon(data=z, aes(x,y,group=base),fill="transparent",color="black")
				
				#res <- z[grDevices::chull(x=z[,1], y=z[,2]),]
				#gp_main + geom_polygon(data=res, aes(x,y,group=base),fill="transparent",color="black")
				
				#res <- alphahull::ashape(z[,1:2], alpha=0.5)
				#plot(res)
				
			})
			df <- do.call(rbind, ls_res)
			#gp_main + geom_polygon(data=subset(df,base==2), aes(x,y,group=base),fill="transparent",color="black")
			
			if(boundary.type=="line"){
				gp <- gp + geom_polygon(data=df, aes(x,y,group=base), fill="transparent", color=boundary.color)			
			}else if(boundary.type=="point"){
				gp <- gp + geom_point(data=df, aes(x,y), fill="transparent", color=boundary.color)
			}
			
			######################	
		}
		
	}

	invisible(gp)
}
