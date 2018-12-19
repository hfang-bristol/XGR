#' Function to draw multiple component planes reorded within a tree using ggplot2
#'
#' \code{xSHtree} is supposed to draw multiple component planes reorded within a tree using ggplot2.
#'
#' @param sMap an object of class "sMap"
#' @param phylo an object of class "phylo" or NULL
#' @param layout the visual layout. It can be "rectangular" (via ggtree) or circular layout via ggraph including "fan" and "partition"
#' @param tip.size the tip/title size
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
#' @seealso \code{\link{xSHtree}}
#' @include xSHtree.r
#' @examples
#' \dontrun{
#' library(XGR)
#' 
#' gp <- xSHtree(sMap, layout="rectangular", width=1)
#' gp
#' ls_gp <- gp$ls_gp
#' 
#' # advanced use
#' tree <- visTreeBootstrap(t(sMap$codebook))
#' gp <- xSHtree(sMap, tree, layout="rectangular", embed.width=1)
#' gp <- xSHtree(sMap, tree, layout="fan", embed.width=0.4)
#' gp <- xSHtree(sMap, tree, layout="partition", embed.width=0.7)
#' }

xSHtree <- function(sMap, phylo=NULL, layout=c("rectangular","fan","partition"), tip.size=10, embed.width=1, colormap="spectral", zlim=NULL, border.color="transparent", barwidth=0.4, barheight=NULL, nbin=64, legend.title='', legend.text.size=6, legend.title.size=8)
{

    if(class(sMap) != "sMap"){
        stop("The funciton must apply to 'sMap' object.\n")
    }
    
    codebook <- sMap$codebook
    cnames <- colnames(codebook)
    if(is.null(cnames)){
        cnames <- seq(1,ncol(codebook))
    }
    
    if(class(phylo)=="phylo"){
    	tree <- phylo
    }else{
    	tree <- visTreeBootstrap(t(codebook))
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
    	xSH(sMap, customised.comp=codebook[,i], colormap=colormap, zlim=zlim, border.color=border.color, barwidth=barwidth, barheight=barheight, nbin=nbin, legend.title=legend.title, legend.text.size=legend.text.size, legend.title.size=legend.title.size) + labs(title=cnames[i]) + theme(plot.title=element_text(hjust=0.5,size=tip.size), plot.margin=unit(rep(0,4),rep("lines",4)))
    	#+ theme(legend.position="none")
    })
    names(ls_gp) <- cnames
	
	if(layout=='rectangular'){
		
		x <- y <- fill <- isTip <- bootstrap <- NULL
		
		gp <- xGT(tree, tip.size=tip.size/3)
		
		##################
		## data: initialise the plot with 3 columns: x, y, fill
		##################
		data <- data.frame(x=c(0,0), y=c(max(gp$data$y),max(gp$data$y)), stringsAsFactors=F)
		data$fill <- zlim
		
		gp <- gp + geom_point(data=data, aes(x,y,fill=fill), alpha=0) + scale_fill_gradientn(colors=xColormap(colormap)(nbin), limits=zlim, guide=guide_colorbar(title=legend.title,title.position="top",barwidth=barwidth,barheight=barheight,nbin=nbin,draw.ulim=FALSE,draw.llim=FALSE)) + theme(legend.position="left",legend.box="vertical")
		
		##################
		# df: geom_subview
		##################
		df_tips <- gp$data %>% dplyr::filter(isTip) %>% dplyr::arrange(-y)
		ind <- match(names(ls_gp), df_tips$label)
		df <- df_tips[ind[!is.na(ind)],]
		df$ls_gp <- lapply(ls_gp[!is.na(ind)], function(gp) gp + theme(legend.position="none", plot.title=element_blank()))
		#df$ls_gp <- lapply(ls_gp[!is.na(ind)], function(gp) gp + theme(legend.position="none", plot.title=element_text(size=5)))
		df$x <- max(df$x) + max(1,max(df$x)*0.1)

		gp <- gp + ggimage::geom_subview(data=df, aes(x=x, y=y, subview=ls_gp), width=embed.width, height=embed.width) + scale_x_continuous(position="bottom", limits=c(0,max(df$x)))
		
		##################
		## internal node labellings
		##################
		data <- as.data.frame(gp$data %>% dplyr::filter(!isTip))
		data$bootstrap <- as.numeric(data$label)
		if(is.numeric(data$bootstrap)){
			gp <- gp + geom_point(data=data, aes(x=x,y=y,size=bootstrap), shape=21, color='grey', fill="white") + geom_text(data=data, aes(label=bootstrap,size=bootstrap/3), show.legend=F)
		}
	
	}else if(layout %in% c("fan","partition")){
	
		x <- y <- isTip <- label <- leaf <- ..index.. <- depth <- NULL
		
		## if internal node labels are duplicated, remove all (and add by as.igraph)
		if(any(duplicated(tree$node.label))){
			tree$node.label<-NULL
		}
		ig <- as.igraph(tree, directed=T, use.labels=T)
		
		# append 'tipid': NA for internal node, sequential order for tips
		ind <- match(V(ig)$name, tree$tip.label)
		V(ig)$tipid <- NA
		V(ig)$tipid[!is.na(ind)] <- 1:sum(!is.na(ind))
		
		if(layout=="fan"){
			g <- ig
			
			#####################################
			name <- angle <- hjust <- NULL
			## orientation: angle and hjust			
			V(g)$angle <- 90 - 360 * V(g)$tipid / sum(!is.na(V(g)$tipid))
			V(g)$hjust <- ifelse(V(g)$angle < -90, 1, 0)
			V(g)$angle <- ifelse(V(g)$angle < -90, V(g)$angle+180, V(g)$angle)
			#####################################

			gp <- ggraph::ggraph(g, layout='dendrogram', circular=TRUE)
			#gp <- gp + ggraph::geom_edge_diagonal(aes(color=1 - ..index..), color='black')
			gp <- gp + ggraph::geom_edge_diagonal(color='grey')
			gp <- gp + ggraph::geom_node_point(aes(filter=!leaf), color='grey')
			gp <- gp + ggraph::geom_node_text(aes(x=x*1.1, y=y*1.1, filter=leaf, label=name, angle=angle, hjust=hjust),show.legend=F, size=2, alpha=1)
			gp <- gp + theme_void() + coord_fixed() + theme(legend.position="none",plot.margin=unit(c(1,1,1,1),"cm")) + expand_limits(x=c(-1.2,1.2), y=c(-1.2,1.2))
			
		}else if(layout=="partition"){
			
			gp <- ggraph::ggraph(ig, layout='partition', circular=TRUE)
			#gp <- gp + ggraph::geom_edge_diagonal(aes(alpha=1-..index..), color='black')
			gp <- gp + ggraph::geom_edge_diagonal(color='grey')
			gp <- gp + ggraph::geom_node_point(aes(filter=!leaf,alpha=1-depth), color='black')
			#gp <- gp + ggraph::geom_node_text(aes(filter=leaf, label=name))
			gp <- gp + theme_void() + coord_fixed() + theme(legend.position="none",plot.margin=unit(c(1,1,1,1),"cm"))
			
		}
		
		##################
		## data: initialise the plot with 3 columns: x, y, fill
		##################
		data <- data.frame(x=c(0,0), y=c(max(gp$data$y),max(gp$data$y)), stringsAsFactors=F)
		data$fill <- zlim
		
		gp <- gp + geom_point(data=data, aes(x,y,fill=fill), alpha=0) + scale_fill_gradientn(colors=xColormap(colormap)(nbin), limits=zlim, guide=guide_colorbar(title=legend.title,title.position="top",barwidth=barwidth,barheight=barheight,nbin=nbin,draw.ulim=FALSE,draw.llim=FALSE)) + theme(legend.position="right",legend.box="vertical",legend.justification="top")
		gp <- gp + guides(alpha="none")
		
		##################
		# df: geom_subview
		##################
		df_tips <- subset(gp$data, leaf==T)
		ind <- match(names(ls_gp), df_tips$name)
		df <- df_tips[ind[!is.na(ind)],]
		df$ls_gp <- lapply(ls_gp[!is.na(ind)], function(gp){
			gp + theme(legend.position="none")
		})
		
		gp <- gp + ggimage::geom_subview(data=df, aes(x=x, y=y, subview=ls_gp), width=embed.width, height=embed.width)
	
	}
	
	gp$ls_gp <- ls_gp
	
	invisible(gp)
}
