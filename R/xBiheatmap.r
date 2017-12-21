#' Function to visualise bipartitle graph communities using heatmap
#'
#' \code{xBiheatmap} is supposed to visualise bipartitle graph communities using heatmap.
#'
#' @param g an object of class "igraph" for a bipartitel graph with node attributes 'type', 'community' and 'contribution'
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
#' @param intercept.color intercept color
#' @param intercept.size intercept size
#' @return 
#' a gpplot2 object
#' @export
#' @seealso \code{\link{xBigraph}}, \code{\link{xHeatmap}}
#' @include xBiheatmap.r
#' @examples
#' # 1) generate a random bipartite graph
#' set.seed(123)
#' g <- sample_bipartite(100, 50, p=0.1)
#' V(g)$name <- V(g)
#' 
#' \dontrun{
#' # 2) obtain its community
#' ig <- xBigraph(g)
#' gp <- xBiheatmap(ig)
#' }

xBiheatmap <- function(g, colormap="spectral", ncolors=64, zlim=NULL, barwidth=0.3, barheight=NULL, nbin=64, legend.title='', x.rotate=60, x.text.size=3, y.text.size=3, legend.text.size=4, legend.title.size=6, shape=19, size=0.5, plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt"), font.family="sans", na.color='transparent', intercept.color="grey95", intercept.size=0.3)
{
    
    if (class(g) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }else{
    	ig <- g
    }
	
	ls_res <- igraph::bipartite_mapping(ig)
	if(!(ls_res$res)){
		stop("The igraph object is not bipartitle.\n")
	}
	
	if(!all(c("type","community","contribution") %in% igraph::vertex_attr_names(ig))){
		stop("The igraph object must have vertex attributes 'type','community' and 'contribution'.\n")
		
	}else{
	
		message(sprintf("The igraph object has %d nodes and %d edges (%s) ...", vcount(ig), ecount(ig), as.character(Sys.time())), appendLF=T)
	
		adj <- as.matrix(xConverter(ig, from='igraph', to='dgCMatrix', verbose=FALSE))
		
		type <- community <- contribution <- name <- NULL
		df_nodes <- igraph::get.data.frame(ig,what="vertices")
		## sorted by type, community (1,2,3,...), -contribution
		df_nodes <- df_nodes %>% dplyr::arrange(type,community,-contribution)
		
		#############################
		## heatmap
		#############################
		## ynode.memb
		ynode.memb <- subset(df_nodes, type=='ynode')
		ynode.order <- ynode.memb$name
		adj <- adj[ynode.order, ]
		## xnode.memb
		xnode.memb <- subset(df_nodes, type=='xnode')
		xnode.order <- xnode.memb$name
		adj <- adj[, xnode.order]
		## rowsep, colsep
		rowsep <- cumsum(as.vector(table(ynode.memb$community)))
		colsep <- cumsum(as.vector(table(xnode.memb$community)))
		if(0){
			labCol <- as.character(sort(xnode.memb$community))
			labCol[duplicated(labCol)] <- ""
			labRow <- as.character(sort(ynode.memb$community))
			labRow[duplicated(labRow)] <- ""
		}
		data <- adj
		data[data==0] <- NA
		
		gp <- xHeatmap(data, colormap=colormap, ncolors=ncolors, zlim=zlim, barwidth=barwidth, barheight=barheight, nbin=nbin, legend.title=legend.title, x.rotate=x.rotate, x.text.size=x.text.size, y.text.size=y.text.size, legend.text.size=legend.text.size, legend.title.size=legend.title.size, shape=shape, size=size, plot.margin=plot.margin, font.family=font.family, na.color=na.color, data.label=NULL)
		gp <- gp + geom_vline(xintercept=colsep+0.5,color=intercept.color,size=intercept.size) + geom_hline(yintercept=nrow(data)-rowsep+0.5,color=intercept.color,size=intercept.size) + theme(axis.ticks=element_line(size=0.25),axis.ticks.length=unit(0.05,"cm"))
		
    }

    invisible(gp)
}


