#' Function to draw heatmap using ggplot2
#'
#' \code{xHeatmap} is supposed to draw heatmap using ggplot2.
#'
#' @param data a data frame/matrix for coloring. The coloring can be continuous (numeric matrix) or discrete (factor matrix)
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
#' @param ... additional graphic parameters for supraHex::visTreeBootstrap
#' @return 
#' a gpplot2 object
#' @note none
#' @export
#' @seealso \code{\link{xHeatmap}}
#' @include xHeatmap.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' data(mtcars)
#' gp <- xHeatmap(mtcars, reorder="none", colormap='jet.top', x.rotate=45, shape=19, size=3, x.text.size=8,y.text.size=8, legend.title='mtcars')
#' gp + theme(legend.position="bottom",legend.direction="horizontal") + guides(color=guide_colorbar(title="mtcars",title.position="top",barwidth=5,barheight=0.3))
#' gp + theme(legend.position="bottom",legend.direction="horizontal") + guides(color=guide_legend(title="mtcars",title.position="top",barwidth=5,barheight=0.3))
#' gp + geom_text(aes(x, y, label=val),size=1.8,color='black',fontface='bold',na.rm=TRUE,angle=45)
#' }

xHeatmap <- function(data, reorder=c("none","row","col","both"), colormap="spectral", ncolors=64, zlim=NULL, barwidth=0.3, barheight=NULL, nbin=64, legend.title='', x.rotate=60, x.text.size=6, y.text.size=6, legend.text.size=4, legend.title.size=6, shape=19, size=2, plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt"), font.family="sans", na.color='grey80', data.label=NULL, label.size=1, label.color="black", ...)
{

    reorder <- match.arg(reorder)
    
	if(is.null(rownames(data))){
		rownames(data) <- paste('R', 1:nrow(data), sep=' ')
	}
	if(is.null(colnames(data))){
		colnames(data) <- paste('C', 1:ncol(data), sep=' ')
	}
	mat_val <- data
	
	flag_factor <- FALSE
	if(class(unlist(mat_val))=='factor'){
		lvl <- levels(unlist(mat_val))
		indx <- sapply(mat_val, is.factor)
		mat_val[indx] <- sapply(mat_val[,indx], function(x) as.numeric(x))
		
		################
		flag_factor <- F
		################
	}
	
	ind_row <- 1:nrow(mat_val)
	if(ncol(mat_val)>1 & (reorder=="row" | reorder=="both")){
		mat <- mat_val
		colnames(mat) <- 1:ncol(mat)
		rownames(mat) <- 1:nrow(mat)
		####
		mat[is.na(mat)] <- 0
		####
		tree_bs <- visTreeBootstrap(mat, visTree=FALSE, ...)
		ind_row <- match(tree_bs$tip.label, rownames(mat))
	}
	ind_col <- 1:ncol(mat_val)
	if(nrow(mat_val)>1 & (reorder=="col" | reorder=="both")){
		mat <- mat_val
		colnames(mat) <- 1:ncol(mat)
		rownames(mat) <- 1:nrow(mat)
		####
		mat[is.na(mat)] <- 0
		####
		tree_bs <- visTreeBootstrap(t(mat), visTree=FALSE, ...)
		ind_col <- match(tree_bs$tip.label, colnames(mat))
	}
	
	################
	if(0){
		mat_val <- mat_val[ind_row, ind_col]
	}else{
		mat_tmp <- as.matrix(mat_val[ind_row, ind_col], ncol=length(ind_col))
		rownames(mat_tmp) <- rownames(mat_val)[ind_row]
		colnames(mat_tmp) <- colnames(mat_val)[ind_col]
		mat_val <- mat_tmp
	}
	################
		
	if(class(mat_val)=='matrix'){
		mat_val <- as.data.frame(mat_val)
	}
	
	if(class(mat_val)=='data.frame'){
		
		if(is.null(zlim)){
			zlim <- c(ceiling(min(mat_val)*10)/10, ceiling(max(mat_val)*10)/10)
		}
		mat_val[mat_val<=zlim[1]] <- zlim[1]
		mat_val[mat_val>=zlim[2]] <- zlim[2]
		
		gene <- sample <- val <- NULL
		df <- mat_val %>% dplyr::mutate(gene=rownames(mat_val)) %>% tidyr::gather(sample, val, -gene)
		df$gene <- factor(df$gene, levels=rev(rownames(mat_val)))
		df$sample <- factor(df$sample, levels=colnames(mat_val))
		
		df <- df %>% dplyr::mutate(uid=paste(gene,sample,sep=':'))
		
		df$y <- as.numeric(df$gene)
		df$x <- as.numeric(df$sample)
		
		if(flag_factor){
			df$val <- factor(df$val, levels=unique(df$val))
		}
		
		gp <- ggplot(df, aes(x=sample, y=gene, color=val))
		#gp <- ggplot(df, aes(x=x, y=y, color=val))
		gp <- gp + geom_point(size=size, shape=shape)
		
		if(!flag_factor){
			gp <- gp + scale_colour_gradientn(colors=xColormap(colormap)(ncolors), limits=zlim, guide=guide_colorbar(title=legend.title,title.position="top",barwidth=barwidth,barheight=barheight,nbin=nbin,draw.ulim=FALSE,draw.llim=FALSE), na.value=na.color)
		}else{
			gp <- gp + scale_colour_manual(legend.title, values=xColormap(colormap)(length(lvl)), labels=lvl)
			if(is.null(barheight)){
				barheight <- barwidth
			}
			gp <- gp + theme(legend.key.width=unit(barwidth,'pt'), legend.key.height=unit(barheight,'pt')) + guides(col=guide_legend(ncol=1))
		}
		
		gp <- gp + theme_bw() + theme(legend.position="right", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold",color="black",size=x.text.size,angle=x.rotate,hjust=0), axis.text.y=element_text(face="bold",color="black",size=y.text.size,angle=0), panel.background=element_rect(fill="transparent")) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + theme(plot.margin=plot.margin) + theme(legend.title=element_text(face="bold",color="black",size=legend.title.size),legend.text=element_text(face="bold",color="black",size=legend.text.size),legend.title.align=0.5) + theme(legend.background=element_rect(fill="transparent"))
		gp <- gp + theme(text=element_text(family=font.family))
		gp_main <- gp + scale_x_discrete(position="top")
		
		#######################
		if(!is.null(data.label)){
			mat_label <- as.data.frame(data.label[ind_row, ind_col])
			gene <- sample <- val <- NULL
			df_label <- suppressWarnings(mat_label %>% dplyr::mutate(gene=rownames(mat_label)) %>% tidyr::gather(sample, val, -gene))
			df_label$gene <- factor(df_label$gene, levels=rev(rownames(mat_label)))
			df_label$sample <- factor(df_label$sample, levels=colnames(mat_label))
			df_label <- df_label %>% dplyr::mutate(uid=paste(gene,sample,sep=':'))
			df_label$y <- as.numeric(df_label$gene)
			df_label$x <- as.numeric(df_label$sample)
			gp_main <- gp_main + geom_text(data=df_label, aes(x=sample, y=gene, label=val),hjust=0.5,vjust=0.5,size=label.size,color=label.color)
		}
		#######################
		
		invisible(gp_main)
		
	}else{
		return(NULL)	
	}
	
}
