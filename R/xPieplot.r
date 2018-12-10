#' Function to visualise data frame using pie plots
#'
#' \code{xPieplot} is supposed to visualise data frame using pie plots. It returns an object of class "ggplot". 
#'
#' @param df a data frame
#' @param columns a vector containing column names of the input data frame. These columns are used to draw pie charts
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param pie.radius the radius of a pie. If NULL, it equals roughly 1/75
#' @param pie.color the border color of a pie
#' @param pie.color.alpha the 0-1 value specifying transparency of pie border colors
#' @param pie.thick the pie border thickness
#' @param legend.title the legend title
#' @param gp an existing ggplot object or NULL. It is used for overlapping
#' @return
#' a ggplot object.
#' @export
#' @seealso \code{\link{xPieplot}}
#' @include xPieplot.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xPieplot(df, columns=c('dGene','pGene','fGene','nGene','eGene','cGene'), legend.title='Seeds')
#' }

xPieplot <- function(df, columns, colormap="ggplot2", pie.radius=NULL, pie.color='transparent', pie.color.alpha=1, pie.thick=0.1, legend.title='', gp=NULL)
{
	
	if(class(df) != "data.frame" ){
		stop("The function must apply to a data frame.\n")
	}
	
	ind <- colnames(df) %in% c('x','y', columns)
	if(sum(ind) != (length(columns)+2)){
		stop("The input data frame must contain columns 'x', 'y' and those in columns provided.\n")
	}
	df_sub <- df[, ind]
	
	## remove duplicated rows
	df_sub <- df_sub[!duplicated(df_sub),]
	
	## remove rows with all zeros
	df_sub <- df_sub[apply(df_sub[,c(-1,-2)],1,sum)!=0, ]
	
	#######
	geom_pie <- function(mapping=NULL, data, columns, ...) {
		if (is.null(mapping)){
			mapping <- ggplot2::aes_(x=~x, y=~y)
		}
		mapping <- utils::modifyList(mapping, ggplot2::aes_(r0=0, fill=~type, amount=~value))

		if (!('r' %in% names(mapping))) {
			#xvar <- as.character(mapping)[["x"]]
			xvar <- "x"
			tmp <- base::diff(range(data[, xvar]))
			## makde sure range difference is not zero (otherwise 1)
			if(tmp==0){
				tmp <- 1
			}
			size <- tmp/75
			mapping <- utils::modifyList(mapping, aes_(r=size))
		}

		names(mapping)[match(c("x", "y"), names(mapping))] <- c("x0", "y0")

		df <- tidyr::gather(data, "type", "value", columns, factor_key=TRUE)
		#df$type <- factor(df$type, levels=columns)
		ggforce::geom_arc_bar(mapping, data=df, stat='pie', inherit.aes=FALSE, ...)
	}
	#######
	
	x <- y <- NULL
	
	if(all(class(gp) %in% c("gg","ggplot"))){
		gp <- gp
	}else{
		gp <- ggplot()	
	}
	
	if(is.null(pie.radius)){
		gp <- gp + geom_pie(data=df_sub, aes(x=x, y=y), columns=columns, show.legend=TRUE, color=pie.color, alpha=pie.color.alpha, size=pie.thick)
	}else{
		df_sub$pie.radius <- rep(pie.radius, nrow(df_sub))
		gp <- gp + geom_pie(data=df_sub, aes(x=x, y=y, r=pie.radius), columns=columns, show.legend=TRUE, color=pie.color, alpha=pie.color.alpha, size=pie.thick)
	}
	
	gp <- gp + coord_equal() + guides(fill=guide_legend(title=legend.title, title.position="top", keywidth=0.6, keyheight=0.6, override.aes=list(shape=19)))
	
	## ggplot2: Fix colors to factor levels
	my_colors <- xColormap(colormap)(length(columns))
	names(my_colors) <- columns
	gp <- gp + scale_fill_manual(values=my_colors)
	
    invisible(gp)
}
