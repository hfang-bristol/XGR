#' Function to visualise a data frame using a polar dotplot
#'
#' \code{xPolarDot} is supposed to visualise a data frame using a polar dotplot. It returns an object of class "ggplot".
#'
#' @param df a data frame with two columns ('name' and 'value')
#' @param colormap either NULL or color names ('blue-yellow-red' by default) for points according to the value column
#' @param shape an integer specifying point shape
#' @param size an integer specifying the point size. By default, it sets to 2
#' @param parallel logical to indicate whether the label is parallel to polar coordinate. By default, it sets FALSE
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xPolarDot}}
#' @include xPolarDot.r
#' @examples
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' \dontrun{
#' # a) provide the seed nodes/genes with the weight info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get genes within 500kb away from AS GWAS lead SNPs
#' seeds.genes <- ImmunoBase$AS$genes_variants
#' ## seeds weighted according to distance away from lead SNPs
#' data <- 1- seeds.genes/500000
#'
#' # b) prepare a data frame
#' df <- data.frame(name=names(data), value=data, stringsAsFactors=FALSE)
#' 
#' # c) do correlation
#' gp <- xPolarDot(df[1:50,])
#' gp
#' }

xPolarDot <- function(df, colormap='blue-yellow-red', shape=19, size=2, parallel=FALSE, font.family="sans", signature=TRUE) 
{
    
    if(class(df) == "data.frame"){
    	df <- df[,c(1:2)]
    	yname <- colnames(df)[2]
    	colnames(df) <- c("name","value")
    }else{
    	stop("The function must apply to a 'data frame' object.\n")
    }
	
	if(class(df$value)=='factor'){
		df$value <- as.numeric(as.character(df$value))
	}
	
	name <- value <- rank <- NULL
	
	df$rank <- rank(-1*df$value,ties.method="first")
	df <- df %>% dplyr::arrange(rank)
	df$name <- factor(df$name, levels=df$name)
	
	color <- unlist(strsplit(colormap, "-"))
	
	if(parallel){
		angle <- 90 - 360/length(df$name) * seq_along(df$name)
		angle[angle < -90] <- -180 + angle[angle < -90]
	}else{
		angle <- 0
	}
	
	## polar plot
	gp <- ggplot(df, aes(x=name, y=value)) 
	gp <- gp + geom_point(aes(color=value), shape=shape, size=size)
	if(length(color)==2){
		gp <- gp + scale_color_gradient(low=color[1],high=color[2]) 
	}else if(length(color)==3){
		gp <- gp + scale_color_gradient2(low=color[1],mid=color[2],high=color[3]) 
	}else{
		gp <- gp + scale_color_gradient(low='black',high='black')
	}
	gp <- gp + theme_bw() + theme(legend.position="none") + labs(x='',y=yname)
	gp <- gp + theme(axis.title.y=element_text(size=12,color="black"), axis.text.y=element_text(size=8,color="black"), axis.title.x=element_text(size=12,color="black"), axis.text.x=element_text(size=6,color="black", angle=angle), panel.background=element_rect(fill="transparent"))
	#gp <- gp + geom_polygon(color='grey', fill=NA) 
	gp <- gp + coord_polar(start=0)
	
	## caption
    if(signature){
    	caption <- paste("Created by xPolarDot from XGR version", utils ::packageVersion("XGR"))
    	gp <- gp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
	
	## change font family to 'Arial'
	gp <- gp + theme(text=element_text(family=font.family))
		
	invisible(gp)
}
