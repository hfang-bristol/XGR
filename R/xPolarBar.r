#' Function to visualise a data frame using a polar barplot
#'
#' \code{xPolarBar} is supposed to visualise a data frame using a polar dotplot. It returns an object of class "ggplot".
#'
#' @param df a data frame with two columns ('name' and 'value')
#' @param colormap either NULL or color names ('spectral' by default) for bars according to the name column
#' @param size.name an integer specifying the text size for the name column. By default, it sets to 10
#' @param size.value an integer specifying the text size for the value column. By default, it sets to 3
#' @param parallel logical to indicate whether the label is parallel to polar coordinate. By default, it sets FALSE
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xPolarBar}}
#' @include xPolarBar.r
#' @examples
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
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
#' gp <- xPolarBar(df[1:20,], parallel=TRUE)
#' gp
#' }

xPolarBar <- function(df, colormap='spectral', size.name=10, size.value=3, parallel=FALSE, font.family="sans", signature=TRUE)
{
    
    if(class(df) == "data.frame"){
    	df <- df[,c(1:2)]
    	colnames(df) <- c("name","value")
    }else{
    	stop("The function must apply to a 'data frame' object.\n")
    }
	
	if(class(df$value)=='factor'){
		df$value <- as.numeric(as.character(df$value))
	}
	
	name <- value <- label <- NULL
	
	my_colors <- xColormap(colormap)(length(df$name))
	names(my_colors) <- df$name
	
	df <- as.data.frame(df %>% dplyr::arrange(-value))
	df$name <- factor(df$name, levels=df$name)
	df$label <- signif(df$value, digits=3)
	
	if(parallel){
		angle <- 90 - 360/length(df$name) * seq_along(df$name)
		angle[angle < -90] <- -180 + angle[angle < -90]
	}else{
		angle <- 0
	}
	
	## polar plot
	gp <- ggplot(df, aes(x=name, y=value)) 
	gp <- gp + geom_col(aes(fill=name)) + scale_fill_manual(values=my_colors)
	gp <- gp + theme_bw() + theme(legend.position="none",axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(size=size.name,color="black",angle=angle), panel.border=element_blank()) 
	gp <- gp + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())
	
	if(max(df$value) * min(df$value) < 0){
		gp <- gp + geom_text(aes(x=name,label=label),y=0,size=size.value,hjust=0.5)
	}else{
		if(min(df$value) > 0){
			gp <- gp + geom_text(aes(x=name,y=value,label=label),size=size.value,hjust=0.5)
		}else{
			gp <- gp + geom_text(aes(x=name,label=label),y=0,size=size.value,hjust=0.5)
		}
	}
	
	gp <- gp + coord_polar(start=0)
	
	## caption
    if(signature){
    	caption <- paste("Created by xPolarBar from XGR version", utils ::packageVersion("XGR"))
    	gp <- gp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
	
	## change font family to 'Arial'
	gp <- gp + theme(text=element_text(family=font.family))
		
	invisible(gp)
}
