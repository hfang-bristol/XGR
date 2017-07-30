#' Function to draw a volcano plot
#'
#' \code{xVolcano} is supposed to draw a volcano plot
#'
#' @param data a data frame
#' @param column.lfc a character specifying 'lfc' column (log2-transformed fold change)
#' @param column.fdr a character specifying 'fdr' column
#' @param cutoff.lfc a numeric defining 'lfc' cutoff. By default, it is 1 (at least 2-fold changes)
#' @param cutoff.fdr a numeric defining 'fdr' cutoff. By default, it is 0.05
#' @param colors a 4-element vector for color-coded points. By default, it is c("#EEEEEE","darkgrey","pink","red")
#' @param column.label a character specifying 'label' column
#' @param top an integer specifying the number of the top points for labellings
#' @param top.direction the direction (up- and down-regulated) of the top points. It can be one of 'both' (up- and down-regulated), 'up' (up-regulated only) and 'down' (down-regulated only)
#' @param label.size the label size
#' @param label.color the label color
#' @param label.alpha the 0-1 value specifying transparency of labelling
#' @param label.padding the padding around the labeled
#' @param label.arrow the arrow pointing to the labeled
#' @param label.force the repelling force between overlapping labels
#' @param xlim the limits in the x-axis
#' @param ylim the limits in the y-axis
#' @param y.scale how to transform the y scale. It can be "normal" for no transformation, and "log" for log-based transformation
#' @param xlab the x labelling. By default, it is expression(log[2]("fold change"))
#' @param ylab the y labelling. By default, it is expression(-log[10]("FDR"))
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE
#' @return
#' a ggplot object
#' @note none
#' @export
#' @seealso \code{\link{xVolcano}}
#' @include xVolcano.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' }

xVolcano <- function(data, column.lfc='lfc', column.fdr='fdr', cutoff.lfc=1, cutoff.fdr=5e-2, colors=c("#EEEEEE","darkgrey","pink","red"), column.label=NULL, top=10, top.direction=c('both','up','down'), label.size=2, label.color='black', label.alpha=0.8, label.padding=0.5, label.arrow=0.01, label.force=0.5, xlim=NULL, ylim=NULL, y.scale=c("normal","log"), xlab=expression(log[2]("fold change")), ylab=expression(-log[10]("FDR")), font.family="sans", signature=TRUE)
{
	
	top.direction <- match.arg(top.direction)
	y.scale <- match.arg(y.scale)
	
	if(all(c(column.lfc, column.fdr) %in% colnames(data))){
		df <- data.frame(LFC=data[,column.lfc], FDR=data[,column.fdr], stringsAsFactors=FALSE)
	}else{
		return(NULL)
	}
	
	if(!is.null(top)){
		if(all(column.label %in% colnames(data))){
			df$label <- data[,column.label]
		}else{
			top <- NULL
		}
	}
	
	LFC <- FDR <- label <- flag <- NULL
	
	## only keep finite values
	df <- df %>% dplyr::filter(is.finite(LFC), is.finite(FDR))

	## add a column 'flag'
	df <- df %>% dplyr::mutate(flag=ifelse(FDR>=cutoff.fdr & abs(LFC)>=cutoff.lfc, 2,
						ifelse(FDR<cutoff.fdr & abs(LFC)<cutoff.lfc, 3, 
						ifelse(FDR<cutoff.fdr & abs(LFC)>=cutoff.lfc, 4, 1)
						)))
	df$flag <- factor(df$flag, levels=1:4)
	
	## plot
	gp <- ggplot(df, aes(x=LFC, y=-log10(FDR)))
	gp <- gp + geom_point(aes(color=flag), size=0.9, na.rm=TRUE)
	names(colors) <- 1:4
	gp <- gp + scale_colour_manual(values=colors)
	gp <- gp + theme_classic() + xlab(xlab) + ylab(ylab) 
	#gp <- gp + geom_vline(xintercept=c(-cutoff.lfc,cutoff.lfc),colour="black") + geom_hline(yintercept=-log10(cutoff.fdr), colour="black")
	
	gp <- gp + theme(legend.position="none",legend.title=element_blank(),legend.key=element_rect(colour="transparent"), axis.title.y=element_text(size=12), axis.title.x=element_text(size=12))
	if(!is.null(xlim)){
		gp <- gp + xlim(xlim)
	}
	if(!is.null(ylim)){
		gp <- gp + ylim(ylim)
	}
	
	## y scale
    if(y.scale=="log"){
    	gp <- gp + scale_y_continuous(trans="log1p", limits=ylim)
    	#gp <- gp + annotation_logticks(sides='l')
    }
    

	if(!is.null(top)){
	
		## up
		if(top.direction %in% c('both','up')){
			df_sub <- df %>% dplyr::filter(FDR<cutoff.fdr, LFC>=cutoff.lfc) %>% dplyr::arrange(FDR, desc(LFC))
			df_sub <- df_sub[1:min(top,nrow(df_sub)), ]	
			gp <- gp + ggrepel::geom_text_repel(data=df_sub, aes(x=LFC,y=-log10(FDR),label=label), size=label.size, color=label.color, fontface="bold", alpha=label.alpha, box.padding=unit(0.5,"lines"), point.padding=unit(label.padding,"lines"), segment.alpha=0.5, segment.color="grey50", segment.size=0.5, arrow=arrow(length=unit(label.arrow,'npc')), force=label.force)
		}
		
		## down
		if(top.direction %in% c('both','down')){
			df_sub <- df %>% dplyr::filter(FDR<cutoff.fdr, LFC<=-cutoff.lfc) %>% dplyr::arrange(FDR, LFC)
			df_sub <- df_sub[1:min(top,nrow(df_sub)), ]	
			gp <- gp + ggrepel::geom_text_repel(data=df_sub, aes(x=LFC,y=-log10(FDR),label=label), size=label.size, color=label.color, fontface="bold", alpha=label.alpha, box.padding=unit(0.5,"lines"), point.padding=unit(label.padding,"lines"), segment.alpha=0.5, segment.color="grey50", segment.size=0.5, arrow=arrow(length=unit(label.arrow,'npc')), force=label.force)
		}
	}

	## caption
    if(signature){
    	caption <- paste("Created by xVolcano from XGR version", utils::packageVersion("XGR"))
    	gp <- gp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
	
	## change font family to 'Arial'
	gp <- gp + theme(text=element_text(family=font.family))
	
	## put arrows on x-axis
	gp <- gp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"),type="open")), axis.line.y=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"),type="open")))
	
	return(gp)
}
