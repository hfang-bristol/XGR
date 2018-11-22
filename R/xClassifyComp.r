#' Function to compare prediction performance results
#'
#' \code{xClassifyComp} is supposed to compare prediction performance results. It returns an object of class "ggplot".
#'
#' @param list_pPerf a list of "pPerf" objects or a "pPerf" object 
#' @param displayBy which performance will be used for comparison. It can be "ROC" for ROC curve (by default), "PR" for PR curve
#' @param type the type of plot to draw. It can be "bar" for bar plot (by default) or "curve" for curve plot
#' @param sort logical to indicate whether to sort methods according to performance. By default, it sets FALSE
#' @param detail logical to indicate whether to label performance and score direction (together with methods). By default, it sets TRUE. Only works for the curve
#' @param facet logical to indicate whether to facet/wrap a 1d of panels into 2d. By default, it sets FALSE. Only works for the curve
#' @return an object of class "ggplot" or NULL (if all input pPerf objects are NULL)
#' @note none
#' @export
#' @seealso \code{\link{xClassifyPerf}}
#' @include xClassifyComp.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xClassifyComp(list_pPerf, displayBy="ROC")
#' print(gp)
#' ## modify legend position
#' gp + theme(legend.position=c(0.75,0.25))
#' }

xClassifyComp <- function(list_pPerf, displayBy=c("ROC","PR"), type=c("bar","curve"), sort=FALSE, detail=TRUE, facet=FALSE)
{

    displayBy <- match.arg(displayBy)
    type <- match.arg(type)
    
    if(is.null(list_pPerf)){
        warnings("Input is NULL.\n")
        return(NULL)
    }
    
   	if(class(list_pPerf)=="pPerf"){
		list_pPerf <- list(list_pPerf)
	}else if(class(list_pPerf)=="list"){
		## Remove null elements in a list
		list_pPerf <- base::Filter(base::Negate(is.null), list_pPerf)
		if(length(list_pPerf)==0){
			warnings("All pPerf objects are NULL!")
			return(NULL)
		}
	}else{
		return(NULL)
	}
    
	## Combine into a data frame called 'df_data'
	list_names <- names(list_pPerf)
	if(is.null(list_names)){
		list_names <- paste('Method', 1:length(list_pPerf), sep=' ')
	}
	ls_df <- lapply(1:length(list_pPerf), function(i){
		df_data <- list_pPerf[[i]]$data[,c('Accuracy','Precision','Recall','Specificity')]
		
		auroc <- signif(list_pPerf[[i]]$auroc, digits=3)
		fmax <- signif(list_pPerf[[i]]$fmax, digits=3)
		amax <- signif(list_pPerf[[i]]$amax, digits=3)		
		direction <- list_pPerf[[i]]$direction
		method <- list_names[i]
		
		#label <- paste(method, ' (AUC=', auroc, '; Fmax=', fmax,')', sep='')
		if(displayBy=='ROC'){
			label <- method
			if(detail){
				label <- paste0(method, ' (AUC=', auroc, ') ',direction)
			}
		}else if(displayBy=='PR'){
			label <- method
			if(detail){
				label <- paste0(method, ' (Fmax=', fmax, ') ',direction)
			}
		}
		
		data.frame(df_data, methods=method, auroc=auroc, fmax=fmax, amax=amax, direction=direction, labels=label, stringsAsFactors=FALSE)
	})
	df_data <- do.call(rbind, ls_df)

	## Method factor
	df_data$methods <- factor(df_data$methods, levels=list_names)
	
	if(type=='curve'){
		## draw curves
		Recall <- Precision <- Specificity <- methods <- labels <- auroc <- fmax <- amax <- direction <- NULL
		if(displayBy=='ROC'){
			## sort by: auroc
			if(sort){
				df_data <- df_data[with(df_data,order(-auroc)), ]
				## define levels
				if(detail){
					df_data$labels <- factor(df_data$labels, levels=unique(df_data$labels))
				}else{
					df_data$methods <- factor(df_data$methods, levels=unique(df_data$methods))
				}
			}
		
			## ggplot
			p <- ggplot(df_data, aes(x=1-Specificity,y=Recall))
		
			if(detail){
				p <- p + geom_line(aes(colour=factor(labels)))
			}else{
				p <- p + geom_line(aes(colour=factor(methods)))
			}
		
			p <- p + ylab("Sensitivity = TP/(TP+FN)") + xlab("1-Specificity = FP/(FP+TN)") + ylim(0,max(df_data$Recall)) + xlim(0,max(1-df_data$Specificity))
		
		}else if(displayBy=='PR'){
			## sort by: fmax
			if(sort){
				df_data <- df_data[with(df_data, order(-fmax)), ]
				## define levels
				if(detail){
					df_data$labels <- factor(df_data$labels, levels=unique(df_data$labels))
				}else{
					df_data$methods <- factor(df_data$methods, levels=unique(df_data$methods))
				}
			}
			## ggplot
			p <- ggplot(df_data, aes(x=Recall,y=Precision)) 

			if(detail){
				p <- p + geom_line(aes(colour=factor(labels)))
			}else{
				p <- p + geom_line(aes(colour=factor(methods)))
			}

			p <- p + ylab("Precision = TP/(TP+FP)") + xlab("Recall = TP/(TP+FN)") + ylim(0,max(df_data$Precision)) + xlim(0,max(df_data$Recall))
		
		}
	
		p <- p + theme_bw() + theme(axis.title.y=element_text(size=12,color="black"), axis.title.x=element_text(size=12,color="black"))
	
		if(facet){
			if(detail){
				p <- p + facet_wrap(~labels)
			}else{
				p <- p + facet_wrap(~methods)
			}
		
			## strip
			p <- p + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(face="italic"))
		
			p <- p + theme(legend.position="none", legend.title=element_blank())
		}else{
			p <- p + theme(legend.title=element_blank(), legend.key=element_rect(colour="transparent"))
		
			#p + theme(legend.position=c(0.75,0.25))
		}
	
	}else if(type=='bar'){
		
		df <- unique(df_data[,c("methods","auroc","fmax","amax","direction")])
		
		## draw bar
		methods <- auroc <- fmax <- NULL
		if(displayBy=='ROC'){
			## sort by: auroc
			if(sort){
				df <- df[with(df, order(-auroc)), ]
			}
			
			## define levels
			df$methods <- factor(df$methods, levels=rev(unique(df$methods)))
				
			## ggplot
			p <- ggplot(df, aes(x=methods,y=auroc))
			p <- p + geom_col(aes(fill=factor(methods)))
			p <- p + ylab("AUC\n(a measure of ROC)")
			p <- p + geom_text(aes(label=auroc), hjust=1)
			
			if(0){
				ylim_low <- ifelse(min(df$auroc)>0.5, 0.5, min(df$auroc))
				p <- p + coord_cartesian(ylim=c(ylim_low,1))
			}
			
		}else if(displayBy=='PR'){
			## sort by: fmax
			if(sort){
				df <- df[with(df, order(-fmax)), ]
			}
		
			## define levels
			df$methods <- factor(df$methods, levels=rev(unique(df$methods)))
		
			## ggplot
			p <- ggplot(df, aes(x=methods,y=fmax))
			p <- p + geom_col(aes(fill=factor(labels)))
			p <- p + ylab("F-max\n(a measure of Precision-Recall curve)")
			p <- p + geom_text(aes(label=fmax), hjust=1)
		}
		
		p <- p + theme_bw() + theme(legend.position="none",axis.title.y=element_blank(), axis.text.y=element_text(size=12,color="black"), axis.title.x=element_text(size=14,color="black")) + coord_flip()
	
	}
	
	p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	## put arrows on both axes
	p <- p + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))

	## y-axis position
	p <- p + scale_y_continuous(position="right")
	
	invisible(p)
}
