#' Function to visualise predictors
#'
#' \code{xClassifyPred} is supposed to visualise predictors from a 'sClass' object. It returns an object of class "ggplot".
#'
#' @param sClass an object of class "sClass"
#' @param displayBy which statistics will be used for displaying. It can be "importance" for predictor importance, "ROC" for AUC in ROC (shaped by direction), "Fmax" for F-max in Precision-Recall curve, and "Amax" for maximum accuracy
#' @param sort logical to indicate whether to sort methods according to performance. By default, it sets TRUE
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xClassifyPred}}
#' @include xClassifyPred.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xClassifyPred(sClass, displayBy="ROC")
#' gp <- xClassifyPred(sClass, displayBy="importance")
#' }

xClassifyPred <- function(sClass, displayBy=c("ROC","Fmax","importance","Amax"), sort=TRUE) 
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    displayBy <- match.arg(displayBy)

    if(class(sClass) != "sClass"){
    	warnings("The function must apply to a 'sClass' object.\n")
        return(NULL)
    }
    
    if(displayBy=='importance'){
    	df <- data.frame(Predictor=rownames(sClass$importance), Val=sClass$importance[,1], stringsAsFactors=FALSE)
    	xlab <- "Predictor importance"
    }else if(displayBy=='ROC'){
    	df <- data.frame(Predictor=rownames(sClass$performance)[-1], Val=sClass$performance[-1,'auroc'], direction=sClass$performance[-1,'direction'], stringsAsFactors=FALSE)
    	xlab <- "AUC (a measure of ROC)"
    }else if(displayBy=='Fmax'){
    	df <- data.frame(Predictor=rownames(sClass$performance)[-1], Val=sClass$performance[-1,'fmax'], stringsAsFactors=FALSE)
    	xlab <- "F-max (a measure of PR curve)"
    }else if(displayBy=='Amax'){
    	df <- data.frame(Predictor=rownames(sClass$performance)[-1], Val=sClass$performance[-1,'amax'], stringsAsFactors=FALSE)
    	xlab <- "Maximum accuracy (along with ROC)"
    }
    
    Predictor <- Val <- direction <- NULL
	if(sort){
		df <- df %>% dplyr::arrange(-Val)
	}
	df$Predictor <- factor(df$Predictor, levels=unique(df$Predictor))
	
	p <- ggplot(df, aes(x=Val, y=Predictor))
	p <- p + geom_point(fill='transparent',shape=21)
	p <- p  + theme_bw() + theme(legend.position="none")
	p <- p + xlab(xlab)
	p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	## put arrows on both axes
	p <- p + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	
	## x-axis position
    if(displayBy == "ROC"){
    	p <- p + geom_point(aes(shape=direction)) + scale_shape_manual(values=c(95,43))
    	p <- p + scale_x_continuous(position="top", limits=c(0.5,1))
    }else if(displayBy == "importance"){
    	p <- p + geom_vline(xintercept=0, color='black', linetype='dashed')
    	p <- p + scale_x_continuous(position="top")
    }else{
		p <- p + scale_x_continuous(position="top")
	}
	
	invisible(p)
}
