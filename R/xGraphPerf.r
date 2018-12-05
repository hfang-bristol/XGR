#' Function to evaluate the graph node similarity prediction performance
#'
#' \code{xGraphPerf} is supposed to evaluate the graph node similarity prediction performance. It returns an object of class "pPerf".
#'
#' @param prediction a symmetric sparse matrix containing predictive scores
#' @param g an object of class "igraph" (or "graphNEL") for a graph. If either GSP or GSN is NULL, its adjacency matrix is used to define GSP (all edges) and the rest as GSN
#' @param GSP a vector containing Gold Standard Positives (GSP)
#' @param GSN a vector containing Gold Standard Negatives (GSN)
#' @param plot the way to plot performance curve. It can be 'none' for no curve returned, 'ROC' for ROC curve, and 'PR' for PR curve. 
#' @param highlight logical to indicate whether a dot is highlighted. It only works when plot is drawn. When true, the maximum accuracy highlighted in ROC curve, and the Fmax highlighted in PR curve. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' an object of class "pPerf", a list with following components:
#' \itemize{
#'  \item{\code{data}: a data frame with 8 columns, including 4 performance measures ('Accuracy', 'Precision', 'Recall' and 'Specificity'), 'name' (subjects), 'pred' (predictive scores), 'label' (1 for GSP and 0 for GSN), 'corrected' (corrected/transformed predictiv scores, always the higher the better)}
#'  \item{\code{auroc}: a scalar value for ROC AUC}
#'  \item{\code{fmax}: a scalar value for maximum F-measure}
#'  \item{\code{amax}: a scalar value for maximum accuracy}
#'  \item{\code{direction}: '+' (the higher score the better prediction) and '-' (the higher score the worse prediction)}
#'  \item{\code{gp}: a ggplot object (if plotted) or NULL}
#'  \item{\code{Pred_obj}: a ROCR prediction-class object (potentially used for calculating other performance measures)}
#' }
#' @note
#' AUC: the area under ROC
#' F-measure: the maximum of a harmonic mean between precision and recall along PR curve
#' @export
#' @include xGraphPerf.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' g <- make_graph("Zachary")
#' # Common Neighbors (CN)
#' prediction <- xGraphSim(g, measure="CN")
#' pPerf <- xGraphPerf(prediction, g)
#' gp <- xClassifyComp(pPerf, displayBy="PR")
#' }

xGraphPerf <- function(prediction, g=NULL, GSP=NULL, GSN=NULL, plot=c("none","ROC","PR"), highlight=F, verbose=TRUE)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    plot <- match.arg(plot)

    if(!(!is.null(g) | (!is.null(GSP) & !is.null(GSN)))){
    	warnings("'g' should be provided or both 'GSP' and 'GSN' is not null!")
    	return(NULL)
    }
	
	if(is.null(GSP) | is.null(GSN)){
		A  <- igraph::get.adjacency(g)
		df_A <- xSM2DF(A+1, verbose=F)
		df_A <- df_A[df_A[,1]<df_A[,2],]
		df_GSP <- df_A[df_A[,3]==2,]
		df_GSN <- df_A[df_A[,3]==1,]
		GSP <- paste0(df_GSP[,1],'--',df_GSP[,2])
		GSN <- paste0(df_GSN[,1],'--',df_GSN[,2])
	}
	
	df <- xSM2DF(prediction, verbose=F)
	df <- df[df[,1]<df[,2],]
	prediction <- data.frame(name=paste0(df[,1],'--',df[,2]), score=df[,3], stringsAsFactors=F)
	
	pPerf <- xClassifyPerf(prediction, GSP, GSN, plot=plot, highlight=highlight, verbose=verbose)

	invisible(pPerf)
}


