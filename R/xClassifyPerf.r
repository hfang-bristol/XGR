#' Function to evaluate the prediction performance via ROC and Precision-Recall (PR) analysis
#'
#' \code{xClassifyPerf} is supposed to assess the prediction performance via Receiver Operating Characteristic (ROC) and Precision-Recall (PR) analysis. It requires three inputs: 1) Gold Standard Positive (GSP) targets; 2) Gold Standard Negative (GSN) targets; 3) prediction containing predictive scores on subjects. It returns an object of class "pPerf".
#'
#' @param prediction a data frame containing predictions along with predictive scores. It has two columns: 1st column for subjects, 2nd column for predictive scores on subjects
#' @param GSP a vector containing Gold Standard Positives (GSP)
#' @param GSN a vector containing Gold Standard Negatives (GSN)
#' @param rescale logical to indicate whether to linearly rescale predictive scores for GSP/GSN to the range [0,1]. By default, it sets to false
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
#' @include xClassifyPerf.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' pPerf <- xClassifyPerf(prediction, GSP, GSN)
#' }

xClassifyPerf <- function(prediction, GSP, GSN, rescale=F, plot=c("none","ROC","PR"), highlight=F, verbose=TRUE)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    plot <- match.arg(plot)

    res_ls <- split(x=as.numeric(prediction[,2]), f=prediction[,1])
    pred <- unlist(lapply(res_ls, base::max))
    ####
    pred <- pred[!is.na(pred)]
    ####
    if(verbose){
        message(sprintf("There are %d subjects in predictions (%s).", length(pred), as.character(Sys.time())), appendLF=TRUE)
    }
	
	## GSP
	gsp <- unique(GSP)
	### GSP but only predicted
    ind <- match(gsp, names(pred))
    gsp_predicted <- gsp[!is.na(ind)]
    if(verbose){
        message(sprintf("Of %d GSP, %d also predicted for evaluation (%s).", length(gsp), length(gsp_predicted), as.character(Sys.time())), appendLF=TRUE)
    }
    
	## GSN
	gsn <- unique(GSN)
	### GSN but only predicted
    ind <- match(gsn, names(pred))
    gsn_predicted <- gsn[!is.na(ind)]
    if(verbose){
        message(sprintf("Of %d GSN, %d also predicted for evaluation (%s).", length(gsn), length(gsn_predicted), as.character(Sys.time())), appendLF=TRUE)
    }
    
    ########################################
    # NULL if no GSP and/or GSN is predicted
    ########################################
    if(length(gsp_predicted)==0 | length(gsn_predicted)==0){
    	warnings("No GSP and/or GSN is predicted!")
    	return(NULL)
    }
    ########################################
    
    ######################################
	## prepare input for ROCR
	ind <- match(gsp_predicted, names(pred))
	gsp_pred_label <- data.frame(pred=pred[ind], label=1, stringsAsFactors=FALSE)
	ind <- match(gsn_predicted, names(pred))
	gsn_pred_label <- data.frame(pred=pred[ind], label=0, stringsAsFactors=FALSE)	
	pred_label <- rbind(gsp_pred_label, gsn_pred_label)
			
	## whether prediction scores are rescaled to the range [0,1]
	pred_label$scaled <- pred_label$pred
	if(rescale){
		x <- pred_label$scaled
		x_scaled <- (x - min(x)) / (max(x) - min(x))
		pred_label$scaled <- x_scaled
	}
	
	## ROCR
	direction <- '+'
	if(length(suppressWarnings(tryCatch(pred_obj <- ROCR::prediction(predictions=pred_label$scaled, labels=pred_label$label), error=function(e) e, warning=function(w) w)))==2){
		warnings("Failed to evaluate!")
		return(NULL)
	}else{
		## check whether AUC < 0.5. If so, reverse predictive scores
		## auc: Area under the ROC curve
		auroc <- unlist(ROCR::performance(pred_obj, measure="auc")@y.values)
		if(auroc < 0.5){
			direction <- '-'
			auroc <- 1 - auroc
			pred_obj <- ROCR::prediction(predictions=-pred_label$scaled, labels=pred_label$label)
		}
	}

	## acc: Accuracy=(TP+TN)/(P+N)
	res <- ROCR::performance(pred_obj, measure="acc")
	acc <- unlist(res@y.values)
    
	## ROC curves
	perf_roc <- ROCR::performance(pred_obj, measure="tpr", x.measure="fpr")
    tpr <- unlist(perf_roc@y.values)
    fpr <- unlist(perf_roc@x.values)
    
    ## PR curves
	perf_pr <- ROCR::performance(pred_obj, measure="prec", x.measure="rec")
    prec <- unlist(perf_pr@y.values)
    rec <- unlist(perf_pr@x.values)
    
    ## TPR vs pcfall (Prediction-conditioned fallout: FP/(TP+FP))
	#perf_tprfdr <- ROCR::performance(pred_obj, measure="tpr", x.measure="pcfall")    
    #tpr <- unlist(perf_tprfdr@y.values)
    #fdr <- unlist(perf_tprfdr@x.values)
    
    ###################
    # df_data
    ###################
	## df_APRS
	df_APRS <- data.frame(Accuracy=acc, Precision=prec, Recall=rec, Specificity=1-fpr, stringsAsFactors=F)
    
    corrected <- NULL
    
	if(nrow(df_APRS) == nrow(pred_label) + 1){
		## construct df_pred_corrected_label
		df_pred_corrected_label <- data.frame(name=rownames(pred_label), pred=pred_label$pred, label=pred_label$label, corrected=pred_label$scaled, stringsAsFactors=F)
		if(direction=='+'){
			df_pred_corrected_label <- df_pred_corrected_label %>% dplyr::arrange(-corrected)
		}else{
			df_pred_corrected_label <- df_pred_corrected_label %>% dplyr::arrange(corrected)
			df_pred_corrected_label$corrected <- rev(df_pred_corrected_label$corrected)
		}
		## combine df_pred_corrected_label and df_PRSA
		df_data <- cbind(df_APRS[-1,], df_pred_corrected_label)
		rownames(df_data) <- 1:nrow(df_data)
		
	}else{
		df_data <- data.frame(df_APRS, name=NA, pred=NA, label=NA, corrected=NA, stringsAsFactors=F)
		rownames(df_data) <- 1:nrow(df_data)
	}
    
    ####################################
    
	## fmax: Precision-recall F measure
	if(1){
    	vec <- 2 * df_data$Precision * df_data$Recall / (df_data$Precision + df_data$Recall)
    	fmax <- base::max(vec, na.rm=TRUE)
    	fmax_i <- which(vec==fmax)[1]
    }else{
		res <- ROCR::performance(pred_obj, measure="f", x.measure="rec")
		fmax <- base::max(unlist(res@y.values), na.rm=TRUE)
		fmax_i <- which(unlist(res@y.values)==fmax)[1]
    }
    
	## amax: maximum Accuracy
	amax <- base::max(df_data$Accuracy)
    amax_i <- which(df_data$Accuracy==amax)[1]
    
    if(verbose){
        message(sprintf("In summary, Area under ROC: %.3f", auroc), appendLF=TRUE)
        message(sprintf("\tPR F-max: %.3f (Precision=%.3f and Recall=%.3f).", fmax, df_data$Precision[fmax_i], df_data$Recall[fmax_i]), appendLF=TRUE)
        message(sprintf("\tAccuracy (maximum): %.3f (Sensitivity=%.3f and Specificity=%.3f).", amax, df_data$Recall[amax_i], df_data$Specificity[amax_i]), appendLF=TRUE)
    }
    
    ##############################################################################################  
    
    pPerf <- list(
    		data=df_data,
    		auroc=auroc,
    		fmax=fmax,
    		amax=amax,
    		direction=direction,
    		gp=NULL,
    		Pred_obj=pred_obj
    )
    class(pPerf) <- "pPerf"
    
    Recall <- Precision <- Specificity <- NULL
    if(plot=='PR'){
		p <- ggplot(df_data, aes(x=Recall,y=Precision)) 
		p <- p + geom_line() + theme_bw() + ylab("Precision = TP/(TP+FP)") + xlab("Recall = TP/(TP+FN)") + ylim(0,max(df_data$Precision)) + xlim(0,max(df_data$Recall)) 
		p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		
		## title
		title <- paste0('PR curve (Fmax = ',signif(fmax,digits=3),')')
		p <- p + labs(title=title) + theme(plot.title=element_text(hjust=0.5))
		
		if(highlight){
			label <- paste0('Fmax = ',signif(fmax,digits=3), '\n', 'Precision = ',signif(df_data$Precision[fmax_i],digits=3), '\n', 'Specificity = ',signif(df_data$Recall[fmax_i],digits=3))
			## Add the point with the maximum F
			p <- p + geom_point(data=df_data[fmax_i,], colour="red")
			p <- p + ggrepel::geom_label_repel(data=df_data[fmax_i,], label=label, alpha=0.2, seed=825, fill='red') + ggrepel::geom_label_repel(data=df_data[fmax_i,], label=label, alpha=1, seed=825, fill=NA, color='red')
		}
		
		## put arrows on both axes
		gp <- p + theme(axis.line=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))

		pPerf$gp <- gp
		
    }else if(plot=='ROC'){
		p <- ggplot(df_data, aes(x=1-Specificity,y=Recall)) 
		p <- p + geom_line() + theme_bw() + ylab("Sensitivity = TP/(TP+FN)") + xlab("1 - Specificity = FP/(FP+TN)") + ylim(0,max(df_data$Recall)) + xlim(0,max(1-df_data$Specificity))
		p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		
		## title
		title <- paste0('ROC (AUC = ',signif(auroc,digits=3),')')
		p <- p + labs(title=title) + theme(plot.title=element_text(hjust=0.5))
		
		if(highlight){
			label <- paste0('Accuracy = ',signif(amax,digits=3), '\n', 'Sensitivity = ',signif(df_data$Recall[amax_i],digits=3), '\n', 'Specificity = ',signif(df_data$Specificity[amax_i],digits=3))
			## Add the point with the maximum accuracy
			p <- p + geom_point(data=df_data[amax_i,], colour="red")
			p <-  p + ggrepel::geom_label_repel(data=df_data[amax_i,], label=label, alpha=0.2, seed=825, fill='red') + ggrepel::geom_label_repel(data=df_data[amax_i,], label=label, alpha=1, seed=825, fill=NA, color='red')
		}
		
		## put arrows on both axes
		gp <- p + theme(axis.line=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
		
		pPerf$gp <- gp
		
	}

	invisible(pPerf)
}


