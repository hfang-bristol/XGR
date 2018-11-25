#' Function to combine predictor matrix in a supervised manner via machine learning algorithm random forest.
#'
#' \code{xClassifyRF} is supposed to integrate predictor matrix in a supervised manner via machine learning algorithm random forest. It requires three inputs: 1) Gold Standard Positive (GSP) targets; 2) Gold Standard Negative (GSN) targets; 3) a predictor matrix containing genes in rows and predictors in columns, with their predictive scores inside it. It returns an object of class 'sClass'.
#'
#' @param df_predictor a data frame with columns as predictors, with their predictive scores inside it. This data frame must has row names
#' @param GSP a vector containing Gold Standard Positive (GSP)
#' @param GSN a vector containing Gold Standard Negative (GSN)
#' @param nfold an integer specifying the number of folds for cross validataion. Per fold creates balanced splits of the data preserving the overall distribution for each class (GSP and GSN), therefore generating balanced cross-vallidation train sets and testing sets. By default, it is 3 meaning 3-fold cross validation
#' @param nrepeat an integer specifying the number of repeats for cross validataion. By default, it is 10 indicating the cross-validation repeated 10 times
#' @param seed an integer specifying the seed
#' @param mtry an integer specifying the number of predictors randomly sampled as candidates at each split. If NULL, it will be tuned by `randomForest::tuneRF`, with starting value as sqrt(p) where p is the number of predictors. The minimum value is 3
#' @param ntree an integer specifying the number of trees to grow. By default, it sets to 2000
#' @param cv.aggregateBy the aggregate method used to aggregate results from k-fold cross validataion. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param ... additional parameters. Please refer to 'randomForest::randomForest' for the complete list.
#' @return 
#' an object of class "sClass", a list with following components:
#' \itemize{
#'  \item{\code{prediction}: a data frame of nrow X 5 containing priority information, where nrow is the number of rows in the input data frame, and the 5 columns are "GS" (either 'GSP', or 'GSN', or 'NEW'), "name" (row names of the input data frame), "rank" (ranks of the priority scores), "pvalue" (the cross-fold aggregated p-value of being GSP, per-fold p-value converted from empirical cumulative distribution of the probability of being GSP), "priority" (sqrt(-log10(pvalue)) rescaled into [0,100]])}
#'  \item{\code{predictor}: a data frame, which is the same as the input data frame but inserting two additional columns ('GS' and 'name')}
#'  \item{\code{performance}: a data frame of 1+nPredictor X 4 containing the supervised/predictor performance info, where nPredictor is the number of predictors, and the 4 columns are "auroc" (AUC values), "fmax" (F-max values), "amx" (maximum accuracy), and "direction" ('+' indicating the higher score the better prediction; '-' indicating the higher score the worse prediction)}
#'  \item{\code{importance}: a data frame of nPredictor X 2 containing the predictor importance info, where nPredictor is the number of predictors, two columns are predictor importance measures ("MeanDecreaseAccuracy" and "MeanDecreaseGini") . "MeanDecreaseAccuracy" sees how worse the model performs without each predictor (a high decrease in accuracy would be expected for very informative predictors), while "MeanDecreaseGini" measures how pure the nodes are at the end of the tree (a high score means the predictor was important if each predictor is taken out)}
#'  \item{\code{parameter}: NULL or a data frame containing tuning parameters and their associated AUC, Sensitivity and Specificity.}
#'  \item{\code{model}: the best model.}
#'  \item{\code{algorithm}: the classifying algorithm.}
#'  \item{\code{cv_model}: a list of models, results from per-fold train set}
#'  \item{\code{cv_prob}: a data frame of nrow X 2+nfold containing the probability of being GSP, where nrow is the number of rows in the input data frame, nfold is the number of folds for cross validataion, and the first two columns are "GS" (either 'GSP', or 'GSN', or 'NEW'), "name" (gene names), and the rest columns storing the per-fold probability of being GSP}
#'  \item{\code{cv_auroc}: a data frame of 1+nPredictor X 4+nfold containing the supervised/predictor ROC info (AUC values), where nPredictor is the number of predictors, nfold is the number of folds for cross validataion, and the first 4 columns are "median" (the median of the AUC values across folds), "mad" (the median of absolute deviation of the AUC values across folds), "min" (the minimum of the AUC values across folds), "max" (the maximum of the AUC values across folds), and the rest columns storing the per-fold AUC values}
#'  \item{\code{cv_fmax}: a data frame of 1+nPredictor X 4+nfold containing the supervised/predictor PR info (F-max values), where nPredictor is the number of predictors, nfold is the number of folds for cross validataion, and the first 4 columns are "median" (the median of the F-max values across folds), "mad" (the median of absolute deviation of the F-max values across folds), "min" (the minimum of the F-max values across folds), "max" (the maximum of the F-max values across folds), and the rest columns storing the per-fold F-max values}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note none
#' @export
#' @include xClassifyRF.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' sClass <- xClassifyRF(df_predictor, GSP, GSN)
#' }

xClassifyRF <- function(df_predictor, GSP, GSN, nfold=3, nrepeat=10, seed=825, mtry=NULL, ntree=1000, cv.aggregateBy=c("none","logistic","Ztransform","fishers","orderStatistic"), verbose=TRUE, ...)
{
	
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    cv.aggregateBy <- match.arg(cv.aggregateBy)
	
    if(!is.data.frame(df_predictor)){
    	warnings("The function must apply to a 'data.frame'.\n")
        return(NULL)
    }else{
    	if(is.null(rownames(df_predictor))){
			warnings("The input 'data.frame' must have rownames.\n")
			return(NULL)
    	}
    }
	
	# convert characters into numeric
	vec_class <- sapply(df_predictor,class)
	for(i in 1:ncol(df_predictor)){
		if(vec_class[i]=='character'){
			df_predictor[,i] <- as.numeric(as.factor(df_predictor[,i]))
			if(verbose){
				message(sprintf("Characters in column %d ('%s') converted into the numeric.", i, names(vec_class)[i], as.character(Sys.time())), appendLF=TRUE)
			}
		}
	}
		
	## pre-process GSP and GSN
	gsp <- unique(GSP)
	gsn <- unique(GSN)
	gsn <- setdiff(gsn, gsp)
	df_gsp <- data.frame(name=gsp, gs=1, stringsAsFactors=F)
	df_gsn <- data.frame(name=gsn, gs=0, stringsAsFactors=F)
	df_gs <- rbind(df_gsp, df_gsn)
	
	## predictors + class
	ind <- match(df_gs$name, rownames(df_predictor))
	df_predictor_class <- df_predictor[ind[!is.na(ind)],]
	## class as factor ("GSN","GSP")
	class <- df_gs$gs[!is.na(ind)]
	class <- ifelse(class==1, 'GSP', 'GSN')
	class <- factor(class, levels=c("GSN","GSP"))
	df_predictor_class$class <- class
	
    if(verbose){
        message(sprintf("1. Gold standards (%d in GSP, %d in GSN) are used for supervised integration of %d predictors (%s).", sum(df_predictor_class$class==1), sum(df_predictor_class$class==0), ncol(df_predictor), as.character(Sys.time())), appendLF=TRUE)
    }
	
	#####################################

	nfold <- as.integer(nfold)
	nrepeat <- as.integer(nrepeat)
    if(verbose){
        if(nfold==1){
			message(sprintf("2. GS matrix of %d rows X %d columns (predictors+class) are used as train set (%s) ...", nrow(df_predictor_class), ncol(df_predictor_class), as.character(Sys.time())), appendLF=TRUE)
        }else{
        	message(sprintf("2. GS matrix of %d rows X %d columns (predictors+class) are split (by rows) into %d non-redundant sets: each fold '%d/%d' are used as train set and the remaining '1/%d' as test set. These splits are repeated over %d times (%s) ...", nrow(df_predictor_class), ncol(df_predictor_class), nfold, nfold-1, nfold, nfold, nrepeat, as.character(Sys.time())), appendLF=TRUE)
        }
    }
    
    ## non-redundant sets
    # balanced splits of the data preserving the overall class distribution to generate balanced cross-validation groupings
	if(!is.null(seed)){
		set.seed(seed)
	}
	index_sets <- caret::createMultiFolds(y=df_predictor_class$class, k=nfold, times=nrepeat)
	
	#####################################
	
    if(verbose){
        if(nfold==1){
			message(sprintf("3. Do random forest: '1/%d' as train set (%s) ...", nfold, as.character(Sys.time())), appendLF=TRUE)
		}else{
			message(sprintf("3. Do random forest: %d-fold cross-validation (%s) ...", nfold, as.character(Sys.time())), appendLF=TRUE)
		}
    }
	
	if(nfold==1 & nrepeat==1){
		ls_model <- lapply(1:length(index_sets), function(i){
			trainset <- df_predictor_class[index_sets[[i]],]
			
			if(verbose){
				message(sprintf("\tFold %d: %d GSP + %d GSN", i, table(trainset$class)[2], table(trainset$class)[1]), appendLF=TRUE)
			}
			
			if(is.null(mtry)){
				mtry <- as.integer(sqrt(ncol(trainset)-1))
				suppressMessages(df_mtry <- randomForest::tuneRF(x=trainset[,-ncol(trainset)], y=trainset[, ncol(trainset)], ntreeTry=ntree, mtryStart=mtry, stepFactor=2, trace=FALSE, plot=FALSE))
				ind <- which(df_mtry[,2] == min(df_mtry[,2]))[1]
				mtry <- df_mtry[ind,1]
				if(mtry<3){
					mtry <- 3
				}
			}
			set.seed(i)
			suppressMessages(rf.model <- randomForest::randomForest(class ~ ., data=trainset, importance=TRUE, ntree=ntree, mtry=mtry, ...))
		})
		
	}else{
		ls_model <- lapply(1:length(index_sets), function(i){
			trainset <- df_predictor_class[index_sets[[i]],]
		
			if(verbose){
				message(sprintf("\tRepeatFold %d: %d GSP + %d GSN", i, table(trainset$class)[2], table(trainset$class)[1]), appendLF=TRUE)
			}
		
			if(is.null(mtry)){
				mtry <- as.integer(sqrt(ncol(trainset)))
				
				invisible(utils::capture.output(df_mtry <- randomForest::tuneRF(x=trainset[,-ncol(trainset)], y=trainset[, ncol(trainset)], ntreeTry=ntree, mtryStart=mtry, stepFactor=2, trace=FALSE, plot=FALSE, na.action=na.omit)))
				
				ind <- which(df_mtry[,2] == min(df_mtry[,2]))[1]
				mtry <- df_mtry[ind,1]
				if(mtry<3){
					mtry <- 3
				}
			}

			set.seed(i)
			suppressMessages(rf.model <- randomForest::randomForest(class ~ ., data=trainset, importance=TRUE, ntree=ntree, mtry=mtry, ...))
		})
		names(ls_model) <- names(index_sets)
	
	}
	
	#####################################

    if(verbose){
        message(sprintf("4. Performance evaluation using test sets (%s).", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("Extract the ROC matrix of %d rows (Supervised + predictors) X %d columns/repeats*folds (%s).", ncol(df_predictor_class), nfold*nrepeat, as.character(Sys.time())), appendLF=TRUE)
    }
	
	######################
	## evaluation per fold
	######################
	lsls_predictors <- lapply(1:length(ls_model), function(i){
		rf.model <- ls_model[[i]]
		## prediction for testset: ?predict.randomForest
		testset <- df_predictor_class[-index_sets[[i]],]
		vec_predict_test <- predict(rf.model, newdata=testset[,-ncol(testset)], type='prob')[,'GSP']
		### do preparation
		ind <- match(rownames(testset), rownames(df_predictor))
		df_predictor_test <- cbind(Supervised=as.numeric(vec_predict_test), df_predictor[ind,])
		rownames(df_predictor_test) <- rownames(df_predictor[ind,])
		df_pred <- df_predictor_test
		ls_predictors <- lapply(colnames(df_pred), function(x){
			data.frame(rownames(df_pred), df_pred[,x], stringsAsFactors=FALSE)
		})
		names(ls_predictors) <- colnames(df_pred)
		return(ls_predictors)
	})
	names(lsls_predictors) <- names(ls_model)
	ls_res <- lapply(1:length(lsls_predictors), function(i){
		ls_predictors <- lsls_predictors[[i]]
		# do evaluation
		ls_pPerf <- lapply(ls_predictors, function(x){
			pPerf <- xClassifyPerf(prediction=x, GSP=GSP, GSN=GSN, verbose=FALSE)
		})
		
		# do plotting
		bp <- xClassifyComp(ls_pPerf)
		
		if(is.null(bp)){
			df <- NULL
		}else{
			df <- unique(bp$data[,c('methods','auroc','fmax','amax','direction')])
			df$model <- names(ls_model)[i]
		}
		
		return(df)
	})
	df_res <- do.call(rbind, ls_res)
	if(!is.null(df_res)){
		## df_auroc
		df_res <- as.matrix(xSparseMatrix(df_res[,c('methods','model','auroc')], verbose=FALSE))
		ind <- match(c("Supervised",colnames(df_predictor)), rownames(df_res))
		if(nfold==1 & nrepeat==1){
			df_res <- as.matrix(df_res[ind,], ncol=nfold)
			colnames(df_res) <- 'Fold1'
		}else{
			df_res <- df_res[ind,]
		}
		vec_median <- apply(df_res, 1, stats::median)
		vec_mad <- apply(df_res, 1, stats::mad)
		vec_min <- apply(df_res, 1, base::min)
		vec_max <- apply(df_res, 1, base::max)
		df_auroc <- cbind(median=vec_median, mad=vec_mad, min=vec_min, max=vec_max, df_res)
		
		## df_fmax
		df_res <- do.call(rbind, ls_res)
		df_res <- as.matrix(xSparseMatrix(df_res[,c('methods','model','fmax')], verbose=FALSE))
		ind <- match(c("Supervised",colnames(df_predictor)), rownames(df_res))
		if(nfold==1 & nrepeat==1){
			df_res <- as.matrix(df_res[ind,], ncol=nfold)
			colnames(df_res) <- 'Fold1'
		}else{
			df_res <- df_res[ind,]
		}
		vec_median <- apply(df_res, 1, stats::median)
		vec_mad <- apply(df_res, 1, stats::mad)
		vec_min <- apply(df_res, 1, base::min)
		vec_max <- apply(df_res, 1, base::max)
		df_fmax <- cbind(median=vec_median, mad=vec_mad, min=vec_min, max=vec_max, df_res)
    
    }else{
    	df_auroc <- NULL
    	df_fmax <- NULL
    }
    
	#####################################
    if(verbose){
        message(sprintf("5. Do prediction for fullset (%s).", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("Extract the full prediction matrix of %d rows/genes X %d columns/repeats*folds, aggregated via '%s' (%s) ...", nrow(df_predictor_class), nfold*nrepeat, cv.aggregateBy, as.character(Sys.time())), appendLF=TRUE)
    }
	
	if(cv.aggregateBy=="none"){
		# rf.model
		if(is.null(mtry)){
			mtry <- as.integer(sqrt(ncol(df_predictor_class)-1))
			invisible(utils::capture.output(df_mtry <- randomForest::tuneRF(x=df_predictor_class[,-ncol(df_predictor_class)], y=df_predictor_class[, ncol(df_predictor_class)], ntreeTry=ntree, mtryStart=mtry, stepFactor=2, trace=FALSE, plot=FALSE)))
			ind <- which(df_mtry[,2] == min(df_mtry[,2]))[1]
			mtry <- df_mtry[ind,1]
			if(mtry<3){
				mtry <- 3
			}
		}
		if(!is.null(seed)){
			set.seed(seed)
		}
		suppressMessages(rf.model <- randomForest::randomForest(class ~ ., data=df_predictor_class, importance=TRUE, ntree=ntree, mtry=mtry, ...))
	
		## prediction for fullset
		fullset <- df_predictor_class
		vec_predict_full <- predict(rf.model, newdata=fullset, type='prob')[,'GSP']
		names(vec_predict_full) <- rownames(df_predictor)
		vec_full <- sort(vec_predict_full, decreasing=TRUE)
	
		## get rank
		vec_rank <- rank(-1*vec_full, ties.method="min")
	
		## priority: being rescaled into the [0,100] range
		priority <- vec_full
		priority <- 100 * (priority - min(priority))/(max(priority) - min(priority))
		priority <- signif(priority, digits=2)
		
	}else{
		ls_full <- lapply(1:length(ls_model), function(i){
			rf.model <- ls_model[[i]]
			## prediction for fullset: ?predict.randomForest
			vec_predict_full <- predict(rf.model, newdata=df_predictor, type='prob')[,'GSP']
			# output
			df <- data.frame(name=names(vec_predict_full), model=names(ls_model)[i], prob=vec_predict_full, stringsAsFactors=FALSE)
		})
		df_full <- do.call(rbind, ls_full)
		df_full <- as.matrix(xSparseMatrix(df_full, verbose=FALSE))

		## Convert into p-values by computing an empirical cumulative distribution function
		ls_pval <- lapply(1:ncol(df_full), function(j){
			x <- df_full[,j]
			my.CDF <- stats::ecdf(x)
			pval <- 1 - my.CDF(x)
		})
		df_pval <- do.call(cbind, ls_pval)
		rownames(df_pval) <- rownames(df_full)
	
		## aggregate p values
		df_ap <- dnet::dPvalAggregate(pmatrix=df_pval, method=cv.aggregateBy)
		df_ap <- sort(df_ap, decreasing=FALSE)
	
		## get rank
		df_rank <- rank(df_ap, ties.method="min")

		## priority: first log10-transformed ap and then being rescaled into the [0,100] range
		df_ap[df_ap==0] <- min(df_ap[df_ap!=0])
		priority <- sqrt(-log10(df_ap))
		priority <- 100* (priority - min(priority))/(max(priority) - min(priority))
		priority <- signif(priority, digits=2)
	
	}
	
	#########################################
	## output
	### df_priority
	output_gs <- rep('NEW', length(df_ap))
	names(output_gs) <- names(df_ap)
	ind <- match(names(output_gs), df_gs$name)
	output_gs[!is.na(ind)] <- df_gs$gs[ind[!is.na(ind)]]
	output_gs[output_gs=='0'] <- 'GSN'
	output_gs[output_gs=='1'] <- 'GSP'
	df_priority <- data.frame(GS=output_gs, name=names(df_ap), rank=df_rank, pvalue=df_ap, priority=priority, stringsAsFactors=FALSE)

	### df_prob
	ind <- match(names(df_ap), rownames(df_full))
	if(nfold==1 & nrepeat==1){
		output_df_full <- as.matrix(df_full[ind,], ncol=nfold)
		colnames(output_df_full) <- 'Fold1'
	}else{
		output_df_full <- df_full[ind,]
	}
	df_prob <- data.frame(GS=output_gs, name=names(df_ap), output_df_full, stringsAsFactors=FALSE)

	### df_predictor_gs
	ind <- match(names(df_ap), rownames(df_predictor))
	output_df_predictor <- df_predictor[ind,]
	df_predictor_gs <- data.frame(GS=output_gs, name=names(df_ap), output_df_predictor, stringsAsFactors=FALSE)
	
    ####################################################################################

	######################
	## overall importance
	######################
	if(is.null(mtry)){
		mtry <- as.integer(sqrt(ncol(df_predictor_class)-1))
		invisible(utils::capture.output(df_mtry <- randomForest::tuneRF(x=df_predictor_class[,-ncol(df_predictor_class)], y=df_predictor_class[, ncol(df_predictor_class)], ntreeTry=ntree, mtryStart=mtry, stepFactor=2, trace=FALSE, plot=FALSE, na.action=na.omit)))
		ind <- which(df_mtry[,2] == min(df_mtry[,2]))[1]
		mtry <- df_mtry[ind,1]
		if(mtry<3){
			mtry <- 3
		}
	}
	if(!is.null(seed)){
		set.seed(seed)
	}
	suppressMessages(rf.model.overall <- randomForest::randomForest(class ~ ., data=df_predictor_class, importance=TRUE, ntree=ntree, mtry=mtry, ...))
	df_importance <- as.data.frame(randomForest::importance(rf.model.overall, type=NULL, class=NULL, scale=TRUE)[,3:4])
	#randomForest::varImpPlot(rf.model.overall)
	# For classification, the nclass columns are the class-specific measures computed as mean descrease in accuracy. The nclass + 1st column is the mean descrease in accuracy over all classes. The last column is the mean decrease in Gini index
	
	######################
	## overall evaluation
	######################
	### do preparation
	df_predictor_overall <- cbind(Supervised=df_priority$priority, df_predictor_gs[,-c(1,2)])
	rownames(df_predictor_overall) <- rownames(df_priority)
	df_pred <- df_predictor_overall
	ls_predictors <- lapply(colnames(df_pred), function(x){
		data.frame(rownames(df_pred), df_pred[,x], stringsAsFactors=FALSE)
	})
	names(ls_predictors) <- colnames(df_pred)
	# do evaluation
	ls_pPerf <- lapply(ls_predictors, function(x){
		pPerf <- xClassifyPerf(prediction=x, GSP=GSP, GSN=GSN, verbose=FALSE)
	})
	# do plotting
	bp <- xClassifyComp(ls_pPerf)
	df <- unique(bp$data[,c('methods','auroc','fmax','amax','direction')])
	df_evaluation <- df[,-1]
	rownames(df_evaluation) <- df[,1]
	#####################
	#####################

    res <- list(
    				prediction = df_priority,
    				predictor = df_predictor_gs,
    				performance = df_evaluation,
    				importance = df_importance,
    				parameter = NULL,
    				model = rf.model.overall, 
    				algorithm = "RandomForest",
    				cv_model = ls_model,				
    				cv_prob = df_prob,
    				cv_auroc = df_auroc,
    				cv_fmax = df_fmax,
                  	call  = match.call()
                 )
    class(res) <- "sClass"    
  ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(res)
}


