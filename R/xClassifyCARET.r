#' Function to combine predictor matrix in a supervised manner via classifying algorithms.
#'
#' \code{xClassifyCARET} is supposed to integrate predictor matrix in a supervised manner via classifying algorithm. It requires three inputs: 1) Gold Standard Positive (GSP) targets; 2) Gold Standard Negative (GSN) targets; 3) a predictor matrix containing genes in rows and predictors in columns, with their predictive scores inside it. It returns an object of class 'sClass'.
#'
#' @param df_predictor a data frame with columns as predictors, with their predictive scores inside it. This data frame must has row names
#' @param GSP a vector containing Gold Standard Positive (GSP)
#' @param GSN a vector containing Gold Standard Negative (GSN)
#' @param algorithm classifying algorithm. It can be one of "gbm" for Gradient Boosting Machine (GBM), "svmRadial" for Support Vector Machines with Radial Basis Function Kernel (SVM), "rda" for Regularized Discriminant Analysis (RDA), "knn" for k-nearest neighbor (KNN), "pls" for Partial Least Squares (PLS), "nnet" for Neural Network (NNET), "rf" for Random Forest (RF), "myrf" for customised Random Forest (RF), "cforest" for Conditional Inference Random Forest, "glmnet" for glmnet, "glm" for Generalized Linear Model (GLM), "bayesglm" for Bayesian Generalized Linear Model (BGLM), "LogitBoost" for Boosted Logistic Regression (BLR), "xgbLinear" for eXtreme Gradient Boosting as linear booster (XGBL), "xgbTree" for eXtreme Gradient Boosting as tree booster (XGBT)
#' @param nfold an integer specifying the number of folds for cross validataion. Per fold creates balanced splits of the data preserving the overall distribution for each class (GSP and GSN), therefore generating balanced cross-vallidation train sets and testing sets. By default, it is 3 meaning 3-fold cross validation
#' @param nrepeat an integer specifying the number of repeats for cross validataion. By default, it is 10 indicating the cross-validation repeated 10 times
#' @param seed an integer specifying the seed
#' @param cv.aggregateBy the aggregate method used to aggregate results from k-fold cross validataion. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return 
#' an object of class "sClass", a list with following components:
#' \itemize{
#'  \item{\code{prediction}: a data frame of nrow X 5 containing priority information, where nrow is the number of rows in the input data frame, and the 5 columns are "GS" (either 'GSP', or 'GSN', or 'NEW'), "name" (row names of the input data frame), "rank" (ranks of the priority scores), "pvalue" (the cross-fold aggregated p-value of being GSP, per-fold p-value converted from empirical cumulative distribution of the probability of being GSP), "priority" (sqrt(-log10(pvalue)) rescaled into [0,100]])}
#'  \item{\code{predictor}: a data frame, which is the same as the input data frame but inserting two additional columns ('GS' and 'name')}
#'  \item{\code{performance}: a data frame of 1+nPredictor X 4 containing the supervised/predictor performance info, where nPredictor is the number of predictors, and the 4 columns are "auroc" (AUC values), "fmax" (F-max values), "amx" (maximum accuracy), and "direction" ('+' indicating the higher score the better prediction; '-' indicating the higher score the worse prediction)}
#'  \item{\code{importance}: a data frame of nPredictor X 2 containing the predictor importance info, where nPredictor is the number of predictors, two columns are predictor importance measures ("MeanDecreaseAccuracy" and "MeanDecreaseGini") . "MeanDecreaseAccuracy" sees how worse the model performs without each predictor (a high decrease in accuracy would be expected for very informative predictors), while "MeanDecreaseGini" measures how pure the nodes are at the end of the tree (a high score means the predictor was important if each predictor is taken out)}
#'  \item{\code{parameter}: NULL or a data frame detailing tuning parameters and their associated AUC, Sensitivity and Specificity.}
#'  \item{\code{model}: the best model.}
#'  \item{\code{algorithm}: the classifying algorithm.}
#'  \item{\code{cv_model}: a list of models, results from per-fold train set}
#'  \item{\code{cv_prob}: a data frame of nrow X 2+nfold containing the probability of being GSP, where nrow is the number of rows in the input data frame, nfold is the number of folds for cross validataion, and the first two columns are "GS" (either 'GSP', or 'GSN', or 'NEW'), "name" (gene names), and the rest columns storing the per-fold probability of being GSP}
#'  \item{\code{cv_auroc}: a data frame of 1+nPredictor X 4+nfold containing the supervised/predictor ROC info (AUC values), where nPredictor is the number of predictors, nfold is the number of folds for cross validataion, and the first 4 columns are "median" (the median of the AUC values across folds), "mad" (the median of absolute deviation of the AUC values across folds), "min" (the minimum of the AUC values across folds), "max" (the maximum of the AUC values across folds), and the rest columns storing the per-fold AUC values}
#'  \item{\code{cv_fmax}: a data frame of 1+nPredictor X 4+nfold containing the supervised/predictor PR info (F-max values), where nPredictor is the number of predictors, nfold is the number of folds for cross validataion, and the first 4 columns are "median" (the median of the F-max values across folds), "mad" (the median of absolute deviation of the F-max values across folds), "min" (the minimum of the F-max values across folds), "max" (the maximum of the F-max values across folds), and the rest columns storing the per-fold F-max values}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note It will depend on whether a package "caret" and its suggested packages have been installed. It can be installed via: \code{BiocManager::install(c("caret","e1071","gbm","kernlab","klaR","pls","nnet","randomForest","party","glmnet","arm","caTools","xgboost"))}.
#' @export
#' @include xClassifyCARET.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' sClass <- xClassifyCARET(df_predictor, GSP, GSN)
#' }

xClassifyCARET <- function(df_predictor, GSP, GSN, algorithm=c("glm","bayesglm","glmnet","lda","rda","pls","svmRadial","knn","nnet", "gbm","xgbLinear","xgbTree","LogitBoost","cforest","rf","myrf"), nfold=3, nrepeat=10, seed=825, cv.aggregateBy=c("logistic","Ztransform","fishers","orderStatistic","none"), verbose=TRUE)
{
	
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    algorithm <- match.arg(algorithm)
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
				message(sprintf("0. Characters in column %d ('%s') converted into the numeric.", i, names(vec_class)[i], as.character(Sys.time())), appendLF=TRUE)
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
        message(sprintf("1. Gold standards (%d in GSP, %d in GSN) are used for supervised integration of %d predictors (%s).", sum(df_predictor_class$class=='GSP'), sum(df_predictor_class$class=='GSN'), ncol(df_predictor), as.character(Sys.time())), appendLF=TRUE)
    }
	
	#####################################
	# function: varImp.train
	#####################################
	varImp.train <- function(train_fit, useModel=TRUE, nonpara=TRUE, scale=FALSE, ...){
		code <- train_fit$modelInfo
	  	if(is.null(code$varImp)){
	  		useModel <- FALSE
	  	}
	  	if(useModel){
	  		# Model-specific metrics: estimating the contribution of each variable to the model 
			imp <- code$varImp(train_fit$finalModel, ...)
	  	}else{
	  		# Model-independent metrics
			isX <- which(!(colnames(train_fit$trainingData) %in% ".outcome"))
		  	x_dat <- train_fit$trainingData[, isX, drop=FALSE]
		  	y_dat <- train_fit$trainingData[, -isX]
			# For classification, AUC for each predictor
			# For regression, the relationship between each predictor and the outcome is evaluated
			# When nonpara = FALSE, a linear model is fit and the absolute value of the t-value for the slope of the predictor is used. Otherwise, a loess smoother is fit between the outcome and the predictor. The R2 statistic is calculated for this model against the intercept only null model
			imp <- caret::filterVarImp(x_dat, y_dat, nonpara=nonpara, ...)
	  	}
	  	if(scale){
			if(class(train_fit$finalModel)[1] == "pamrtrained"){
				imp <- abs(imp)
			}
			imp <- imp - min(imp, na.rm=TRUE)
			imp <- imp/max(imp, na.rm = TRUE)*100
	  	}
	  	
	  	return(imp)
	}
	###############################################################################################################
	nfold <- as.integer(nfold)
	nrepeat <- as.integer(nrepeat)
    if(verbose){
        message(sprintf("2. Optimise tuning parameters of classifying algorithm '%s' using %d-fold cross validation repeated %d times (%s) ....", algorithm, nfold, nrepeat, as.character(Sys.time())), appendLF=TRUE)
    }
	
	fitControl <- caret::trainControl(method=c("repeatedcv","cv","oob")[1], number=nfold, repeats=nrepeat, classProbs=TRUE, summaryFunction=caret::twoClassSummary, allowParallel=FALSE)
	fitControl_withoutParameters <- caret::trainControl(method="none", classProbs=TRUE, allowParallel=FALSE)
	
	if(algorithm=="glm"){
		#####################################
		# Generalized Linear Model (GLM)
		# library(stats)
		# caret::getModelInfo("glm")$glm$grid
		#####################################
		
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "glm", 
								family = "binomial",
								trControl = fitControl,
								metric = "ROC"
								)
								
		# get fit parameters
		df_fit <- train_fit$results[,2:4]
		
		# get predictor importance 
		# Model-specific metrics: the absolute value of the t-statistic for each predictor
		df_importance <- varImp.train(train_fit, useModel=T)
		
		# define prob function
		func_prob <- caret::getModelInfo("glm")$glm$prob

		# get a list of models from repeated cross-validation
		suppressWarnings(
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "glm", 
								family = "binomial",
								trControl = fitControl_withoutParameters,
								metric = "ROC"
								)
			fit$finalModel
		})
		)
	
    }else if(algorithm=="bayesglm"){
		#####################################
		# Bayesian Generalized Linear Model (BGLM)
		# library(arm)
		# caret::getModelInfo("bayesglm")$bayesglm$grid
		#####################################
		
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "bayesglm", 
								family = "binomial",
								trControl = fitControl,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,2:4]
		
		# get predictor importance 
		# Model-independent metrics: class-specific AUC (maximum across multi-class)
		df_importance <- varImp.train(train_fit, useModel=F)
		
		# get prob function
		func_prob <- caret::getModelInfo("bayesglm")$bayesglm$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "bayesglm",
								family = "binomial",
								trControl = fitControl_withoutParameters,
								metric = "ROC"
								)
			fit$finalModel
		})
	
    }else if(algorithm=="glmnet"){
		#####################################
		# Lasso and Elastic-Net Regularized Generalized Linear Models (GLMNET)
		# library(glmnet)
		# caret::getModelInfo("glmnet")$glmnet$grid
		#####################################
		
		grid_glmnet <- caret::getModelInfo("glmnet")$glmnet$grid(x=df_predictor_class[,-ncol(df_predictor_class)], y=df_predictor_class$class, len=10)				
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "glmnet", 
								family = "binomial",
								trControl = fitControl,
								tuneGrid = grid_glmnet,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(3:5,1,2)]

		# get predictor importance (based on finalModel)
		# Model-specific metrics: the absolute value of the coefficients corresponding the the tuned model
		df_importance <- varImp.train(train_fit, useModel=T)

		# get prob function
		func_prob <- caret::getModelInfo("glmnet")$glmnet$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "glmnet",
								family = "binomial",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
	
	}else if(algorithm=="lda"){
		#####################################
		# Linear Discriminant Analysis
		# library(MASS)
		# caret::getModelInfo("lda")$lda$grid
		#####################################
		
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "lda", 
								trControl = fitControl,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,2:4]
		
		# get predictor importance 
		# Model-independent metrics: class-specific AUC (maximum across multi-class)
		df_importance <- varImp.train(train_fit, useModel=F)
		
		# get prob function
		func_prob <- caret::getModelInfo("lda")$lda$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "lda", 
								trControl = fitControl_withoutParameters,
								metric = "ROC"
								)
			fit$finalModel
		})
	
	}else if(algorithm=="rda"){
		#####################################
		# Regularized Discriminant Analysis
		# library(klaR)
		# caret::getModelInfo("rda")$rda$grid
		#####################################
		
		grid_rda <-  base::expand.grid(
						gamma = seq(0.1,1,0.1), 
						lambda = seq(0.1,1,0.1)
						)
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "rda", 
								trControl = fitControl,
								tuneGrid = grid_rda,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(3:5,1,2)]
		
		# get predictor importance 
		# Model-independent metrics: class-specific AUC (maximum across multi-class)
		df_importance <- varImp.train(train_fit, useModel=F)
		
		# get prob function
		func_prob <- caret::getModelInfo("rda")$rda$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "rda", 
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(algorithm=="pls"){
		#####################################
		# Partial Least Squares (PLS)
		# library(pls)
		# caret::getModelInfo("pls")$pls$grid
		#####################################
		
		grid_pls <-  base::expand.grid(
						ncomp = seq(1,ncol(df_predictor_class)-2,1)
						)
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "pls", 
								trControl = fitControl,
								tuneGrid = grid_pls,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(2:4,1)]
		
		# get predictor importance (based on finalModel)
		# Model-specific metrics: based on weighted sums of the absolute regression coefficients. The weights are a function of the reduction of the sums of squares across the number of PLS components and are computed separately for each outcome
		df_importance <- varImp.train(train_fit, useModel=T)
		
		# get prob function
		func_prob <- caret::getModelInfo("pls")$pls$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "pls",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(algorithm=="svmRadial"){
		#####################################
		# Support Vector Machines with Radial Basis Function Kernel (SVM)
		# library(kernlab)
		# caret::getModelInfo("svmRadial")$svmRadial$grid
		#####################################
		
		grid_svm <- caret::getModelInfo("svmRadial")$svmRadial$grid(x=df_predictor_class[,-ncol(df_predictor_class)], y=df_predictor_class$class, len=10)
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "svmRadial", 
								trControl = fitControl,
								tuneGrid = grid_svm,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(3:5,2)]
		
		# get predictor importance 
		# Model-independent metrics: class-specific AUC (maximum across multi-class)
		df_importance <- varImp.train(train_fit, useModel=F)
		
		# get prob function
		func_prob <- caret::getModelInfo("svmRadial")$svmRadial$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "svmRadial",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
	
    }else if(algorithm=="knn"){
		#####################################
		# Support Vector Machines with Radial Basis Function Kernel (SVM)
		# library(caret)
		# caret::getModelInfo("knn")$knn$grid
		#####################################
		
		grid_knn <- caret::getModelInfo("knn")$knn$grid(len=ncol(df_predictor_class)-1)
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "knn", 
								trControl = fitControl,
								tuneGrid = grid_knn,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(2:4,1)]
		
		# get predictor importance 
		# Model-independent metrics: class-specific AUC (maximum across multi-class)
		df_importance <- varImp.train(train_fit, useModel=F)
		
		# get prob function
		func_prob <- caret::getModelInfo("knn")$knn$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "knn",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(algorithm=="nnet"){
		#####################################
		# Neural Network (NNET)
		# library(nnet)
		# caret::getModelInfo("nnet")$nnet$grid
		#####################################
		
		len <- 10
		grid_nnet <- base::expand.grid(
										size = ((1:len)*2) - 1,
										decay = c(0, 10^seq(-1,-4,length=len-1))
									   )
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "nnet", 
								trControl = fitControl,
								tuneGrid = grid_nnet,
								trace = FALSE,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(3:5,1,2)]
		
		# get predictor importance (based on finalModel)
		# Model-specific metrics: combinations of the absolute values of the weights
		df_importance <- varImp.train(train_fit, useModel=T)
		
		# get prob function
		func_prob <- caret::getModelInfo("nnet")$nnet$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "nnet",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								trace = FALSE,
								metric = "ROC"
								)
			fit$finalModel
		})
		
   }else if(algorithm=="gbm"){
		#####################################
		# Gradient Boosting Machine (GBM) model
		# library(gbm)
		# caret::getModelInfo("gbm")$gbm$grid
		#####################################
		
		len <- 10
		grid_gbm <- base::expand.grid(
										n.trees = (1:len)*50, 
										interaction.depth = 1:len,
										shrinkage = 0.1,
										n.minobsinnode = 10
									   )
		#grid_gbm <- caret::getModelInfo("gbm")$gbm$grid(len=len)			   
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "gbm", 
								trControl = fitControl,
								tuneGrid = grid_gbm,
								verbose = FALSE,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(5:7,2,4)]
		
		# get predictor importance (based on finalModel)
		# Model-specific metrics: sums the importances over each boosting iteration
		df_importance <- varImp.train(train_fit, useModel=T)
		
		# get prob function
		func_prob <- caret::getModelInfo("gbm")$gbm$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "gbm",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								verbose = FALSE,
								metric = "ROC"
								)
			fit$finalModel
		})
		
   }else if(algorithm=="xgbLinear"){
		#####################################
		# eXtreme Gradient Boosting linear (XGBL)
		# library(xgboost)
		# caret::getModelInfo("xgbLinear")$xgbLinear$grid
		#####################################
		
		len <- 10
		grid_xgbl <- base::expand.grid(
						  				lambda = 0,
										alpha = c(0, 10^seq(-1, -4, length=len-1)),
										nrounds = floor((1:len)*50),
										eta = 0.3
						  				)
		#len <- 5
		#grid_xgbl <- caret::getModelInfo("xgbLinear")$xgbLinear$grid(len=len)	   
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "xgbLinear", 
								trControl = fitControl,
								tuneGrid = grid_xgbl,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(5:7,2,3)]

		# get predictor importance (based on finalModel)
		# Model-specific metrics: the absolute magnitude of linear coefficients
		df_importance <- varImp.train(train_fit, useModel=T)

		# get prob function
		func_prob <- caret::getModelInfo("xgbLinear")$xgbLinear$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "xgbLinear",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
   }else if(algorithm=="xgbTree"){
		#####################################
		# eXtreme Gradient Boosting tree (XGBT)
		# library(xgboost)
		# caret::getModelInfo("xgbTree")$xgbTree$grid
		#####################################
		
		len <- 10
		grid_xgbt <- base::expand.grid(
										max_depth = seq(1, len),
										nrounds = floor((1:len) * 50),
										eta = 0.3,
										gamma = 0,
										colsample_bytree = 1,
										min_child_weight = 1,
										subsample = 0.75
									   )
		#len <- 3
		#grid_xgbt <- caret::getModelInfo("xgbTree")$xgbTree$grid(len=len)	   
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "xgbTree", 
								trControl = fitControl,
								tuneGrid = grid_xgbt,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(8:10,2,7)]

		# get predictor importance (based on finalModel)
		# Model-specific metrics: fractional contribution of each feature to the model based on the total gain of this feature's splits. Higher percentage means a more important predictive feature
		df_importance <- varImp.train(train_fit, useModel=T)
		
		# get prob function
		func_prob <- caret::getModelInfo("xgbTree")$xgbTree$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "xgbTree",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})

   }else if(algorithm=="LogitBoost"){
		#####################################
		# Boosted Logistic Regression (BLR)
		# library(caTools)
		# caret::getModelInfo("LogitBoost")$LogitBoost$grid
		#####################################
		
		len <- 30
		grid_blr <- caret::getModelInfo("LogitBoost")$LogitBoost$grid(len=len)
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "LogitBoost", 
								trControl = fitControl,
								tuneGrid = grid_blr,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(2:4,1)]

		# get predictor importance (based on finalModel)
		# Model-independent metrics: class-specific AUC (maximum across multi-class)
		df_importance <- varImp.train(train_fit, useModel=F)
		
		# get prob function
		func_prob <- caret::getModelInfo("LogitBoost")$LogitBoost$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "LogitBoost",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
   }else if(algorithm=="cforest"){
		#####################################
		# Conditional Inference Random Forest (CRF)
		# library(party)
		# caret::getModelInfo("cforest")$cforest$grid
		#####################################
		
		grid_crf <- caret::getModelInfo("cforest")$cforest$grid(x=df_predictor_class[,-ncol(df_predictor_class)], y=df_predictor_class$class, len=10)
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "cforest", 
								trControl = fitControl,
								tuneGrid = grid_crf,
								metric = "ROC"
								)
								
		# get fit parameters
		df_fit <- train_fit$results[,c(2:4,1)]
		
		# get predictor importance (based on finalModel)
		# Model-specific metrics: standard and conditional variable importance for 'cforest', following the permutation principle of the mean decrease in accuracy importance in randomForest
		df_importance <- varImp.train(train_fit, useModel=T)
		
		# get prob function
		func_prob <- caret::getModelInfo("cforest")$cforest$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "cforest",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
   }else if(algorithm=="rf"){
		#####################################
		# Random Forest (RF)
		# library(randomForest)
		# caret::getModelInfo("rf")$rf$grid
		#####################################
		
		ntree <- 1000
		grid_rf <- caret::getModelInfo("rf")$rf$grid(x=df_predictor_class[,-ncol(df_predictor_class)], y=df_predictor_class$class, len=10)
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "rf", 
								trControl = fitControl,
								tuneGrid = grid_rf,
								ntree = ntree,
								importance = TRUE,
								metric = "ROC"
								)

		# get fit parameters
		df_fit <- train_fit$results[,c(2:4,1)]
		
		# get predictor importance (based on finalModel)
		df_importance <- varImp.train(train_fit, useModel=T)
		
		# get prob function
		func_prob <- caret::getModelInfo("rf")$rf$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = "rf",
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								ntree = ntree,
								importance = TRUE,
								metric = "ROC"
								)
			fit$finalModel
		})
		
   }else if(algorithm=="myrf"){
		#####################################
		# my customised Random Forest
		# library(randomForest)
		# caret::getModelInfo("rf")$rf
		#####################################
		
		##############################
		my_rf <- list(label = "my_RF",
               library = "randomForest",
               type = "Classification",
               parameters = data.frame(parameter=c('ntree','mtry'),
                                       class=c("numeric",'numeric'),
                                       label=c('#Number of Trees','#Randomly Selected Predictors')
                                       ),
               grid = function(x, y, len=NULL, search="grid"){
                	if(search=="grid"){
                   		grid <- expand.grid(ntree=seq(100,200*len,by=200),
                                       		mtry=caret::var_seq(p=ncol(x),len=len)
                                       		)
                   	}
               },
               loop = NULL,
               fit = function(x, y, wts, param, lev, last, classProbs, ...) { 
            		randomForest::randomForest(x, y, mtry=param$mtry, ntree=param$ntree, ...)
               },
               predict = function(modelFit, newdata, submodels = NULL) {  
               		if(!is.null(newdata)) predict(modelFit, newdata) else predict(modelFit)
               },
               prob = function(modelFit, newdata, submodels = NULL){
               		if(!is.null(newdata)) predict(modelFit, newdata, type="prob") else predict(modelFit, type = "prob")
               },
               varImp = function(object, ...){
                    varImp <- randomForest::importance(object, ...)
                    if(object$type=="regression"){
                      varImp <- data.frame(Overall=varImp[,"%IncMSE"], stringsAsFactors=F)
                    }else{
                      retainNames <- levels(object$y)
                      varImp <- varImp[, (1+length(retainNames)):ncol(varImp)]
                    }
                    out <- as.data.frame(varImp)
                },
               predictors = function(x, ...) {
                    ## After doing some testing, it looks like randomForest
                    ## will only try to split on plain main effects (instead
                    ## of interactions or terms like I(x^2).
                    varIndex <- as.numeric(names(table(x$forest$bestvar)))
                    varIndex <- varIndex[varIndex > 0]
                    varsUsed <- names(x$forest$ncat)[varIndex]
                    varsUsed
                  },
               levels = function(x) x$classes,
               sort = function(x) x[order(x[,1]),]
            )
		##############################
		
		len <- 10
		grid_myrf <-  my_rf$grid(x=df_predictor_class[,-ncol(df_predictor_class)], len=len)
		if(!is.null(seed)) set.seed(seed)
		train_fit <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = my_rf, 
								trControl = fitControl,
								tuneGrid = grid_myrf,
								importance = TRUE,
								metric = "ROC"
								)
		
		# get fit parameters
		df_fit <- train_fit$results[,c(3:5,1,2)]
		
		# get predictor importance (based on finalModel)
		df_importance <- varImp.train(train_fit, useModel=T)
		
		# get prob function
		func_prob <- my_rf$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(train_fit$control$index, function(k){
			x <- train_fit$trainingData[k, -1]
			y <- train_fit$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = my_rf,
								trControl = fitControl_withoutParameters,
								tuneGrid = train_fit$bestTune,
								importance = TRUE,
								metric = "ROC"
								)
			fit$finalModel
		})
		
	}
	
	if(verbose){
		if(length(train_fit$bestTune)==1 && train_fit$bestTune=='none'){
			message(sprintf("\tThe best model with %s = %.3f", train_fit$metric, df_fit$ROC), appendLF=TRUE)
		}else{
			message(sprintf("\tThe best model with optimised parameters (%s) to achieve the highest %s = %.3f", paste(names(train_fit$bestTune),"=",format(train_fit$bestTune),collapse=' and '), train_fit$metric, df_fit$ROC[as.integer(rownames(train_fit$bestTune))]), appendLF=TRUE)
		}
	}
	###############################################################################################################

    if(verbose){
        message(sprintf("3. Performance evaluation using test sets (%s).", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("\tExtract the ROC matrix of %d rows (Supervised + predictors) X %d columns/repeats*folds (%s).", ncol(df_predictor_class), nfold*nrepeat, as.character(Sys.time())), appendLF=TRUE)
    }
	
	######################
	## evaluation per fold
	######################
	lsls_predictors <- lapply(1:length(ls_model), function(i){
		## prediction for testset
		testset <- train_fit$trainingData[train_fit$control$indexOut[[i]], -1]
		vec_predict_test <- func_prob(ls_model[[i]], newdata=testset)[,'GSP']
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
        message(sprintf("4. Do prediction for fullset (%s).", as.character(Sys.time())), appendLF=TRUE)
        message(sprintf("\tExtract the full prediction matrix of %d rows/genes X %d columns/repeats*folds, aggregated via '%s' (%s) ...", nrow(df_predictor_class), nfold*nrepeat, cv.aggregateBy, as.character(Sys.time())), appendLF=TRUE)
    }
	
	if(cv.aggregateBy=="none"){
		## prediction for fullset
		fullset <- df_predictor
		vec_predict_full <- func_prob(train_fit$finalModel, newdata=fullset)[,'GSP']
		names(vec_predict_full) <- rownames(fullset)
		vec_full <- sort(vec_predict_full, decreasing=TRUE)
	
		## get rank
		vec_rank <- rank(-1*vec_full, ties.method="min")
	
		## priority: being rescaled into the [0,100] range
		priority <- vec_full
		priority <- 100 * (priority - min(priority))/(max(priority) - min(priority))
		priority <- signif(priority, digits=2)
		
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
		suppressMessages(rf.model.overall <- randomForest::randomForest(class ~ ., data=df_predictor_class, importance=TRUE, ntree=ntree, mtry=mtry))
		df_importance <- as.data.frame(randomForest::importance(rf.model.overall, type=NULL, class=NULL, scale=TRUE)[,3:4])
		#randomForest::varImpPlot(rf.model.overall)
		# For classification, the nclass columns are the class-specific measures computed as mean descrease in accuracy. The nclass + 1st column is the mean descrease in accuracy over all classes. The last column is the mean decrease in Gini index
		
		res <- list(importance = df_importance)
		
		
	}else{
	
		ls_full <- lapply(1:length(ls_model), function(i){
			## prediction for fullset
			fullset <- df_predictor
			vec_predict_full <- func_prob(ls_model[[i]], newdata=fullset)[,'GSP']
			names(vec_predict_full) <- rownames(fullset)
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
		#df_importance <- caret::varImp(train_fit, scale=F)

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
	
		######################
		## final model
		######################
		model <- train_fit$finalModel

		res <- list(
						prediction = df_priority,
						predictor = df_predictor_gs,
						performance = df_evaluation,
						importance = df_importance,
						parameter = df_fit,
						model = model,
						algorithm = algorithm,
						cv_model = ls_model,
						cv_prob = df_prob,
						cv_auroc = df_auroc,
						cv_fmax = df_fmax,
						call  = match.call()
					 )
		class(res) <- "sClass"
	
   	}
	
  	####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(res)
}


