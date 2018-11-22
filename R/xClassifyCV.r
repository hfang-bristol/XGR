#' Function to visualise cross-validation performance
#'
#' \code{xClassifyCV} is supposed to visualise cross-validation performance. It returns an object of class "ggplot".
#'
#' @param list_sClass a list of "sClass" objects or an object of class "sClass"
#' @param displayBy which statistics will be used for displaying. It can be "ROC" for AUC in ROC (shaped by direction), "Fmax" for F-max in Precision-Recall curve
#' @param type the plot type. It can be "bar" and "ridge"
#' @param sort logical to indicate whether to sort methods according to performance. By default, it sets TRUE
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param ... additional graphic parameters for ggridges::geom_density_ridges_gradient (such as scale=1)
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xClassifyCV}}
#' @include xClassifyCV.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xClassifyCV(list_sClass, displayBy="ROC", type="bar")
#' gp <- xClassifyCV(list_sClass, displayBy="ROC", type="ridge")
#' gp <- xClassifyCV(sClass)
#' }

xClassifyCV <-function(list_sClass, displayBy=c("ROC","Fmax"), type=c("bar","ridge"), sort=TRUE, colormap="spectral", ncolors=64, ...)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    displayBy <- match.arg(displayBy)
    type <- match.arg(type)

    if(class(list_sClass) %in% "sClass"){
    	list_sClass <- list(list_sClass)
    }else if(class(list_sClass)=="list"){
		## Remove null elements in a list
		list_sClass <- base::Filter(base::Negate(is.null), list_sClass)
		if(length(list_sClass)==0){
			warnings("All sClass objects are NULL!")
			return(NULL)
		}
    }else{
    	return(NULL)
    }
	
	## check list names
	list_names <- names(list_sClass)
	if(is.null(list_names)){
		list_names <- paste('Method', 1:length(list_sClass), sep=' ')
	}
	
	if(displayBy=="ROC"){
		ls_df <- lapply(1:length(list_sClass), function(i){
			x <- list_sClass[[i]]$cv_auroc[1,-c(1:4)]
			data.frame(Method=list_names[i], Val=x, stringsAsFactors=FALSE)
		})
		df_data <- do.call(rbind, ls_df)
		rownames(df_data) <- NULL
		if(is.na(xlab)){
			xlab <- "AUC: 95% Confidence Interval\n(Repeated Cross-Validation)"
		}
	}else if(displayBy=="Fmax"){
		ls_df <- lapply(1:length(list_sClass), function(i){
			x <- list_sClass[[i]]$cv_fmax[1,-c(1:4)]
			data.frame(Method=list_names[i], Val=x, stringsAsFactors=FALSE)
		})
		df_data <- do.call(rbind, ls_df)
		rownames(df_data) <- NULL
		if(is.na(xlab)){
			xlab <- "F-max: 95% Confidence Interval\n(Repeated Cross-Validation)"
		}
	}
	
	if(type=='bar'){
		ls_vec <- split(x=df_data$Val, f=df_data$Method)
		res_ttest <- lapply(ls_vec, stats::t.test)
		ls_res <- lapply(1:length(res_ttest), function(i){
			x <- res_ttest[[i]]
			data.frame(Method=names(res_ttest)[i], mean=x$estimate, conf_lower=x$conf.int[1], conf_upper=x$conf.int[2], stringsAsFactors=F)
		})
		df <- do.call(rbind, ls_res)
		rownames(df) <- NULL
		
		Method <- mean <- conf_upper <- conf_lower <- NULL
		if(sort){
			df <- df %>% dplyr::arrange(-mean)
		}
		df$Method <- factor(df$Method, levels=rev(unique(df$Method)))
		
		p <- ggplot(df, aes(x=mean, y=Method))
		p <- p + geom_point() + geom_errorbarh(aes(xmax=conf_upper, xmin=conf_lower, height=.2))		
		p <- p  + theme_bw() + theme(legend.position="none", axis.title.y=element_blank(), panel.background=element_rect(fill=rgb(0.95,0.95,0.95,1)))
		p <- p + xlab(xlab)
		p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		## put arrows on x-axis
		p <- p + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))

		## x-axis position
		if(displayBy == "ROC"){
			p <- p + scale_x_continuous(position="top", limits=c(0.5,1))
		}else{
			p <- p + scale_x_continuous(position="top")
		}
		
	}else if(type=="ridge"){
		
		Method <- Val <- mean <- ..density.. <- NULL
		if(sort){
			methods <- as.data.frame(df_data %>% dplyr::group_by(Method) %>% dplyr::summarise(mean=mean(Val)) %>% dplyr::arrange(mean))$Method
		}else{
			methods <- unique(df_data$Method)
		}
		df_data$Method <- factor(df_data$Method, levels=methods)
		
		p <- ggplot(df_data) + ggridges::geom_density_ridges_gradient(aes(x=Val, y=Method, group=Method, fill=..density..), lwd=0.2, jittered_points=T, position=ggridges::position_points_jitter(width=0.05,height=0), point_shape='|', point_size=1, point_alpha=0.5) + scale_fill_gradientn(colors=xColormap(colormap)(ncolors), guide=F)
		p <- p  + theme_bw() + theme(legend.position="none", axis.title.y=element_blank(), panel.background=element_rect(fill=rgb(0.95,0.95,0.95,1)))
		p <- p + xlab(xlab)
		p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		## put arrows on x-axis
		p <- p + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
		
		## x-axis position
		if(displayBy == "ROC"){
			p <- p + scale_x_continuous(position="top", limits=c(0.5,1))
		}else{
			p <- p + scale_x_continuous(position="top")
		}
	}
	
	invisible(p)
}

    
