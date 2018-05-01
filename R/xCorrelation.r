#' Function to calculate and visualise correlation
#'
#' \code{xCorrelation} is supposed to calculate and visualise correlation between a data frame and a named vector (or a list of named vectors). 
#
#' @param df a data frame with two columns ('name' and 'value')
#' @param list_vec a named vector containing numeric values. Alternatively it can be a list of named vectors
#' @param method the method used to calcualte correlation. It can be 'pearson' for Pearson's correlation or 'spearman' for Spearman rank correlation
#' @param p.type the type of the p-value calcualted. It can be 'nominal' for nominal p-value or 'empirical' for empirical p-value
#' @param seed an integer specifying the seed
#' @param nperm the number of random permutations
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param plot logical to indicate whether scatter plot is drawn
#' @param plot.smooth the smooth method for the scatter plot. It can be NA (depending on correlation type), "lm" for the linear line or 'loess' for the loess curve
#' @return 
#' a list with three componets:
#' \itemize{
#'  \item{\code{df_summary}: a data frame of n x 5, where n is the number of named vectors, and the 5 columns are "name", "num" (i.e. number of data points used for calculation), "cor" (i.e. correlation), "pval" (i.e. p-value), "fdr"}
#'  \item{\code{ls_gp_curve}: NULL if the plot is not drawn; otherwise, a list of 'ggplot' objects for scatter plot together with an estimated curve}
#'  \item{\code{ls_gp_pdf}: NULL if the plot is not drawn; otherwise, a list of 'ggplot' objects for pdf plot for null distribution of correlation together with a vertical line for observed correlation}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xCorrelation}}
#' @include xCorrelation.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
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
#' ls_res <- xCorrelation(df, data, method="pearson", p.type="empirical", nperm=2000, plot=TRUE)
#' }

xCorrelation <- function(df, list_vec, method=c("pearson","spearman"), p.type=c("nominal","empirical"), seed=825, nperm=2000, p.adjust.method=c("BH","BY","bonferroni","holm","hochberg","hommel"), plot=FALSE, plot.smooth=c(NA, "lm","loess"))
{
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    method <- match.arg(method)
    p.type <- match.arg(p.type)
    p.adjust.method <- match.arg(p.adjust.method)
    plot.smooth <- match.arg(plot.smooth)
    
    if(class(df) == "data.frame"){
    	df_priority <- df[,c(1:2)]
    	colnames(df_priority) <- c("name","priority")
    }else{
    	stop("The function must apply to a 'data frame' object.\n")
    }
    
    ############
    if(length(list_vec)==0){
    	return(NULL)
    }
    ############
    if(is.vector(list_vec) & class(list_vec)!="list"){
    	# assume a vector
		if(is.null(names(list_vec))){
			stop("The input vector must have names.\n")
		}else{
			list_vec <- list(list_vec)
		}
	}else if(class(list_vec)=="list"){
		## Remove null elements in a list
		list_vec <- base::Filter(base::Negate(is.null), list_vec)
		if(length(list_vec)==0){
			return(NULL)
		}
    }else{
        stop("The input data must be a named vector or a list of named vectors.\n")
    }
    
	list_names <- names(list_vec)
	if(is.null(list_names)){
		list_names <- paste0('V', 1:length(list_vec))
		names(list_vec) <- list_names
	}
    
    ls_df_gp <- lapply(1:length(list_vec), function(i){
		
		data <- list_vec[[i]]
		df_priority_data <- df_priority
		ind <- match(df_priority_data$name, names(data))
		df_priority_data$data <- data[ind]
		df <- subset(df_priority_data, !is.na(data))
		
		################
		if(nrow(df)<=5){
			return(NULL)
		}
		################
				
		##############
		res <- stats::cor.test(x=df$priority, y=as.numeric(df$data), method=method, exact=FALSE)
		cor_obs <- signif(res$estimate, 3)
		##############
	
		if(p.type == 'nominal'){
			pval_obs <- res$p.value
			
			if(pval_obs < 0.05){
				pval_obs <- as.numeric(format(signif(res$p.value,2),scientific=TRUE))
			}else{
				pval_obs <- signif(res$p.value,3)
			}
			
			gp_pdf <- NULL
			
		}else if(p.type == 'empirical'){
			B <- nperm
			set.seed(seed)
			vec_exp <- sapply(1:B, function(i){
				df$priority <- sample(df_priority_data$priority, nrow(df))
				cor_exp <- stats::cor(x=df$priority, y=df$data, method=method)
			})
			pval_obs <- sum(abs(vec_exp) > abs(cor_obs))/B
			
			if(pval_obs < 0.05){
				pval_obs <- as.numeric(format(signif(res$p.value,2),scientific=TRUE))
			}else{
				pval_obs <- signif(res$p.value,3)
			}
			
			##################################
			## gp_pdf
			cor <- NULL
			gp <- ggplot(data.frame(cor=vec_exp), aes(cor))
			gp <- gp + geom_density(adjust=0.8,fill='grey',colour='grey',alpha=0.5)
			gp <- gp + geom_vline(xintercept=cor_obs,color='red')

			gp <- gp + scale_x_continuous(limits=c(-1,1))
			gp <- gp + theme_bw() + theme(axis.title.y=element_text(size=10,color="black"), axis.text.y=element_text(size=8,color="black"), axis.title.x=element_text(size=10,color="black"), axis.text.x=element_text(size=8,color="black"), panel.background=element_rect(fill="transparent"))
			gp <- gp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
			gp_pdf <- gp + labs(x=paste0("Correlation (",method,")"), y=paste0("Probability density (null distribution)\nestimated based on ",B," permutations"), title=names(list_vec)[i], subtitle=paste0("correlation: ",cor_obs,' (empirical p: ',pval_obs,')')) + theme(plot.title=element_text(hjust=0.5, size=10), plot.subtitle=element_text(hjust=0.5, size=8))
			##################################
			
		}

		df <- data.frame(name=names(list_vec)[i], num=nrow(df), cor=cor_obs, pval=pval_obs, stringsAsFactors=FALSE)
		list(df=df, gp=gp_pdf)
    })
    
    ## ls_gp_pdf
    ls_gp_pdf <- lapply(ls_df_gp, function(x){
    	if(!is.null(x)){
    		x$gp
    	}
    })
    names(ls_gp_pdf) <- names(list_vec)
    ## dff
    ls_df <- lapply(ls_df_gp, function(x){
    	if(!is.null(x)){
    		x$df
    	}
    })
    dff <- do.call(rbind, ls_df)
    
    ###########
    if(is.null(dff)){
    	return(NULL)
    }
    ###########    
    fdr <- stats::p.adjust(dff$pval, method=p.adjust.method)
    dff$fdr <- ifelse(fdr<0.05, as.numeric(format(signif(fdr,2),scientific=TRUE)), signif(fdr,3))
    rownames(dff) <- 1:nrow(dff)
    
    if(plot){
		ls_gp_curve <- lapply(1:length(list_vec), function(i){
		
			data <- list_vec[[i]]
			df_priority_data <- df_priority
			ind <- match(df_priority_data$name, names(data))
			df_priority_data$data <- data[ind]
			df <- subset(df_priority_data, !is.na(data))
			
			ind <- match(names(list_vec)[i], dff$name)
			if(is.na(ind)){
				return(NULL)
			}else{
				name_obs <- dff$name[ind]
				cor_obs <- dff$cor[ind]
				pval_obs <- dff$pval[ind]
				fdr_obs <- dff$fdr[ind]
			}
			
			priority <- name <- NULL
			m <- ggplot(df, aes(x=priority, y=data))
			m <- m + geom_point()
			if(is.na(plot.smooth)){
				if(method=="pearson"){
					plot.smooth <- "lm"
				}else if(method=="spearman"){
					plot.smooth <- "loess"
				}
			}
			m <- m + geom_smooth(method=plot.smooth, se=TRUE, span=4)
			m <- m + theme_bw() + theme(legend.position="top", axis.title.y=element_text(size=10,color="black"), axis.text.y=element_text(size=8,color="black"), axis.title.x=element_text(size=10,color="black"), axis.text.x=element_text(size=8,color="black"), panel.background=element_rect(fill="transparent"))
			m <- m + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
			subtitle <- paste0("correlation: ",cor_obs,' (',p.type,' p: ',pval_obs,', fdr: ',fdr_obs,')')
			if(length(list_vec)==1){
				subtitle <- paste0("correlation: ",cor_obs,', (',p.type,' p: ',pval_obs,')')
			}
			gp_curve <- m + labs(x="data frame", y=name_obs, title=paste0("Correlation (",method,"; n=",nrow(df),")"), subtitle=subtitle) + theme(plot.title=element_text(hjust=0.5, size=10), plot.subtitle=element_text(hjust=0.5, size=8))
		
			if(0){
			gp_curve <- gp_curve + scale_x_continuous(limits=c(0,ceiling(max(df$priority)*10)/10)) 
			gp_curve <- gp_curve + scale_y_reverse(limits=c(ceiling(max(df$data)*10)/10, floor(min(df$data)*10)/10))
			gp_curve <- gp_curve + ggrepel::geom_text_repel(aes(label=name), size=2, fontface='bold', box.padding=unit(0.2,"lines"), point.padding=unit(0.2,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')), force=0.1)
			}
			return(gp_curve)
    	})
    	names(ls_gp_curve) <- names(list_vec)
    }else{
    	ls_gp_curve <- NULL
    }
    
    ls_res <- list(df_summary = dff,
    			  ls_gp_curve = ls_gp_curve,
    			  ls_gp_pdf	= ls_gp_pdf
                 )
    
    invisible(ls_res)
}
