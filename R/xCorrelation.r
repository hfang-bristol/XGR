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
#' @return 
#' a list with two componets:
#' \itemize{
#'  \item{\code{df_summary}: a data frame of n x 4, where n is the number of named vectors, and the 4 columns are "name", "cor" (i.e. "correlation"), "pval" (i.e. p-value), "fdr"}
#'  \item{\code{ls_gp}: NULL if the plot is not drawn; otherwise, a list of 'ggplot' objects}
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
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' # a) provide the seed nodes/genes with the weight info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get genes within 500kb away from AS GWAS lead SNPs
#' seeds.genes <- ImmunoBase$AS$genes_variants
#' ## seeds weighted according to distance away from lead SNPs
#' data <- 1- seeds.genes/500000
#'
#' # b) perform priority analysis
#' df <- data.frame(name=names(data), value=data, stringsAsFactors=FALSE)
#' 
#' # c) do correlation
#' ls_res <- xCorrelation(df, data, method="pearson", p.type="empirical", nperm=2000, plot=TRUE)
#' }

xCorrelation <- function(df, list_vec, method=c("pearson","spearman"), p.type=c("nominal","empirical"), seed=825, nperm=2000, p.adjust.method=c("BH","BY","bonferroni","holm","hochberg","hommel"), plot=FALSE)
{
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    method <- match.arg(method)
    p.type <- match.arg(p.type)
    p.adjust.method <- match.arg(p.adjust.method)
    
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
    if(is.vector(list_vec)){
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
    
	## Combine into a data frame called 'df_disease'
	list_names <- names(list_vec)
	if(is.null(list_names)){
		list_names <- paste0('V', 1:length(list_vec))
		names(list_vec) <- list_names
	}
    
    ls_df <- lapply(1:length(list_vec), function(i){
		
		data <- list_vec[[i]]
		df_priority_data <- df_priority
		ind <- match(df_priority_data$name, names(data))
		df_priority_data$data <- data[ind]
		df <- subset(df_priority_data, !is.na(data))
	
		##############
		res <- stats::cor.test(x=df$priority, y=as.numeric(df$data), method=method, exact=FALSE)
		cor_obs <- signif(res$estimate, 3)
		##############
	
		if(p.type == 'nominal'){
			pval_obs <- res$p.value
		}else if(p.type == 'empirical'){
			B <- nperm
			set.seed(seed)
			vec_p <- sapply(1:B, function(i){
				df$priority <- sample(df_priority_data$priority, nrow(df))
				cor_exp <- stats::cor(x=df$priority, y=df$data, method=method)
			})
			pval_obs <- sum(abs(vec_p) > abs(cor_obs))/B
		}
	
		if(pval_obs < 0.05){
			pval_obs <- format(signif(res$p.value,2),scientific=TRUE)
		}else{
			pval_obs <- signif(res$p.value,3)
		}
	
		data.frame(name=names(list_vec)[i], cor=cor_obs, pval=pval_obs, stringsAsFactors=FALSE)
    })
    dff <- do.call(rbind, ls_df)
    fdr <- stats::p.adjust(dff$pval, method=p.adjust.method)
    dff$fdr <- ifelse(fdr<0.05, format(signif(fdr,2),scientific=TRUE), signif(fdr,3))
    rownames(dff) <- 1:nrow(dff)
    
    if(plot){
		ls_gp_curve <- lapply(1:length(list_vec), function(i){
		
			data <- list_vec[[i]]
			df_priority_data <- df_priority
			ind <- match(df_priority_data$name, names(data))
			df_priority_data$data <- data[ind]
			df <- subset(df_priority_data, !is.na(data))
			
			name_obs <- dff$name[i]
			cor_obs <- dff$cor[i]
			pval_obs <- dff$pval[i]
			fdr_obs <- dff$fdr[i]
			
			priority <- name <- NULL
			m <- ggplot(df, aes(x=priority, y=data))
			m <- m + geom_point()
			m <- m + geom_smooth(method=c("lm","loess")[1], se=TRUE, span=4)
			m <- m + theme_bw() + theme(legend.position="top", axis.title.y=element_text(size=12,color="black"), axis.text.y=element_text(size=8,color="black"), axis.title.x=element_text(size=12,color="black"), axis.text.x=element_text(size=8,color="black"), panel.background=element_rect(fill="transparent"))
			m <- m + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
			gp_curve <- m + labs(x="data frame", y=name_obs, title=paste0("Correlation (",method,")"), subtitle=paste0("correlation: ",cor_obs,', ',p.type,' p-value: ',pval_obs,', fdr: ',fdr_obs)) + theme(plot.title=element_text(hjust=0.5, size=12), plot.subtitle=element_text(hjust=0.5, size=10))
		
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
    			  ls_gp = ls_gp_curve
                 )
    
    invisible(ls_res)
}
