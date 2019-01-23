#' Function to visualise eQTL analysis using violin plot
#'
#' \code{xMEplot} is supposed to visualise eQTL analysis using violin plot.
#
#' @param gene an integer specifying ArrayAddress or a character specifying a gene symbol (possibly several ArrayAddress queried returning a list of ggplot objects)
#' @param snp a dbSNP
#' @param summary a data frame storing eQTL summary statistics. It must contain columns 'context', 'gene', 'snps', 'FDR', 'statistic'
#' @param expression a data matrix/frame with genes in rows and samples in columns
#' @param genotype a data (sparse) matrix/frame with SNPs in rows and samples in columns
#' @param contexts a vector specifying contexts included
#' @param colormap the colormap ("sci_jco")
#' @param ncolumns an integer specifying the number of columns for contexts. By defaul, it is NULL (decided on according to the number of contexts that will be visualised)
#' @return 
#' a ggplot object (or a list of ggplot objects)
#' @note none
#' @export
#' @seealso \code{\link{xME}}
#' @include xMEplot.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' summary <- xRDataLoader('JK_cohort_xMEdb', RData.location=RData.location)
#' genotype <- xRDataLoader('JK_cohort_genotype', RData.location=RData.location)
#' expression <- xRDataLoader('JK_cohort_expression_cis', RData.location=RData.location)
#' gp <- xMEplot(gene="MED7", snp='rs13182119', summary, expression, genotype, contexts=c("CD14","LPS2","LPS24","IFN","Bcell","NK","Neutrophil","Monocyte"))
#' }

xMEplot <- function(gene=580484, snp='rs13182119', summary, expression, genotype, contexts=c("CD14","LPS2","LPS24","IFN","Bcell","NK","Neutrophil","Monocyte","shared_CD14","shared_IFN","shared_LPS2","shared_LPS24"), colormap="sci_jco", ncolumns=NULL)
{
	
	## df_summary_all
	if(class(summary)=='data.frame'){
		
		if(class(gene)=='numeric'){
			ind <- which(!is.na(match(summary$gene, gene)) & !is.na(match(summary$snps, snp)))
			df_summary_all <- summary[ind,]
		}else{
			ind <- which(!is.na(match(summary$Symbol, gene)) & !is.na(match(summary$snps, snp)))
			df_summary_all <- summary[ind,]
		}
		
		if(nrow(df_summary_all)==0){
			return(NULL)
		}
	}else{
		return(NULL)
	}
	## vec_genotype
	vec_genotype <- c(0,1,2)
	effect <- unique(df_summary_all$effect_allele)
	other <- unique(df_summary_all$other_allele)
	names(vec_genotype) <- c(paste0(rep(other,2),collapse=''), paste0(c(other,effect),collapse=''), paste0(rep(effect,2),collapse=''))
	
	## df_genotype
	ind <- match(rownames(genotype), snp)
	vec <- genotype[!is.na(ind), ]
	df_genotype <- data.frame(individual=names(vec), dosage=vec, stringsAsFactors=F)
	### append genotype
	ind <- match(df_genotype$dosage, vec_genotype)
	df_genotype$genotype <- names(vec_genotype)[ind]
	
	##################################################
	ls_gene <- split(x=df_summary_all, f=df_summary_all$gene)
	
	ls_gp <- lapply(1:length(ls_gene), function(i){
		gene <- as.numeric(names(ls_gene)[i])
		df_summary <- ls_gene[[i]]
	
		## df_expression
		ind <- match(expression$gene, gene)
		df_expression <- expression[!is.na(ind),]
	
		## df_output
		df_output <- df_expression
		### dosage and genotype
		ind <- match(df_expression$individual, df_genotype$individual)
		df_output$dosage <- df_genotype$dosage[ind]
		df_output$genotype <- df_genotype$genotype[ind]
		### remove dosage='NA'
		ind <- !is.na(df_output$dosage)
		df_output <- df_output[ind,]
		### FDR
		df_output$FDR <- 'ns'
		ind <- match(df_output$context, df_summary$context)
		df_output$FDR[!is.na(ind)] <- format(signif(df_summary$FDR[ind[!is.na(ind)]],digits=2), scientific=T)
		### statistic
		df_output$statistic <- 0
		ind <- match(df_output$context, df_summary$context)
		df_output$statistic[!is.na(ind)] <- signif(df_summary$statistic[ind[!is.na(ind)]], digits=3)
		df_output$direction <- ifelse(df_output$FDR=='ns', 'ns', ifelse(df_output$statistic>0, 'increasing', 'decreasing'))
		df_output$direction <- factor(df_output$direction, c('decreasing','increasing','ns'))
	
		### strip
		if(is.null(contexts)){
			vec_context <- unique(df_output$context)
		}else{
			vec_context <- contexts
		}
		#### keep only contexts specified
		ind <- match(df_output$context, vec_context)
		df_output <- df_output[!is.na(ind),]
		#### order by contexts specified
		id <- NULL
		df_output$id <- as.numeric(factor(df_output$context,vec_context))
		df_output <- df_output %>% dplyr::arrange(id)
		#### strip labelling
		df_output$strip <- ifelse(df_output$FDR=='ns', df_output$context, paste0(df_output$context, ' (FDR ', df_output$FDR, ')'))
		df_output$strip <- factor(df_output$strip, unique(df_output$strip))

		## ggviolin plot
		strip <- NULL
		df_output$genotype <- factor(df_output$genotype, names(vec_genotype))
		ylab <- paste0("Expression of ", unique(df_summary$Symbol), " (", unique(df_summary$gene), ")")
		xlab <- unique(df_summary$snps)
		#gp <- ggpubr::ggboxplot(df_output, x="genotype", y="val", xlab=xlab, ylab=ylab, color="genotype", palette ="jco", add=c("jitter","dotplot")[2], add.params=list(binwidth=0.1, dotsize=0.2))
		gp <- ggpubr::ggviolin(df_output, x="genotype", y="val", xlab=xlab, ylab=ylab, color="direction", palette="jco", add=c("median_iqr","boxplot")[1], add.params=list(size=0.2))
		
		if(is.null(ncolumns)){
			ncolumns <- ceiling(sqrt(length(levels(df_output$strip))))
		}		
		gp <- gp + facet_wrap(~strip,ncol=ncolumns)
		
		values <- xColormap(colormap)(length(levels(df_output$direction)))
		names(values) <- levels(df_output$direction)
		gp <- gp + scale_colour_manual(values=values, guide=guide_legend(title="Direction",keywidth=0.7,keyheight=0.7))
		gp <- gp + theme(axis.text.y=element_text(size=7,color="black"), 
	axis.title.y=element_text(size=9,color="black"), 
	axis.text.x=element_text(size=7,color="black"), 
	axis.title.x=element_text(size=9,color="black"), panel.background=element_rect(fill="grey99")) + theme(strip.background=element_rect(fill="grey95",color="black"), strip.text=element_text(size=8))
	})
	
	if(length(ls_gp)==1){
		gp <- ls_gp[[1]]
	}else{
		gp <- ls_gp
	}
	
    invisible(gp)
}
