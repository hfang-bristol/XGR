#' Function to visualise enrichment results using a barplot
#'
#' \code{xEnrichBarplot} is supposed to visualise enrichment results using a barplot. It returns an object of class "ggplot".
#'
#' @param eTerm an object of class "eTerm"
#' @param top_num the number of the top terms (sorted according to FDR or adjusted p-values). If it is 'auto', only the significant terms (FDR < 0.05) will be displayed
#' @param displayBy which statistics will be used for displaying. It can be "fc" for enrichment fold change (by default), "adjp" or "fdr" for adjusted p value (or FDR), "pvalue" for p value, "zscore" for enrichment z-score
#' @param bar.label logical to indicate whether to label each bar with FDR. By default, it sets to true for bar labelling
#' @param bar.label.size an integer specifying the bar labelling text size. By default, it sets to 3
#' @param wrap.width a positive integer specifying wrap width of name
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichBarplot.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location="~/Sites/SVN/github/bigdata"
#' 
#' # 1) load eQTL mapping results: cis-eQTLs significantly induced by IFN
#' cis <- xRDataLoader(RData.customised='JKscience_TS2A', RData.location=RData.location)
#' ind <- which(cis$IFN_t > 0 & cis$IFN_fdr < 0.05)
#' df_cis <- cis[ind, c('variant','Symbol','IFN_t','IFN_fdr')]
#' data <- df_cis$variant
#' 
#' # 2) Enrichment analysis using Experimental Factor Ontology (EFO)
#' # Considering LD SNPs and respecting ontology tree
#' eTerm <- xEnricherSNPs(data, ontology="EF", include.LD="EUR", LD.r2=0.8, ontology.algorithm="lea", RData.location=RData.location)
#'
#' # 3) Barplot of enrichment results
#' bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="fc")
#' #pdf(file="enrichment_barplot.pdf", height=6, width=12, compress=TRUE)
#' print(bp)
#' #dev.off()
#' ## modify y axis text
#' bp + theme(axis.text.y=element_text(size=10,color="blue"))
#' ## modify x axis title
#' bp + theme(axis.title.x=element_text(color="blue"))
#' }

xEnrichBarplot <- function(eTerm, top_num=10, displayBy=c("fc","adjp","fdr","zscore","pvalue"), bar.label=TRUE, bar.label.size=3, wrap.width=NULL) 
{
    
    displayBy <- match.arg(displayBy)
    
    if(is.logical(eTerm)){
        stop("There is no enrichment in the 'eTerm' object.\n")
    }
    
    ## when 'auto', will keep the significant terms
	df <- xEnrichViewer(eTerm, top_num="all")
	if(top_num=='auto'){
		top_num <- sum(df$adjp<0.05)
		if(top_num<=1){
			top_num <- sum(df$adjp<0.1)
		}
		if(top_num <= 1){
			top_num <- 10
		}
	}
	df <- xEnrichViewer(eTerm, top_num=top_num, sortBy="adjp")
	
	## text wrap
	if(!is.null(wrap.width)){
		width <- as.integer(wrap.width)
		res_list <- lapply(df$name, function(x){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=width)
			if(length(y)>1){
				paste0(y[1], '...')
			}else{
				y
			}
		})
		df$name <- unlist(res_list)
	}
	
	if(displayBy=='adjp' | displayBy=='fdr'){
		df <- df[with(df,order(-adjp,zscore)),]
		df$name <- factor(df$name, levels=df$name)
		p <- ggplot(df, aes(x=df$name, y=-1*log10(df$adjp))) 
		p <- p + ylab("Enrichment significance: -log10(FDR)")
	}else if(displayBy=='fc'){
		df <- df[with(df,order(fc,-adjp)),]
		df$name <- factor(df$name, levels=df$name)
		p <- ggplot(df, aes(x=df$name, y=df$fc))
		p <- p + ylab("Enrichment changes")
	}else if(displayBy=='pvalue'){
		df <- df[with(df,order(-pvalue,zscore)),]
		df$name <- factor(df$name, levels=df$name)
		p <- ggplot(df, aes(x=df$name, y=-1*log10(df$pvalue)))
		p <- p + ylab("Enrichment significance: -log10(p-value)")
	}else if(displayBy=='zscore'){
		df <- df[with(df,order(zscore,-adjp)),]
		df$name <- factor(df$name, levels=df$name)
		p <- ggplot(df, aes(x=df$name, y=df$zscore))
		p <- p + ylab("Enrichment z-scores")
	}
	
	bp <- p + geom_bar(stat="identity", fill="orange") + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_text(size=12,color="black"), axis.title.x=element_text(size=14,color="black")) + coord_flip()
	
	if(bar.label){
		## get text label
		to_scientific_notation <- function(x) {
			res <- format(x, scientific=T)
			res <- sub("\\+0?", "", res)
			sub("-0?", "-", res)
		}
		label <- to_scientific_notation(df$adjp)
		label <- paste('FDR', as.character(label), sep='=')
		
		bp <- bp + geom_text(aes(label=label),hjust=1,size=bar.label.size)
	}
	
	invisible(bp)
}
