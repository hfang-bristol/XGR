#' Function to visualise enrichment results using a barplot
#'
#' \code{xEnrichBarplot} is supposed to visualise enrichment results using a barplot. It returns an object of class "ggplot".
#'
#' @param eTerm an object of class "eTerm"
#' @param top_num the number of the top terms (sorted according to 'displayBy' below). If it is 'auto', only the signficant terms will be displayed
#' @param displayBy which statistics will be used for displaying. It can be "adjp" for adjusted p value, "pvalue" for p value, "zscore" for enrichment z-score, "fc" for enrichment fold change
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichBarplot.r
#' @examples
#' \dontrun{
#' xEnrichBarplot(eTerm)
#' }

xEnrichBarplot <- function(eTerm, top_num=10, displayBy=c("fc","adjp","zscore","pvalue")) 
{
    
    displayBy <- match.arg(displayBy)
    
    if(is.logical(eTerm)){
        stop("There is no enrichment in the 'eTerm' object.\n")
    }
    
    ## when 'auto', will keep the significant terms
	df <- xEnrichViewer(eTerm, top_num="all", sortBy=displayBy)
	if(top_num=='auto'){
		top_num <- sum(df$adjp<0.05)
		if(top_num<=1){
			top_num <- sum(df$adjp<0.1)
		}
		if(top_num <= 1){
			top_num <- 10
		}
	}
	df <- xEnrichViewer(eTerm, top_num=top_num, sortBy=displayBy)
	
	## get text label
	to_scientific_notation <- function(x) {
	  	res <- format(x, scientific=T)
	  	res <- sub("\\+0?", "", res)
	  	sub("-0?", "-", res)
	}
    label <- to_scientific_notation(df$adjp)
	label <- paste('FDR', as.character(label), sep='=')
	
	df$name <- factor(df$name, levels=df$name)
	if(displayBy=='adjp'){
		p <- ggplot(df, aes(x=df$name, y=-1*log10(df$adjp))) 
		p <- p + ylab("Enrichment signficance: -log10(FDR)") + geom_hline(yintercept=-log10(0.05),color="grey")
	}else if(displayBy=='fc'){
		p <- ggplot(df, aes(x=df$name, y=df$fc))
		p <- p + ylab("Enrichment changes") + geom_hline(yintercept=1,color="grey")
	}else if(displayBy=='pvalue'){
		p <- ggplot(df, aes(x=df$name, y=-1*log10(df$pvalue)))
		p <- p + ylab("Enrichment signficant: -log10(p-value)")
	}else if(displayBy=='zscore'){
		p <- ggplot(df, aes(x=df$name, y=df$zscore))
		p <- p + ylab("Enrichment z-scores") + geom_hline(yintercept=0,color="grey")
	}
	
	bp <- p + geom_bar(stat="identity", fill="deepskyblue") + geom_text(aes(label=label),hjust=1,size=3) + theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_text(size=12,color="blue"), axis.title.x=element_text(size=14,color="blue")) + coord_flip()
	
	invisible(bp)
}
