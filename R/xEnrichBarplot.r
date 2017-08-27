#' Function to visualise enrichment results using a barplot
#'
#' \code{xEnrichBarplot} is supposed to visualise enrichment results using a barplot. It returns an object of class "ggplot".
#'
#' @param eTerm an object of class "eTerm"
#' @param top_num the number of the top terms (sorted according to FDR or adjusted p-values). If it is 'auto', only the significant terms (see below FDR.cutoff) will be displayed
#' @param displayBy which statistics will be used for displaying. It can be "fc" for enrichment fold change (by default), "adjp" or "fdr" for adjusted p value (or FDR), "pvalue" for p value, "zscore" for enrichment z-score
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05. This option only works when setting top_num (see above) is 'auto'
#' @param bar.label logical to indicate whether to label each bar with FDR. By default, it sets to true for bar labelling
#' @param bar.label.size an integer specifying the bar labelling text size. By default, it sets to 3
#' @param bar.color either NULL or fill color names ('lightyellow-orange' by default)
#' @param bar.width bar width. By default, 80% of the resolution of the data
#' @param wrap.width a positive integer specifying wrap width of name
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichBarplot.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
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
#' }
#' 
#' # 4) use font family (Arial)
#' \dontrun{
#' source("http://bioconductor.org/biocLite.R"); biocLite("extrafont")
#' library(extrafont)
#' font_import()
#' fonttable()
#' ## creating PDF files with fonts
#' library(extrafont)
#' loadfonts()
#' bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="fc", font.family="Arial Black")
#' pdf(file="enrichment_barplot_fonts.pdf", height=6, width=12, family="Arial Black")
#' print(bp)
#' dev.off()
#' }

xEnrichBarplot <- function(eTerm, top_num=10, displayBy=c("fc","adjp","fdr","zscore","pvalue"), FDR.cutoff=0.05, bar.label=TRUE, bar.label.size=3, bar.color='lightyellow-orange', bar.width=0.8, wrap.width=NULL, font.family="sans", signature=TRUE) 
{
    
    displayBy <- match.arg(displayBy)
    
    if(is.null(eTerm)){
        warnings("There is no enrichment in the 'eTerm' object.\n")
        return(NULL)
    }
    
    ## when 'auto', will keep the significant terms
	df <- xEnrichViewer(eTerm, top_num="all")
	if(top_num=='auto'){
		top_num <- sum(df$adjp<FDR.cutoff)
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
	
	name <- height <- direction <- hjust <- NULL
	adjp <- zscore <- pvalue <- fc <- NULL
	
	##########
	## consider the direction of z-score
	df <- df %>% dplyr::mutate(direction=ifelse(zscore>0,1,-1))
	##########	
	
	if(displayBy=='adjp' | displayBy=='fdr'){
		df <- df %>% dplyr::arrange(direction, desc(adjp), zscore) %>% dplyr::mutate(height=-1*log10(adjp)) %>% dplyr::mutate(hjust=1)
		df$name <- factor(df$name, levels=df$name)
		####
		df$height[is.infinite(df$height)] <- max(df$height[!is.infinite(df$height)])
		####
		p <- ggplot(df, aes(x=name, y=height))
		p <- p + ylab(expression(paste("Enrichment significance: ", -log[10]("FDR"))))
		
	}else if(displayBy=='fc'){
		df <- df %>% dplyr::arrange(direction, fc, desc(adjp)) %>% dplyr::mutate(height=log2(fc)) %>% dplyr::mutate(hjust=ifelse(height>=0,1,0))
		df$name <- factor(df$name, levels=df$name)
		p <- ggplot(df, aes(x=name, y=height))
		p <- p + ylab(expression(paste("Enrichment changes: ", log[2]("FC"))))
		
	}else if(displayBy=='pvalue'){
		df <- df %>% dplyr::arrange(direction, desc(pvalue), zscore) %>%  dplyr::mutate(height=-1*log10(pvalue)) %>% dplyr::mutate(hjust=1)
		df$name <- factor(df$name, levels=df$name)
		####
		df$height[is.infinite(df$height)] <- max(df$height[!is.infinite(df$height)])
		####
		p <- ggplot(df, aes(x=name, y=height))
		p <- p + ylab(expression(paste("Enrichment significance: ", -log[10]("p-value"))))

	}else if(displayBy=='zscore'){
		df <- df %>% dplyr::arrange(direction, zscore, desc(adjp)) %>% dplyr::mutate(height=zscore) %>% dplyr::mutate(hjust=ifelse(height>=0,1,0))
		df$name <- factor(df$name, levels=df$name)
		p <- ggplot(df, aes(x=name, y=height))
		p <- p + ylab("Enrichment z-scores")

	}
	
	if(is.null(bar.color)){
		bp <- p + geom_col(color='grey80',fill='transparent', width=bar.width)
	}else{
		bar.color <- unlist(strsplit(bar.color, "-"))
		if(length(bar.color)==2){
			bp <- p + geom_col(aes(fill=height), width=bar.width) + scale_fill_gradient(low=bar.color[1],high=bar.color[2]) 
		}else{
			bp <- p + geom_col(color='grey80',fill='transparent', width=bar.width)
		}
	}
	#bp <- p + geom_col(aes(fill=height)) + scale_fill_gradient2(low="cyan", mid="grey", high="orange", midpoint=0)
	bp <- bp + theme_bw() + theme(legend.position="none",axis.title.y=element_blank(), axis.text.y=element_text(size=10,color="black"), axis.title.x=element_text(size=12,color="black"))
	
	bp <- bp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	if(bar.label){
		## get text label
		to_scientific_notation <- function(x) {
			res <- format(x, scientific=T)
			res <- sub("\\+0?", "", res)
			sub("-0?", "-", res)
		}
		label <- to_scientific_notation(df$adjp)
		df$label <- paste('FDR', as.character(label), sep='=')
		
		## hjust==1
		df_sub <- df %>% dplyr::filter(hjust==1)
		bp <- bp + geom_text(data=df_sub, aes(label=label),hjust=1,size=bar.label.size,family=font.family)
		## hjust==0
		df_sub <- df %>% dplyr::filter(hjust==0)
		bp <- bp + geom_text(data=df_sub, aes(label=label),hjust=0,size=bar.label.size,family=font.family)
	}
	
	## caption
    if(signature){
    	caption <- paste("Created by xEnrichBarplot from XGR version", utils ::packageVersion("XGR"))
    	bp <- bp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
	
	## change font family to 'Arial'
	bp <- bp + theme(text=element_text(family=font.family))
	
	## put arrows on x-axis
	bp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	
	## x-axis (actually y-axis) position
	bp <- bp + scale_y_continuous(position="top")
		
	invisible(bp)
}
