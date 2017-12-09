#' Function to visualise enrichment results using a forest plot
#'
#' \code{xEnrichForest} is supposed to visualise enrichment results using a forest plot. A point is colored by the significance level, and a horizontal line for the 95% confidence interval (CI) of odds ratio (OR; the wider the CI, the less reliable). It returns an object of class "ggplot".
#'
#' @param eTerm an object of class "eTerm" or "ls_eTerm". Alterntively, it can be a data frame having all these columns (named as 'group','ontology','name','adjp','or','CIl','CIu')
#' @param top_num the number of the top terms (sorted according to FDR or adjusted p-values). If it is 'auto', only the significant terms (see below FDR.cutoff) will be displayed
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05. This option only works when setting top_num (see above) is 'auto'
#' @param CI.one logical to indicate whether to allow the inclusion of one in CI. By default, it is TURE (allowed)
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to the range of the -log10(FDR)
#' @param barwidth the width of the colorbar. Default value is 'legend.key.width' or 'legend.key.size' in 'theme' or theme
#' @param barheight the height of the colorbar. Default value is 'legend.key.height' or 'legend.key.size' in 'theme' or theme
#' @param wrap.width a positive integer specifying wrap width of name
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichForest.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' 
#' # provide the input Genes of interest (eg 100 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' set.seed(825)
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 100))
#' data
#' 
#' # optionally, provide the test background (if not provided, all human genes)
#' #background <- as.character(org.Hs.eg$gene_info$Symbol)
#' 
#' # 1) Gene-based enrichment analysis using REACTOME pathways
#' # perform enrichment analysis
#' eTerm <- xEnricherGenes(data, ontology="REACTOME", RData.location=RData.location)
#' ## forest plot of enrichment results
#' gp <- xEnrichForest(eTerm, top_num="auto", FDR.cutoff=0.05)
#' 
#' # 2) Gene-based enrichment analysis using ontologies (REACTOME and GOMF)
#' # perform enrichment analysis
#' ls_eTerm <- xEnricherGenesAdv(data, ontologies=c("REACTOME","GOMF"), RData.location=RData.location)
#' ## forest plot of enrichment results
#' gp <- xEnrichForest(ls_eTerm, FDR.cutoff=0.1)
#' }

xEnrichForest <- function(eTerm, top_num=10, FDR.cutoff=0.05, CI.one=T, colormap="ggplot2.top", ncolors=64, zlim=NULL, barwidth=0.5, barheight=NULL, wrap.width=NULL, font.family="sans", signature=TRUE)
{
    
    if(is.null(eTerm)){
        warnings("There is no enrichment in the 'eTerm' object.\n")
        return(NULL)
    }
    
    if(class(eTerm)=='eTerm'){
		## when 'auto', will keep the significant terms
		df <- xEnrichViewer(eTerm, top_num="all")
		if(top_num=='auto'){
			top_num <- sum(df$adjp<FDR.cutoff)
			if(top_num <= 1){
				top_num <- 10
			}
		}
		df <- xEnrichViewer(eTerm, top_num=top_num, sortBy="or")
		df$group <- 'group'
		df$ontology <- 'ontology'
		
	}else if(class(eTerm)=='ls_eTerm' | class(eTerm)=='data.frame'){
	
		if(class(eTerm)=='ls_eTerm'){
			## when 'auto', will keep the significant terms
			df <- eTerm$df
			
		}else if(class(eTerm)=='data.frame'){
			if(all(c('group','ontology','name','adjp','or','CIl','CIu') %in% colnames(eTerm))){
				df <- eTerm[,c('group','ontology','name','adjp','or','CIl','CIu')]
			
			}else if(all(c('group','name','adjp','or','CIl','CIu') %in% colnames(eTerm))){
				df <- eTerm[,c('group','name','adjp','or','CIl','CIu')]
				df$ontology <- 'ontology'
			
			}else if(all(c('ontology','name','adjp','or','CIl','CIu') %in% colnames(eTerm))){
				df <- eTerm[,c('ontology','name','adjp','or','CIl','CIu')]
				df$group <- 'group'
			
			}else if(all(c('name','adjp','or','CIl','CIu') %in% colnames(eTerm))){
				df <- eTerm[,c('name','adjp','or','CIl','CIu')]
				df$group <- 'group'
				df$ontology <- 'ontology'
			
			}else{
				warnings("The input data.frame does not contain required columns: c('group','ontology','name','adjp','or','CIl','CIu').\n")
				return(NULL)
			}
			
		}
		
		or <- group <- ontology <- rank <- NULL
		df <- df %>% dplyr::arrange(-or)
		if(top_num=='auto'){
			df <- subset(df, df$adjp<FDR.cutoff)
		}else{
			top_num <- as.integer(top_num)
			df <- as.data.frame(df %>% dplyr::group_by(group,ontology) %>% dplyr::group_by(rank=order(or,decreasing=T),add=TRUE) %>% dplyr::filter(rank<=top_num))
		}
		
	}

	##########################
	##########################
	if(!CI.one){
		ind <- which(df$CIl>1 | df$CIu<1)
		df <- df[ind,]
	}
	
	if(nrow(df)==0){
		return(NULL)
	}
	##########################
	##########################

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
	
	name <- fdr <- or <- CIl <- CIu <- NULL
	group <- ontology <- NULL
	
	df$fdr <- -log10(df$adjp)
	if(is.null(zlim)){
		tmp <- df$fdr
		zlim <- c(floor(min(tmp)), ceiling(max(tmp[!is.infinite(tmp)])))
	}
	df$fdr[df$fdr<=zlim[1]] <- zlim[1]
	df$fdr[df$fdr>=zlim[2]] <- zlim[2]
	
	## order by 'or', 'adjp'
	df$group <- factor(df$group, levels=unique(df$group))
	df <- df[with(df,order(group, ontology, or, fdr)),]
	df$name <- factor(df$name, levels=unique(df$name))
	
	###########################################
	bp <- ggplot(df, aes(x=name, y=log2(or), ymin=log2(CIl), ymax=log2(CIu), color=fdr))
	bp <- bp + geom_pointrange() + ylab(expression(log[2]("odds ratio")))
	bp <- bp + geom_hline(yintercept=0, color='black', linetype='dashed') + coord_flip()
	bp <- bp + theme_bw() + theme(legend.position="right",axis.title.y=element_blank(), axis.text.y=element_text(size=10,color="black"), axis.title.x=element_text(size=10,color="black")) + coord_flip()
	bp <- bp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	bp <- bp + scale_colour_gradientn(colors=xColormap(colormap)(ncolors), limits=zlim, guide=guide_colorbar(title=expression(-log[10]("FDR")),title.position="top",barwidth=barwidth,barheight=barheight,draw.ulim=FALSE,draw.llim=FALSE))
	
	## caption
    if(signature){
    	caption <- paste("Created by xEnrichForest from XGR version", utils ::packageVersion("XGR"))
    	bp <- bp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
	
	## change font family to 'Arial'
	bp <- bp + theme(text=element_text(family=font.family))
	
	## put arrows on x-axis
	bp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	
	## x-axis (actually y-axis) position
	#bp <- bp + scale_y_continuous(position="bottom")
	
	# facet_grid: partitions a plot into a matrix of panels
	## group (columns), ontology (rows)
	ngroup <- length(unique(df$group))
	nonto <- length(unique(df$ontology))
	if(ngroup!=1 | nonto!=1){
		scales <- "free_y"
		space <- "free_y"
		
		if(ngroup==1){
			bp <- bp + facet_grid(ontology~., scales=scales, space=space)
		}else if(nonto==1){
			bp <- bp + facet_grid(.~group, scales=scales, space=space)
		}else{
			bp <- bp + facet_grid(ontology~group, scales=scales, space=space)
		}
		
		## strip
		bp <- bp + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=10,face="bold.italic"))
	}
	
	invisible(bp)
}

