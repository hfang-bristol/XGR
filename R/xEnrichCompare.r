#' Function to compare enrichment results using side-by-side barplots
#'
#' \code{xEnrichCompare} is supposed to compare enrichment results using side-by-side barplots. It returns an object of class "ggplot".
#'
#' @param list_eTerm a list of "eTerm" objects
#' @param displayBy which statistics will be used for comparison. It can be "fc" for enrichment fold change (by default), "adjp" or "fdr" for adjusted p value (or FDR), "pvalue" for p value, "zscore" for enrichment z-score
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05
#' @param bar.label logical to indicate whether to label each bar with FDR. By default, it sets to true for bar labelling
#' @param bar.label.size an integer specifying the bar labelling text size. By default, it sets to 3
#' @param wrap.width a positive integer specifying wrap width of name
#' @param sharings a numeric vector specifying whether only shared terms will be displayed. For example, when comparing three groups of enrichment results, it can be set into c(2,3) to display only shared terms by any two or all three. By default, it is NULL meaning no such restriction
#' @param font.family the font family for texts
#' @param facet logical to indicate whether to facet/wrap a 1d of panels into 2d. By default, it sets TRUE
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot", but appended with a 'g' (an igraph object to represent DAG after being unionised)
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichDAGplotAdv}}
#' @include xEnrichCompare.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' # 1) load eQTL mapping results: cis-eQTLs significantly induced by IFN
#' cis <- xRDataLoader(RData.customised='JKscience_TS2A', RData.location=RData.location)
#' ind <- which(cis$IFN_t > 0 & cis$IFN_fdr < 0.05)
#' df_cis <- cis[ind, c('variant','Symbol','IFN_t','IFN_fdr')]
#' data <- df_cis$variant
#' 
#' # 2) Enrichment analysis using Experimental Factor Ontology (EFO)
#' # 2a) Without considering LD SNPs and without respecting ontology tree
#' eTerm_noLD_noTree <- xEnricherSNPs(data, ontology="EF_disease", include.LD=NA, ontology.algorithm="none", RData.location=RData.location)
#' # 2b) Without considering LD SNPs but respecting ontology tree
#' eTerm_noLD_Tree <- xEnricherSNPs(data, ontology="EF_disease", include.LD=NA, ontology.algorithm="lea", RData.location=RData.location)
#' # 2c) Considering LD SNPs but without respecting ontology tree
#' eTerm_LD_noTree <- xEnricherSNPs(data, ontology="EF_disease", include.LD="EUR", LD.r2=0.8, ontology.algorithm="none", RData.location=RData.location)
#' # 2d) Considering LD SNPs and respecting ontology tree
#' eTerm_LD_Tree <- xEnricherSNPs(data, ontology="EF_disease", include.LD="EUR", LD.r2=0.8, ontology.algorithm="lea", RData.location=RData.location)
#'
#' # 3) Compare enrichment results
#' list_eTerm <- list(eTerm_noLD_noTree, eTerm_noLD_Tree, eTerm_LD_noTree, eTerm_LD_Tree)
#' names(list_eTerm) <- c('LD(-) & Tree(-)','LD(-) & Tree(+)','LD(+) & Tree(-)','LD(+) & Tree(+)')
#' ## side-by-side comparisons 
#' bp <- xEnrichCompare(list_eTerm, displayBy="fc")
#' #pdf(file="enrichment_compared.pdf", height=6, width=12, compress=TRUE)
#' print(bp)
#' #dev.off()
#' ## modify y axis text
#' bp + theme(axis.text.y=element_text(size=10,color="black"))
#' ## modify x axis title
#' bp + theme(axis.title.x=element_text(color="black"))
#' ## modify fill colors
#' bp + scale_fill_manual(values=c("black","#888888"))
#' ## show legend title but hide strip
#' bp + theme(legend.position="right", strip.text=element_blank())
#'
#' # 4) DAGplot of comparative enrichment results in the context of ontology tree
#' xEnrichDAGplotAdv(bp, graph.node.attrs=list(fontsize=100))
#' }

xEnrichCompare <- function(list_eTerm, displayBy=c("fc","adjp","fdr","zscore","pvalue"), FDR.cutoff=0.05, bar.label=TRUE, bar.label.size=2.5, wrap.width=NULL, sharings=NULL, font.family="sans", facet=TRUE, signature=TRUE) 
{
    
    displayBy <- match.arg(displayBy)
    
    ## Remove null elements in a list
	list_eTerm <- base::Filter(base::Negate(is.null), list_eTerm)
    if(length(list_eTerm)==0){
    	return(NULL)
    }
    
	## Combine into a data frame called 'df_all'
	list_names <- names(list_eTerm)
	if(is.null(list_names)){
		list_names <- paste('Enrichment', 1:length(list_eTerm), sep=' ')
	}
	res_ls <- lapply(1:length(list_eTerm), function(i){
		df <- xEnrichViewer(list_eTerm[[i]], top_num="all", sortBy="none")
		if(is.null(df)){
			return(NULL)
		}else{
			cbind(group=rep(list_names[i],nrow(df)), id=rownames(df), df, stringsAsFactors=F)
		}
	})
	df_all <- do.call(rbind, res_ls)
	rownames(df_all) <- NULL

	## extract the columns: name fc adjp group
	ind <- which(df_all$adjp < FDR.cutoff)
	d <- df_all[ind, c("id","name","fc","adjp","zscore","pvalue","group")]
	
	## group factor
	d$group <- factor(d$group, levels=rev(list_names))

	## text wrap
	if(!is.null(wrap.width)){
		width <- as.integer(wrap.width)
		res_list <- lapply(d$name, function(x){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=width)
			if(length(y)>1){
				paste0(y[1], '...')
			}else{
				y
			}
		})
		d$name <- unlist(res_list)
	}

	## append 'nSig' and 'code' to the data frame 'd'
	### nSig: the number of sharings per significant term
	### code: indicative of being present/absent for each eTerm (the same order as the input)
	id_ls <- split(x=d$group, f=d$id)
	ind <- match(d$id, names(id_ls))
	id_full_ls <- id_ls[ind]
	#### for nSig
	nSig <- unlist(lapply(id_full_ls, length))
	d$nSig <- nSig
	#### for code
	code <- lapply(id_full_ls, function(x){
		res <- rep(0, length(levels(x)))
		ind <- match(x, levels(x))
		res[ind] <- 1
		paste(res, collapse='-')
	})
	d$code <- unlist(code)
	
	## restrict to sharings?
	if(!is.null(sharings)){
		sharings <- as.numeric(sharings)
		ind <- match(sharings, unique(d$nSig))
		found <- sharings[!is.na(ind)]
		if(length(found)>0){
			flag <- match(d$nSig, found)
			d <- d[!is.na(flag), ]
			
			nSig <- nSig[!is.na(flag)]
		}
	}
	
	## draw side-by-side barplot
	name <- fc <- group <- adjp <- zscore <- pvalue <- label <- direction <- height <- hjust <- NULL
	
	##########
	## consider the direction of z-score
	d <- d %>% dplyr::mutate(direction=ifelse(zscore>0,1,-1))
	##########
	
	if(displayBy=='fc'){
		## sort by: nSig group fc (adjp)
		d <- d %>% dplyr::arrange(nSig,group, direction, fc, desc(adjp)) %>% dplyr::mutate(height=log2(fc)) %>% dplyr::mutate(hjust=ifelse(height>=0,1,0))
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$id), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		p <- ggplot(d, aes(x=name,y=height,fill=group))
		p <- p + ylab(expression(paste("Enrichment changes: ", log[2]("FC"))))
		
	}else if(displayBy=='adjp' | displayBy=='fdr'){
		## sort by: nSig group adjp (zscore)
		d <- d %>% dplyr::arrange(nSig,group, direction, desc(adjp), zscore) %>% dplyr::mutate(height=-1*log10(adjp)) %>% dplyr::mutate(hjust=1)
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$id), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		p <- ggplot(d, aes(x=name,y=height,fill=group))
		p <- p + ylab(expression(paste("Enrichment significance: ", -log[10]("FDR"))))
		
	}else if(displayBy=='zscore'){
		## sort by: nSig group zcore (adjp)
		d <- d %>% dplyr::arrange(nSig,group, direction, zscore, desc(adjp)) %>% dplyr::mutate(height=zscore) %>% dplyr::mutate(hjust=ifelse(height>=0,1,0))
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$id), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		p <- ggplot(d, aes(x=name,y=height,fill=group))
		p <- p + ylab("Enrichment z-scores")
		
	}else if(displayBy=='pvalue'){
		## sort by: nSig group pvalue (zscore)
		d <- d %>% dplyr::arrange(nSig,group, direction, desc(pvalue), zscore) %>% dplyr::mutate(height=-1*log10(pvalue)) %>% dplyr::mutate(hjust=1)
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$id), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		p <- ggplot(d, aes(x=name,y=height,fill=group))
		p <- p + ylab(expression(paste("Enrichment significance: ", -log[10]("p-value"))))
		
	}
	
	bp <- p + geom_bar(stat="identity")+ theme_bw() + theme(legend.position="none",legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=10,color="black"), axis.title.x=element_text(size=14,color="black")) + geom_vline(xintercept=xintercept-0.5,color="black",linetype="dotdash") + coord_flip()
	bp <- bp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	## strip
	bp <- bp + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=12,face="italic"))
	
	if(bar.label){
		## get text label
		to_scientific_notation <- function(x) {
			res <- format(x, scientific=T)
			res <- sub("\\+0?", "", res)
			sub("-0?", "-", res)
		}
		label <- to_scientific_notation(d$adjp)
		d$label <- paste('FDR', as.character(label), sep='=')
	
		## hjust==1
		df_sub <- d %>% dplyr::filter(hjust==1)
		bp <- bp + geom_text(data=df_sub, aes(label=label),hjust=1,size=bar.label.size,family=font.family)
		## hjust==0
		df_sub <- d %>% dplyr::filter(hjust==0)
		bp <- bp + geom_text(data=df_sub, aes(label=label),hjust=0,size=bar.label.size,family=font.family)
	}
	
	## title
	title <- paste0('Enrichments (FDR < ', FDR.cutoff, ')')
	bp <- bp + labs(title=title) + theme(plot.title=element_text(hjust=0.5))
	## caption
    if(signature){
    	caption <- paste("Created by xEnrichCompare from XGR version", utils ::packageVersion("XGR"))
    	bp <- bp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
	
	## change font family to 'Arial'
	bp <- bp + theme(text=element_text(family=font.family))
	
	## put arrows on x-axis
	bp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	
	## group
	if(facet){
		bp <- bp + facet_grid(~group)
	}else{
		#bp <- bp + theme(legend.position="bottom",legend.title=element_blank()) + guides(fill=guide_legend(reverse=TRUE,nrow=1,byrow=TRUE))
		bp <- bp + theme(legend.position="bottom",legend.title=element_text(size=10,color="black",face="bold")) + guides(fill=guide_legend(title="",title.position="left",reverse=TRUE,nrow=1,byrow=TRUE))
	}
	
	##############################
	## append 'g' to 'bp' (if DAG)
	flag <- sapply(list_eTerm, function(x) !is.null(x$g))
	if(all(flag)){
		### edges
		ls_edges <- lapply(list_eTerm, function(x){
			df_edge <- igraph::get.data.frame(x$g, what="edges")
		})
		relations <- do.call(rbind, ls_edges)
		relations <- relations[!duplicated(relations), ]
		### nodes
		ls_nodes <- lapply(list_eTerm, function(x){
			df_nodes <- igraph::get.data.frame(x$g, what="vertices")
		})
		nodes <- do.call(rbind, ls_nodes)[,1:4]
		nodes <- nodes[!duplicated(nodes), ]
		### igraph
		ig <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
	
		bp$g <- ig
	}
	##############################

		
	invisible(bp)
}
