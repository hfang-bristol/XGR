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
#' @return an object of class "ggplot", but appended with a 'g' (an igraph object to represent DAG after being unionised)
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichDAGplotAdv}}
#' @include xEnrichCompare.r
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
#' names(list_eTerm) <- c('LD (-) & Tree (-)', 'LD (-) & Tree (+)', 'LD (+) & Tree (-)', 'LD (+) & Tree (+)')
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
#' bp + theme(legend.position="right", strip.text = element_blank())
#'
#' # 4) DAGplot of comparative enrichment results in the context of ontology tree
#' xEnrichDAGplotAdv(bp, graph.node.attrs=list(fontsize=100))
#' }

xEnrichCompare <- function(list_eTerm, displayBy=c("fc","adjp","fdr","zscore","pvalue"), FDR.cutoff=0.05, bar.label=TRUE, bar.label.size=3, wrap.width=NULL, sharings=NULL) 
{
    
    displayBy <- match.arg(displayBy)
    
	## Combine into a data frame called 'df_all'
	list_names <- names(list_eTerm)
	if(is.null(list_names)){
		list_names <- paste('Enrichment', 1:length(list_eTerm), sep=' ')
	}
	res_ls <- lapply(1:length(list_eTerm), function(i){
		df <- xEnrichViewer(list_eTerm[[i]], top_num="all", sortBy="none")
		cbind(group=rep(list_names[i],nrow(df)), id=rownames(df), df, stringsAsFactors=F)
	})
	df_all <- do.call(rbind, res_ls)
	rownames(df_all) <- NULL

	## extract the columns: name fc adjp group
	ind <- which(df_all$adjp < FDR.cutoff)
	d <- df_all[ind, c("id","name","fc","adjp","zscore","pvalue","group")]
	
	## group factor
	d$group <- factor(d$group, levels=list_names)

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
	
	## append 'label' (FDR) to the data frame 'd'
	d$label <- rep('',nrow(d))
	if(bar.label){
		## get text label
		to_scientific_notation <- function(x) {
			res <- format(x, scientific=T)
			res <- sub("\\+0?", "", res)
			sub("-0?", "-", res)
		}
		label <- to_scientific_notation(d$adjp)
		label <- paste('FDR', as.character(label), sep='=')
		d$label <- label
	}
	
	## draw side-by-side barplot
	if(displayBy=='fc'){
		## sort by: nSig group fc (adjp)
		d <- d[with(d, order(nSig,group,fc,-adjp)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$id), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		#p <- ggplot(d, aes(x=name,y=fc,fill=group))
		p <- ggplot(d, eval(parse(text=paste("aes(x=name,y=fc,fill=group)",sep=""))))
		p <- p + ylab("Enrichment changes")
	}else if(displayBy=='adjp' | displayBy=='fdr'){
		## sort by: nSig group adjp (zscore)
		d <- d[with(d, order(nSig,group,-adjp,zscore)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$id), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		#p <- ggplot(d, aes(x=name,y=-1*log10(adjp),fill=group))
		p <- ggplot(d, eval(parse(text=paste("aes(x=name,y=-1*log10(adjp),fill=group)",sep=""))))
		p <- p + ylab("Enrichment significance: -log10(FDR)")
	}else if(displayBy=='zscore'){
		## sort by: nSig group zcore (adjp)
		d <- d[with(d, order(nSig,group,zscore,-adjp)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$id), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		#p <- ggplot(d, aes(x=name,y=zscore,fill=group))
		p <- ggplot(d, eval(parse(text=paste("aes(x=name,y=zscore,fill=group)",sep=""))))
		p <- p + ylab("Enrichment z-scores")
	}else if(displayBy=='pvalue'){
		## sort by: nSig group pvalue (zscore)
		d <- d[with(d, order(nSig,group,-pvalue,zscore)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$id), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		#p <- ggplot(d, aes(x=name,y=-1*log10(pvalue),fill=group))
		p <- ggplot(d, eval(parse(text=paste("aes(x=name,y=-1*log10(pvalue),fill=group)",sep=""))))
		p <- p + ylab("Enrichment significance: -log10(p-value)")
	}
	
	p <- p + geom_bar(stat="identity")+ theme_bw() + theme(legend.position="none",legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=8,color="blue"), axis.title.x=element_text(size=14,color="blue")) + geom_vline(xintercept=xintercept-0.5,color="black",linetype="dotdash") + coord_flip()
	
	## title
	p <- p + ggtitle(paste0('Enrichments under FDR < ', FDR.cutoff))
	
	## strip
	p <- p + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=12,face="italic"))
	
	if(bar.label){
		#p <- p + geom_text(aes(label=label),hjust=1,size=bar.label.size)
		p <- p + geom_text(eval(parse(text=paste("aes(label=label)",sep=""))) ,hjust=1,size=bar.label.size)
	}
	
	## group
	#bp <- p + facet_grid(~group)
	bp <- p + eval(parse(text=paste("facet_grid(~group)",sep="")))
	
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
