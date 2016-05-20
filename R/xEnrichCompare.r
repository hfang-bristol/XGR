#' Function to compare enrichment results using side-by-side barplots
#'
#' \code{xEnrichCompare} is supposed to compare enrichment results using side-by-side barplots. It returns an object of class "ggplot".
#'
#' @param list_eTerm a list of "eTerm" objects
#' @param displayBy which statistics will be used for comparison. It can be "fc" for enrichment fold change (by default), "adjp" for adjusted p value, "pvalue" for p value, "zscore" for enrichment z-score
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
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
#' eTerm_noLD_noTree <- xEnricherSNPs(data, ontology="EF", include.LD=NA, ontology.algorithm="none", RData.location=RData.location)
#' # 2b) Without considering LD SNPs but respecting ontology tree
#' eTerm_noLD_Tree <- xEnricherSNPs(data, ontology="EF", include.LD=NA, ontology.algorithm="lea", RData.location=RData.location)
#' # 2c) Considering LD SNPs but without respecting ontology tree
#' eTerm_LD_noTree <- xEnricherSNPs(data, ontology="EF", include.LD="EUR", LD.r2=0.8, ontology.algorithm="none", RData.location=RData.location)
#' # 2d) Considering LD SNPs and respecting ontology tree
#' eTerm_LD_Tree <- xEnricherSNPs(data, ontology="EF", include.LD="EUR", LD.r2=0.8, ontology.algorithm="lea", RData.location=RData.location)
#'
#' # 3) Compare enrichment results
#' list_eTerm <- list(eTerm_noLD_noTree, eTerm_noLD_Tree, eTerm_LD_noTree, eTerm_LD_Tree)
#' names(list_eTerm) <- c('LD (-) & Tree (-)', 'LD (-) & Tree (+)', 'LD (+) & Tree (-)', 'LD (+) & Tree (+)')
#' ## side-by-side comparisons 
#' bp <- xEnrichCompare(list_eTerm, displayBy="fc")
#' #pdf(file="enrichment_compared.pdf", height=6, width=12, compress=TRUE)
#' print(bp)
#' #dev.off()
#' }

xEnrichCompare <- function(list_eTerm, displayBy=c("fc","adjp","zscore","pvalue")) 
{
    
    displayBy <- match.arg(displayBy)
    
	## Combine into a data frame called 'df_all'
	list_names <- names(list_eTerm)
	res_ls <- lapply(1:length(list_eTerm), function(i){
		df <- xEnrichViewer(list_eTerm[[i]], top_num="all", sortBy="none")
		cbind(group=rep(list_names[i],nrow(df)), df)
	})
	df_all <- do.call(rbind, res_ls)

	## extract the columns: name fc adjp group
	ind <- which(df_all$adjp<0.05)
	d <- df_all[ind, c("name","fc","adjp","zscore","pvalue","group")]

	## append the number of sharings per significant term: nSig
	nSig <- base::table(d$name)
	ind <- match(d$name, names(nSig))
	d$nSig <- nSig[ind]

	## append the flag to indicate the enrichment significant level: 
	## * (<0.05)
	## ** (<0.01)
	## ** (<0.001)
	d$flag <- rep('', nrow(d))
	d$flag[d$adjp < 0.05] <- '*'
	d$flag[d$adjp < 0.01] <- '**'
	d$flag[d$adjp < 0.001] <- '***'
	
	## draw side-by-side barplot
	if(displayBy=='fc'){
		## sort by: nSig group fc (adjp)
		d <- d[with(d, order(nSig,group,fc,-adjp)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$name), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		#p <- ggplot(d, aes(x=name,y=d$fc,fill=group))
		p <- ggplot(d, eval(parse(text=paste("aes(x=name,y=d$fc,fill=group)",sep=""))))
		p <- p + ylab("Enrichment changes")
	}else if(displayBy=='adjp'){
		## sort by: nSig group adjp (zscore)
		d <- d[with(d, order(nSig,group,-adjp,zscore)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$name), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		#p <- ggplot(d, aes(x=name,y=-1*log10(d$adjp),fill=group))
		p <- ggplot(d, eval(parse(text=paste("aes(x=name,y=-1*log10(d$adjp),fill=group)",sep=""))))
		p <- p + ylab("Enrichment significance: -log10(FDR)")
	}else if(displayBy=='zscore'){
		## sort by: nSig group zcore (adjp)
		d <- d[with(d, order(nSig,group,zscore,-adjp)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$name), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		#p <- ggplot(d, aes(x=name,y=d$zscore,fill=group))
		p <- ggplot(d, eval(parse(text=paste("aes(x=name,y=d$zscore,fill=group)",sep=""))))
		p <- p + ylab("Enrichment z-scores")
	}else if(displayBy=='pvalue'){
		## sort by: nSig group pvalue (zscore)
		d <- d[with(d, order(nSig,group,-pvalue,zscore)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		## define horizontal lines: separating terms according to their sharings
		ind <- match(unique(d$name), names(nSig))
		xintercept <- which(!duplicated(nSig[ind]))[-1]
		## ggplot
		#p <- ggplot(d, aes(x=name,y=-1*log10(d$pvalue),fill=group))
		p <- ggplot(d, eval(parse(text=paste("aes(x=name,y=-1*log10(d$pvalue),fill=group)",sep=""))))
		p <- ggplot(d, aes(x=name, y=-1*log10(d$pvalue), fill=group))
		p <- p + ylab("Enrichment significance: -log10(p-value)")
	}
	
	bp <- p + geom_bar(stat="identity")+ theme_bw() + theme(legend.position="none", axis.title.y=element_blank(), axis.text.y=element_text(size=8,color="black"), axis.title.x=element_text(size=14,color="black")) + geom_vline(xintercept=xintercept-0.5,color="black",linetype="dotdash") + coord_flip() + facet_grid(~group) + geom_text(aes(label=d$flag),hjust=1,vjust=1,size=4)
	
	invisible(bp)
}
