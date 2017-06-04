#' Function to compare enrichment results using matrix plots
#'
#' \code{xEnrichMatrix} is supposed to compare enrichment results using matrix plots.
#'
#' @param list_eTerm a list of "eTerm" objects, or a data frame (with at least 3 columns "group", "name" and "adjp")
#' @param method which method will be used for plotting. It can be "circle" (by default), "square", "color" and "pie"
#' @param displayBy which statistics will be used for comparison. It can be "fc" for enrichment fold change (by default), "adjp" for adjusted p value (or FDR), "pvalue" for p value, "zscore" for enrichment z-score
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05
#' @param wrap.width a positive integer specifying wrap width of name
#' @param sharings a numeric vector specifying whether only shared terms will be displayed. For example, when comparing three groups of enrichment results, it can be set into c(2,3) to display only shared terms by any two or all three. By default, it is NULL meaning no such restriction
#' @param reorder how to reorder rows and columns. It can be "none" for no reordering, "row" for reordering rows according to number of sharings (by default), "col" for reordering columns, and "both" for reordering rows and columns
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to the range of the finite values of displayed matrix
#' @param ... additional graphic parameters for corrplot::corrplot
#' @return 
#' a data frame
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichMatrix.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' xEnrichMatrix(list_eTerm, method="color", displayBy="zscore", reorder="row", colormap="jet", zlim=c(0,10), cl.pos="b", cl.ratio=0.1, cl.align.text="c", mar=c(0.1,0.1,0.5,0.1), tl.col="black", tl.cex=0.6, tl.srt=90)
#' xEnrichMatrix(list_eTerm, method="pie", displayBy="fdr", FDR.cutoff=0.05, wrap.width=NULL, sharings=NULL, reorder="row", colormap="green-white-blue", ncolors=20, zlim=c(0,10), cl.pos="n", mar=c(0.1,0.1,0.5,0.1), tl.col="black", tl.cex=0.6, tl.srt=90)
#' }

xEnrichMatrix <- function(list_eTerm, method=c("circle","square","color","pie"), displayBy=c("zscore","fc","adjp","pvalue"), FDR.cutoff=0.05, wrap.width=NULL, sharings=NULL, reorder=c("row","none","col","both"), colormap="jet", ncolors=20, zlim=NULL, ...)
{
    method <- match.arg(method)    
    displayBy <- match.arg(displayBy)
    reorder <- match.arg(reorder)
    
    if(class(list_eTerm)=="list"){
		## Remove null elements in a list
		list_eTerm <- base::Filter(base::Negate(is.null), list_eTerm)
	
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
		
	}else if(class(list_eTerm)=="data.frame"){
		ind <- colnames(list_eTerm) %in% c("group","name","adjp",displayBy)
		df_all <- list_eTerm[,ind]
		## extract the columns: name fc adjp group
		ind <- which(df_all$adjp < FDR.cutoff)
		d <- df_all[ind,]
		
		list_names <- unique(d$group)
	}
		
	## group factor
	d$group <- factor(d$group, levels=rev(list_names))
	
	## append 'nSig' and 'code' to the data frame 'd'
	### nSig: the number of sharings per significant term
	### code: indicative of being present/absent for each eTerm (the same order as the input)
	id_ls <- split(x=d$group, f=d$name)
	ind <- match(d$name, names(id_ls))
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
	
	if(displayBy=='fc'){
		## sort by: nSig group fc (adjp)
		d <- d[with(d, order(nSig,group,fc,-adjp)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		d$val <- d$fc
		mat_val <- as.matrix(xSparseMatrix(d[,c('name','group','val')]))
		mat_val[mat_val==0] <- 1
	}else if(displayBy=='adjp'){
		########
		d$adjp[d$adjp==0] <- min(d$adjp[d$adjp!=0])
		########
		## sort by: nSig group adjp
		d <- d[with(d, order(nSig,group,-adjp)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		d$val <- -1*log10(d$adjp)
		mat_val <- as.matrix(xSparseMatrix(d[,c('name','group','val')]))
	}else if(displayBy=='zscore'){
		## sort by: nSig group zcore (adjp)
		d <- d[with(d, order(nSig,group,zscore,-adjp)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		d$val <- d$zscore
		mat_val <- as.matrix(xSparseMatrix(d[,c('name','group','val')]))
	}else if(displayBy=='pvalue'){
		########
		d$pvalue[d$pvalue==0] <- min(d$pvalue[d$pvalue!=0])
		########
		## sort by: nSig group pvalue
		d <- d[with(d, order(nSig,group,-pvalue)), ]
		## define levels
		d$name <- factor(d$name, levels=unique(d$name))
		d$val <- -1*log10(d$pvalue)
		mat_val <- as.matrix(xSparseMatrix(d[,c('name','group','val')]))
	}
	
	mat_fdr <- as.matrix(xSparseMatrix(d[,c('name','group','adjp')]))
	mat_fdr[mat_fdr==0] <- 1
	
	ind_row <- 1:nrow(mat_val)
	if(reorder=="row" | reorder=="both"){
		ind_row <- match(levels(d$name), rownames(mat_val))
	}
	ind_row <- rev(ind_row)
	ind_col <- 1:ncol(mat_val)
	if(reorder=="col" | reorder=="both"){
		mat <- mat_val
		colnames(mat) <- 1:ncol(mat)
		rownames(mat) <- 1:nrow(mat)
		tree_bs <- visTreeBootstrap(t(mat), visTree=FALSE)
		ind_col <- match(tree_bs$tip.label, colnames(mat))
	}
	mat_val <- mat_val[ind_row, ind_col]
	mat_fdr <- mat_fdr[ind_row, ind_col]
	
	if(is.null(zlim)){
		zlim <- c(min(mat_val), max(mat_val))
	}
	mat_val[mat_val<=zlim[1]] <- zlim[1]
	mat_val[mat_val>=zlim[2]] <- zlim[2]
	
	corrplot::corrplot(mat_val, method=method, is.cor=FALSE, col=xColormap(colormap)(ncolors), cl.lim=c(zlim[1],zlim[2]), p.mat=mat_fdr, sig.level=FDR.cutoff, insig="blank", addgrid.col="transparent", mar=c(0.1,0.1,1,0.1), ...)

	invisible(d)
}
