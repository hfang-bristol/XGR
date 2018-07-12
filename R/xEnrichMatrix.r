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
#' @param slim the minimum and maximum displaying values for which sizes should be plotted
#' @param title the title of the plot. By default, it is NULL
#' @param flip logical to indicate whether to flip the coordiate. By default, it sets to false
#' @param y.rotate the angle to rotate the y tick labelings. By default, it is 45
#' @param shape the number specifying the shape. By default, it is 19
#' @param font.family the font family for texts
#' @param ... additional graphic parameters for corrplot::corrplot
#' @return 
#' If the method is 'ggplot2', it returns a ggplot object. Otherwise, it is a data frame
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichMatrix.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' xEnrichMatrix(list_eTerm, method="circle", displayBy="adjp", FDR.cutoff=0.05, wrap.width=50, sharings=NULL, reorder="row", colormap="black-yellow-red", ncolors=16, zlim=c(0,8), cl.pos="b", cl.ratio=0.1, cl.align.text="c", tl.col="black", tl.cex=0.7, tl.srt=90, title=paste0(ontology,": log10(FDR)"))
#' xEnrichMatrix(list_eTerm, method="pie", displayBy="adjp", FDR.cutoff=0.05, wrap.width=50, sharings=NULL, reorder="row", colormap="grey-grey", ncolors=1, zlim=c(0,8), cl.pos="n", cl.ratio=0.1, cl.align.text="c", tl.col="black", tl.cex=0.7, tl.srt=90, title=paste0(ontology,": log10(FDR)"))
#' gp <- xEnrichMatrix(list_eTerm, method="ggplot2", displayBy="zscore", FDR.cutoff=0.05, wrap.width=40, sharings=NULL, reorder="row", colormap="yellow-red", flip=T, y.rotate=45, font.family=font.family)
#' }

xEnrichMatrix <- function(list_eTerm, method=c("ggplot2","circle","square","color","pie"), displayBy=c("zscore","fc","adjp","pvalue"), FDR.cutoff=0.05, wrap.width=NULL, sharings=NULL, reorder=c("row","none","col","both"), colormap="jet", ncolors=20, zlim=NULL, slim=NULL, title=NULL, flip=FALSE, y.rotate=45, shape=19, font.family="sans", ...)
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
		if(sum(ind)==4){
			#df_all <- list_eTerm[,ind]
			df_all <- list_eTerm
			## extract the columns: name fc adjp group
			ind <- which(df_all$adjp < FDR.cutoff)
			d <- df_all[ind,]
			if(is.factor(d$group)){
				list_names <- levels(d$group)
			}else{
				list_names <- unique(d$group)
			}
		}else{
			return(NULL)
		}
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
		if(!is.factor(d$name)){
			d$name <- factor(d$name, levels=unique(d$name))
		}
		d$val <- d$fc
		mat_val <- as.matrix(xSparseMatrix(d[,c('name','group','val')], rows=levels(d$name), columns=levels(d$group)))
		mat_val[mat_val==0] <- 1
		title_size <- "FC"
	}else if(displayBy=='adjp'){
		########
		d$adjp[d$adjp==0] <- min(d$adjp[d$adjp!=0])
		########
		## sort by: nSig group adjp
		d <- d[with(d, order(nSig,group,-adjp)), ]
		## define levels
		if(!is.factor(d$name)){
			d$name <- factor(d$name, levels=unique(d$name))
		}
		d$val <- -1*log10(d$adjp)
		mat_val <- as.matrix(xSparseMatrix(d[,c('name','group','val')], rows=levels(d$name), columns=levels(d$group)))
		title_size <- "-log10(FDR)"
		
	}else if(displayBy=='zscore'){
		## sort by: nSig group zcore (adjp)
		d <- d[with(d, order(nSig,group,zscore,-adjp)), ]
		## define levels
		if(!is.factor(d$name)){
			d$name <- factor(d$name, levels=unique(d$name))
		}
		d$val <- d$zscore
		mat_val <- as.matrix(xSparseMatrix(d[,c('name','group','val')], rows=levels(d$name), columns=levels(d$group)))
		title_size <- "Z-score"
	}else if(displayBy=='pvalue'){
		########
		d$pvalue[d$pvalue==0] <- min(d$pvalue[d$pvalue!=0])
		########
		## sort by: nSig group pvalue
		d <- d[with(d, order(nSig,group,-pvalue)), ]
		## define levels
		if(!is.factor(d$name)){
			d$name <- factor(d$name, levels=unique(d$name))
		}
		d$val <- -1*log10(d$pvalue)
		mat_val <- as.matrix(xSparseMatrix(d[,c('name','group','val')], rows=levels(d$name), columns=levels(d$group)))
		title_size <- "-log10(p-value)"
	}
	
	d$bycol <- -1*log10(d$adjp)
	
	mat_fdr <- as.matrix(xSparseMatrix(d[,c('name','group','adjp')], rows=levels(d$name), columns=levels(d$group)))
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
	
	if(is.null(slim)){
		slim <- c(min(mat_val), max(mat_val))
	}
	mat_val[mat_val<=slim[1]] <- slim[1]
	mat_val[mat_val>=slim[2]] <- slim[2]
	
	if(method!="ggplot2"){
		if(is.null(title)){
			corrplot::corrplot(mat_val, method=method, is.cor=FALSE, col=xColormap(colormap)(ncolors), cl.lim=c(zlim[1],zlim[2]), p.mat=mat_fdr, sig.level=FDR.cutoff, insig="blank", addgrid.col="transparent", mar=c(0.1,0.1,1,0.1), ...)
		}else{
			corrplot::corrplot(mat_val, method=method, is.cor=FALSE, col=xColormap(colormap)(ncolors), cl.lim=c(zlim[1],zlim[2]), p.mat=mat_fdr, sig.level=FDR.cutoff, insig="blank", addgrid.col="transparent", mar=c(0.1,0.1,1,0.1), title=title, ...)
		}
		invisible(d)
		
	}else{
		####
		if(flip){
			d$group <- factor(d$group, levels=rev(colnames(mat_val)))
			d$name <- factor(d$name, levels=rownames(mat_val))
		}else{
			d$group <- factor(d$group, levels=colnames(mat_val))
			d$name <- factor(d$name, levels=rev(rownames(mat_val)))
		}
		d$val[d$val<=slim[1]] <- slim[1]
		d$val[d$val>=slim[2]] <- slim[2]
		####
		
		if(is.null(zlim)){
			zlim <- c(0,ceiling(max(d$bycol)*10)/10)
		}
		d$bycol[d$bycol>=zlim[2]] <- zlim[2]
		
		group <- name <- val <- zscore <- bycol <- NULL
		
		gp <- ggplot(d, aes(x=group, y=name, color=bycol))
		gp <- gp + geom_point(aes(size=val),shape=shape)
		gp <- gp + scale_colour_gradientn(colors=xColormap(colormap)(ncolors), limits=zlim, guide=guide_colorbar(title=expression(-log[10]("FDR")),title.position="top",barwidth=0.5,nbin=5,draw.ulim=FALSE,draw.llim=FALSE))
		gp <- gp + scale_size_continuous(limits=c(floor(min(d$val)*10)/10, ceiling(max(d$val)*10)/10), range=c(1,4), guide=guide_legend(title_size,title.position="top",ncol=1))
		
		gp <- gp + theme_bw() + theme(legend.position="right", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", color="black", size=8, angle=y.rotate), axis.text.y=element_text(face="bold", color="black", size=8, angle=0), panel.background=element_rect(fill="transparent")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		gp <- gp + labs(title=title)
		gp <- gp + theme(text=element_text(family=font.family))
		gp <- gp + scale_x_discrete(position="top")
		
		if(flip){
			gp <- gp + coord_flip() + scale_x_discrete(position="bottom") + scale_y_discrete(position="right")
			gp <- gp + theme(axis.text.x=element_text(face="bold", color="black", size=8, angle=y.rotate, hjust=0))
		}
		
		invisible(gp)
	}
	
}
