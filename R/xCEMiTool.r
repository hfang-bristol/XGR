#' Function to post-process results from CEMiTool
#'
#' \code{xCEMiTool} is supposed to post-process results from CEMiTool (an object of class "CEMiTool").
#'
#' @param cem an object of class "CEMiTool"
#' @param filename the without-extension part of the name of the output file. By default, it is 'xCEMiTool'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a data frame with three columns 'genes' (ordered by intramodular connectivity), 'modules' (M1, M2, ...) and 'rank' (integer with a module).
#' @note none
#' @export
#' @seealso \code{\link{xCEMiTool}}
#' @include xCEMiTool.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' cem <- cemitool(expr, filter=FALSE, merge_similar=TRUE, min_ngen=70, diss_thresh=0.8, verbose=T)
#' df_res <- xCEMiTool(cem)
#' }

xCEMiTool <- function(cem, filename='xCEMiTool', verbose=T)
{
    
    if (class(cem) != "CEMiTool" & is.na(cem@parameters$beta)){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
	# df_genes_modules_rank: containing genes (grouped by modules and ranked by sum of connectivity to other genes within each module)
	if(verbose){
		message(sprintf("Extract module membership data info (%s) ...", as.character(Sys.time())), appendLF=T)
	}
	modules <- NULL
	df_genes_modules <- cem@module %>% dplyr::filter(modules != 'Not.Correlated')
	ls_modules <- split(x=df_genes_modules$genes, f=df_genes_modules$modules)
	ls_res <- lapply(1:length(ls_modules), function(i){
		x <- ls_modules[[i]]
		if(length(x)>1){
			mod_adj <- cem@adjacency[x, x]			
			genes <- names(sort(apply(mod_adj, 1, sum), decreasing=TRUE))
			data.frame(genes=genes, modules=names(ls_modules)[i], rank=1:length(genes), stringsAsFactors=F)
		}else{
			NULL
		}
	})
	df_genes_modules_rank <- do.call(rbind, ls_res)
	
	#######################################################
	# TOM Heatmap plot
	#######################################################
    if(!is.null(filename)){
    
		if(verbose){
			message(sprintf("Calculate Topological Overlap Matrix (TOM) (%s) ...", as.character(Sys.time())), appendLF=T)
		}
		# extract 'expr' and 'adj'; genes ordered according to df_genes_modules_rank
		expr <- cem@expression[cem@selected_genes,]
		ind <- match(df_genes_modules_rank$genes, rownames(expr))
		expr <- expr[ind,]
		adj <- cem@adjacency[ind,ind]

		# Topological Overlap Matrix (TOM)
		scor <- sign(WGCNA::cor(t(expr), use="p", method=cem@parameters$cor_method))
		tom <- WGCNA::TOMsimilarity(adj*scor, TOMType="signed")
		diss <- 1 - tom
		tom_pseudo <- 1 - diss^(cem@parameters$beta)
		diag(tom_pseudo) <- 0
	
		# Colors for modules
		modules_size <- table(df_genes_modules_rank$modules)
		modules_color <- xColormap(colormap='ggplot2')(length(modules_size))
		names(modules_color) <- names(modules_size)
		labeltree <- modules_color[df_genes_modules_rank$modules]
		
		if(verbose){
			message(sprintf("Plot TOM heatmap (%s) ...", as.character(Sys.time())), appendLF=T)
		}
		# output the file
		filename <- gsub('.jpg$', '', filename)
		outputfile <- paste0(filename, ".jpg")
		grDevices::jpeg(outputfile)
		ntop <- nrow(tom_pseudo)
		supraHex::visHeatmap(tom_pseudo[1:ntop,1:ntop], row.metric="none", column.metric="none", colormap="grey-orange-darkred", zlim=c(0,1), scale="none", revC=TRUE, ColSideColors=labeltree[1:ntop], RowSideColors=labeltree[1:ntop], labRow=FALSE, labCol=FALSE)
		dev.off()
		message(sprintf("Congratulations! A file '%s' (in the directory %s) has been created!", outputfile, getwd()), appendLF=T)
		############################
    }
    
    invisible(df_genes_modules_rank)
}






