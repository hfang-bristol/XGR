#' Function to draw TOM plot
#'
#' \code{xWGCNAtom} is supposed to draw TOM plot from an object of class "cModule".
#'
#' @param cModule an object of class "cModule"
#' @param beta_pseudo NULL or an integer used to enhance TOM plot. If NULL, beta will be the one determined internally. It can be specified explicitly by the beta value (an integer)
#' @param which.modules an integer vector specifying which modules are visualised. If NULL (by default), all modules will be used
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to [0,1]
#' @param filename the without-extension part of the name of the output file. By default, it is 'xWGCNAtom'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a named vector specifying module colors
#' @note none
#' @export
#' @seealso \code{\link{xWGCNAtom}}
#' @include xWGCNAtom.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' modules_color <- xWGCNAtom(cModule)
#' modules_color <- xWGCNAtom(cModule, beta_pseudo=4, which.modules=c(1,3,10))
#' }

xWGCNAtom <- function(cModule, beta_pseudo=NULL, which.modules=NULL, colormap="grey-orange-red", zlim=c(0,1), filename='xWGCNAtom', verbose=T)
{
    
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    if(class(cModule)!= "cModule" & is.na(cModule$io$beta)){
        stop("The function must apply to a 'cModule' object.\n")
    }
    
	df_nodes_modules_rank <- cModule$modules
	
	#######################################################
	# TOM Heatmap plot
	#######################################################
	modules_color <- NULL
    if(!is.null(filename)){
    
		if(verbose){
			message(sprintf("Calculate Topological Overlap Matrix (TOM) (%s) ...", as.character(Sys.time())), appendLF=T)
		}
		expr <- cModule$expr
		adj <- cModule$adj
		TOMType <- cModule$io$TOMType
		num_modules <- cModule$io$num_modules
		
		if(is.null(beta_pseudo)){
			beta_pseudo <- cModule$io$beta
		}
		
		# Topological Overlap Matrix (TOM)
		adjMat <- adj
		if(TOMType=='signed'){
			adjMat <- adjMat * sign(WGCNA::cor(t(expr), use="p", method="pearson"))
		}
		tom <- WGCNA::TOMsimilarity(adjMat, TOMType=TOMType, verbose=0)
		diss <- 1 - tom
		tom_pseudo <- 1 - diss^beta_pseudo
		diag(tom_pseudo) <- 0
	
		# Colors for modules
		modules_size <- table(df_nodes_modules_rank$modules)
		modules_color <- xColormap(colormap='ggplot2')(length(modules_size))
		names(modules_color) <- names(modules_size)
		labeltree <- modules_color[df_nodes_modules_rank$modules]
		
		if(verbose){
			message(sprintf("Plot TOM heatmap powered by %d (%s) ...", beta_pseudo, as.character(Sys.time())), appendLF=T)
		}
		
		if(!is.null(which.modules)){
			which.modules <- as.integer(which.modules)
			if(all(which.modules %in% 1:num_modules)){
				if(verbose){
					message(sprintf("\tRestricted to modules '%s'(%s)", paste(which.modules,collapse=','), as.character(Sys.time())), appendLF=T)
				}
				ind <- match(df_nodes_modules_rank$modules, which.modules)
				ind <- which(!is.na(ind))
				tom_pseudo <- tom_pseudo[ind, ind]
				labeltree <- labeltree[ind]
			}
			
		}
		
		# output the file
		filename <- gsub('.jpg$', '', filename)
		outputfile <- paste0(filename, ".jpg")
		grDevices::jpeg(outputfile)
		supraHex::visHeatmap(tom_pseudo, row.metric="none", column.metric="none", colormap=colormap, zlim=zlim, scale="none", revC=TRUE, ColSideColors=labeltree, RowSideColors=labeltree, labRow=FALSE, labCol=FALSE)
		dev.off()
		
		if(verbose){
			message(sprintf("Congratulations! A file '%s' (in the directory %s) has been created! (%s)", outputfile, getwd(), as.character(Sys.time())), appendLF=T)
		}
    }
    
####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
	############################
	invisible(modules_color)
}