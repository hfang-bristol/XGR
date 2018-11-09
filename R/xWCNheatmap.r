#' Function to draw heatmap from a cModule object
#'
#' \code{xWCNheatmap} is supposed to draw heatmap from an object of class "cModule".
#'
#' @param cModule an object of class "cModule"
#' @param displayBy which matrix used for display. It can be "tom" for the topological overlap matrix, "adj" for the adjacency matrix
#' @param power_pseudo NULL or an integer used to enhance the heatmap blocks. If NULL, no such enhancement. It can be specified explicitly by an integer
#' @param which.modules an integer vector specifying which modules are visualised. If NULL (by default), all modules will be used
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to [0,1]
#' @param filename the without-extension part of the name of the output file. By default, it is 'xWCNheatmap'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a named vector specifying module colors
#' @note none
#' @export
#' @seealso \code{\link{xWCNheatmap}}
#' @include xWCNheatmap.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' modules_color <- xWCNheatmap(cModule)
#' modules_color <- xWCNheatmap(cModule, displayBy="tom", power_pseudo=4, which.modules=c(1,3,10))
#' modules_color <- xWCNheatmap(cModule, displayBy="adj")
#' }

xWCNheatmap <- function(cModule, displayBy=c("tom","adj"), power_pseudo=NULL, which.modules=NULL, colormap="lightgrey-orange-red", zlim=c(0,1), filename='xWCNheatmap', verbose=T)
{
    
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    displayBy <- match.arg(displayBy)
    
    if(class(cModule)!= "cModule" & is.na(cModule$io$beta)){
        stop("The function must apply to a 'cModule' object.\n")
    }
    
	data <- cModule$data
	adj <- cModule$adj
	TOMType <- cModule$io$TOMType
	num_modules <- cModule$io$num_modules
	dataIsExpr <- T
	if(cModule$io$inputType=="similarity"){
		dataIsExpr <- F
	}
	df_nodes_modules_rank <- cModule$mem
	
	# Colors for modules
	modules_color <- xColormap(colormap='ggplot2')(num_modules)
	names(modules_color) <- 1:num_modules
	labeltree <- modules_color[df_nodes_modules_rank$modules]
	
	#######################################################
	# Heatmap plot
	#######################################################
    if(!is.null(filename)){
		
		# enhance the heatmap blocks
		if(is.null(power_pseudo)){
			power_pseudo <- 1
		}else{
			power_pseudo <- as.integer(power_pseudo)
		}
		
		if(verbose){
			message(sprintf("Heatmap of '%s' further powered by %d (%s) ...", displayBy, power_pseudo, as.character(Sys.time())), appendLF=T)
		}
		
		if(displayBy=="adj"){
			diss <- 1 - adj
		}else if(displayBy=="tom"){
			if(verbose){
				message(sprintf("Calculate Topological Overlap Matrix (TOM) typed as '%s' (%s) ...", TOMType, as.character(Sys.time())), appendLF=T)
			}
			adjMat <- adj
			if(TOMType=='signed'){
				if(dataIsExpr){
					adjMat <- adjMat * sign(WGCNA::cor(t(data), use="p", method="pearson", verbose=0))
				}else{
					adjMat <- adjMat * sign(data)
				}
			}
			diss <- 1 - WGCNA::TOMsimilarity(adjMat, TOMType=TOMType, verbose=F)
		}
		dataDisplay <- 1 - diss^power_pseudo
		diag(dataDisplay) <- 0
		
		# module restriction
		if(!is.null(which.modules)){
			which.modules <- as.integer(which.modules)
			if(all(which.modules %in% 1:num_modules)){
				if(verbose){
					message(sprintf("\tRestricted to modules '%s'(%s)", paste(which.modules,collapse=','), as.character(Sys.time())), appendLF=T)
				}
				ind <- match(df_nodes_modules_rank$modules, which.modules)
				ind <- which(!is.na(ind))
				dataDisplay <- dataDisplay[ind, ind]
				labeltree <- labeltree[ind]
			}
		}
		
		# output the file
		filename <- gsub('.jpg$', '', filename)
		outputfile <- paste0(filename, ".jpg")
		grDevices::jpeg(outputfile)
		supraHex::visHeatmap(dataDisplay, row.metric="none", column.metric="none", colormap=colormap, zlim=zlim, scale="none", revC=TRUE, ColSideColors=labeltree, RowSideColors=labeltree, labRow=FALSE, labCol=FALSE)
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