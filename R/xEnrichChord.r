#' Function to visualise enrichment results using a chord plot
#'
#' \code{xEnrichChord} is supposed to visualise enrichment results using a chord plot. The thickness of links is proportional to the enrichment Z-scores. Particularly useful for multiple groups and/or ontologies. The left-half part sorted by the input groups (anti-clockwise), and the right-half part sorted first by the input ontolgies and then by the number of links within an ontology (clockwise).
#'
#' @param eTerm an object of class "eTerm" or "ls_eTerm". Alterntively, it can be a data frame having all these columns (named as 'group','ontology','name','adjp','zscore')
#' @param top_num the number of the top terms (sorted according to adjp). For the eTerm object, if it is 'auto' (for eTerm), only the significant terms (see below FDR.cutoff) will be displayed
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05
#' @param colormap.group short name for the group sector colormap
#' @param colormap.ontology short name for the ontology sector colormap
#' @param wrap.width a positive integer specifying wrap width of labellings
#' @param text.size the text size of the labelings. By default, it is 0.6
#' @param legend logical to indicate whether to show the legend. If NULL, automatically determined. For the group sector, the legends shown on the bottom-left corner. For the ontology sector, the legends shown on the bottom-right corner
#' @param vline logical to indicate whether to vertically put a line seperating two symmetric parts of sectors
#' @param ... additional graphic parameters (such as big.gap=15, small.gap=1.5) used in circlize::chordDiagram
#' @return a data frame used for visualisation
#' @note none
#' @export
#' @seealso \code{\link{xEnrichChord}}
#' @include xEnrichChord.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' xEnrichChord(eTerm)
#' }

xEnrichChord <- function(eTerm, top_num=5, FDR.cutoff=0.05, colormap.group="ggplot2", colormap.ontology=NULL, wrap.width=NULL, text.size=0.6, legend=NULL, vline=F, ...)
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
		df <- xEnrichViewer(eTerm, top_num=top_num, sortBy="adjp")
		df$group <- 'group'
		df$ontology <- 'ontology'
		
	}else if(class(eTerm)=='ls_eTerm' | class(eTerm)=='data.frame'){
	
		if(class(eTerm)=='ls_eTerm'){
			## when 'auto', will keep the significant terms
			df <- eTerm$df[,c('group','ontology','name','adjp','zscore')]
			
		}else if(class(eTerm)=='data.frame'){
			if(all(c('group','ontology','name','adjp','zscore') %in% colnames(eTerm))){
				df <- eTerm[,c('group','ontology','name','adjp','zscore')]
			
			}else if(all(c('group','name','adjp','zscore') %in% colnames(eTerm))){
				df <- eTerm[,c('group','name','adjp','zscore')]
				df$ontology <- 'ontology'
			
			}else if(all(c('ontology','name','adjp','zscore') %in% colnames(eTerm))){
				df <- eTerm[,c('ontology','name','adjp','zscore')]
				df$group <- 'group'
			
			}else if(all(c('name','adjp','zscore') %in% colnames(eTerm))){
				df <- eTerm[,c('name','adjp','zscore')]
				df$group <- 'group'
				df$ontology <- 'ontology'
			
			}else{
				warnings("The input data.frame does not contain required columns: c('group','ontology','name','adjp','zscore').\n")
				return(NULL)
			}
			
		}
		
		group <- ontology <- rank <- adjp <- NULL
		df <- df %>% dplyr::arrange(adjp)
		if(top_num=='auto'){
			df <- subset(df, df$adjp<FDR.cutoff)
		}else{
			top_num <- as.integer(top_num)
			df <- as.data.frame(df %>% dplyr::group_by(group,ontology) %>% dplyr::group_by(rank=rank(adjp),add=TRUE) %>% dplyr::filter(rank<=top_num & adjp<FDR.cutoff))
		}
		
	}

	## text wrap
	if(!is.null(wrap.width)){
		width <- as.integer(wrap.width)
		res_list <- lapply(df$name, function(x){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=width)
			if(length(y)>1){
				paste(y,collapse='\n')
			}else{
				y
			}
		})
		df$name <- unlist(res_list)
	}

	## define order in Chord graph (the left half for group, the right half for ontology)
	name <- n <- ontology <- zscore <- group <- NULL
	df <- as.data.frame(df %>% dplyr::group_by(name) %>% dplyr::group_by(n=n(),add=T) %>% dplyr::arrange(ontology, n,zscore))
	order <- c(sort(unique(df$group)), unique(df$name))

	#######################################
	# how to deal with colormap.ontology when NULL
	if(is.null(colormap.ontology)){
		if(length(table(df$ontology))==1){
			## always grey
			colormap.ontology <- "grey-grey"
		}else{
			## always spectral
			colormap.ontology <- "spectral"
		}
	}
	#######################################
	
	# chordDiagram
	circlize::circos.clear()
	circlize::circos.par(start.degree=90, clock.wise=FALSE)
	
	## group colors
	tmp <- table(df$group)
	color_group <- xColormap(colormap.group)(length(tmp))
	names(color_group) <- names(tmp)
	## ontology colors
	tmp <- table(df$ontology)
	color_ontology <- xColormap(colormap.ontology)(length(tmp))
	names(color_ontology) <- names(tmp)
	## ontology colors -> name colors
	color_name <- color_ontology[df$ontology]
	names(color_name) <- df$name
	## grid.col
	grid.col <- c(color_group,color_name)

	## df_data to plot
	df_data <- df[,c('group','name','zscore')]
	#circlize::chordDiagram(df_data, order=order, grid.col=grid.col)

	circlize::chordDiagram(df_data, order=order, annotationTrack="grid", grid.col=grid.col, preAllocateTracks=1, ...)
	
	circlize::circos.trackPlotRegion(track.index=1, panel.fun = function(x, y) {
  		xlim = circlize::get.cell.meta.data("xlim")
  		ylim = circlize::get.cell.meta.data("ylim")
  		sector.name = circlize::get.cell.meta.data("sector.index")
  		circlize::circos.text(mean(xlim), ylim[1] + .1, sector.name, facing="clockwise", niceFacing=TRUE, adj=c(0,0.5), cex=text.size)
  		#circlize::circos.axis(h="top", labels.cex=0.5, major.tick.percentage=0.2, sector.index=sector.name, track.index=2)
	}, bg.border=NA)

	# symmetric line
	if(vline){
		abline(v=0, lty=2, col="#00000080")
	}
	
	# legend for groups
	if(any((is.null(legend) & length(color_group)>1), legend)){
		legend("topleft", title="Groups", legend=names(color_group), border=color_group, fill=color_group, horiz=F, box.col="transparent", cex=text.size)	
	}
	# legend for ontology
	if(any((is.null(legend) & length(color_ontology)>1), legend)){
		legend("topright", title="Ontologies", legend=names(color_ontology), border=color_ontology, fill=color_ontology, horiz=F, box.col="transparent", cex=text.size)
	}
	
	invisible(df)
}

