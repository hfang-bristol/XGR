#' Function to visualise enrichment results using radial-like plot
#'
#' \code{xEnrichRadial} is supposed to visualise enrichment results using radial-like plot. It returns three ggplot2 objects, the first for visualing the network with nodes lablelled by codes, the second for listing code meaning in a table, and the third for the network with nodes colored/sized with enrichment results.
#'
#' @param eTerm an object of class "eTerm" or "ls_eTerm". Alterntively, it can be a data frame having all these columns ('name','adjp','or','zscore'; 'group' optionally). Be aware that multiple ontologies are not supported here
#' @param ig the igraph object. If provided, only those terms within it will be visualised. By default, it is NULL meaning no surch restriction
#' @param fixed logical to indicate whether all terms in ig will be visualised. By default, it is TURE; otherwise only overlapped terms from eTerm will be visualised
#' @param node.color which statistics will be used for node coloring. It can be "or" for the odds ratio, "adjp" for adjusted p value (FDR) and "zscore" for enrichment z-score
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param zlim the minimum and maximum values for which colors should be plotted
#' @param node.size which statistics will be used for node size. It can be "or" for the odds ratio, "adjp" for adjusted p value (FDR) and "zscore" for enrichment z-score
#' @param slim the minimum and maximum values for which sizes should be plotted
#' @param node.size.range the range of actual node size
#' @param edge.color a character specifying which edge attribute defining the the edge colors
#' @param edge.color.alpha the 0-1 value specifying transparency of edge colors
#' @param edge.curve a numeric value specifying the edge curve. 0 for the straight line
#' @param edge.arrow.gap a gap between the arrow and the node
#' @param ... additional graphic parameters used in xGGnetwork
#' @return
#' a list with 3 components, three ggplot objects (code, table, data) and an igraph object (ig appended with node attributes 'zscore', 'adjp' and 'or')
#' @note none
#' @export
#' @seealso \code{\link{xEnrichViewer}}, \code{\link{xOBOcode}}, \code{\link{xGGnetwork}}
#' @include xEnrichRadial.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' ls_res <- xEnrichRadial(eTerm, ig, fixed=T, node.color="or", colormap="grey-orange-darkred", zlim=c(0,7), node.size="adjp", slim=c(0,30), node.size.range=c(1,3))
#' pdf("xEnrichRadial.pdf", width=6.5, height=6.5)
#' print(ls_res$data + coord_equal(ratio=1.3))
#' print(ls_res$code + coord_equal(ratio=1.3))
#' print(ls_res$table)
#' dev.off()
#' 
#' # advanced use: customise layout
#' ig.PhasedTargets <- xRDataLoader('ig.PhasedTargets', RData.location=RData.location)
#' ig <- xLayout(ig.PhasedTargets, layout="gplot.layout.fruchtermanreingold")
#' ls_res <- xEnrichRadial(df, ig=ig, fixed=F, node.color="or", node.size="adjp", node.xcoord="xcoord", node.ycoord="ycoord")
#' pdf("xEnrichRadial.pdf", width=6.5, height=6.5)
#' print(ls_res$data + coord_equal())
#' gridExtra::grid.arrange(grobs=c(list(ls_res$code+coord_equal()),ls_res$table), ncol=2)
#' dev.off()
#' }

xEnrichRadial <- function(eTerm, ig=NULL, fixed=T, node.color=c("or","adjp","zscore"), colormap="grey-orange-darkred", zlim=NULL, node.size=c("adjp","or","zscore"), slim=NULL, node.size.range=c(0.5,3.5), edge.color="skyblue", edge.color.alpha=0.5, edge.curve=0.1, edge.arrow.gap=0.02, ...)
{
	## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    node.color <- match.arg(node.color)
    node.size <- match.arg(node.size)
    
    gp_data <- NULL
    gp_code <- NULL
    gp_table <- NULL
    
    if(is.null(eTerm)){
        warnings("There is no enrichment in the 'eTerm' object.\n")
        return(NULL)
    }
    
    if(is(eTerm,'eTerm')){
		df_enrichment_group <- xEnrichViewer(eTerm, top_num="all", sortBy="or")
		df_enrichment_group$group <- 'group'
	}else if(is(eTerm,'ls_eTerm') | is(eTerm,'data.frame')){
	
		if(is(eTerm,'ls_eTerm')){
			df_enrichment_group <- eTerm$df
			
		}else if(is(eTerm,'data.frame')){
			if(all(c('group','name','adjp','or','zscore') %in% colnames(eTerm))){
				df_enrichment_group <- eTerm[,c('group','name','adjp','or','zscore')]
			}else if(all(c('name','adjp','or','zscore') %in% colnames(eTerm))){
				df_enrichment_group <- eTerm[,c('name','adjp','or','zscore')]
				df_enrichment_group$group <- 'group'
			}
		}
	}
	
	if(is(df_enrichment_group$group,'factor')){
		if(length(unique(df_enrichment_group$group)) != length(levels(df_enrichment_group$group))){
			df_enrichment_group$group <- factor(df_enrichment_group$group, levels=sort(unique(df_enrichment_group$group)))
		}
	}
	
	##########################################################
	# restrict those nodes provided in 'ig'
	if(!is(ig,"igraph")){
		if(is(eTerm,'eTerm')){
			ig <- eTerm$g
		}else{
			return(NULL)
		}
	}
	
	if(!fixed){
		ind <- match(V(ig)$term_name, df_enrichment_group$name)
		nodes_query <- V(ig)$name[!is.na(ind)]
		if(is(suppressWarnings(try(subg <- dnet::dDAGinduce(ig, nodes_query, path.mode="all_paths"), TRUE)),"try-error")){
			subg <- NULL
		}
	}else{
		subg <- ig
	}
	##########################################################
	
	if(is(subg,"igraph")){
		gp_code_table <- xOBOcode(g=subg, node.level='term_distance', node.level.value=2, node.label.color='black', node.shape=21, node.size.range=4, edge.color.alpha=0.2, table.base.size=7, table.row.space=2, table.nrow=min(vcount(subg),55), ...)
		gp_code <- gp_code_table$code
		gp_table <- gp_code_table$table
		
		#########################
		## replace those infinite
		df_enrichment_group$or[is.infinite(df_enrichment_group$or)] <- max(df_enrichment_group$or[!is.infinite(df_enrichment_group$or)])
		#########################
		
		ls_df <- split(x=df_enrichment_group[,c("name","zscore","adjp","or")], f=df_enrichment_group$group)
		ls_ig <- lapply(ls_df, function(df_enrichment){
			V(subg)$zscore <- 0
			ind <- match(V(subg)$term_name, df_enrichment$name)
			V(subg)$zscore[!is.na(ind)] <- df_enrichment$zscore[ind[!is.na(ind)]]

			V(subg)$adjp <- 0
			ind <- match(V(subg)$term_name, df_enrichment$name)
			V(subg)$adjp[!is.na(ind)] <- -1*log10(df_enrichment$adjp)[ind[!is.na(ind)]]
	
			V(subg)$or <- 0
			ind <- match(V(subg)$term_name, df_enrichment$name)
			V(subg)$or[!is.na(ind)] <- log2(df_enrichment$or)[ind[!is.na(ind)]]
			subg
		})
		
		if(node.color=="or"){
			node.color.title <- expression(log[2](OR))
			if(is.null(zlim)){
				zlim <- c(0, ceiling(max(log2(df_enrichment_group$or))))
			}
		}else if(node.color=="adjp"){
			node.color.title <- expression(-log[10](FDR))
			if(is.null(zlim)){
				zlim <- c(0, ceiling(max(-1*log10(df_enrichment_group$adjp))))
			}
		}else if(node.color=="zscore"){
			node.color.title <- "Z-score"
			if(is.null(zlim)){
				zlim <- c(0, ceiling(max(df_enrichment_group$zscore)))
			}
		}
		
		if(node.size=="or"){
			node.size.title <- expression(log[2](OR))
			if(is.null(slim)){
				slim <- c(0, ceiling(max(log2(df_enrichment_group$or))))
			}
		}else if(node.size=="adjp"){
			node.size.title <- expression(-log[10](FDR))
			if(is.null(slim)){
				slim <- c(0, ceiling(max(-1*log10(df_enrichment_group$adjp))))
			}
		}else if(node.size=="zscore"){
			node.size.title <- "Z-score"
			if(is.null(slim)){
				slim <- c(0, ceiling(max(df_enrichment_group$zscore)))
			}
		}
		
		gp_data <- xGGnetwork(g=ls_ig, node.color=node.color, node.color.title=node.color.title, colormap=colormap, zlim=zlim, node.size=node.size, node.size.title=node.size.title, slim=slim, node.size.range=node.size.range, edge.color=edge.color, edge.color.alpha=edge.color.alpha, edge.curve=edge.curve, edge.arrow.gap=edge.arrow.gap, ...)
	}
	
	ls_gp <- list(code=gp_code, table=gp_table, data=gp_data, ig=ls_ig)
	
    return(ls_gp)
}
