#' Function to visualise enrichment results using a tree-like circular plot
#'
#' \code{xEnrichCtree} is supposed to visualise enrichment results using a tree-like circular plot.
#'
#' @param eTerm an object of class "eTerm" or "ls_eTerm". Alterntively, it can be a data frame having all these columns ('name','adjp','or','zscore'; 'group' optionally)
#' @param ig an object of class "igraph" with node attribute 'name'. It could be a 'phylo' object converted to. Note: the leave labels would be the node attribute 'name' unless the node attribute 'label' is explicitely provided
#' @param node.color which statistics will be used for node coloring. It can be "or" for the odds ratio, "adjp" for adjusted p value (FDR) and "zscore" for enrichment z-score
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param zlim the minimum and maximum values for which colors should be plotted
#' @param node.size which statistics will be used for node size. It can be "or" for the odds ratio, "adjp" for adjusted p value (FDR) and "zscore" for enrichment z-score
#' @param slim the minimum and maximum values for which sizes should be plotted
#' @param node.size.range the range of actual node size
#' @param group.gap the gap between group circles. Only works when multiple groups provided
#' @param group.color the color of group circles. Only works when multiple groups provided
#' @param group.label.size the size of group circle labelling. Always a sequential integer located at the top middle. Only works when multiple groups provided
#' @param group.label.color the color of group circle labelling. Only works when multiple groups provided
#' @param ... additional graphic parameters used in xCtree
#' @return
#' a ggplot2 object appended with 'ig', 'data' which should contain columns 'x','y', 'leaf' (T/F), 'name' (the same as V(ig)$name), 'tipid' (tip id), 'label' (if not given in ig, a 'name' varient), 'angle' and 'hjust' (assist in leave label orientation), and 'data_enrichment' (enrichment results for tips)
#' @note none
#' @export
#' @seealso \code{\link{xCtree}}
#' @include xEnrichCtree.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' 
#' # load the atlas of AA pathways
#' AA.template <- xRDataLoader("AA.template", RData.location=RData.location)
#' # consensus tree
#' ig <- AA.template$consensus$ig
#' 
#' # enrichment analysis using AA pathways
#' input <- xRDataLoader('Haploid_regulators_all', RData.location=RData.location)
#' data <- subset(input, Phenotype=="AKT")
#' genes <- data$Gene[data$FDR<0.05]
#' background <- data$Gene
#' eTerm <- xEnricherGenes(genes, background=background, ontology="AA", min.overlap=5, test="fisher", RData.location=RData.location)
#' 
#' # circular visualisation of enriched AA pathways
#' gp <- xEnrichCtree(eTerm, ig)
#' 
#' ###############################
#' # advanced use: multiple groups
#' # enrichment analysis using AA pathways
#' Haploid <- subset(input, FDR<0.05)
#' ls_group <- split(x=Haploid$Gene, f=Haploid$Phenotype)
#' background <- unique(input$Gene)
#' ls_eTerm <- xEnricherGenesAdv(ls_group, background=background, ontologies="AA", test="fisher", min.overlap=5, RData.location=RData.location)
#' 
#' # circular visualisation of enriched AA pathways
#' gp <- xEnrichCtree(ls_eTerm, ig)
#' }

xEnrichCtree <- function(eTerm, ig, node.color=c("zscore","adjp","or"), colormap="grey-orange-darkred", zlim=NULL, node.size=c("adjp","zscore","or"), slim=NULL, node.size.range=c(0.5,4.5), group.gap=0.08, group.color="lightblue", group.label.size=2, group.label.color="black", ...)
{
	## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    node.color <- match.arg(node.color)
    node.size <- match.arg(node.size)
    
    if(is.null(eTerm)){
        warnings("There is no enrichment in the 'eTerm' object.\n")
        return(NULL)
    }
    
	if(class(ig)!="igraph"){
		warnings("The 'ig' object must be provided.\n")
		return(NULL)
	}else{
		if(!("name" %in% igraph::vertex_attr_names(ig))){
			warnings("The 'ig' object must contain a node attribute 'name'.\n")
			return(NULL)
		}
	}
    
    if(class(eTerm)=='eTerm'){
		df_enrichment_group <- xEnrichViewer(eTerm, top_num="all")
		df_enrichment_group$group <- 'group'
	}else if(class(eTerm)=='ls_eTerm' | class(eTerm)=='data.frame'){
	
		if(class(eTerm)=='ls_eTerm'){
			df_enrichment_group <- eTerm$df
			
		}else if(class(eTerm)=='data.frame'){
			if(all(c('group','name','adjp','or','zscore') %in% colnames(eTerm))){
				df_enrichment_group <- eTerm[,c('group','name','adjp','or','zscore')]
			}else if(all(c('name','adjp','or','zscore') %in% colnames(eTerm))){
				df_enrichment_group <- eTerm[,c('name','adjp','or','zscore')]
				df_enrichment_group$group <- 'group'
			}
		}
	}
	
	if(class(df_enrichment_group$group)=='factor'){
		if(length(unique(df_enrichment_group$group)) != length(levels(df_enrichment_group$group))){
			df_enrichment_group$group <- factor(df_enrichment_group$group, levels=sort(unique(df_enrichment_group$group)))
		}
	}
	
	gp <- NULL
	
	if(class(ig)=="igraph"){
		
		#########################
		## replace those infinite
		df_enrichment_group$or[is.infinite(df_enrichment_group$or)] <- max(df_enrichment_group$or[!is.infinite(df_enrichment_group$or)])
		#########################
		
		# convert 'ig' into 'ls_igg' by group
		ls_df <- split(x=df_enrichment_group[,c("name","zscore","adjp","or")], f=df_enrichment_group$group)
		ls_igg <- lapply(ls_df, function(df_enrichment){
			igg <- ig
			V(igg)$zscore <- 0
			ind <- match(V(igg)$name, df_enrichment$name)
			V(igg)$zscore[!is.na(ind)] <- df_enrichment$zscore[ind[!is.na(ind)]]

			V(igg)$adjp <- 0
			ind <- match(V(igg)$name, df_enrichment$name)
			V(igg)$adjp[!is.na(ind)] <- -1*log10(df_enrichment$adjp)[ind[!is.na(ind)]]
	
			V(igg)$or <- 0
			ind <- match(V(igg)$name, df_enrichment$name)
			V(igg)$or[!is.na(ind)] <- log2(df_enrichment$or)[ind[!is.na(ind)]]
			igg
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
		
		x <- y <- leaf <- x_group <- y_group <- group <- group_id <- color <- size <- NULL
		
		gp <- xCtree(ig, leave.label.orientation='inwards', ...)
		df_tips <- subset(gp$data, leaf==T)
		
		ls_df <- lapply(1:length(ls_igg), function(i){
			g <- ls_igg[[i]]

			group <- paste0(i, ": ", names(ls_igg)[i])			
			x_group <- 1.05+(i-1)*group.gap
			y_group <- 1.05+(i-1)*group.gap

			ind <- match(V(g)$name, df_tips$name)
			if(node.color=="or"){
				color <- V(g)$or[!is.na(ind)]
			}else if(node.color=="adjp"){
				color <- V(g)$adjp[!is.na(ind)]
			}else if(node.color=="zscore"){
				color <- V(g)$zscore[!is.na(ind)]
			}
			color[color<=zlim[1]] <- zlim[1]
			color[color>=zlim[2]] <- zlim[2]
		
			if(node.size=="or"){
				size <- V(g)$or[!is.na(ind)]
			}else if(node.size=="adjp"){
				size <- V(g)$adjp[!is.na(ind)]
			}else if(node.size=="zscore"){
				size <- V(g)$zscore[!is.na(ind)]
			}
			size[size<=slim[1]] <- slim[1]
			size[size>=slim[2]] <- slim[2]
		
			df <- data.frame(df_tips[ind[!is.na(ind)],], color=color, size=size, group=group, group_id=i, x_group=x_group, y_group=y_group, stringsAsFactors=F)
		})
		df <- do.call(rbind, ls_df)
		
		# if multiple groups
		if(length(ls_igg)>1){
			## add polygon
			gp <- gp + geom_polygon(data=df, aes(x=x*x_group, y=y*y_group, group=group, linetype=df$group),color=group.color,fill="transparent") + guides(linetype=guide_legend(title="Circles\n(outwards)",keywidth=0.7, keyheight=0.7))
			
			## add polygon labelling
			tipid <- NULL
			gp <- gp + geom_text(data=subset(df,tipid==1), aes(x=x*x_group+0.05, y=y*y_group, label=group_id),color=group.label.color,size=group.label.size)
		}		
		
		gp <- gp + geom_point(data=df, aes(x=x*x_group, y=y*y_group, color=color, size=size), alpha=1)
		
		gp <- gp + scale_size_continuous(limits=slim, range=node.size.range, guide=guide_legend(node.size.title,title.position="top",ncol=1)) + scale_colour_gradientn(colors=xColormap(colormap)(64), limits=zlim, guide=guide_colorbar(title=node.color.title,title.position="top",barwidth=0.5))
		
		gp$data_enrichment <- df
	}
	
    return(gp)
}
