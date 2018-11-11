#' Function to visualise enrichment results using a forest plot
#'
#' \code{xGT} is supposed to visualise enrichment results using a forest plot. A point is colored by the significance level, and a horizontal line for the 95% confidence interval (CI) of odds ratio (OR; the wider the CI, the less reliable). It returns an object of class "ggplot".
#'
#' @param phylo an object of class "phylo"
#' @param membership a named input vector containing the membership for tips. For this named vector, the element names are clusters, the element values for tip labels. It can be NULL
#' @param layout the visual layout. It can be "rectangular" or "fan"
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param extendto the maximum value of x. By default it is NULL. Only works for the rectangular visual
#' @param offset the offset to the tips. By default it is NULL. Only works for the rectangular visual
#' @param open.angle the angle for the fan visual
#' @return an object of class "ggplot" appended with a phylo object "tree" (appended with a list "cluster")
#' @note none
#' @export
#' @seealso \code{\link{xGT}}
#' @include xGT.r
#' @examples
#' \dontrun{
#' set.seed(825)
#' tree <- ape::rtree(50)
#' gp <- xGT(tree)
#' gp + ggtree::geom_tiplab(size=2,align=T)
#' gp$data
#' gp$tree
#' }

xGT <- function(phylo, membership=NULL, layout=c("rectangular","fan"), colormap="ggplot2", extendto=NULL, offset=NULL, open.angle=0)
{
    layout <- match.arg(layout)

    if(class(phylo)== "dendrogram"){
        #phylo <- ape::as.phylo(stats::as.hclust(phylo))
    }    
    if(class(phylo)== "hclust"){
        #phylo <- ape::as.phylo(phylo)
    }
    
    if(class(phylo) != "phylo"){
    	warnings("The function must apply to a 'phylo', 'hclust' or 'dendrogram' object.\n")
        return(NULL)
    }else{
    	tree <- phylo
    }
    
    flag <- F
    if(is.vector(membership)){
    	if(!is.null(names(membership))){
    		ind <- match(names(membership), tree$tip.label)
    		if(length(!is.na(ind))>1){
    			#warnings("The input vector 'membership' must be a named vector (with names for memberships).\n")
    			flag <- T
    		}
    	}
    }
    
    if(!flag){
		if(layout=='rectangular'){
			p <- ggtree::ggtree(tree, layout="rectangular")
			gp <- p + ggtree::theme_tree2()
		}else if(layout=='fan'){
			p <- ggtree::ggtree(tree, layout="fan", open.angle=open.angle)
			gp <- p + ggtree::geom_tiplab2(align=T,size=2)
		}
    	
		## append 'cluster'
		tree$cluster <- NULL
    	
    }else{
    
		# clusters
		vec_clusters <- sort(unique(membership))
		# color
		col_clusters <- xColormap("ggplot2")(length(vec_clusters))
		names(col_clusters) <- vec_clusters
		# mrca
		mrca_clusters <- col_clusters
    
		if(layout=='rectangular'){
			p <- ggtree::ggtree(tree, layout="rectangular")
		
			if(is.null(extendto)){
				extendto <- max(p$data$x) * 1.1
			}
			if(is.null(offset)){
				offset <- max(p$data$x) * 0.1
			}

			for(i in 1:length(col_clusters)){
				x <- col_clusters[i]
				ind <- match(membership, names(x))
				tip <- as.character(names(membership[!is.na(ind)]))
				node <- ggtree::MRCA(tree, tip=tip)
				mrca_clusters[i] <- node
				p <- p + ggtree::geom_hilight(node, fill=x, extendto=extendto) + ggtree::geom_cladelabel(node, label=names(x), color=x, align=T, offset=offset)
			}
			gp <- p + ggtree::theme_tree2()
	
		}else if(layout=='fan'){
			p <- ggtree::ggtree(tree, layout="fan", open.angle=open.angle)
			for(i in 1:length(col_clusters)){
				x <- col_clusters[i]
				ind <- match(membership, names(x))
				tip <- as.character(names(membership[!is.na(ind)]))
				node <- ggtree::MRCA(tree, tip=tip)
				mrca_clusters[i] <- node
				p <- p + ggtree::geom_hilight(node, fill=x) + ggtree::geom_cladelabel(node, label=names(x), color=x, barsize=0, geom='label')
			}
			gp <- p + ggtree::geom_tiplab2(align=T,size=2)
		}
		
		## append 'cluster'
		tree$cluster <- list(layout=layout, color=col_clusters, mrca=mrca_clusters)
	
	}
	
	# append 'tree'
	gp$tree <- tree
	
	invisible(gp)
}

