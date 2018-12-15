#' Function to visualise enrichment results using a D3 plot
#'
#' \code{xEnrichD3} is supposed to visualise enrichment results using a D3 plot. It returns an object of class "htmlwidget".
#'
#' @param eTerm an object of class "eTerm" or "ls_eTerm". Alterntively, it can be a data frame having all these columns (named as 'group','ontology','name','adjp','zscore')
#' @param top_num the number of the top terms (sorted according to adjp). For the eTerm object, if it is 'auto' (for eTerm), only the significant terms (see below FDR.cutoff) will be displayed
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05
#' @param type the D3 type of the plot. It can be "sankey" for sankey network, "force" for force directed network graph, "radial" for radial network and "diagonal" for diagonal network
#' @param colormap short name for the group/ontology colormap
#' @param filename the without-extension part of the name of the output html file. By default, it is 'xEnrichD3'
#' @param ... additional graphic parameters used in networkD3::sankeyNetwork, networkD3::forceNetwork, networkD3::radialNetwork and networkD3::diagonalNetwork
#' @return an object of class "htmlwidget", appended with an "igraph" object
#' @note none
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichD3.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' res <- xEnrichD3(eTerm, type="sankey", width=500, height=500)
#' res <- xEnrichD3(eTerm,type="radial", fontSize=12, nodeColour="steelblue", nodeStroke="fff")
#' res
#' res$ig
#' ig <- xConverter(res$ig, from='igraph', to='igraph_tree')
#'
#' # BiocManager::install('webshot')
#' # webshot::install_phantomjs()
#' # BiocManager::install('r2d3')
#' # r2d3::save_d3_png(res, file='xEnrichD3.png', zoom=2)
#' }

xEnrichD3 <- function(eTerm, top_num=10, FDR.cutoff=0.05, type=c("sankey","force","radial","diagonal"), colormap="ggplot2", filename='xEnrichD3', ...)
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
	type <- match.arg(type)
    
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
			df <- eTerm$df
			
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
	
	## df_nodes
	df_nodes_group <- data.frame(name=unique(df$group), category='group', stringsAsFactors=FALSE)
	df_nodes_ontology <- unique(data.frame(name=df$name, category=df$ontology, stringsAsFactors=FALSE))
	df_nodes <- rbind(df_nodes_group, df_nodes_ontology)
	## category colors
	tmp <- table(df_nodes$category)
	color_category <- xColormap(colormap)(length(tmp))
	names(color_category) <- names(tmp)
	df_nodes$color <- color_category[df_nodes$category]
	
	## df_edges
	df_edges <- df[,c('group','name','zscore')]
	g <- igraph::graph.data.frame(d=df_edges, directed=TRUE, vertices=df_nodes)
	
	if(type=='sankey'){
		d3 <- networkD3::igraph_to_networkD3(g, group=df_nodes$category)
	
		## category colors used by d3
		colors <- paste(color_category, collapse = '", "')
		colourScale <- paste('d3.scaleOrdinal(["', colors, '"])')

		#res <- networkD3::sankeyNetwork(Links=d3$links, Nodes=d3$nodes, Source="source", Target="target", Value="value", NodeID="name", NodeGroup='group', colourScale=colourScale, units="TWh", fontSize=12, nodeWidth=20, sinksRight=F, width=500, height=500)
		res <- networkD3::sankeyNetwork(Links=d3$links, Nodes=d3$nodes, Source="source", Target="target", Value="value", NodeID="name", NodeGroup='group', colourScale=colourScale, fontSize=12, sinksRight=F, ...)

	}else if(type=='force'){
		# D3 JavaScript Force Directed Network Graph
		d3 <- networkD3::igraph_to_networkD3(g, group=df_nodes$category)
		
		# rescale to [1,5]
		d3 <- networkD3::igraph_to_networkD3(g, group=df_nodes$category)
		d3$links$value <- 1 + 4 *(d3$links$value - min(d3$links$value)) / (max(d3$links$value) - min(d3$links$value))
		#linkWidth <- networkD3::JS("function(d) { return Math.sqrt(d.value); }")
		linkWidth <- networkD3::JS("function(d) { return d.value; }")

		## category colors used by d3
		colors <- paste(color_category, collapse = '", "')
		colourScale <- paste('d3.scaleOrdinal(["', colors, '"])')

		#res <- networkD3::forceNetwork(Links=d3$links, Nodes=d3$nodes, Source='source', Target='target', Value="value", NodeID='name', Group='group', colourScale=colourScale, zoom=T, linkDistance=200, linkWidth=linkWidth, legend=T, linkColour = "#666", opacity=0.9, fontSize=15, opacityNoHover=c(0,F)[1], charge=-180)
		res <- networkD3::forceNetwork(Links=d3$links, Nodes=d3$nodes, Source='source', Target='target', Value="value", NodeID='name', Group='group', colourScale=colourScale, zoom=T, linkDistance=200, linkWidth=linkWidth, ...)

	}else if(type %in% c('radial','diagonal')){
		
		if(0){
			#y <- data.frame(group='root', name=unique(df_edges[,1]), zscore=1, stringsAsFactors=F)
			#x <- rbind(df_edges, y)
			#x$name <- gsub('\\W+','_',x$name)
			#phylo <- treeio::get.tree(x[,c(1:3)])
			#phylo$tip.label <- gsub('_',' ',phylo$tip.label)
			#dend <- stats::as.dendrogram(phylo)
			#data <- networkD3::as.radialNetwork(dend)
			#data <- jsonlite::read_json("flare.json", simplifyDataFrame=F)
		}
		
		data <- xConverter(g, from="igraph", to="lol", verbose=F)
		if(type=='radial'){
			#res <- networkD3::radialNetwork(data, fontSize=8, linkColour="#ccc", nodeColour="#fff", nodeStroke="steelblue", textColour="#111", opacity=0.9)
			res <- networkD3::radialNetwork(data, ...)
		}else{
			#res <- networkD3::diagonalNetwork(data, fontSize=8, linkColour="#ccc", nodeColour="#fff", nodeStroke="steelblue", textColour="#111", opacity=0.9)
			res <- networkD3::diagonalNetwork(data, ...)
		}
	}
	
	## save into a html file 
	if(all(!is.null(filename),!is.na(filename),filename!='')){
		filename <- gsub('.html$', '', filename)
		filename <- paste0(filename, ".html")
		res %>% networkD3::saveNetwork(file=filename) 
	}
	
	#library(r2d3)
	#g <- make_graph("Zachary")
	#json <- xConverter(g, from='igraph', to='json')
	#res <- r2d3(json, d3_version=4, script="radialtree.js"
	
	## appended with an "igraph" object
	res$ig <- g
	
	invisible(res)
}

