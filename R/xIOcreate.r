#' Function to construct an immune ontology
#'
#' \code{xIOcreate} is supposed to construct an immune ontology.
#'
#' @param cutoff.anno an integer specifying the annotation cutoff obove which ontology terms are considered
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param filename the without-extension part of the name of the output file. By default, it is 'xA2GraphML'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a list with 3 components, two ggplot objects (code and table) and an igraph object (ig) having attributes:
#' \itemize{
#'  \item{\code{name}: a node attribute for name (ie node.code)}
#'  \item{\code{term_id}: a node attribute for term_id (the primary source)}
#'  \item{\code{term_name}: a node attribute for term_name (the primary source)}
#'  \item{\code{term_distance}: a node attribute for term_distance (starting with 0, 1, ...)}
#'  \item{\code{term_namespace}: a node attribute for term_namespace (ie IO, Disease, Function, Phenotype)}
#'  \item{\code{anno}: a node attribute as a list storing annotated gene symbols}
#'  \item{\code{net}: a node attribute as a list storing connectivity between genes}
#'  \item{\code{path}: a node attribute as a list storing pathways enriched}
#'  \item{\code{node.label}: a node attribute for the node label in the form of '[nAnno->nNode@nEdge->nPath] term_name'}
#'  \item{\code{node.label.size}: a node attribute for the node label size (decreasing radiately)}
#'  \item{\code{node.label.color}: a node attribute for the node label color (grouped by the subontologies)}
#'  \item{\code{xcoord}: a node attribute for the node x coordinate}
#'  \item{\code{ycoord}: a node attribute for the node y coordinate}
#'  \item{\code{node.code}: a node attribute for the 3-letter node code, 1st letter for the subontology, 2nd for the distance, 3rd for the sequential order at the same distance}
#'  \item{\code{node.table}: a node attribute for the node label in the lookup table}
#'  \item{\code{edge.color}: an edge attribute for the edge color (grouped by the subontologies)}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xIOcreate}}
#' @include xIOcreate.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' 
#' ls_IO <- xIOcreate(RData.location=RData.location)
#' ls_IO$ig
#' ls_IO$code
#' ls_IO$table
#' 
#' lapply(V(ls_IO$ig)[['b5c']]$path, function(x) xA2GraphML(query=x, filename=x, RData.location=RData.location))
#' }

xIOcreate <- function(cutoff.anno=15, colormap='ggplot2', filename='xIOcreate', verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    #########################################
	# DO: immune system disease ('DOID:2914')
    #########################################
	g <- xRDataLoader(RData.customised='ig.DO', RData.location=RData.location)
	node.root <- "DOID:2914"
	neighs.out <- igraph::neighborhood(g, order=igraph::vcount(g), nodes=node.root, mode="out")
	nodeClan <- names(unlist(neighs.out))
	## subg under the node
	subg <- igraph::induced.subgraph(g, vids=nodeClan)
	## dag with anno
	anno <- xRDataLoader(RData.customised='org.Hs.egDO', RData.location=RData.location)
	dag <- XGR::xDAGanno(subg, annotation=anno, path.mode="all_paths")
	vec <- sapply(V(dag)$anno, length)
	nodeAnno <- V(dag)$name[vec >= cutoff.anno]
	## ig under the node (with anno)
	vids <- intersect(nodeClan, nodeAnno)
	ig <- igraph::induced.subgraph(g, vids=vids)
	## ig being simplified
	ig <- XGR::xSimplifyNet(ig)
	## recalculate term_distance
	V(ig)$term_distance <- igraph::distances(ig, v=node.root, to=V(ig), mode="out")[1,] + 1
	## ig appended with term_namespace
	V(ig)$term_namespace <- 'Disease'	
	## ig appended with anno
	ind <- match(V(ig)$name, V(dag)$name)
	V(ig)$anno <- lapply(V(dag)$anno[ind], as.numeric)
	ig_do <- ig
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Disease: a graph of %d nodes and %d edges with %d depth, constructed from the immune part (%d nodes and %d edges) of the whole (%d nodes and %d edges)", igraph::vcount(ig),igraph::ecount(ig),max(V(ig)$term_distance), igraph::vcount(subg),igraph::ecount(subg), igraph::vcount(g),igraph::ecount(g)), appendLF=T)
	}
    
    #########################################
	# GOBP: immune system process ('GO:0002376')
    #########################################
    node.root <- "GO:0002376"
	g <- xRDataLoader(RData.customised='ig.GOBP', RData.location=RData.location)
	neighs.out <- igraph::neighborhood(g, order=igraph::vcount(g), nodes=node.root, mode="out")
	nodeClan <- names(unlist(neighs.out))
	## subg under the node
	subg <- igraph::induced.subgraph(g, vids=nodeClan)
	## dag with anno
	anno <- xRDataLoader(RData.customised='org.Hs.egGOBP', RData.location=RData.location)
	dag <- XGR::xDAGanno(subg, annotation=anno, path.mode="all_paths")
	vec <- sapply(V(dag)$anno, length)
	nodeAnno <- V(dag)$name[vec >= cutoff.anno]
	## ig under the node (with anno)
	vids <- intersect(nodeClan, nodeAnno)
	ig <- igraph::induced.subgraph(g, vids=vids)
	## ig being simplified
	ig <- XGR::xSimplifyNet(ig)
	## recalculate term_distance
	V(ig)$term_distance <- igraph::distances(ig, v=node.root, to=V(ig), mode="out")[1,] + 1
	## ig appended with term_namespace
	V(ig)$term_namespace <- 'Function'	
	## ig appended with anno
	ind <- match(V(ig)$name, V(dag)$name)
	V(ig)$anno <- lapply(V(dag)$anno[ind], as.numeric)
	ig_gobp <- ig
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Function: a graph of %d nodes and %d edges with %d depth, constructed from the immune part (%d nodes and %d edges) of the whole (%d nodes and %d edges)", igraph::vcount(ig),igraph::ecount(ig),max(V(ig)$term_distance), igraph::vcount(subg),igraph::ecount(subg), igraph::vcount(g),igraph::ecount(g)), appendLF=T)
	}
    
    #########################################
	# HPPA: Abnormality of the immune system ('HP:0002715')
    #########################################
    node.root <- "HP:0002715"
	g <- xRDataLoader(RData.customised='ig.HPPA', RData.location=RData.location)
	neighs.out <- igraph::neighborhood(g, order=igraph::vcount(g), nodes=node.root, mode="out")
	nodeClan <- names(unlist(neighs.out))
	## subg under the node
	subg <- igraph::induced.subgraph(g, vids=nodeClan)
	## dag with anno
	anno <- xRDataLoader(RData.customised='org.Hs.egHPPA', RData.location=RData.location)
	dag <- XGR::xDAGanno(subg, annotation=anno, path.mode="all_paths")
	vec <- sapply(V(dag)$anno, length)
	nodeAnno <- V(dag)$name[vec >= cutoff.anno]
	## ig under the node (with anno)
	vids <- intersect(nodeClan, nodeAnno)
	ig <- igraph::induced.subgraph(g, vids=vids)
	## ig being simplified
	ig <- XGR::xSimplifyNet(ig)
	## recalculate term_distance
	V(ig)$term_distance <- igraph::distances(ig, v=node.root, to=V(ig), mode="out")[1,] + 1
	## ig appended with term_namespace
	V(ig)$term_namespace <- 'Phenotype'	
	## ig appended with anno
	ind <- match(V(ig)$name, V(dag)$name)
	V(ig)$anno <- lapply(V(dag)$anno[ind], as.numeric)
	ig_hppa <- ig
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Phenotype: a graph of %d nodes and %d edges with %d depth, constructed from the immune part (%d nodes and %d edges) of the whole (%d nodes and %d edges)", igraph::vcount(ig),igraph::ecount(ig),max(V(ig)$term_distance), igraph::vcount(subg),igraph::ecount(subg), igraph::vcount(g),igraph::ecount(g)), appendLF=T)
	}
    
    #########################################
	# MP: abnormal immune system physiology ('MP:0001790')
    #########################################
    node.root <- "MP:0001790"
	g <- xRDataLoader(RData.customised='ig.MP', RData.location=RData.location)
	neighs.out <- igraph::neighborhood(g, order=igraph::vcount(g), nodes=node.root, mode="out")
	nodeClan <- names(unlist(neighs.out))
	## subg under the node
	subg <- igraph::induced.subgraph(g, vids=nodeClan)
	## dag with anno
	anno <- xRDataLoader(RData.customised='org.Hs.egMP', RData.location=RData.location)
	dag <- XGR::xDAGanno(subg, annotation=anno, path.mode="all_paths")
	vec <- sapply(V(dag)$anno, length)
	nodeAnno <- V(dag)$name[vec >= cutoff.anno]
	## ig under the node (with anno)
	vids <- intersect(nodeClan, nodeAnno)
	ig <- igraph::induced.subgraph(g, vids=vids)
	## ig being simplified
	ig <- XGR::xSimplifyNet(ig)
	## recalculate term_distance
	V(ig)$term_distance <- igraph::distances(ig, v=node.root, to=V(ig), mode="out")[1,] + 1
	## ig appended with term_namespace
	V(ig)$term_namespace <- 'MPhenotype'	
	## ig appended with anno
	ind <- match(V(ig)$name, V(dag)$name)
	V(ig)$anno <- lapply(V(dag)$anno[ind], as.numeric)
	ig_mp <- ig
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("MPhenotype: a graph of %d nodes and %d edges with %d depth, constructed from the immune part (%d nodes and %d edges) of the whole (%d nodes and %d edges)", igraph::vcount(ig),igraph::ecount(ig),max(V(ig)$term_distance), igraph::vcount(subg),igraph::ecount(subg), igraph::vcount(g),igraph::ecount(g)), appendLF=T)
	}
    
    #########################################
	# REACTOME: Immune System ('R-HSA-168256')
    #########################################
    node.root <- "R-HSA-168256"
	g <- xRDataLoader(RData.customised='ig.REACTOME', RData.location=RData.location)
	neighs.out <- igraph::neighborhood(g, order=igraph::vcount(g), nodes=node.root, mode="out")
	nodeClan <- names(unlist(neighs.out))
	## subg under the node
	subg <- igraph::induced.subgraph(g, vids=nodeClan)
	## dag with anno
	anno <- xRDataLoader(RData.customised='org.Hs.egREACTOME', RData.location=RData.location)
	dag <- XGR::xDAGanno(subg, annotation=anno, path.mode="all_paths")
	vec <- sapply(V(dag)$anno, length)
	nodeAnno <- V(dag)$name[vec >= cutoff.anno]
	## ig under the node (with anno)
	vids <- intersect(nodeClan, nodeAnno)
	ig <- igraph::induced.subgraph(g, vids=vids)
	## ig being simplified
	ig <- XGR::xSimplifyNet(ig)
	## recalculate term_distance
	V(ig)$term_distance <- igraph::distances(ig, v=node.root, to=V(ig), mode="out")[1,] + 1
	## ig appended with term_namespace
	V(ig)$term_namespace <- 'Pathway'	
	## ig appended with anno
	ind <- match(V(ig)$name, V(dag)$name)
	V(ig)$anno <- lapply(V(dag)$anno[ind], as.numeric)
	ig_reactome <- ig
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Pathway: a graph of %d nodes and %d edges with %d depth, constructed from the immune part (%d nodes and %d edges) of the whole (%d nodes and %d edges)", igraph::vcount(ig),igraph::ecount(ig),max(V(ig)$term_distance), igraph::vcount(subg),igraph::ecount(subg), igraph::vcount(g),igraph::ecount(g)), appendLF=T)
	}
    
    
    #########################################
    #########################################
    #ls_ig <- list(ig_do, ig_gobp, ig_reactome, ig_hppa, ig_mp)	
    #ls_ig <- list(ig_do, ig_gobp, ig_reactome, ig_hppa)
    ls_ig <- list(ig_do, ig_gobp, ig_hppa)
	ig_io <- XGR::xCombineNet(ls_ig, combineBy='union', attrBy="intersect", keep.all.vertices=TRUE)
	
	## anno: replace EntrezGenes with gene symbols	
	EG <- xRDataLoader(RData.customised=paste('org.Hs.eg', sep=''), RData.location=RData.location, verbose=verbose)
	allGeneID <- EG$gene_info$GeneID
	allSymbol <- as.vector(EG$gene_info$Symbol)
	V(ig_io)$anno <- lapply(V(ig_io)$anno,function(x){
		ind <- match(x, allGeneID)
		allSymbol[ind]
	})
	
	## extract the names of the subontology nodes
	subonto.names <- V(ig_io)$name[V(ig_io)$term_distance==1]
	
	## add the root node 'IO'
	ig_io <- ig_io %>% igraph::add_vertices(1, attr=list(name='IO', term_id='IO', term_name='Immune Ontology', term_distance=0, anno=list(unique(unlist(V(ig_io)$anno)))), term_namespace='IO') 
	## add the edges coming from the root node 'IO'
	for (k in 1:length(subonto.names)){
		ig_io <- ig_io %>% igraph::add_edges(c('IO',subonto.names[k]))
	}
	
	###########
	## add net
	###########
	networks <- c("KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal")[1:5]
	ls_networks <- lapply(networks, function(network){
		g <- xDefineNet(network=network, RData.location=RData.location)	
	})
	kegg <- xCombineNet(ls_networks, combineBy='union', attrBy="intersect", keep.all.vertices=TRUE, verbose=verbose)
	if(0){
		g <- ls_networks[[2]]
		glayout <- igraph::layout_with_kk(g)
		V(g)$xcoord <- glayout[,1]
		V(g)$ycoord <- glayout[,2]
		gp <- xA2Net(g, node.size.range=0.5, node.xcoord='xcoord', node.ycoord='ycoord', edge.color="black", edge.color.alpha=0.1, edge.curve=0, edge.arrow.gap=0.005)
	}
	V(ig_io)$net <- lapply(V(ig_io)$anno,function(x){
		ind <- match(V(kegg)$name, x)
		nodes_query <- V(kegg)$name[!is.na(ind)]
		if(0){
			## including incoming neighbors
			neighs.out <- igraph::neighborhood(kegg, order=1, nodes=nodes_query, mode="in")
			nodes_query <- unique(names(unlist(neighs.out)))
		}
		subg <- dnet::dNetInduce(kegg, nodes_query=nodes_query, knn=0, remove.loops=TRUE, largest.comp=FALSE, min.comp.size=2)
	})
	###########

	###########
	## add path
	###########
	anno <- V(ig_io)$anno
	names(anno) <- V(ig_io)$name
	ind <- which(V(ig_io)$term_distance>1)
	anno <- anno[ind]
	ls_df <- lapply(1:length(anno), function(i){
		x <- anno[[i]]
		data.frame(Symbol=x, Code=rep(names(anno)[i],length(x)), stringsAsFactors=FALSE)
	})
	annotation.file <- do.call(rbind, ls_df)
	## analyse AA.path
	adjp <- fc <- NULL
 	AA.path <- xRDataLoader(RData.customised="AA.path", RData.location=RData.location)
	ls_detail <- AA.path$detail
	ls_df <- lapply(1:length(ls_detail), function(i){
		x <- ls_detail[[i]]$ig
		data.file <- V(x)$Symbol
		# perform enrichment analysis
		eTerm <- xEnricherYours(data.file=data.file, annotation.file=annotation.file, size.range=c(15,3000), min.overlap=15, test="fisher", p.adjust.method=c("BH","bonferroni")[2], verbose=FALSE)
		df_eTerm <- xEnrichViewer(eTerm, top_num='all')
		df_eTerm <- df_eTerm %>% dplyr::filter(adjp<1e-2 & fc>=2)
		data.frame(path=rep(names(ls_detail)[i],nrow(df_eTerm)), code=df_eTerm$name, stringsAsFactors=FALSE)
	})
	df_path <- do.call(rbind, ls_df)
    ls_path <- split(x=df_path$path, f=df_path$code)
    ## append path
    V(ig_io)$path <- NA
    ind <- match(V(ig_io)$name, names(ls_path))
    V(ig_io)$path[!is.na(ind)] <- ls_path[ind[!is.na(ind)]]
    ###########
    
	## adjust term_name
	V(ig_io)$term_name <- sapply(V(ig_io)$term_name, function(x) {
		paste0(toupper(substr(x,1,1)), substr(x,2,nchar(x)))
	})
	
	## add node.label
	#V(ig_io)$node.label <- paste0('[',V(ig_io)$term_id, '] ', V(ig_io)$term_name)
	nAnno <- sapply(V(ig_io)$anno,length)
	nNode <- sapply(V(ig_io)$net,igraph::vcount)
	nEdge <- sapply(V(ig_io)$net,igraph::ecount)
	nPath <- sapply(V(ig_io)$path, function(x) sum(!is.na(x)))
	#V(ig_io)$node.label <- paste0('[',nAnno, '] ', V(ig_io)$term_name)
	V(ig_io)$node.label <- paste0('[',nAnno, '->', nNode, '@', nEdge, '->', nPath, '] ', V(ig_io)$term_name)
	
	## add node.label.size
	V(ig_io)$node.label.size <- ifelse(V(ig_io)$term_distance==0, 8,
								ifelse(V(ig_io)$term_distance==1, 8, 
								ifelse(V(ig_io)$term_distance==2, 2, 
								1)))
	### adjust for nodes with distance>=3 to the size ranged from 1.2 to 2 (decreased by 0.25 at each step)
	ind <- which(V(ig_io)$term_distance >2)
	tmp <- (2 - (V(ig_io)$term_distance[ind] -2) * 0.25)
	tmp[tmp<=1.25] <- 1.25
	V(ig_io)$node.label.size[ind] <- tmp
	
	## add node.label.color and edge.color
	V(ig_io)$node.label.color <- 'black'
	E(ig_io)$edge.color <- 'black'
	for (k in 1:length(subonto.names)){
		neighs.out <- igraph::neighborhood(ig_io, order=igraph::vcount(ig_io), nodes=subonto.names[k], mode="out")
		vids <- names(unlist(neighs.out))
		tmp_color <- xColormap(colormap)(3)[k]
		V(ig_io)[vids]$node.label.color <- tmp_color
		E(ig_io)[vids %--% vids]$edge.color <- tmp_color
		### adjust alpha to the edge color
		#E(ig_io)[vids %--% vids]$edge.color <- grDevices::adjustcolor(tmp_color, alpha.f=0)
	}
	
	## add xcoord and ycoord
	vec_scale <- c(0.6, 1.05, 1, 0.95)
	layouts <- lapply(1:length(ls_ig), function(i){
		x <- ls_ig[[i]]
		glayout <- igraph::layout_as_tree(x,root=dnet::dDAGroot(x),circular=TRUE,flip.y=TRUE)
		node.xcoord <- glayout[,1]
		node.ycoord <- glayout[,2]
		
		## scale into [-1,1]
		if(max(node.xcoord) != min(node.xcoord)){
			node.xcoord <- (node.xcoord - min(node.xcoord)) / (max(node.xcoord) - min(node.xcoord)) * 2 - 1
		}
		if(max(node.ycoord) != min(node.ycoord)){
			node.ycoord <- (node.ycoord - min(node.ycoord)) / (max(node.ycoord) - min(node.ycoord)) * 2 - 1
		}
		
		## shift
		if(i==1){
			node.xcoord <- node.xcoord * vec_scale[i]
			node.ycoord <- node.ycoord * vec_scale[i] +1			
		}else if(i==2){
			node.xcoord <- node.xcoord * vec_scale[i] -1
			node.ycoord <- node.ycoord * 1
		}else if(i==3){
			node.xcoord <- node.xcoord * vec_scale[i] +1
			node.ycoord <- node.ycoord * 1 -0.25
		}
		
		if(0){
			if(i==1){
				node.xcoord <- node.xcoord * vec_scale[i] -1
				node.ycoord <- node.ycoord * vec_scale[i] +1
			}else if(i==2){
				node.xcoord <- node.xcoord * vec_scale[i] +1
				node.ycoord <- node.ycoord * vec_scale[i] +1
			}else if(i==3){
				node.xcoord <- node.xcoord * vec_scale[i] -1
				node.ycoord <- node.ycoord * vec_scale[i] -1
			}else if(i==4){
				node.xcoord <- node.xcoord * vec_scale[i] +1
				node.ycoord <- node.ycoord * vec_scale[i] -1
			}
		}
		
		return(cbind(node.xcoord, node.ycoord))
	})	
	glayout <- do.call(rbind, layouts)
	## add the root node coordinates
	glayout <- rbind(glayout, c(0.1,0.25))
	V(ig_io)$xcoord <- glayout[,1]
	V(ig_io)$ycoord <- glayout[,2]
	
	## create node code and table
	ls_IO <- xA2NetCode(g=ig_io, node.level='term_distance', node.level.value=1, node.label.size='node.label.size', node.label.color='node.label.color', node.label.alpha=1, node.label.padding=0, node.label.arrow=0, node.label.force=0, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.size.range=2, colormap='white-white', edge.size=0.25, edge.color="edge.color", edge.color.alpha=0.4, edge.curve=0.1, edge.arrow=0, edge.arrow.gap=0.01, node.table='node.label', node.table.wrap=45, table.nrow=66, root.code='IO')
	
    if(!is.null(filename)){
		############################
		outputfile <- paste0(filename, ".pdf")
		pdf(outputfile, useDingbats=FALSE, width=8, height=8)
		print(ls_IO$code + coord_equal(ratio=1))
		print(ls_IO$table)
		dev.off()
    }
    
    V(ls_IO$ig)$name <- V(ls_IO$ig)$node.code
    
    if(1){
		## output: org.Hs.egIO.RData
		### set_info
		node_attr <- c('name','term_name','term_namespace','term_distance')
		set_info <- igraph::get.data.frame(ls_IO$ig, what="vertices")[,node_attr]
		colnames(set_info) <- c("setID","name","namespace","distance")
		### gs
		gs <- V(ls_IO$ig)$anno
		names(gs) <- V(ls_IO$ig)$name
		EG <- xRDataLoader(RData.customised=paste('org.Hs.eg', sep=''), RData.location=RData.location, verbose=verbose)
		allGeneID <- EG$gene_info$GeneID
		allSymbol <- as.vector(EG$gene_info$Symbol)
		gs <- lapply(gs, function(x){
			ind <- match(x, allSymbol)
			allGeneID[ind]
		})
		### GS
		GS <- list(set_info=set_info, gs=gs)
		class(GS) <- 'GS'
		save(list=c("GS"), file='org.Hs.egIO.RData', compress="xz")
		
		## output: IO.RData
		save(list=c("ls_IO"), file='IO.RData', compress="xz")
    }
  
    invisible(ls_IO)
}






