#' Function to visualise enrichment results using a direct acyclic graph (DAG)
#'
#' \code{xEnrichDAGplot} is supposed to visualise enrichment results using a direct acyclic graph (DAG) with node colorings. By default, significant terms (of interest) are highlighted by box-shaped nodes, the others by ellipse nodes. It returns an object of class 'Ragraph'.
#'
#' @param eTerm an object of class "eTerm"
#' @param top_num the number of the top terms (sorted according to FDR or adjusted p-values). If it is 'auto', only the significant terms (FDR < 0.05) will be displayed
#' @param ig the igraph object. If provided, only those terms within it will be visualised. By default, it is NULL meaning no surch restriction
#' @param displayBy which statistics will be used for displaying. It can be "fc" for enrichment fold change (by default), "adjp" or "fdr" for adjusted p value (or FDR), "pvalue" for p value, "zscore" for enrichment z-score
#' @param path.mode the mode of paths induced by nodes in query. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param height a numeric value specifying the height of device
#' @param width a numeric value specifying the width of device
#' @param margin margins as units of length 4 or 1
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/data values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param colorbar logical to indicate whether to append a colorbar. If data is null, it always sets to false
#' @param colorbar.fraction the relative fraction of colorbar block against the device size
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param layout.orientation the orientation of the DAG layout. It can be one of "left_right" for the left-right layout (viewed from the DAG root point), "top_bottom" for the top-bottom layout, "bottom_top" for the bottom-top layout, and "right_left" for the right-left layout
#' @param node.info tells the ontology term information used to label nodes. It can be one of "none" for no node labeling, "term_id" for using Term ID, "term_name" for using Term Name (the first 15 characters), "both" for using both of Term ID and Name (the first 15 characters), and "full_term_name" for using the full Term Name
#' @param wrap.width a positive integer specifying wrap width of Term Name
#' @param graph.node.attrs a list of global node attributes. These node attributes will be changed globally. See 'Note' below for details on the attributes
#' @param graph.edge.attrs a list of global edge attributes. These edge attributes will be changed globally. See 'Note' below for details on the attributes
#' @param node.attrs a list of local edge attributes. These node attributes will be changed locally; as such, for each attribute, the input value must be a named vector (i.e. using Term ID as names). See 'Note' below for details on the attributes
#' @return 
#' An object of class 'Ragraph'
#' @note
#' A list of global node attributes used in "graph.node.attrs":
#' \itemize{
#' \item{"shape": the shape of the node: "circle", "rectangle", "rect", "box" and "ellipse"}
#' \item{"fixedsize": the logical to use only width and height attributes. By default, it sets to true for not expanding for the width of the label}
#' \item{"fillcolor": the background color of the node}
#' \item{"color": the color for the node, corresponding to the outside edge of the node}
#' \item{"fontcolor": the color for the node text/labelings}
#' \item{"fontsize": the font size for the node text/labelings}
#' \item{"height": the height (in inches) of the node: 0.5 by default}
#' \item{"width": the width (in inches) of the node: 0.75 by default}
#' \item{"style": the line style for the node: "solid", "dashed", "dotted", "invis" and "bold"}
#' }
#' A list of global edge attributes used in "graph.edge.attrs":
#' \itemize{
#' \item{"color": the color of the edge: gray by default}
#' \item{"weight": the weight of the edge: 1 by default}
#' \item{"style": the line style for the edge: "solid", "dashed", "dotted", "invis" and "bold"}
#' }
#' A list of local node attributes used in "node.attrs" (only those named Term IDs will be changed locally!):
#' \itemize{
#' \item{"label": a named vector specifying the node text/labelings}
#' \item{"shape": a named vector specifying the shape of the node: "circle", "rectangle", "rect", "box" and "ellipse"}
#' \item{"fixedsize": a named vector specifying whether it sets to true for not expanding for the width of the label}
#' \item{"fillcolor": a named vector specifying the background color of the node}
#' \item{"color": a named vector specifying the color for the node, corresponding to the outside edge of the node}
#' \item{"fontcolor": a named vector specifying the color for the node text/labelings}
#' \item{"fontsize": a named vector specifying the font size for the node text/labelings}
#' \item{"height": a named vector specifying the height (in inches) of the node: 0.5 by default}
#' \item{"width": a named vector specifying the width (in inches) of the node: 0.75 by default}
#' \item{"style": a named vector specifying the line style for the node: "solid", "dashed", "dotted", "invis" and "bold"}
#' }
#' @export
#' @seealso \code{\link{xEnricherGenes}}, \code{\link{xEnricherSNPs}}, \code{\link{xEnrichViewer}}
#' @include xEnrichDAGplot.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' 
#' # 1) load eQTL mapping results: cis-eQTLs significantly induced by IFN
#' cis <- xRDataLoader(RData.customised='JKscience_TS2A', RData.location=RData.location)
#' ind <- which(cis$IFN_t > 0 & cis$IFN_fdr < 0.05)
#' df_cis <- cis[ind, c('variant','Symbol','IFN_t','IFN_fdr')]
#' data <- df_cis$variant
#' 
#' # 2) Enrichment analysis using Experimental Factor Ontology (EFO)
#' # Considering LD SNPs and respecting ontology tree
#' eTerm <- xEnricherSNPs(data, ontology="EF", include.LD="EUR", LD.r2=0.8, ontology.algorithm="lea", RData.location=RData.location)
#'
#' # 3) DAG plot of enrichment results
#' agDAG <- xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc", node.info=c("full_term_name"))
#' ## modify node labels
#' xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc", node.info=c("full_term_name"), graph.node.attrs=list(fontsize=25,fontcolor="blue",color="transparent"))
#' ## modify node shapes
#' xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc", node.info=c("full_term_name"), graph.node.attrs=list(fixedsize=FALSE,shape=c("ellipse","box","circle","plaintext")[2]))
#' ## further modify edge color
#' xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc", node.info=c("full_term_name"), graph.node.attrs=list(fontsize=25), graph.edge.attrs=list(color="black"))
#'
#' # 4) hide labels for ellipse nodes
#' library(Rgraphviz)
#' name_nodes <- sapply(AgNode(agDAG), name)
#' shape_nodes <- sapply(AgNode(agDAG), shape)
#' names(shape_nodes) <- name_nodes
#' ind <- which(shape_nodes=='ellipse')
#' label_nodes <- rep('', length(ind))
#' names(label_nodes) <- name_nodes[ind]
#' xEnrichDAGplot(eTerm, top_num="auto", displayBy="fc", node.info=c("full_term_name"), node.attrs=list(label=label_nodes, shape=shape_nodes))
#' 
#' }

xEnrichDAGplot <- function(eTerm, top_num=10, ig=NULL, displayBy=c("fc","adjp","fdr","zscore","pvalue"), path.mode=c("all_paths","shortest_paths","all_shortest_paths"), height=7, width=7, margin=rep(0.1,4), colormap=c("yr","bwr","jet","gbr","wyr","br","rainbow","wb","lightyellow-orange"), ncolors=40, zlim=NULL, colorbar=T, colorbar.fraction=0.1, newpage=T, layout.orientation=c("top_bottom","left_right","bottom_top","right_left"), node.info=c("none", "term_id", "term_name", "both", "full_term_name"), wrap.width=NULL, graph.node.attrs=NULL, graph.edge.attrs=NULL, node.attrs=NULL)
{
    
    displayBy <- match.arg(displayBy)
    path.mode <- match.arg(path.mode)
    layout.orientation <- match.arg(layout.orientation)
    node.info<- match.arg(node.info)
    
    if(is.null(eTerm)){
        warnings("There is no enrichment in the 'eTerm' object.\n")
        return(NULL)
    }
    
    if(class(eTerm)[1]=="eTerm"){
	
		## when 'auto', will keep the significant terms
		df <- xEnrichViewer(eTerm, top_num="all")
		if(top_num=='auto'){
			top_num <- sum(df$adjp<0.05)
			if(top_num<=1){
				top_num <- sum(df$adjp<0.1)
			}
			if(top_num <= 1){
				top_num <- 10
			}
		}
		df <- xEnrichViewer(eTerm, top_num=top_num, sortBy="adjp")
		nodes_query <- rownames(df)
		
		##########################################################
		# restrict those nodes provided in 'ig'
		if(class(ig)=="igraph"){
			nodes_query_tmp <- intersect(nodes_query, V(ig)$name)
			if(length(nodes_query_tmp)>0){
				nodes_query <- nodes_query_tmp
			}
		}
		##########################################################
				
		g <- eTerm$g
		subg <- dnet::dDAGinduce(g, nodes_query, path.mode=path.mode)
	
		if(displayBy=='adjp' | displayBy=='fdr'){
			data <- -1*log10(df$adjp)
		}else if(displayBy=='fc'){
			data <- df$fc
		}else if(displayBy=='pvalue'){
			data <- -1*log10(df$pvalue)
		}else if(displayBy=='zscore'){
			data <- df$zscore
		}
	
		## for data
		names(data) <- rownames(df)
		data <- data[V(subg)$name]
	}else{
		stop("Cannot find an 'eTerm' object.\n")
	}
	
	
	## for globally graph.edge.attrs
	if(is.null(graph.edge.attrs)){
		#graph.edge.attrs <- list(color="black",weight=1,style="solid")
	}
	
	## for locally node.attrs
	if(is.null(node.attrs)){
		hightlighted.shape <- rep("box", length(nodes_query))
		names(hightlighted.shape) <- V(subg)[nodes_query]$name
		node.attrs <- list(shape=hightlighted.shape)
	}
	
	agDAG <- visDAG(g=subg, data=data, height=height, width=width, margin=margin, colormap=colormap, ncolors=ncolors, zlim=zlim, colorbar=colorbar, colorbar.fraction=colorbar.fraction, newpage=newpage, layout.orientation=layout.orientation, node.info=node.info, numChar=wrap.width, graph.node.attrs=graph.node.attrs, graph.edge.attrs=graph.edge.attrs, node.attrs=node.attrs)
	
	#slotNames(agDAG)
	#str(agDAG@AgNode[[1]])
	#name_nodes <- sapply(agDAG@AgNode,name)
	#shape_nodes <- sapply(agDAG@AgNode,shape)
	#txtLabelText_nodes <- sapply(lapply(agDAG@AgNode,txtLabel), labelText)
	
	invisible(agDAG)
}
