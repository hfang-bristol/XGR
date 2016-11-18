#' Function to visualise comparative enrichment results using a direct acyclic graph (DAG)
#'
#' \code{xEnrichDAGplotAdv} is supposed to visualise the comparative enrichment results (see the function \code{\link{xEnrichCompare}})) using a direct acyclic graph (DAG). Nodes/terms can be colored according to how many times being called significant. If two enrichment results are compared, node names are prefixed with the form of 'x1-x2', where x1 is for result 1 and x2 for result 2 (the value for x1 or x2 can be '0' for being insignificant, and '1' for being significant). It takes input an 'ggplot' object (with two componets alreadly appended 'g' and 'data'), and returns an object of class 'Ragraph'.
#'
#' @param ggplot an object "ggplot" (resulting from \code{\link{xEnrichCompare}})
#' @param displayBy which statistics will be used for displaying. It can be "nSig" for how many times being called significant (by default), "none" for no color-coding on nodes/terms
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
#' @param node.info tells the ontology term information used to label nodes. It can be "term_id" for using Term ID, "term_name" for using Term Name, 'none' for no labellings
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
#' @seealso \code{\link{xEnrichCompare}}
#' @include xEnrichDAGplotAdv.r
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
#' # 2a) Without considering LD SNPs and without respecting ontology tree
#' eTerm_noLD_noTree <- xEnricherSNPs(data, ontology="EF_disease", include.LD=NA, ontology.algorithm="none", RData.location=RData.location)
#' # 2b) Without considering LD SNPs but respecting ontology tree
#' eTerm_noLD_Tree <- xEnricherSNPs(data, ontology="EF_disease", include.LD=NA, ontology.algorithm="lea", RData.location=RData.location)
#' # 2c) Considering LD SNPs but without respecting ontology tree
#' eTerm_LD_noTree <- xEnricherSNPs(data, ontology="EF_disease", include.LD="EUR", LD.r2=0.8, ontology.algorithm="none", RData.location=RData.location)
#' # 2d) Considering LD SNPs and respecting ontology tree
#' eTerm_LD_Tree <- xEnricherSNPs(data, ontology="EF_disease", include.LD="EUR", LD.r2=0.8, ontology.algorithm="lea", RData.location=RData.location)
#'
#' # 3) Compare enrichment results
#' list_eTerm <- list(eTerm_noLD_noTree, eTerm_noLD_Tree, eTerm_LD_noTree, eTerm_LD_Tree)
#' names(list_eTerm) <- c('LD (-) & Tree (-)', 'LD (-) & Tree (+)', 'LD (+) & Tree (-)', 'LD (+) & Tree (+)')
#' ## side-by-side comparisons 
#' bp <- xEnrichCompare(list_eTerm, displayBy="fc")
#' #pdf(file="enrichment_compared.pdf", height=6, width=12, compress=TRUE)
#' print(bp)
#' #dev.off()
#'
#' # 4) DAGplot of comparative enrichment results in the context of ontology tree
#' xEnrichDAGplotAdv(bp, graph.node.attrs=list(fontsize=100))
#' }

xEnrichDAGplotAdv <- function(ggplot, displayBy=c("nSig","none"), path.mode=c("all_paths","shortest_paths","all_shortest_paths"), height=7, width=7, margin=rep(0.1,4), colormap=c("white-lightcyan-cyan","yr","bwr","jet","gbr","wyr","br","rainbow","wb","lightyellow-orange"), ncolors=40, zlim=NULL, colorbar=T, colorbar.fraction=0.1, newpage=T, layout.orientation=c("left_right","top_bottom","bottom_top","right_left"), node.info=c("term_name","term_id","none"), wrap.width=NULL, graph.node.attrs=NULL, graph.edge.attrs=NULL, node.attrs=NULL)
{
    
    displayBy <- match.arg(displayBy)
    path.mode <- match.arg(path.mode)
    layout.orientation <- match.arg(layout.orientation)
    node.info<- match.arg(node.info)
    
   	if(any(class(ggplot) %in% c("gg","ggplot"))){
		bp <- ggplot
		if(!is.null(bp$g)){
			bp_ig <- bp$g
			bp_data <- bp$data
			nodes_query <- unique(bp_data$id)
			subg <- dDAGinduce(g=bp_ig, nodes_query=nodes_query, path.mode=path.mode)
		}else{
			stop("No 'igraph' object can not be found in the input 'ggplot' object.\n")
		}
	}else{
		stop("Cannot find either an 'eTerm' object or an 'ggplot' object.\n")
	}
	
	## add 'code' to subg
	code <- bp_data$code
	names(code) <- bp_data$id
	ind <- match(V(subg)$name, names(code))
	code_in_subg <- code[ind]
	names(code_in_subg) <- V(subg)$name
	V(subg)$code <- code_in_subg

	## add 'nSig' to subg
	nSig <- bp_data$nSig
	names(nSig) <- bp_data$id
	ind <- match(V(subg)$name, names(nSig))
	nSig_in_subg <- nSig[ind]
	names(nSig_in_subg) <- V(subg)$name
	V(subg)$nSig <- nSig_in_subg
	
	## prefixed code to tern id and term name
	for(i in 1:length(V(subg)$code)){
		if(!is.na(V(subg)$code[i])){
			V(subg)$name[i] <- paste(V(subg)$code[i], V(subg)$name[i], sep="\\\n")
			V(subg)$term_name[i] <- paste(V(subg)$code[i], V(subg)$term_name[i], sep="\\\n")	
		}
	}
	
	####################################################################
	## how to color nodes/terms
	if(displayBy=='nSig'){
		data <- V(subg)$nSig
		names(data) <- V(subg)$name
		data <- data[!is.na(data)]
		if(is.null(zlim)){
			zlim <- c(0, max(data))
		}
	}else{
		data <- NULL
	}
	
	agDAG <- visDAG(g=subg, data=data, height=height, width=width, margin=margin, colormap=colormap, ncolors=ncolors, zlim=zlim, colorbar=colorbar, colorbar.fraction=colorbar.fraction, newpage=newpage, layout.orientation=layout.orientation, node.info=node.info, numChar=wrap.width, graph.node.attrs=graph.node.attrs, graph.edge.attrs=graph.edge.attrs, node.attrs=node.attrs)
	
	#slotNames(agDAG)
	#str(agDAG@AgNode[[1]])
	#name_nodes <- sapply(agDAG@AgNode,name)
	#shape_nodes <- sapply(agDAG@AgNode,shape)
	#txtLabelText_nodes <- sapply(lapply(agDAG@AgNode,txtLabel), labelText)
	
	invisible(agDAG)
}
