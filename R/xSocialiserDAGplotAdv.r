#' Function to draw DAG plot for comparing two sets of terms used to annotate two SNPs or genes in query
#'
#' \code{xSocialiserDAGplotAdv} is supposed to use DAG plot for comparing two sets of terms used to annotate two SNPs or genes in query. Per term, comparative results are coded in the form of 'x1-x2', where x1 is for query 1 and x2 for query 2 (the value for x1 or x2 can be '0' encoding for no anntation, '1' for inherited annotation, '2' for direct annotation). It returns an object of class 'Ragraph' or class 'igraph'.
#'
#' @param g an object of class "igraph" (resulting from similarity analysis)
#' @param query1 the first object in query (for example, an SNP or Gene)
#' @param query2 the second object in query (for example, an SNP or Gene)
#' @param displayBy which statistics will be used for displaying. It can be "IC" for information content (by default), "none" for no color-coding on nodes/terms
#' @param path.mode the mode of paths induced by nodes/terms. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
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
#' @param node.info tells the ontology term information used to label nodes. It can be "term_id" for using Term ID, "term_name" for using Term Name
#' @param wrap.width a positive integer specifying wrap width of Term Name
#' @param graph.node.attrs a list of global node attributes. These node attributes will be changed globally. See 'Note' below for details on the attributes
#' @param graph.edge.attrs a list of global edge attributes. These edge attributes will be changed globally. See 'Note' below for details on the attributes
#' @param node.attrs a list of local edge attributes. These node attributes will be changed locally; as such, for each attribute, the input value must be a named vector (i.e. using Term ID as names). See 'Note' below for details on the attributes
#' @param output.format the format specifying the return value. It can be "Ragraph" (by default) or "igraph"
#' @return 
#' An object of class 'Ragraph' or 'igraph'.  If the returned is an Ragraph object, an image will be shown. If the returned is an igraph object, no image will be shown; in this case, the returned igraph object stores ontology terms used to annotate the query, including a new node attribute 'code' indicative of how terms are shared or unique to two queries (in the form of 'x1-x2', x1 for query 1 and x2 for query 2, x1 or x2 can be '0' for no anntation, '1' for inherited annotation, '2' for direct annotation).
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
#' @seealso \code{\link{xSocialiserGenes}}, \code{\link{xSocialiserSNPs}}, \code{\link{xSocialiserDAGplot}}
#' @include xSocialiserDAGplotAdv.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location="~/Sites/SVN/github/bigdata"
#' 
#' # 1) SNP-based similarity analysis using GWAS Catalog traits (mapped to EF)
#' # provide genes and SNPs reported in AS GWAS studies
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' ## get lead SNPs reported in AS GWAS
#' example.snps <- names(ImmunoBase$AS$variants)
#' SNP.g <- xSocialiserSNPs(example.snps, include.LD=NA, RData.location=RData.location)
#' 
#' # 2) Circos plot involving nodes 'rs6871626'
#' xCircos(g=SNP.g, entity="SNP", nodes.query="rs6871626", RData.location=RData.location)
#'
#' # 3) DAG plot visualising terms used to annotate an SNP
#' ## 3a) for 'rs6871626'
#' xSocialiserDAGplot(g=SNP.g, query='rs6871626', displayBy="IC", node.info=c("term_name"), graph.node.attrs=list(fontsize=20,fontcolor="blue",color="transparent"))
#' ## 3b) for 'rs1250550'
#' xSocialiserDAGplot(g=SNP.g, query='rs1250550', displayBy="IC", node.info=c("term_name"), graph.node.attrs=list(fontsize=20,fontcolor="blue",color="transparent"))
#'
#' # 4) DAG plot comparing two sets of terms used to annotate two queried SNPs
#' xSocialiserDAGplotAdv(g=SNP.g, query1='rs6871626', query2='rs1250550', node.info=c("term_name"), graph.node.attrs=list(fontsize=25,fontcolor="blue",color="transparent"))
#'
#' # 5) Return an igraph object storing ontology terms used to annotate an SNP 'rs6871626'
#' dag <- xSocialiserDAGplotAdv(g=SNP.g, query1='rs6871626', query2='rs1250550', output.format="igraph")
#' }

xSocialiserDAGplotAdv <- function(g, query1, query2, displayBy=c("IC","none"), path.mode=c("all_paths","shortest_paths","all_shortest_paths"), height=7, width=7, margin=rep(0.1,4), colormap=c("wyr","bwr","jet","gbr","yr","br","rainbow","wb","lightyellow-orange"), ncolors=40, zlim=NULL, colorbar=T, colorbar.fraction=0.1, newpage=T, layout.orientation=c("top_bottom","left_right","bottom_top","right_left"), node.info=c("term_name","term_id"), wrap.width=NULL, graph.node.attrs=NULL, graph.edge.attrs=NULL, node.attrs=NULL, output.format=c("Ragraph","igraph")) 
{

    displayBy <- match.arg(displayBy)
    path.mode <- match.arg(path.mode)
    layout.orientation <- match.arg(layout.orientation)
    node.info <- match.arg(node.info)
    output.format <- match.arg(output.format)
    
	ig_SNP1 <- xSocialiserDAGplot(g, query=query1, output.format="igraph")
	ig_SNP2 <- xSocialiserDAGplot(g, query=query2, output.format="igraph")
	
	if(is.null(ig_SNP1) | is.null(ig_SNP2)){
		warnings("No terms are found to annotate one or two entities in query!")
		return(NULL)
	}
	
	list_igraph <- list(ig_SNP1, ig_SNP2)
	names(list_igraph) <- c(query1, query2)

	## Combine into an igraph object called "subg"
	### edges
	ls_edges <- lapply(list_igraph, function(x){
		df_edge <- igraph::get.data.frame(x, what="edges")
	})
	relations <- do.call(rbind, ls_edges)
	relations <- relations[!duplicated(relations), ]
	### nodes
	ls_nodes <- lapply(list_igraph, function(x){
		df_nodes <- igraph::get.data.frame(x, what="vertices")
	})
	nodes <- do.call(rbind, ls_nodes)
	nodes <- nodes[,-7] ## remove 'inherited'
	nodes <- nodes[!duplicated(nodes), ]
	### igraph
	ig <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
	
	### code for inherited
	ls_res <- lapply(list_igraph, function(x){
		ind <- match(V(ig)$name, V(x)$name)
		res <- V(x)$inherited[ind] * (-1) + 2
		res[is.na(res)] <- 0
		res
	})
	df_res <- do.call(cbind, ls_res)
	code <- apply(df_res, 1, paste,collapse='-')
	V(ig)$code <- code
	
	######################################################################################
	
	subg <- ig
	
	V(subg)$name <- paste(V(subg)$code, V(subg)$name, sep="\\\n")
	V(subg)$term_name <- paste(V(subg)$code, V(subg)$term_name, sep="\\\n")
	
	## how to color nodes/terms
	if(displayBy=='IC'){
		data <- V(subg)$IC
		names(data) <- V(subg)$name
	}else{
		data <- NULL
	}
	
	if(output.format=="igraph"){
		invisible(ig)
	}else if(output.format=="Ragraph"){
	
		## for globally graph.edge.attrs
		if(is.null(graph.edge.attrs)){
			#graph.edge.attrs <- list(color="black",weight=1,style="solid")
		}
	
		agDAG <- visDAG(g=subg, data=data, height=height, width=width, margin=margin, colormap=colormap, ncolors=ncolors, zlim=zlim, colorbar=colorbar, colorbar.fraction=colorbar.fraction, newpage=newpage, layout.orientation=layout.orientation, node.info=node.info, numChar=wrap.width, graph.node.attrs=graph.node.attrs, graph.edge.attrs=graph.edge.attrs, node.attrs=node.attrs)
	
		#slotNames(agDAG)
		#str(agDAG@AgNode[[1]])
		#name_nodes <- sapply(agDAG@AgNode,name)
		#shape_nodes <- sapply(agDAG@AgNode,shape)
		#txtLabelText_nodes <- sapply(lapply(agDAG@AgNode,txtLabel), labelText)
		
		invisible(agDAG)
	}
}
