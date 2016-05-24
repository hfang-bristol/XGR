#' Function to draw DAG plot for visualising terms used to annotate an input SNP or gene
#'
#' \code{xSocialiserDAGplot} is supposed to draw DAG plot for visualising terms used to annotate an input SNP or gene. By default, original terms are highlighted by box-shaped nodes, inherited terms by ellipse nodes. It returns an object of class 'Ragraph' or class 'igraph'.
#'
#' @param g an object of class "igraph" (resulting from similarity analysis)
#' @param query an object in query (for example, an SNP or Gene)
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
#' @param node.info tells the ontology term information used to label nodes. It can be one of "none" for no node labeling, "term_id" for using Term ID, "term_name" for using Term Name (the first 15 characters), "both" for using both of Term ID and Name (the first 15 characters), and "full_term_name" for using the full Term Name
#' @param wrap.width a positive integer specifying wrap width of Term Name
#' @param graph.node.attrs a list of global node attributes. These node attributes will be changed globally. See 'Note' below for details on the attributes
#' @param graph.edge.attrs a list of global edge attributes. These edge attributes will be changed globally. See 'Note' below for details on the attributes
#' @param node.attrs a list of local edge attributes. These node attributes will be changed locally; as such, for each attribute, the input value must be a named vector (i.e. using Term ID as names). See 'Note' below for details on the attributes
#' @param output.format the format specifying the return value. It can be "Ragraph" (by default) or "igraph"
#' @return 
#' An object of class 'Ragraph' or 'igraph'.  If the returned is an Ragraph object, an image will be shown. If the returned is an igraph object, no image will be shown; in this case, the returned igraph object stores ontology terms used to annotate the query, including a new node attribute 'inherited' indicative of whether terms are inherited or not.
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
#' @seealso \code{\link{xSocialiserGenes}}, \code{\link{xSocialiserSNPs}}
#' @include xSocialiserDAGplot.r
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
#' # 3) DAG plot visualising terms used to annotate an SNP 'rs6871626'
#' agDAG <- xSocialiserDAGplot(g=SNP.g, query='rs6871626', displayBy="IC", node.info=c("full_term_name"))
#' ## modify node labels
#' xSocialiserDAGplot(g=SNP.g, query='rs6871626', displayBy="IC", node.info=c("full_term_name"), graph.node.attrs=list(fontsize=20,fontcolor="blue",color="transparent"))
#'
#' # 4) Return an igraph object storing ontology terms used to annotate an SNP 'rs6871626'
#' dag <- xSocialiserDAGplot(g=SNP.g, query='rs6871626', displayBy="IC", output.format="igraph")
#' }

xSocialiserDAGplot <- function(g, query, displayBy=c("IC","none"), path.mode=c("all_paths","shortest_paths","all_shortest_paths"), height=7, width=7, margin=rep(0.1,4), colormap=c("yr","bwr","jet","gbr","wyr","br","rainbow","wb","lightyellow-orange"), ncolors=40, zlim=NULL, colorbar=T, colorbar.fraction=0.1, newpage=T, layout.orientation=c("top_bottom","left_right","bottom_top","right_left"), node.info=c("none", "term_id", "term_name", "both", "full_term_name"), wrap.width=NULL, graph.node.attrs=NULL, graph.edge.attrs=NULL, node.attrs=NULL, output.format=c("Ragraph","igraph"))
{
    
    displayBy <- match.arg(displayBy)
    path.mode <- match.arg(path.mode)
    layout.orientation <- match.arg(layout.orientation)
    node.info <- match.arg(node.info)
    output.format <- match.arg(output.format)
    
    if(is.logical(g)){
        stop("There is no similarity in the 'igraph' object.\n")
    }
    
    if (class(g) != "igraph"){
        stop("The function must apply to the 'igraph' object.\n")
    }
    
    if(is.null(g$dag)){
    	dag <- g
    }else{
    	dag <- g$dag
    }
    
    if(is.null(V(dag)$anno) | is.null(V(dag)$IC)){
        stop("The function requires that input graph has already contained annotation data and also information content (IC).\n")
    }
    
	## get terms used to annotate an object in query
	flag <- sapply(V(dag)$anno, function(x){
		ind <- match(query, x)
		if(is.na(ind)){
			0
		}else{
			if(is.null(names(x[ind]))){
				1
			}else if(names(x[ind])=="o"){
				1
			}else if(names(x[ind])=="i"){
				2
			}
		}
	})
	terms <- V(dag)$name[flag>0]
	terms_origin <- V(dag)$name[flag==1]
	
	## return NULL if no terms are found
	if(length(terms)==0){
		warnings("No terms are found to annotate the entity in query!")
		return(NULL)
	}
	
	## get a subgraph induced by terms
	subg <- dnet::dDAGinduce(g=dag, nodes_query=terms, path.mode=path.mode)
	
	## append a node attribute 'inherited' to subg
	inherited <- rep(1, length(V(subg)$name))
	names(inherited) <- V(subg)$name
	if(length(terms_origin)>0){
		inherited[terms_origin] <- 0
	}
	V(subg)$inherited <- inherited	
	
	## how to color nodes/terms
	if(displayBy=='IC'){
		data <- V(subg)$IC
		names(data) <- V(subg)$name
	}else{
		data <- NULL
	}
	
	if(output.format=="igraph"){
		invisible(subg)
	}else if(output.format=="Ragraph"){
	
		## for globally graph.edge.attrs
		if(is.null(graph.edge.attrs)){
			#graph.edge.attrs <- list(color="black",weight=1,style="solid")
		}
	
		## for locally node.attrs
		if(is.null(node.attrs)){
			## by default, box-shaped nodes for terms originally used for annotation
			if(length(terms_origin)>0){
				hightlighted.shape <- rep("box", length(terms_origin))
				names(hightlighted.shape) <- V(subg)[terms_origin]$name
				node.attrs <- list(shape=hightlighted.shape)
			}
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
