#' Function to generate a graphml file from a graph object of class "igraph"
#'
#' \code{xGraphML} is supposed to generate a graphml file from a graph object of class "igraph".
#'
#' @param g an object of class "igraph"
#' @param node.label either a vector labelling nodes or a character specifying which node attribute used for the labelling. If NULL (by default), no node labelling. If provided as a vector, a node with 'NA' will be not labelled
#' @param label.wrap.width a positive integer specifying wrap width of name
#' @param node.tooltip either a vector used for node tooltips or a character specifying which node attribute used for the tooltips. If NULL (by default), node attribute 'name' will be used node lab
#' @param node.link a string specifying hyperlink address. By default, it is NULL meaning no hyperlink
#' @param node.color a character specifying which node attribute used for node coloring. If NULL (by default), it is '#BFFFBF'
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param node.size either a vector specifying node size or a character specifying which node attribute used for the node size. If NULL (by default), it will be 30
#' @param filename the without-extension part of the name of the output file. By default, it is 'xGraphML'
#' @return
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{xGraphML}}
#' @include xGraphML.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' 
#' # load REACTOME
#' # restricted to Immune System ('R-HSA-168256') or Signal Transduction ('R-HSA-162582')
#' g <- xRDataLoader(RData.customised='ig.REACTOME', RData.location=RData.location)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="R-HSA-168256", mode="out")
#' nodeInduced <- V(g)[unique(unlist(neighs.out))]$name
#' ig <- igraph::induced.subgraph(g, vids=nodeInduced)
#'
#' # 4) visualise the graph with vertices being color-coded by the pattern
#' pattern <- runif(vcount(ig))
#' names(pattern) <- V(ig)$name
#' xGraphML(g=ig, colormap="wyr")
#' }

xGraphML <- function(g, node.label=NULL, label.wrap.width=NULL, node.tooltip=NULL, node.link=NULL, node.color=NULL, colormap='wyr', ncolors=64, zlim=NULL, node.size=NULL, filename='xGraphML')
{
    
    if (class(g) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }else{
    	ig <- g
    }
    
    ######################################################################################
    ######################################################################################
    ## node.color (by default, "#BFFFBF")
    if (!is.null(node.color)){
    	pattern <- igraph::vertex_attr(ig, node.color)
     
        if(!is.null(pattern)){
        	
        	pattern <- as.numeric(pattern)
        	pattern_nona <- pattern[!is.na(pattern)]
        	pattern_nona <- as.numeric(pattern_nona)
        	
            if(is.null(zlim)){
                vmin <- floor(stats::quantile(pattern_nona, 0.05))
                vmax <- ceiling(stats::quantile(pattern_nona, 0.95))
                if(vmin < 0 & vmax > 0){
                    vsym <- abs(min(vmin, vmax))
                    vmin <- -1*vsym
                    vmax <- vsym
                }
                zlim <- c(vmin,vmax)
            }
            
            ## A function to map a vector to colors
            vec2color <- function(vec, colormap=colormap, ncolors=ncolors, zlim=zlim){
                palette.name <- xColormap(colormap=colormap)
                colors <- palette.name(ncolors)
                scale <- length(colors)/(max(zlim)-min(zlim))
                sapply(1:length(vec), function(x){
                	if(is.na(vec[x])){
                		'transparent'
                	}else{
						ind <- floor(1+(vec[x]-min(zlim))*scale)
						colors[max(1,min(ncolors,ind))]
					}
                })
            }
            node.color <- vec2color(pattern, colormap=colormap, ncolors=ncolors, zlim=zlim)
        }else{
            warning("The input 'pattern' is ignored. Please check the help for enabling your input")
            node.color <- rep("#BFFFBF", vcount(ig))
        }
    }else{
        node.color <- rep("#BFFFBF", vcount(ig))
    }
    ######################################################################################
    ######################################################################################
    
    #############
    ## head
    #############
    output.head <- '<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:java="http://www.yworks.com/xml/yfiles-common/1.0/java" xmlns:sys="http://www.yworks.com/xml/yfiles-common/markup/primitives/2.0" xmlns:x="http://www.yworks.com/xml/yfiles-common/markup/2.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:y="http://www.yworks.com/xml/graphml" xmlns:yed="http://www.yworks.com/xml/yed/3" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd">
  <!--Created by yEd 3.17-->
  <key for="node" id="d1" attr.name="url" attr.type="string"/>
  <key for="node" id="d2" attr.name="description" attr.type="string"/>
  <key for="node" id="d3" yfiles.type="nodegraphics"/>
  <key for="edge" id="d4" attr.name="description" attr.type="string"/>
  <key for="edge" id="d5" yfiles.type="edgegraphics"/>
  <key for="graphml" id="d6" yfiles.type="resources"/>
  <graph edgedefault="directed" id="G">';
    
    #############
    ## nodes
    #############
    nnode <- igraph::vcount(ig)
    
    ## node.label (by default, NULL)
    if(length(node.label)!=nnode){
		if(!is.null(node.label)){
			node.label <- igraph::vertex_attr(ig, node.label)
		}
		if(is.null(node.label)){
			#node.label <- igraph::vertex_attr(ig, 'name')
			node.label <- rep(NA, nnode)
		}
	}
	node.label <- unlist(lapply(node.label, function(x) gsub('/','-',x)))
	node.label <- unlist(lapply(node.label, function(x) gsub('&','-',x)))
	
	## text wrap
	if(!is.null(label.wrap.width)){
		width <- as.integer(label.wrap.width)
		res_list <- lapply(node.label, function(x){
			if(!is.na(x)){
				x <- gsub('_', ' ', x)
				y <- strwrap(x, width=width)
				if(length(y)>1){
					paste0(y[1], '...')
				}else{
					y
				}
			}else{
				x
			}
		})
		node.label <- unlist(res_list)
	}
	
    ## node.tooltip (by default, the 'name' node attribute)
    if(length(node.tooltip)!=nnode){
		if(!is.null(node.tooltip)){
			node.tooltip <- igraph::vertex_attr(ig, node.tooltip)
		}
		if(is.null(node.tooltip)){
			node.tooltip <- igraph::vertex_attr(ig, 'name')
		}
    }
    node.tooltip <- unlist(lapply(node.tooltip, function(x) gsub('/','-',x)))
    node.tooltip <- unlist(lapply(node.tooltip, function(x) gsub('&','-',x)))
    
    ## node.size (by default, 30)
    if(length(node.size)!=nnode){
		if(!is.null(node.size)){
			node.size <- igraph::vertex_attr(ig, node.size)
			node.size <- 5 + 24 * (node.size - min(node.size)) / (max(node.size) - min(node.size))
		}
		if(is.null(node.size)){
			node.size <- rep(30, nnode)
		}
    }
    
    ## artificially create 'name'
    V(ig)$label <- V(ig)$name
    V(ig)$name <- paste0('n', 1:vcount(ig))
    ## layout
    glayout <- igraph::layout_with_kk(ig)
    
    ## do loop
    df_nodes <- igraph::get.data.frame(ig, what="vertices")
    df_nodes$node.label <- node.label
    df_nodes$node.tooltip <- node.tooltip
    #df_nodes$node.color <- node.color
    df_nodes$node.color <- paste0(node.color, 'cc')
    df_nodes$node.size <- node.size
    
    ### sort: NA comes last
    #df_nodes <- df_nodes[with(df_nodes,order(node.label)), ]
    ###
    
    ls_nodes <- lapply(1:nrow(df_nodes), function(i){
    	
    	k <- 0
    	vec <- vector()
    	
    	k <- k+1
    	vec[k] <- paste0('<node id="', df_nodes$name[i], '">')
    	
    	#######
    	if(!is.null(node.link)){
			k <- k+1
			vec[k] <- paste0('<data key="d1"><![CDATA[', node.link, '/', df_nodes$node.label[i], ']]></data>')
    	}
    	#######
    	
    	k <- k+1
    	vec[k] <- paste0('<data key="d2"><![CDATA[', df_nodes$node.tooltip[i], ']]></data>')
    	k <- k+1
    	vec[k] <- paste0('<data key="d3">')
    	k <- k+1
    	vec[k] <- paste0('<y:ShapeNode>')
    	k <- k+1
    	vec[k] <- paste0('<y:Geometry height="', df_nodes$node.size[i], '" width="', df_nodes$node.size[i] ,'" x="0" y="0"/>')
    	#vec[k] <- paste0('<y:Geometry height="1" width="2" x="', glayout[i,1], '" y="', glayout[i,2], '"/>')
    	k <- k+1
    	vec[k] <- paste0('<y:Fill color="', df_nodes$node.color[i], '" transparent="false"/>')
    	k <- k+1
    	vec[k] <- paste0('<y:BorderStyle color="#dddddd" raised="false" type="line" width="1"/>')
    	
    	########
    	k <- k+1
    	if(!is.na(df_nodes$node.label[i])){
    		vec[k] <- paste0('<y:NodeLabel alignment="center" autoSizePolicy="content" borderDistance="0.0" fontFamily="Arial" fontSize="12" fontStyle="bold" hasBackgroundColor="false" hasLineColor="false" height="30" horizontalTextPosition="center" iconTextGap="4" modelName="sides" modelPosition="n" textColor="#000000" verticalTextPosition="bottom" visible="true" width="30" x="0" y="0">', df_nodes$node.label[i], '</y:NodeLabel>')
    	}else{
    		vec[k] <- paste0('<y:NodeLabel alignment="center" autoSizePolicy="content" borderDistance="0.0" fontFamily="Arial" fontSize="12" fontStyle="bold" hasBackgroundColor="false" hasLineColor="false" height="30" horizontalTextPosition="center" iconTextGap="4" modelName="sides" modelPosition="n" textColor="#000000" verticalTextPosition="bottom" visible="false" width="30" x="0" y="0">', df_nodes$node.tooltip[i], '</y:NodeLabel>')
    	}
    	########
    	
    	k <- k+1
    	vec[k] <- paste0('<y:Shape type="ellipse"/>')
    	k <- k+1
    	vec[k] <- paste0('</y:ShapeNode>')
    	k <- k+1
    	vec[k] <- paste0('</data>')
    	k <- k+1
    	vec[k] <- paste0('</node>')
    	
    	paste(vec, collapse='\n')
    })
    vec_nodes <- unlist(ls_nodes)
    output.nodes <- paste(vec_nodes, collapse='\n')
    
   	############
    ## edges
    #############
    df_edges <- igraph::get.data.frame(ig, what="edges")
    ls_edges <- lapply(1:nrow(df_edges), function(i){
    
    	source <- df_edges$from[i]
    	target <- df_edges$to[i]
    	
    	k <- 0
    	vec <- vector()    	
    	
    	k <- k+1
    	vec[k] <- paste0('<edge id="', 'e', i, '" source="', source, '" target="', target, '">')
    	k <- k+1
    	vec[k] <- paste0('<data key="d5">')
    	k <- k+1
    	vec[k] <- paste0('<y:GenericEdge configuration="DEFAULT">')
    	k <- k+1
    	
    	if(0){
    	vec[k] <- paste0('<y:Path sx="0" sy="0" tx="0" ty="0"/>')
    	k <- k+1
    	}
    	
    	vec[k] <- paste0('<y:LineStyle color="#aaaaaaaa" type="line" width="1.0"/>')
    	k <- k+1
    	vec[k] <- paste0('<y:Arrows source="none" target="standard"/>')
    	k <- k+1
    	vec[k] <- paste0('</y:GenericEdge>')
    	k <- k+1
    	vec[k] <- paste0('</data>')
    	k <- k+1
    	vec[k] <- paste0('</edge>')
    	
    	paste(vec, collapse='\n')
    })
    vec_edges <- unlist(ls_edges)
    output.edges <- paste(vec_edges, collapse='\n')
    
    #############
    ## tail
    #############
    output.tail <- '</graph>
  <data key="d6">
    <y:Resources/>
  </data>
</graphml>'
    
    output <- paste0(output.head, '\n', output.nodes, '\n', output.edges, '\n', output.tail, '\n')
    if(!is.null(filename)){
		############################
		outputfile <- paste0(filename, ".graphml")
		fileConn <- base::file(outputfile)
		base::writeLines(output, fileConn)
		base::close(fileConn)
		message(sprintf("Congratulations! A file '%s' (in the directory %s) has been created!", outputfile, getwd()), appendLF=T)
		############################
    }
    
    invisible(output)
}
