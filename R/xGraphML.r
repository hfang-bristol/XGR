#' Function to visualise a graph object of class "igraph"
#'
#' \code{xGraphML} is supposed to visualise a graph object of class "igraph". It also allows vertices/nodes color-coded according to the input pattern. 
#'
#' @param g an object of class "igraph"
#' @param pattern a numeric vector used to color-code vertices/nodes. Notably, if the input vector contains names, then these names should include all node names of input graph, i.e. V(g)$name, since there is a mapping operation. After mapping, the length of the patern vector should be the same as the number of nodes of input graph; otherwise, this input pattern will be ignored. The way of how to color-code is to map values in the pattern onto the whole colormap (see the next arguments: colormap, ncolors, zlim and colorbar)
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
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
#' xGraphML(g=ig, pattern=pattern, colormap="bwr")
#' }

xGraphML <- function(g, pattern=NULL, colormap='wyr', ncolors=64, zlim=NULL, filename='xGraphML')
{
    
    if (class(g) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }else{
    	ig <- g
    }
    
    ######################################################################################
    ######################################################################################
    nsize <- vcount(ig)
    if (!is.null(pattern)){
    
        flag <- 0
        if(!is.null(names(pattern))){
            pattern <- pattern[V(ig)$name]
        }
        if(length(pattern)==nsize){
            flag <- 1
        }
                
        if(flag==1){
        	
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
                palette.name <- visColormap(colormap=colormap)
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
            vertex.color <- vec2color(pattern, colormap=colormap, ncolors=ncolors, zlim=zlim)
        }else{
            warning("The input 'pattern' is ignored. Please check the help for enabling your input")
            vertex.color <- rep("#BFFFBF", nsize)
        }
    }else{
        vertex.color <- rep("#BFFFBF", nsize)
    }
    ######################################################################################
    ######################################################################################
    
	V(ig)$label <- V(ig)$name
    V(ig)$name <- paste0('n', 1:vcount(ig))
    
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
    df_nodes <- igraph::get.data.frame(ig, what="vertices")
    ls_nodes <- lapply(1:nrow(df_nodes), function(i){
    	
    	k <- 0
    	vec <- vector()
    	
    	k <- k+1
    	vec[k] <- paste0('<node id="', df_nodes$name[i], '">')
    	k <- k+1
    	vec[k] <- paste0('<data key="d1"><![CDATA[http://galahad.well.ox.ac.uk:3030/pidb/target/', df_nodes$label[i], ']]></data>')
    	k <- k+1
    	vec[k] <- paste0('<data key="d2"><![CDATA[', df_nodes$term_name[i], ']]></data>')
    	k <- k+1
    	vec[k] <- paste0('<data key="d3">')
    	k <- k+1
    	vec[k] <- paste0('<y:ShapeNode>')
    	k <- k+1
    	vec[k] <- paste0('<y:Geometry height="150" width="300" x="0" y="0"/>')
    	k <- k+1
    	vec[k] <- paste0('<y:Fill color="', vertex.color[i], '" transparent="false"/>')
    	k <- k+1
    	vec[k] <- paste0('<y:BorderStyle color="#999999" raised="false" type="line" width="2.0"/>')
    	k <- k+1
    	vec[k] <- paste0('<y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="24" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="10" horizontalTextPosition="center" iconTextGap="4" modelName="custom" textColor="#000000" verticalTextPosition="bottom" visible="true" width="10" x="0" y="0">', df_nodes$label[i], '</y:NodeLabel>')
    	k <- k+1
    	vec[k] <- paste0('<y:Shape type="roundrectangle"/>')
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
    
   #############
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
    	vec[k] <- paste0('<y:Path sx="23.0" sy="0.0" tx="-23.589796832899538" ty="45.93944090514677"/>')
    	k <- k+1
    	}
    	
    	vec[k] <- paste0('<y:LineStyle color="#000000" type="line" width="1.0"/>')
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
