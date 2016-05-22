#' Function to visualise terms used to annotate an input SNP or gene using different network layouts
#'
#' \code{xSocialiserNetplot} is supposed to visualise terms used to annotate an input SNP or gene using different network layouts. It returns an object of class 'igraph'.
#'
#' @param g an object of class "igraph" (resulting from similarity analysis)
#' @param query an object in query (for example, an SNP or Gene)
#' @param displayBy which statistics will be used for displaying. It can be "IC" for information content (by default), "none" for no color-coding on nodes/terms
#' @param path.mode the mode of paths induced by nodes in query. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param node.info tells the ontology term information used to label nodes. It can be one of "none" for no node labeling, "term_id" for using Term ID, "term_name" for using Term Name, "both" for using both of Term ID and Name (the first 15 characters), and "full_term_name" for using the full Term Name
#' @param wrap.width a positive integer specifying wrap width of Term Name. By default, first 15 characters
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param colorbar logical to indicate whether to append a colorbar. If pattern is null, it always sets to false
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param glayout either a function or a numeric matrix configuring how the vertices will be placed on the plot. If layout is a function, this function will be called with the graph as the single parameter to determine the actual coordinates. This function can be one of "layout_nicely" (previously "layout.auto"), "layout_randomly" (previously "layout.random"), "layout_in_circle" (previously "layout.circle"), "layout_on_sphere" (previously "layout.sphere"), "layout_with_fr" (previously "layout.fruchterman.reingold"), "layout_with_kk" (previously "layout.kamada.kawai"), "layout_as_tree" (previously "layout.reingold.tilford"), "layout_with_lgl" (previously "layout.lgl"), "layout_with_graphopt" (previously "layout.graphopt"), "layout_with_sugiyama" (previously "layout.kamada.kawai"), "layout_with_dh" (previously "layout.davidson.harel"), "layout_with_drl" (previously "layout.drl"), "layout_with_gem" (previously "layout.gem"), "layout_with_mds". A full explanation of these layouts can be found in \url{http://igraph.org/r/doc/layout_nicely.html}
#' @param vertex.frame.color the color of the frame of the vertices. If it is NA, then there is no frame
#' @param vertex.size the size of each vertex. If it is a vector, each vertex may differ in size
#' @param vertex.color the fill color of the vertices. If it is NA, then there is no fill color. If the pattern is given, this setup will be ignored
#' @param vertex.shape the shape of each vertex. It can be one of "circle", "square", "csquare", "rectangle", "crectangle", "vrectangle", "pie" (\url{http://igraph.org/r/doc/vertex.shape.pie.html}), "sphere", and "none". If it sets to NULL, these vertices with negative will be "csquare" and the rest "circle". 
#' @param vertex.label the label of the vertices. If it is NA, then there is no label. The default vertex labels are the name attribute of the nodes
#' @param vertex.label.cex the font size of vertex labels.
#' @param vertex.label.dist the distance of the label from the center of the vertex. If it is 0 then the label is centered on the vertex. If it is 1 then the label is displayed beside the vertex.
#' @param vertex.label.color the color of vertex labels.
#' @param edge.arrow.size the size of the arrows for the directed edge. The default value is 1.
#' @param ... additional graphic parameters. See \url{http://igraph.org/r/doc/plot.common.html} for the complete list.
#' @return 
#' an igraph object to represent DAG, appended with a node attribute called 'inherited' indicative of whether terms are inherited or not
#' @note none
#' @export
#' @seealso \code{\link{xSocialiserGenes}}, \code{\link{xSocialiserSNPs}}
#' @include xSocialiserNetplot.r
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
#' # 3) Net plot visualising terms used to annotate an SNP 'rs6871626'
#' dag <- xSocialiserNetplot(g=SNP.g, query='rs6871626', displayBy="IC", node.info=c("none"), vertex.label=NA, wrap.width=30)
#' }

xSocialiserNetplot <- function(g, query, displayBy=c("IC","none"), path.mode=c("all_paths","shortest_paths","all_shortest_paths"), node.info=c("none", "term_id", "term_name", "both", "full_term_name"), wrap.width=15, colormap=c("yr","jet","gbr","wyr","br","bwr","rainbow","wb"), ncolors=40, zlim=NULL, colorbar=T, newpage=T, glayout=layout_as_tree, vertex.frame.color=NA, vertex.size=NULL, vertex.color=NULL, vertex.shape=NULL, vertex.label=NULL, vertex.label.cex=NULL, vertex.label.dist=0.3, vertex.label.color="blue", edge.arrow.size=0.3, ...)
{
    
    displayBy <- match.arg(displayBy)
    path.mode <- match.arg(path.mode)
    node.info<- match.arg(node.info)
    
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
		stop("No terms are found to annotate the entity in query!")
		
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
	
	#####################################################################################
	## for vertex.label
    if(is.null(vertex.label)){
		## define node labels
		getTermInfo <- function(g, vids, numChar=15, mulLines=F){
			fullNames <- V(g)[vids]$term_name
			names(fullNames) <- V(g)[vids]$name
	
			if(mulLines==F){
				shortNames <- paste(substr(fullNames,1,numChar), ifelse(nchar(fullNames)>numChar, '...', ''), sep='')
			}else{
				shortNames <- sapply(fullNames,function(x){
					return(paste(strwrap(x, numChar), sep="", collapse = "\n"))
				})
			}
	
			names(shortNames) <- names(fullNames)
			return(shortNames)
		}
		termNames <- getTermInfo(subg, vids=V(subg)$term_id, numChar=15, mulLines=F)
		vertex.label <- switch(node.info,
						   none = NULL,
						   term_id = V(subg)$term_id,
						   term_name = V(subg)$term_name,
						   both = paste(V(subg)$term_id, termNames, sep="\n"),
						   full_term_name = getTermInfo(subg, vids=V(subg)$term_id, numChar=wrap.width, mulLines=T)
		)
    }
	
	######################################################################################
	## for data (pattern)
	pattern <- data
	
    if (!is.null(pattern)){
    
        flag <- 0
        if(!is.null(names(pattern))){
            pattern <- pattern[V(subg)$name]
        }
        if(length(pattern)==vcount(subg)){
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
                palette.name <- supraHex::visColormap(colormap=colormap)
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
            vertex.frame.color <- vec2color(pattern, colormap=colormap, ncolors=ncolors, zlim=zlim)
            vertex.frame.color[vertex.frame.color=="transparent"] <- "grey"
        }else{
            warning("The input 'pattern' is ignored. Please check the help for enabling your input")
            pattern <- NULL
            if(is.null(vertex.color)){
                vertex.color <- "SkyBlue2"
            }
        }
    }else{
        if(is.null(vertex.color)){
            vertex.color <- "SkyBlue2"
        }
    }
    
    ######################################################################################
	
	par_old <- graphics::par()
	
	dnet::visNet(g=subg, pattern=pattern, colormap=colormap, ncolors=ncolors, zlim=zlim, colorbar=colorbar, newpage=newpage, glayout=glayout, vertex.frame.color=vertex.frame.color, vertex.size=vertex.size, vertex.color=vertex.color, vertex.shape=vertex.shape, vertex.label=vertex.label, vertex.label.cex=vertex.label.cex, vertex.label.dist=vertex.label.dist, vertex.label.color=vertex.label.color, edge.arrow.size=edge.arrow.size, ...)
	
	suppressWarnings(graphics::par(par_old))
	
	invisible(subg)
}
