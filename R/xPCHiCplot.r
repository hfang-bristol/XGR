#' Function to visualise promoter capture HiC data using different network layouts
#'
#' \code{xPCHiCplot} is supposed to visualise promoter capture HiC data using different network layouts.
#'
#' @param g an object of both classes "igraph" and "PCHiC" (part of the results from \code{\link{xDefineHIC}})
#' @param node.info tells the information used to label nodes. It can be one of "none" for no node labeling, "GR" for only using genomic regions (GR), "GR_SNP" for using GR and SNP (if any), "GR_SNP_target" for using GR and SNP (if any) and target genes (if any), "SNP_target" for using SNP (if any) and target genes (if any), and "smart" (by default) for only using GR if both SNP and target genes are not available (otherwise GR will be hidden)
#' @param node.colors colors used to flag which nodes contain SNP or not. By default, a node harboring an SNP will be colored in 'skyblue' and the node without an SNP in 'pink'
#' @param nodes.query nodes in query for which edges attached to them will be displayed. By default, it sets to NULL meaning no such restriction
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE
#' @param glayout either a function or a numeric matrix configuring how the vertices will be placed on the plot. If layout is a function, this function will be called with the graph as the single parameter to determine the actual coordinates. This function can be one of "layout_nicely" (previously "layout.auto"), "layout_randomly" (previously "layout.random"), "layout_in_circle" (previously "layout.circle"), "layout_on_sphere" (previously "layout.sphere"), "layout_with_fr" (previously "layout.fruchterman.reingold"), "layout_with_kk" (previously "layout.kamada.kawai"), "layout_as_tree" (previously "layout.reingold.tilford"), "layout_with_lgl" (previously "layout.lgl"), "layout_with_graphopt" (previously "layout.graphopt"), "layout_with_sugiyama" (previously "layout.kamada.kawai"), "layout_with_dh" (previously "layout.davidson.harel"), "layout_with_drl" (previously "layout.drl"), "layout_with_gem" (previously "layout.gem"), "layout_with_mds". A full explanation of these layouts can be found in \url{http://igraph.org/r/doc/layout_nicely.html}
#' @param vertex.frame.color the color of the frame of the vertices. If it is NA, then there is no frame
#' @param vertex.size the size of each vertex. If it is a vector, each vertex may differ in size
#' @param vertex.color the fill color of the vertices. If it is NA, then there is no fill color
#' @param vertex.shape the shape of each vertex. It can be one of "circle", "square", "csquare", "rectangle", "crectangle", "vrectangle", "pie" (\url{http://igraph.org/r/doc/vertex.shape.pie.html}), "sphere", and "none"
#' @param vertex.label the label of the vertices. If it is NA, then there is no label. The default vertex labels are the name attribute of the nodes
#' @param vertex.label.cex the font size of vertex labels.
#' @param vertex.label.font the font of vertex labels. It is interpreted the same way as the the 'font' graphical parameter: 1 is plain text, 2 is bold face, 3 is italic, 4 is bold and italic and 5 specifies the symbol font.
#' @param vertex.label.dist the distance of the label from the center of the vertex. If it is 0 then the label is centered on the vertex. If it is 1 then the label is displayed beside the vertex.
#' @param vertex.label.color the color of vertex labels.
#' @param edge.arrow.size the size of the arrows for the directed edge. The default value is 0.5.
#' @param edge.width the width of the directed edge. If NULL, the width edge is proportional to CHiCAGO scores (quantifying the strength of physical interactions).
#' @param edge.color the color of the directed edge. The default value is 'grey'.
#' @param ... additional graphic parameters. See \url{http://igraph.org/r/doc/plot.common.html} for the complete list.
#' @return an igraph object
#' @note
#' \itemize{
#'  \item{\code{edge arrow}: interactions are represented as a direct graph (bait-prey)}
#'  \item{\code{edge thickness}: the thickness in an edge is proportional to the interaction strength}
#'  \item{\code{node color}: a node is colored in pink if it harbors SNPs in query; otherwise skyblue}
#'  \item{\code{node label}: a node is labelled with three pieces of information (if any): genomic regions, SNPs in query, genes associated (marked by an @ icon)}
#' }
#' @export
#' @seealso \code{\link{xDefineHIC}}
#' @include xPCHiCplot.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' data(ImmunoBase)
#' data <- names(ImmunoBase$AS$variants)
#'
#' # b) extract HiC-gene pairs given a list of AS SNPs
#' PCHiC <- xDefineHIC(data, include.HiC="Monocytes", GR.SNP="dbSNP_GWAS", RData.location=RData.location)
#' head(PCHiC$df)
#' 
#' # c) visualise the interaction (a directed graph: bait->prey)
#' g <- PCHiC$ig
#' ## a node with SNPs colored in 'skyblue' and the one without SNPs in 'pink'
#' ## the width in an edge is proportional to the interaction strength
#' xPCHiCplot(g, vertex.shape="sphere")
#' xPCHiCplot(g, glayout=layout_in_circle, vertex.shape="sphere")
#' 
#' # d) control node labelling info
#' xPCHiCplot(g, node.info="GR_SNP_target")
#' xPCHiCplot(g, node.info="GR_SNP")
#' xPCHiCplot(g, node.info="SNP_target")
#' xPCHiCplot(g, node.info='SNP_target', vertex.label.cex=0.5)
#' }

xPCHiCplot <- function(g, node.info=c("smart", "none", "GR", "GR_SNP", "GR_SNP_target", "SNP_target"), node.colors=c("skyblue","pink1"), nodes.query=NULL, newpage=TRUE, signature=TRUE, glayout=layout_with_kk, vertex.frame.color=NA, vertex.size=NULL, vertex.color=NULL, vertex.shape="sphere", vertex.label=NULL, vertex.label.cex=NULL, vertex.label.font=2, vertex.label.dist=0.3, vertex.label.color="black", edge.arrow.size=0.5, edge.width=NULL, edge.color="grey", ...)
{
    
    node.info<- match.arg(node.info)
    
    if(all(class(g) %in% c("igraph","PCHiC"))){
		subg <- g
		class(subg) <- "igraph"

	}else{
		stop("Cannot find an object of both classes 'igraph' and 'PCHiC'.\n")
	}
	
	## restrict to nodes in query
	if(!is.null(nodes.query)){
		subg <- dnet::dNetInduce(g=subg, nodes_query=nodes.query, knn=1, remove.loops=FALSE, largest.comp=FALSE)
	}

	## for vertex.label
    if(is.null(vertex.label)){
		## define node labels
		vertex.label <- switch(node.info,
						   none = NULL,
						   GR = V(subg)$name,
						   GR_SNP = paste(V(subg)$name, '\n',V(subg)$SNP, sep=""),
						   GR_SNP_target = paste(V(subg)$name, '\n',V(subg)$SNP, '\n@', V(subg)$target, sep=""),
						   SNP_target = paste(V(subg)$SNP, '\n@', V(subg)$target, sep=""),
						   smart = paste(V(subg)$name, '\n',V(subg)$SNP, '\n@', V(subg)$target, sep=""),
		)
		
		for(i in 1:length(vertex.label)){
			vertex.label[i]
		}
		
		if(!is.null(vertex.label)){
			vertex.label <- gsub('\nNA','',vertex.label)
			vertex.label <- gsub('\n@\\.','',vertex.label)
			vertex.label <- gsub('NA\n','',vertex.label)
			
			vertex.label[vertex.label=='NA'] <- ''
			
			if(node.info=='smart'){
				res_ls <- lapply(strsplit(vertex.label,'\n'), function(x){
					if(length(x)>=2){
						x <- x[-1]
					}
					paste(x, collapse='\n')
				})
				vertex.label <- unlist(res_ls)
			}
			
		}
    }
	
	## vertex.color
    if(is.null(vertex.color)){
		vertex.color <- rep(node.colors[1], vcount(subg))
		vertex.color[!is.na(V(subg)$SNP)] <- node.colors[2]
    }
	
	
	## edge.width
	if(is.null(edge.width)){
		## extract edge weight (with 2-digit precision)
		x <- signif(as.numeric(E(subg)$score), digits=2)
		## rescale into an interval [1,4] as edge width
		edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
	}
	
	
	suppressWarnings(par_old <- graphics::par(no.readonly=TRUE))
	
	suppressWarnings(dnet::visNet(g=subg, newpage=newpage, glayout=glayout, vertex.frame.color=vertex.frame.color, vertex.size=vertex.size, vertex.color=vertex.color, vertex.shape=vertex.shape, vertex.label=vertex.label, vertex.label.cex=vertex.label.cex, vertex.label.font=vertex.label.font, vertex.label.dist=vertex.label.dist, vertex.label.color=vertex.label.color, edge.arrow.size=edge.arrow.size, edge.width=edge.width, edge.color=edge.color, ...))
	
	suppressWarnings(graphics::par(par_old))
	
    if(signature){
    	caption <- paste("Created by xPCHiCplot from XGR version", utils::packageVersion("XGR"))
    	graphics::mtext(caption, side=1, line=2, adj=1, cex=.66, font=3)
    }
	
	invisible(subg)
}
