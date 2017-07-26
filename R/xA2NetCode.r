#' Function to create codes annotating nodes in an igraph object
#'
#' \code{xA2NetCode} is supposed to create codes annotating nodes in an igraph object. It returns two ggplot2 objects, one for visualing the network with nodes lablelled by codes, the other for listing code meaning in a table
#'
#' @param g an object of class "igraph"
#' @param node.level a character specifying which node attribute defining the node level. By default, it is 'term_distance'
#' @param node.level.value a positive integer specifying the level value as major branches. By default, it is 2
#' @param node.label.size a vector specifying node size or a character specifying which node attribute used for node label size
#' @param node.label.color the node label color
#' @param node.label.alpha the 0-1 value specifying transparency of node labelling
#' @param node.label.padding the padding around the labeled node
#' @param node.label.arrow the arrow pointing to the labeled node
#' @param node.label.force the repelling force between overlapping labels
#' @param node.shape an integer specifying node shape
#' @param node.xcoord a vector specifying x coordinates. If NULL, it will be created using igraph::layout_with_kk
#' @param node.ycoord a vector specifying y coordinates. If NULL, it will be created using igraph::layout_with_kk
#' @param node.color a character specifying which node attribute used for node coloring
#' @param node.color.title a character specifying the title for node coloring
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum values for which colors should be plotted
#' @param node.size.range the range of actual node size
#' @param title a character specifying the title for the plot
#' @param edge.color a character specifying the edge colors
#' @param edge.color.alpha the 0-1 value specifying transparency of edge colors
#' @param edge.curve a numeric value specifying the edge curve. 0 for the straight line
#' @param edge.arrow.gap a gap between the arrow and the node
#' @param node.table a character specifying which node attribute for coding. By default, it is 'term_name'
#' @param node.table.wrap a positive integer specifying wrap width of coded node labelling
#' @param table.base.size a positive integer specifying font size in the table
#' @param table.row.space a positive numeric value specifying amplying horizental space for a row with wrapped text
#' @param table.nrow a positive integer specifying the number of rows in the table
#' @param table.ncol NULL or a positive integer specifying the number of columns per page. If NULL, it will be 3 or less
#' @param simplify logical to simplify the graph or not. If TRUE, only keep edges from nodes to nodes with higher level
#' @return
#' al list of ggplot objects, with two components: code and table
#' @note none
#' @export
#' @seealso \code{\link{xA2Net}}
#' @include xA2NetCode.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' 
#' # load REACTOME
#' # 1a) restricted to Immune System ('R-HSA-168256') or Signal Transduction ('R-HSA-162582')
#' g <- xRDataLoader(RData.customised='ig.REACTOME', RData.location=RData.location)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="R-HSA-168256", mode="out")
#' vids <- V(g)[unique(unlist(neighs.out))]$name
#' ig <- igraph::induced.subgraph(g, vids=vids)
#'
#' # 1b) visualise the graph with nodes coded
#' ls_gp <- xA2NetCode(g=ig, node.level='term_distance', node.level.value=2, node.shape=19, node.size.range=4, edge.color.alpha=0.2)
#' pdf('xA2NetCode.pdf', useDingbats=FALSE, width=7, height=7)
#' print(ls_gp$code + coord_equal(ratio=1))
#' print(ls_gp$table)
#' dev.off()
#' 
#' # 1c) visualise the graph with nodes coded and colored by information content (IC)
#' V(ig)$IC <- -1*log10(V(ig)$nAnno/max(V(ig)$nAnno))
#' ls_gp <- xA2NetCode(g=ig, node.level='term_distance', node.level.value=2, node.shape=19, node.size.range=4, node.color='IC', node.color.title='IC', colormap='white-cyan-darkcyan')
#' 
#' V(ig)$term_anno <- log10(V(ig)$nAnno)
#' ls_gp <- xA2NetCode(g=ig, node.level='term_distance', node.level.value=2, node.shape=19, node.size.range=4, node.color='term_anno', node.color.title='# genes\n(log10)', colormap='white-cyan-darkcyan', zlim=c(1,4))
#' 
#' 
#' # load EF
#' # 2a) restricted to disease ('EFO:0000408') and annotation (>=10)
#' g <- xRDataLoader(RData.customised='ig.EF', RData.location=RData.location)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="EFO:0000408", mode="out")
#' nodeClan <- V(g)[unique(unlist(neighs.out))]$name
#' anno <- xRDataLoader(RData.customised='org.Hs.egEF', RData.location=RData.location)
#' vec <- sapply(anno$gs, length)
#' nodeAnno <- names(vec[vec>=10])
#' neighs.in <- igraph::neighborhood(g, order=vcount(g), nodes=nodeAnno, mode="in")
#' nodeAnno <- V(g)[unique(unlist(neighs.in))]$name
#' vids <- intersect(nodeClan, nodeAnno)
#' ig <- igraph::induced.subgraph(g, vids=vids)
#'
#' # 2b) visualise the graph with nodes coded
#' ls_gp <- xA2NetCode(g=ig, node.level='term_distance', node.level.value=4, node.shape=19, node.size.range=4, edge.color.alpha=0.2)
#' pdf('xA2NetCode.pdf', useDingbats=FALSE, width=7, height=7)
#' print(ls_gp$code + coord_equal(ratio=1))
#' print(ls_gp$table)
#' dev.off()
#' 
#' 
#' # load GOBP
#' # 3a) restricted to immune system process ('GO:0002376') and annotation (>=10)
#' g <- xRDataLoader(RData.customised='ig.GOBP', RData.location=RData.location)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="GO:0002376", mode="out")
#' nodeClan <- V(g)[unique(unlist(neighs.out))]$name
#' anno <- xRDataLoader(RData.customised='org.Hs.egGOBP', RData.location=RData.location)
#' vec <- sapply(anno$gs, length)
#' nodeAnno <- names(vec[vec>=10])
#' neighs.in <- igraph::neighborhood(g, order=vcount(g), nodes=nodeAnno, mode="in")
#' nodeAnno <- V(g)[unique(unlist(neighs.in))]$name
#' vids <- intersect(nodeClan, nodeAnno)
#' ig <- igraph::induced.subgraph(g, vids=vids)
#'
#' # 3b) visualise the graph with nodes coded
#' ls_gp <- xA2NetCode(g=ig, node.level='term_distance', node.level.value=1, node.shape=19, node.size.range=4, edge.color.alpha=0.2)
#' pdf('xA2NetCode.pdf', useDingbats=FALSE, width=7, height=7)
#' print(ls_gp$code + coord_equal(ratio=1))
#' print(ls_gp$table)
#' dev.off()
#' 
#' 
#' # load GOMF
#' # 4a) restricted to molecular function ('GO:0003674') and annotation (>=50)
#' g <- xRDataLoader(RData.customised='ig.GOMF', RData.location=RData.location)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="GO:0003674", mode="out")
#' nodeClan <- V(g)[unique(unlist(neighs.out))]$name
#' anno <- xRDataLoader(RData.customised='org.Hs.egGOMF', RData.location=RData.location)
#' vec <- sapply(anno$gs, length)
#' nodeAnno <- names(vec[vec>=50])
#' neighs.in <- igraph::neighborhood(g, order=vcount(g), nodes=nodeAnno, mode="in")
#' nodeAnno <- V(g)[unique(unlist(neighs.in))]$name
#' vids <- intersect(nodeClan, nodeAnno)
#' ig <- igraph::induced.subgraph(g, vids=vids)
#'
#' # 4b) visualise the graph with nodes coded
#' ls_gp <- xA2NetCode(g=ig, node.level='term_distance', node.level.value=1, node.shape=19, node.size.range=4, edge.color.alpha=0.2, simplify=T)
#' pdf('xA2NetCode.pdf', useDingbats=FALSE, width=7, height=7)
#' print(ls_gp$code + coord_equal(ratio=1))
#' print(ls_gp$table)
#' dev.off()
#' }

xA2NetCode <- function(g, node.level='term_distance', node.level.value=2, node.label.size=2, node.label.color='darkblue', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0.01, node.label.force=0, node.shape=19, node.xcoord=NULL, node.ycoord=NULL, node.color=NULL, node.color.title=NULL, colormap='grey-grey', ncolors=64, zlim=NULL, node.size.range=4, title='', edge.color="black", edge.color.alpha=0.4, edge.curve=0.1, edge.arrow.gap=0.02, node.table='term_name', node.table.wrap=50, table.base.size=7, table.row.space=2, table.nrow=55, table.ncol=NULL, simplify=FALSE)
{
    
   	if(any(class(g) %in% c("igraph"))){
		ig <- g
	}else{
		stop("The function must apply to a 'igraph' object.\n")
	}

	if(simplify){
		dist_from <- dist_to <- NULL
	
		nodes <- igraph::get.data.frame(ig, what="vertices")
		relations <- igraph::get.data.frame(ig, what="edges")[, c(1,2)]
		colnames(relations) <- c("from","to")
		ind <- match(relations$from, nodes$name)
		relations$dist_from <- nodes[ind, node.level]
		ind <- match(relations$to, nodes$name)
		relations$dist_to <- nodes[ind, node.level]
		relations_subset <- subset(relations, dist_from<=dist_to)
		subg <- igraph::graph.data.frame(d=relations_subset, directed=TRUE, vertices=nodes)
		ig <- subg
	}

	## node.level
	if(is.null(node.level)){
		stop("The node level must be provided!\n")
	}else{
		node.level <- igraph::vertex_attr(ig, node.level)
		if(is.null(node.level)){
			stop("The node level provided does not exist!\n")
		}
	}
	ind <- which(node.level==node.level.value)
	if(length(ind)==0){
		stop("The node level value provided is wrong!\n")
	}
	vec.name <- V(ig)$name[ind]
	
	#################
	## generate code
	V(ig)$code <- 'RT'
	for(k in 1:length(vec.name)){
		x <- vec.name[k]
		neighs.out <- igraph::neighborhood(ig, order=vcount(ig), nodes=x, mode="out")
		neighbors <- names(unlist(neighs.out))
	
		tmp <- igraph::distances(ig, v=x, to=V(ig), mode="out")[1,]
		
		##########
		## remove the rest one in vec.name
		#ind <- match(names(tmp), vec.name[-k])
		## remove the rest one with code
		ind <- match(names(tmp), V(ig)$name[V(ig)$code!='RT'])
		tmp <- tmp[is.na(ind)]
		##########
		
		## only children
		tmp <- tmp[!is.infinite(tmp)]
	
		# self
		V(ig)[names(tmp)[tmp==0]]$code <- LETTERS[k]
		# children
		ind <- which(tmp!=0)
		ttmp <- base::make.unique(as.character(tmp[ind]), sep="-")
		tmp_ls <- lapply(strsplit(ttmp, '-'), function(x){
			x <- as.numeric(x)
			a <- x[1]
			if(a>=10 & a <=35){
				a <- letters[a-9]
			}else if(a>35){
				a <- LETTERS[a-35]
			}
			
			if(length(x)==1){
				c(a,1)
			}else{
				b <- x[2]+1
				if(b>=10 & b <=35){
					b <- letters[b-9]
				}else if(b>35 & b <=61){
					b <- LETTERS[b-35]
				}
				c(a,b)
			}
		})
		tmp_mat <- do.call(rbind, tmp_ls)
		
		if(k<=26){
			code1 <- LETTERS[k]
		}else{
			code1 <- letters[k-26]
		}
		V(ig)[names(tmp)[ind]]$code <- paste0(code1, tmp_mat[,1], tmp_mat[,2])
	}
	#################
	## gp_code
	gp_code <- xA2Net(g=ig, node.label='code', label.wrap.width=30, node.label.size=node.label.size, node.label.color=node.label.color, node.label.alpha=node.label.alpha, node.label.padding=node.label.padding, node.label.arrow=node.label.arrow, node.label.force=node.label.force, node.shape=node.shape, node.xcoord=node.xcoord, node.ycoord=node.ycoord, node.color=node.color, node.color.title=node.color.title, colormap=colormap, ncolors=ncolors, zlim=zlim, node.size.range=node.size.range, title=title, edge.color=edge.color, edge.color.alpha=edge.color.alpha, edge.curve=edge.curve, edge.arrow.gap=edge.arrow.gap)
    #################
    
	## node.table
	if(is.null(node.table)){
		stop("The node name must be provided for the table!\n")
	}else{
		node.table <- igraph::vertex_attr(ig, node.table)
		if(is.null(node.table)){
			stop("The node name for the table provided does not exist!\n")
		}
	}
	
	res_list <- lapply(node.table, function(x){
		if(!is.na(x)){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=node.table.wrap)
			if(length(y)==2){
				paste(y, collapse='\n')
			}else if(length(y)>2){
				paste0(paste(y[1:2],collapse='\n'),'...' )
			}else{
				y
			}
		}else{
			x
		}
	})
	V(ig)$node.table <- unlist(res_list)
    
    df_code <- data.frame(Code=V(ig)$code, Name=V(ig)$node.table, stringsAsFactors=FALSE)
    #tt <- gridExtra::ttheme_default(colhead=list(fg_params=list(parse=TRUE)), base_size=5)
    tt <- gridExtra::ttheme_default(base_size=table.base.size,
    		padding=unit(c(1,1),"mm"), 
    		core=list(fg_params=list(lineheight=0.8, fontfamily='sans'), bg_params=list(fill=c('snow1','snow2'))),
    		colhead=list(bg_params=list(fill="snow3"), fg_params=list(cex=1)),
    		rowhead=list(fg_params=list(hjust=0, x=0, fontface="bold.italic"), bg_params=list(fill=c('snow2','snow1'))),
    )
	
	if(0){
		ind <- base::nchar(df_code$Code)==1 | df_code$Code=='RT'
		df_code_sub <- df_code[ind,]
		df_code_other <- df_code[!ind,]
		ls_gt <- lapply(df_code_sub$Code, function(x){
			ind <- grep(x, df_code_other$Code)
			if(length(ind)>0){
				gridExtra::tableGrob(df_code_other[ind,], theme=tt, rows=NULL)
			}else{
				NULL
			}
		})
		ls_gt <- base::Filter(base::Negate(is.null), ls_gt)
		gridExtra::grid.arrange(grobs=ls_gt, ncol=3, as.table=TRUE)
		gx <- grid::textGrob('aa')
	}else{
		ind <- base::nchar(df_code$Code)==1 | df_code$Code=='RT'
		df_code_sub <- df_code[ind,]
		df_code_other <- df_code[!ind,]
		ls_df <- lapply(df_code_sub$Code, function(x){
			ind <- base::grep(paste0('^',x), df_code_other$Code)
			df_code_other[ind,]
		})
		df_code_order <- rbind(df_code_sub, do.call(rbind, ls_df))
		if(0){
			vec <- ggplot2::cut_interval(1:nrow(df_code_order),length=20)
		}else{
			vec_length <- unlist(lapply(strsplit(df_code_order$Name,'\n'), length))
			vec_length[vec_length==2] <- table.row.space
			vec <- ggplot2::cut_interval(cumsum(vec_length),length=table.nrow)
		}
		ls_gt <- lapply(unique(vec), function(x){
			ind <- which(vec==x)
			y <- df_code_order[ind,]
			
			diff <- table.nrow - sum(vec_length[ind])
			if(diff>0){
				z <- matrix(rep(c('',''),diff),ncol=2)
				colnames(z) <- c('Code','Name')
				y <- rbind(y, z)
			}
			if(0){
				rownames(y) <- y$Code
				y <- subset(y, select='Name')
				gridExtra::tableGrob(y, theme=tt)
			}else{
				gridExtra::tableGrob(y, theme=tt, rows=NULL)
			}
		})
		ls_gt <- base::Filter(base::Negate(is.null), ls_gt)
		#gridExtra::grid.arrange(grobs=ls_gt, ncol=length(unique(vec)), as.table=TRUE)
		
		if(is.null(table.ncol)){
			table.ncol <- 3
		}
		if(length(unique(vec)) < table.ncol){
			table.ncol <- length(unique(vec))
		}
		top <- ''
		if(table.ncol>3){
			top <- quote(paste("page", g, "of", pages))
		}
		gp_table <- gridExtra::marrangeGrob(ls_gt, ncol=table.ncol, nrow=1, as.table=FALSE, top=top)
	}
    
    ls_gp <- list(code=gp_code, table=gp_table)
    invisible(ls_gp)
}
