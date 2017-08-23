#' Function to generate a graphml file from a graph object of class "igraph"
#'
#' \code{xA2GraphML} is supposed to generate a graphml file from a graph object of class "igraph".
#'
#' @param data a data frame
#' @param query the identity of a pathway in query. The full list of pathways in human can be found at \url{http://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=hsa}. For example, 'AA:hsa04672' for 'NOD-like receptor signaling pathway', where the prefix 'AA:' can be ignored. Alternatively, it can be key words describing the pathway
#' @param curation the type of curation. It can be one of "manual" (only the manual one), "automatic" (only the automatic one) or "any" (first the manual one then the automatic one)
#' @param node.label a character specifying which column used for node labelling. By default, it is 'label'
#' @param node.color a character specifying which column used for node coloring. By default, it is 'lfc'
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param nlegend the number of colors specified in the legend. By default, it is 11
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param legend.title the legend title. By default, it is ''
#' @param node.tooltip a character specifying which column used for node tooltip. By default, it is 'tooltip'. If not found, it will be 'Symbol-Name-Color'
#' @param node.highlight a character specifying which column used for node highlighting. By default, it is 'fdr'. If so, those highlighted will have bold and larger labels
#' @param node.highlight.cutoff a numeric specifying the cutoff for node highlighting. By default, it is 0.05 meaninght those less than this cutoff will be highlighted
#' @param edge.color a character specifying the edge colors. By default, it is '#00000033'
#' @param edge.width the edge width. By default, it is 1
#' @param color.gene a character specifying the gene node colors. By default, it is '#dddddd'
#' @param color.thispath a character specifying the color for this pathway node. By default, it is '#bbbbbb'
#' @param color.otherpath a character specifying the color for other pathway nodes. By default, it is '#cccccc'
#' @param filename the without-extension part of the name of the output file. By default, it is 'xA2GraphML'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{xA2GraphML}}
#' @include xA2GraphML.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' 
#' xA2GraphML(query="AA:hsa04630", RData.location=RData.location)
#' 
#' ## load GWAS genes
#' GWAS_Gene <- xRDataLoader(RData.customised='GWAS_Gene', RData.location=RData.location)
#' data <- GWAS_Gene %>% dplyr::filter(Odds_Ratio!='NULL' & Disease_ID=='RA') %>% dplyr::transmute(label=Symbol, lfc=log2(as.numeric(Odds_Ratio)), fdr=Pvalue) %>% dplyr::group_by(label) %>% dplyr::summarise(lfc=max(lfc), fdr=min(fdr))
#' 
#' ## manual one (the same as curation='any')
#' xA2GraphML(data, query="AA:hsa04630", curation='manual', node.label="label", node.color="lfc", node.highlight='fdr', node.highlight.cutoff=5e-8, filename='xA2GraphML', legend.title='log2(Odds ratio)', zlim=c(-1,1), RData.location=RData.location)
#' ## automatic one
#' xA2GraphML(data, query="AA:hsa04630", curation='automatic', node.label="label", node.color="lfc", node.highlight='fdr', node.highlight.cutoff=5e-8, filename='xA2GraphML', legend.title='log2(Odds ratio)', zlim=c(-1,1), RData.location=RData.location)
#' 
#' ## key words 
#' xA2GraphML(data, query="Asthma", curation='any', node.label="label", node.color="lfc", node.highlight='fdr', node.highlight.cutoff=5e-8, filename='xA2GraphML', RData.location=RData.location, legend.title='log2(Odds ratio)', zlim=c(-1,1))
#' }

xA2GraphML <- function(data=NULL, query="AA:hsa04672", curation=c('manual','automatic','any'), node.label='label', node.color='lfc', colormap='deepskyblue-white-darkorange', ncolors=64, nlegend=11, zlim=NULL, legend.title='', node.tooltip='tooltip', node.highlight='fdr', node.highlight.cutoff=0.05, edge.color="#00000033", edge.width=1, color.gene='#dddddd', color.thispath='#bbbbbb', color.otherpath='#cccccc', filename='xA2GraphML', verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    curation <- match.arg(curation)
    
    if(any(class(data) %in% c('tbl_df', 'tbl'))){
    	data <- as.data.frame(data)
    }
    
	if(all(c(node.label, node.color) %in% colnames(data))){
		df <- data.frame(Symbol=data[,node.label], LFC=signif(data[,node.color],digits=2), stringsAsFactors=FALSE)
	}else{
		df <- data.frame(Symbol=NA, LFC=0, stringsAsFactors=FALSE)
		data <- df
		#return(NULL)
	}

    ## FDR (by default, 1)
    df$FDR <- rep(1, nrow(df))
    if(!is.null(node.highlight)){
		if(node.highlight %in% colnames(data)){
			df$FDR <- signif(data[,node.highlight],digits=2)
		}
    }

    ## tooltip
    vec_description <- XGR::xSymbol2GeneID(df$Symbol, details=TRUE, verbose=verbose, RData.location=RData.location)$description
    df$tooltip <- paste0('Symbol: ', df$Symbol, '\nName: ', vec_description, '\nColor: ', df$LFC, '\nHightlight: ', df$FDR)
    if(!is.null(node.tooltip)){
		if(node.tooltip %in% colnames(data)){
			df$tooltip <- data[,node.tooltip]
		}
    }
	
    ######################################################################################
    # df$node.color
    ######################################################################################
    ## node.color (by default, "#BFFFBF")
	df_legends <- NULL
    if (!is.null(node.color)){
    	ind <- match(node.color, colnames(data))

        if(!is.na(ind)){
        	
        	pattern <- as.numeric(df$LFC)
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
            
            #############
            ## df_legends
            colors <- xColormap(colormap=colormap)(ncolors)
            legend_colors <- colors[round(seq(1,ncolors,length.out=nlegend))]
            df_legends <- data.frame(name=paste0('Legend',(nlegend-1):0), colors=legend_colors, labels=signif(seq(min(zlim),max(zlim),length.out=nlegend),digits=2), stringsAsFactors=FALSE)
            ############# 
            
        }else{
            #warning("The input 'pattern' is ignored. Please check the help for enabling your input")
            node.color <- rep("#BFFFBF", nrow(data))
        }
    }else{
        node.color <- rep("#BFFFBF", nrow(data))
    }
    df$node.color <- node.color
    
    ######################################################################################
    ######################################################################################
	
    #####################################
    if(curation %in% c('any','manual')){
		AA.path <- xRDataLoader(RData.customised="AA.path", RData.location=RData.location)
		info <- AA.path$info
		path <- gsub('^AA:', '', info$path)
		query <- gsub('^AA:', '', query)
		query <- gsub('^path:', '', query)
		manual_ind <- match(query, path)
		if(is.na(manual_ind)){
			manual_ind <- grep(query, info$name)
			if(length(manual_ind)>1){
				warning(sprintf("Manual curation: %d found for queried %s: only 1st kept", length(manual_ind), query), appendLF=TRUE)
				manual_ind <- manual_ind[1]
			}else if(length(manual_ind)==0){
				warning(sprintf("Manual curation: no found for queried '%s'", query), appendLF=TRUE)
				if(curation=='manual'){
					return(NULL)
				}
			}
		}
    }
    
    #####################################
    if(length(manual_ind) != 0){
    
		if(verbose){
			now <- Sys.time()
			message(sprintf("Manual curation: visualising '%s: %s' (%s) ...", AA.path$info$path[manual_ind], AA.path$info$name[manual_ind], as.character(now)), appendLF=TRUE)
		}
		detail <- AA.path$detail[[manual_ind]]
		df_nodes <- detail$nodes
		df_edges <- detail$edges

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
		n_total <- 0
		n_matched <- 0
		for( i in 1:nrow(df_nodes)){
			
			########################
			if(grepl('Legend',df_nodes$name[i])){
				if(is.null(df_legends)){
					next
				}
			}
			########################
					
			ind_found <- match(df_nodes$name[i], df$Symbol)
			ind_legend <- match(df_nodes$name[i], df_legends$name)
			#############
			# gene node
			flag_gene <- grepl('#BFFFBF', df_nodes$fill[i])
			# this pathway node
			flag_thispath <- grepl('#00FF00', df_nodes$fill[i])
			# other pathway node
			flag_otherpath <- grepl('#C0C0C0', df_nodes$fill[i])
			#############
			if(flag_gene & is.na(ind_legend)){
				n_total <- n_total+1;
				if(!is.na(ind_found)){
					n_matched <- n_matched+1;
				}
			}
		}
		if(verbose){
			now <- Sys.time()
			message(sprintf("\t %d genes in total and %d genes matched (%s)", n_total, n_matched, as.character(now)), appendLF=TRUE)
		}
		
		############
		ls_nodes <- lapply(1:nrow(df_nodes), function(i){
		
			########################
			if(grepl('Legend',df_nodes$name[i])){
				if(is.null(df_legends)){
					return(NULL)
				}
			}
			########################
		
			ind_found <- match(df_nodes$name[i], df$Symbol)
			ind_legend <- match(df_nodes$name[i], df_legends$name)
		
			#############
			# gene node
			flag_gene <- grepl('#BFFFBF', df_nodes$fill[i])
			# this pathway node
			flag_thispath <- grepl('#00FF00', df_nodes$fill[i])
			# other pathway node
			flag_otherpath <- grepl('#C0C0C0', df_nodes$fill[i])
			#############
			
			if(flag_gene & is.na(ind_legend)){
				n_total <- n_total+1;
				if(!is.na(ind_found)){
					n_matched <- n_matched+1;
				}
			}

			k <- 0
			vec <- vector()
		
			##############
			## node id
			##############
			k <- k+1
			vec[k] <- paste0('<node id="', df_nodes$node_id[i], '">')

			##############
			## data key="d1"
			## data key="d2"
			##############
			if(!is.na(ind_found)){
				k <- k+1
				vec[k] <- paste0('<data key="d1"><![CDATA[', "http://www.genecards.org/cgi-bin/carddisp.pl?gene=", df$Symbol[ind_found], ']]></data>')
			
				k <- k+1
				vec[k] <- paste0('<data key="d2"><![CDATA[', df$tooltip[ind_found], ']]></data>')
			}else{
			
				if(!is.na(ind_legend)){
					k <- k+1
					vec[k] <- paste0('<data key="d2"><![CDATA[', df_legends$labels[ind_legend], ']]></data>')
				
				}else if(flag_gene){
					k <- k+1
					vec[k] <- paste0('<data key="d1"><![CDATA[', "http://www.genecards.org/cgi-bin/carddisp.pl?gene=", df_nodes$name[i], ']]></data>')
			
					k <- k+1
					vec[k] <- paste0('<data key="d2"><![CDATA[', df_nodes$name[i], ']]></data>')
				}
			}
		
			##############
			## data key="d3"
			##############
			k <- k+1
			vec[k] <- paste0('<data key="d3">')
		
			##############
			## y:ShapeNode
			##############
			k <- k+1
			vec[k] <- paste0('<y:ShapeNode>')

			##############
			## y:Geometry
			##############    	
			k <- k+1
			vec[k] <- paste0('<y:Geometry ', df_nodes$geometry[i], '/>')
		
			##############
			## y:Fill
			##############
			fill <- df_nodes$fill[i]
			if(!is.na(ind_found)){
				fill <- gsub('#BFFFBF', df$node.color[ind_found], fill)
			}else{
				if(!is.na(ind_legend)){
					fill <- gsub('#BFFFBF', df_legends$colors[ind_legend], fill)
				}else{
					# gene node
					fill <- gsub('#BFFFBF', color.gene, fill)
					# this pathway node
					fill <- gsub('#00FF00', color.thispath, fill)
					# other pathway node
					fill <- gsub('#C0C0C0', color.otherpath, fill)
				}
			}
			k <- k+1
			vec[k] <- paste0('<y:Fill ', fill, '/>')
		
			##############
			## y:BorderStyle
			##############
			borderstyle <- df_nodes$borderstyle[i]
			if(!is.na(ind_legend)){
				borderstyle <- gsub('color="#\\w+"', 'color="#cccccc"', borderstyle)
				borderstyle <- gsub('width=".*"', 'width="0.8"', borderstyle)
			}
			if(!is.na(ind_found)){
				borderstyle <- gsub('color="#\\w+"', 'color="#ffffff00"', borderstyle)
				borderstyle <- gsub('width=".*"', 'width="0"', borderstyle)
			}else if(flag_gene){
				borderstyle <- gsub('color="#\\w+"', 'color="#ffffff00"', borderstyle)
				borderstyle <- gsub('width=".*"', 'width="0"', borderstyle)
			}
			k <- k+1
			vec[k] <- paste0('<y:BorderStyle ', borderstyle, '/>')
		
			#alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="14" fontStyle="bold" hasBackgroundColor="false" hasLineColor="false" height="20.48828125" horizontalTextPosition="center" iconTextGap="0" modelName="internal" modelPosition="tl" textColor="#000000" verticalTextPosition="bottom" visible="true" width="98.2265625" x="4.0" y="4.0"
			##############
			## y:NodeLabel
			##############
			### nodelabel
			nodelabel <- df_nodes$nodelabel[i]
			nodelabel <- gsub('fontFamily="\\w+"', 'fontFamily="Arial"', nodelabel)
			if(!is.na(ind_found)){
				if(df$FDR[ind_found] < node.highlight.cutoff){
					## found gene nodes highlighted
					nodelabel <- gsub('fontStyle="\\w+"', 'fontStyle="bolditalic"', nodelabel)
					nodelabel <- gsub('fontSize="\\w+"', 'fontSize="13"', nodelabel)
					#nodelabel <- gsub('verticalTextPosition=', 'underlinedText="true" verticalTextPosition=', nodelabel)
				}else{
					## found gene nodes not highlighted
					nodelabel <- gsub('fontStyle="\\w+"', 'fontStyle="italic"', nodelabel)
					nodelabel <- gsub('fontSize="\\w+"', 'fontSize="12"', nodelabel)
				}
			
			}else if(flag_gene & is.na(ind_legend)){
				## gene nodes but not found
				nodelabel <- gsub('fontStyle="\\w+"', 'fontStyle="italic"', nodelabel)
				nodelabel <- gsub('fontSize="\\w+"', 'fontSize="11"', nodelabel)
			}
			### name_display
			name_display <- df_nodes$name[i]
			name_display <- gsub('@@', '\n', name_display)
			if(!is.na(ind_legend)){
				name_display <- df_legends$labels[ind_legend]
			}else{
				if(name_display %in% c("5-star rating", "5-star\nrating")){
					name_display <- legend.title
					nodelabel <- gsub('modelPosition="\\w+"', 'modelPosition="t"', nodelabel)
				}
			}
			k <- k+1
			vec[k] <- paste0('<y:NodeLabel ', nodelabel, '>', name_display, '</y:NodeLabel>')
		
			##############
			## y:Shape
			##############
			k <- k+1
			vec[k] <- paste0('<y:Shape ', df_nodes$shape[i], '/>')
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
		ls_edges <- lapply(1:nrow(df_edges), function(i){
		
			k <- 0
			vec <- vector()    	
		
			k <- k+1
			vec[k] <- paste0('<edge id="', df_edges$edge_id[i], '" source="', df_edges$source[i], '" target="', df_edges$target[i], '">')
			k <- k+1
			vec[k] <- paste0('<data key="d5">')
			k <- k+1
			vec[k] <- paste0('<y:GenericEdge configuration="DEFAULT">')
			k <- k+1
			vec[k] <- paste0('<y:Path ', df_edges$path[i], '/>')
		
			linestyle <- df_edges$linestyle[i]
			linestyle <- gsub('#000000', '#888888', linestyle)
			k <- k+1
			vec[k] <- paste0('<y:LineStyle ', linestyle, '/>')
		
			k <- k+1
			vec[k] <- paste0('<y:Arrows ', df_edges$arrows[i], '/>')
		
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

    }else if(curation %in% c('automatic','any')){
		
		ls_ig <- xRDataLoader(RData.customised="ig.KEGG.list", RData.location=RData.location)
		kegg <- gsub('^path:', '', sapply(ls_ig,function(x) x$path))
		query <- gsub('^AA:', '', query)
		query <- gsub('^path:', '', query)
		automatic_ind <- match(query, kegg)
		if(is.na(automatic_ind)){
			automatic_ind <- grep(query, names(kegg))
			if(length(automatic_ind)==0){
				warning(sprintf("Automated: no found for queried '%s'", query), appendLF=TRUE)
				return(NULL)
				
			}else if(length(automatic_ind)>1){
				warning(sprintf("Automated: %d found for queried %s: only 1st kept", length(automatic_ind), query), appendLF=TRUE)
				automatic_ind <- automatic_ind[1]
			}
		}
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("Automated: visualising '%s: %s' (%s) ...", kegg[automatic_ind], names(kegg[automatic_ind]), as.character(now)), appendLF=TRUE)
		}
		ig <- ls_ig[[automatic_ind]]
		
		if(0){
			ind <- match(V(ig)$name, df$Symbol)
			nodes_query <- V(ig)$name[!is.na(ind)]
			## including incoming neighbors
			order <- vcount(ig)
			order <- 2
			neighs.out <- igraph::neighborhood(ig, order=2, nodes=nodes_query, mode="in")
			neighbors <- unique(names(unlist(neighs.out)))
			#ind <- match(neighbors, df$Symbol)
			#neighbors <- neighbors[!is.na(ind)]
			subg <- dnet::dNetInduce(ig, nodes_query=neighbors, knn=0, remove.loops=TRUE, largest.comp=FALSE)
		}else{
			subg <- ig
		}
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("\t %d genes and %d edges (%s)", vcount(subg), ecount(subg), as.character(now)), appendLF=TRUE)
		}
		
		## calling xGraphML
		glayout <- igraph::layout_with_kk(subg)
		V(subg)$xcoord <- glayout[,1]
		V(subg)$ycoord <- glayout[,2]
		ind <- match(V(subg)$name, df$Symbol)
		V(subg)$LFC <- df$LFC[ind]
		V(subg)$tooltip <- df$tooltip[ind]
		V(subg)$tooltip[is.na(V(subg)$tooltip)] <- V(subg)$description[is.na(V(subg)$tooltip)]
		V(subg)$FDR <- df$FDR[ind]
		V(subg)$node.size <- ifelse(V(subg)$FDR<node.highlight.cutoff & !is.na(V(subg)$FDR), 25, 15)

		output <- xGraphML(g=subg, node.label="name", node.label.size=10, node.tooltip="tooltip", node.xcoord="xcoord", node.ycoord="ycoord", node.color.na=color.gene, node.color="LFC", node.link="http://www.genecards.org/cgi-bin/carddisp.pl?gene=", nlegend=nlegend, node.size='node.size', node.coord.scale=300, zlim=zlim, colormap=colormap, filename=filename)
		invisible(output)
		#######################################################
			
    }
}






