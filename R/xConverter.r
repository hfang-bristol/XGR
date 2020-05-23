#' Function to convert an object between graph classes
#'
#' \code{xConverter} is supposed to convert an object between classes "igraph", "dgCMatrix", "dtree", "lol", and "json".
#'
#' @param obj an object of class "igraph", "dgCMatrix", "dtree", "lol", and "json"
#' @param from a character specifying the class converted from. It can be one of "igraph", "dgCMatrix", "dtree", "lol", "json"
#' @param to a character specifying the class converted to. It can be one of "igraph", "dgCMatrix", "dtree", "lol", "json" and "igraph_tree"
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return an object of class "igraph", "dgCMatrix", "dtree", "lol", or "json".
#' @note Conversion is supported directly: 1) from 'igraph' to "dgCMatrix","dtree","lol","json","igraph_tree"; 2) from 'dgCMatrix' to "igraph"; 3) from 'dtree' to "igraph","lol","json"; 4) from 'lol' to "dtree","json"; 5) from 'json' to "lol","dtree". In summary: "dgCMatrix" -- "igraph" (hub) -- "dtree" (hub) -- "lol" -- "json". Note: 1) igraph --as.igraph-- phylo --as.hclust/as.phylo-- hclust --as.dendrogram/as.hclust-- dendro; 2) igraph --ggraph::den_to_igraph-- dendro
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xConverter.r
#' @examples
#' # generate a ring graph
#' g <- make_ring(10, directed=TRUE)
#' 
#' # convert the object from 'igraph' to 'dgCMatrix' class
#' xConverter(g, from='igraph', to='dgCMatrix')
#' 
#' \dontrun{
#' # Conversion between 'dgCMatrix' and 'igraph'
#' # ig.EF (an object of class "igraph" storing as a directed graph)
#' g <- xRDataLoader('ig.EF')
#' g
#' 
#' # convert the object from 'igraph' to 'dgCMatrix' class
#' s <- xConverter(g, from='igraph', to='dgCMatrix')
#' s[1:10,1:10]
#'
#' # convert the object from 'dgCMatrix' to 'igraph' class
#' ig <- xConverter(s, from="dgCMatrix", to="igraph")
#' ig
#' 
#' ##############
#' g <- make_graph("Zachary")
#' 
#' # from 'igraph' to "dtree","lol","json"
#' dtree <- xConverter(g, from='igraph', to='dtree')
#' lol <- xConverter(g, from='igraph', to='lol')
#' json <- xConverter(g, from='igraph', to='json')
#' 
#' # from "lol","json" to 'dtree' 
#' dtree <- xConverter(lol, from='lol', to='dtree')
#' dtree <- xConverter(json, from='json', to='dtree')
#' 
#' # from 'dtree' to "igraph"
#' g <- xConverter(dtree, from='dtree', to='igraph')
#' 
#' # force 'igraph' to a tree
#' gtree <- xConverter(g, from='igraph', to='igraph_tree')
#' }

xConverter <- function(obj, from=c("igraph","dgCMatrix","dtree","lol","json"), to=c("dgCMatrix","igraph","dtree","lol","json","igraph_tree"), verbose=TRUE)
{
    
    from <- match.arg(from)
    to <- match.arg(to)
    
    #if (!is(obj,from)){
        #stop(sprintf("The class of your input object '%s' is '%s', mismatched as you intended (from='%s').\n", deparse(substitute(obj)), class(obj), from))
    #}
    
    if(from!="igraph" & to!="igraph"){
        #stop(sprintf("Conversion between '%s' and '%s' is not supported.\n", from, to))
    }
    
    if(from==to){
        warnings(sprintf("Since the class '%s' converted from is the same as the class '%s' converted to, it will return exactly what you input.\n", from, to))
        return(obj)
    }
    
    if(from=="igraph"){
        
        if(to=="dgCMatrix"){
			## get adjacency matrix
			if ("weight" %in% list.edge.attributes(obj)){
				E(obj)$weight <- as.numeric(E(obj)$weight)
				objConverted <- igraph::as_adjacency_matrix(obj, type="both", attr="weight", edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
			}else{
				objConverted <- igraph::as_adjacency_matrix(obj, type="both", attr=NULL, edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
			}
			
        }else if(to %in% c("dtree","lol","json","igraph_tree")){
			## get edges data frame
        	df_edges <- igraph::get.data.frame(obj, what="edges")
        	
        	## append 'ROOT' if not a tree
			df <- df_edges
			root <- setdiff(df[,1], df[,2])
			if(length(root)!=1){
				df_root <- data.frame('ROOT', unique(root), stringsAsFactors=F)
				if(ncol(df_edges)>=3){
					for(i in 3:ncol(df_edges)){
						df_root[,i] <- NA
					}
				}
				colnames(df_root) <- colnames(df_edges)
				df <- rbind(df_root, df)
			}
			
			if(to=="igraph_tree"){
				## igraph but a tree
				objConverted <- igraph::graph_from_data_frame(d=df, directed=T)
			}else{
			
				## df -> dtree ('data.tree' object)
				dtree <- data.tree::FromDataFrameNetwork(df)
				if(to %in% c("lol","json")){
					## dtree -> lol (a hierarchical list object with a root node and children)
					lol <- data.tree::ToListExplicit(dtree, unname=TRUE)
				
					if(0){
						# function: traverse next layer and then recurve
						# return: a hierarchical list object with a root node and children
						func.igraph2list <- function(g, thisNode) {
							nm <- igraph::vertex_attr(g, "name", thisNode)
							childNodes <- V(g)[which(igraph::shortest.paths(g, thisNode, mode="out")==1)]
							if(length(childNodes)==0){
								return(list(name=nm))
							}
							list(name=nm, children=unname(lapply(childNodes, func.igraph2list, g=g)))
						}
						ig <- igraph::graph_from_data_frame(d=df, directed=T)
						root <- setdiff(df[,1], df[,2])
						lol <- func.igraph2list(ig, V(ig)[root])
					}
				
					if(to=="json"){
						## lol -> json
						json <- jsonlite::toJSON(lol)
						objConverted <- json
					}else{
						objConverted <- lol
					}
				}else{
					objConverted <- dtree
				}

			}
			
        }
        
    }else if(from=="dgCMatrix"){
        
        if(to=="igraph"){
			## node info
			nodes <- data.frame(name=rownames(obj))
			nodenames <- rownames(obj)
		
			## adjacency matrix
			adjM <- obj
			tmp <- which(as.matrix(adjM!=0), arr.ind=T)
		
			## un-direct graph
			if(from=="dgCMatrix"){
				ind <- which(tmp[,1]<tmp[,2])
				ttmp <- matrix(0, nrow=length(ind), ncol=2)
				ttmp[1:length(ind),] <- tmp[ind,]
				tmp <- ttmp
			}
		
			## weighted or not
			weight_flag <- T
			if(all(adjM[tmp]==1)){
				weight_flag <- F
			}
			if(weight_flag){
				relations <- data.frame(from=nodenames[tmp[,1]], to=nodenames[tmp[,2]], weight=adjM[tmp])
			}else{
				relations <- data.frame(from=nodenames[tmp[,1]], to=nodenames[tmp[,2]])
			}
		
			## convert to "igraph"
			if(from=="dgCMatrix"){
				objConverted <- igraph::graph_from_data_frame(d=relations, directed=F, vertices=nodes)
			}
			
        }else{
        	warnings(sprintf("Conversion between '%s' and '%s' is not supported; instead converted to the 'igraph' object first.\n", from, to))
        	return(obj)
        }
        
     }else if(from=="dtree"){
     	
     	if(to %in% c("lol","json")){
     		## dtree -> lol
			lol <- data.tree::ToListExplicit(obj, unname=TRUE)
			
			if(to=="json"){
				## lol -> json
				json <- jsonlite::toJSON(lol)
				objConverted <- json
			}else{
				objConverted <- lol
			}

     	}else if(to=="igraph"){
     		df_edges <- data.tree::ToDataFrameNetwork(obj)[,c(1,2)]
     		df_edges[,1] <- gsub('.*/','',df_edges[,1])
     		df_edges[,2] <- gsub('.*/','',df_edges[,2])
     		ig <- igraph::graph_from_data_frame(d=df_edges, directed=T)
     		objConverted <- ig
     		
     	}else{
        	warnings(sprintf("Conversion between '%s' and '%s' is not supported; instead converted to the 'igraph' object first.\n", from, to))
        	return(obj)
        }
     	
     }else if(from=="lol"){
     	
     	if(to=="json"){
     		## lol-> json
			json <- jsonlite::toJSON(obj)
			objConverted <- json

     	}else if(to=="dtree"){
     		## lol-> dtree
			dtree <- data.tree::FromListExplicit(obj)
			objConverted <- dtree
     		
     	}else{
        	warnings(sprintf("Conversion between '%s' and '%s' is not supported; instead converted to the 'data.tree' object first.\n", from, to))
        	return(obj)
        }
     	
     }else if(from=="json"){
     	
     	if(to %in% c("lol","dtree")){
     		## lol-> json
			lol <- jsonlite::fromJSON(obj, simplifyDataFrame=F)
			
			if(to=="dtree"){
     			## lol-> dtree
				dtree <- data.tree::FromListExplicit(lol)
				objConverted <- dtree
     		}else{
     			objConverted <- lol
     		}
			
     	}else{
        	warnings(sprintf("Conversion between '%s' and '%s' is not supported; instead converted to the 'data.tree' object first.\n", from, to))
        	return(obj)
        }
     	
     }
     
    
    if(verbose){
        message(sprintf("Your input object of class '%s' has been converted into an object of class '%s'.", from, to), appendLF=T)
    }
    
    return(objConverted)
}
