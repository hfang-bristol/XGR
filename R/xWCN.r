#' Function to perform weighted gene correlated network analysis
#'
#' \code{xWCN} is supposed to perform weighted gene correlated network analysis. It returns an object of class "cModule".
#'
#' @param data a data matrix or data frame. It can be an expresison-like matrix (with rows for genes and columns for samples) or a similarity matrix (symmetric)
#' @param networkType network type. It determines how the adajcency [0,1] is transformed from correlations ([-1,1]), resulting in a weighted network adjacency matrix. It can be one of "unsigned" (|cor|^power), "signed" ((0.5*(1+cor))^power), and "signed hybrid" (cor^power if cor>0 and 0 otherwise)
#' @param powerVector a vector of soft thresholding powers for which the scale-free topology fit indices are to be calculated
#' @param setBeta NULL or "wgcna" or an integer. If NULL, beta will be determined internally; if "wgcna", estimate of an appropriate soft-thresholding power is the lowest power for which the scale free topology fit R^2 exceeds "RsquaredCut". It can be specified explicitly by the beta value (an integer)
#' @param RsquaredCut desired minimum scale free topology fitting index R^2
#' @param TOMType TOM type. It determines what will be used as inputs. It can be one of "unsigned" (taking as inputs the adjacency), "signed" (taking as inputs the adjacency multiplied by the sign of the adjacency)
#' @param minClusterSize minimum cluster size
#' @param deepSplit an integer [0,4] controlling sensitivity to cluster splitting. The higher the value, the more and smaller clusters will be produced
#' @param merge logical to indicate whether to merge close modules (measured by the eigengene correlation). It only works when input data is expression-like matrix
#' @param cutHeight maximum dissimilarity (1-cor) that qualifies modules for merging
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' an object of class "cModule", a list with following components:
#' \itemize{
#'  \item{\code{mem}: a data frame of n X 3 containing module membership information, where n is the number of nodes, and the 3 columns are "nodes", "modules" (an integer), "rank" (intramodule hub rank)}
#'  \item{\code{expr}: an input data for nodes assigned to modules}
#'  \item{\code{adj}: a weighted network adjacency matrix}
#'  \item{\code{io}: a named list containing inputs/outputs associated, including input paramters ("inputType","networkType","TOMType","minClusterSize"), and output results ("beta","r2","fit","num_modules","pattern_modules","gp_fit","tree")}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xWCN}}
#' @include xWCN.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' cModule <- xWCN(data, networkType="unsigned", powerVector=1:25, setBeta="wgcna", RsquaredCut=0.85, TOMType="signed", minClusterSize=30, merge=F, cutHeight=0.2, verbose=T)
#' head(cModule$mem)
#' names(cModule$io)
#' 
#' # Enrichment analysis of modular genes: using REACTOME pathways
#' ls_eTerm <- xEnricherGenesAdv(cModule, background=cModule$mem$nodes, ontologies=c("REACTOME"), size.range=c(10,500), test="fisher", min.overlap=10, p.tail="one-tail", RData.location=RData.location, plot=T, fdr.cutoff=0.05, displayBy="zscore")
#' gp <- xEnrichForest(ls_eTerm, top_num=5, CI.one=F, FDR.cutoff=0.05, signature=F, drop=T, colormap="yellow-red")
#' 
#' # Network analysis of modular genes: intramodular genes in an interaction network
#' # 1) define netowrk to visualise
#' g <- xDefineNet(network="STRING_high", STRING.only=c("experimental_score","database_score"), RData.location=RData.location)
#' # 2) genes/nodes in module 2
#' nodes_query <- (cModule$mem %>% dplyr::filter(modules==2))$nodes
#' # 3) connectivity of intramodular genes in the interaction network
#' ig <- dnet::dNetInduce(g, nodes_query=nodes_query, knn=0, largest.comp=F)
#' # 4) visualisation
#' ig <- xLayout(ig, layout="gplot.layout.fruchtermanreingold")
#' gp <- xGGnetwork(ig, node.xcoord="xcoord", node.ycoord="ycoord", node.color.alpha=0.5, edge.color.alpha=0.2)
#' ## also size by degree
#' V(ig)$degree <- igraph::degree(ig)
#' gp <- xGGnetwork(ig, node.xcoord="xcoord", node.ycoord="ycoord", node.color.alpha=0.5, edge.color.alpha=0.2, node.size="degree", node.size.range=c(1,5))
#' ## also label by intramodule hub genes (top 10)
#' ind <- match(V(ig)$name,cModule$mem$nodes)
#' V(ig)$rank <- cModule$mem$rank[ind]
#' V(ig)$label <- ''
#' V(ig)$label[V(ig)$rank<=10] <- V(ig)$name[V(ig)$rank<=10]
#' gp <- xGGnetwork(ig, node.xcoord="xcoord", node.ycoord="ycoord", node.color.alpha=0.5, edge.color.alpha=0.2, node.size="degree", node.size.range=c(1,5), node.label="label", node.label.size=2, node.label.force=0.01)
#' }

xWCN <- function(data, networkType=c("unsigned","signed","signed hybrid"), powerVector=1:25, setBeta="wgcna", RsquaredCut=0.85, TOMType=c("signed","unsigned"), minClusterSize=30, deepSplit=c(2,1,0,3,4), merge=T, cutHeight=0.2, verbose=T)
{
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    networkType <- match.arg(networkType)
    TOMType <- match.arg(TOMType)
    
    deepSplit <- deepSplit[1]
    
    if(any(class(data) %in% c("matrix","data.frame"))){
    	data <- as.matrix(data)
		if(verbose){
			message(sprintf("Input data: %d rows X %d columns (%s) ...", nrow(data), ncol(data), as.character(Sys.time())), appendLF=T)
		}
    }else{
    	return(NULL)
    }
    
    ## check whether given data is a similarity matrix
    dataIsExpr <- T
    if(!is.null(dim(data)) & length(dim(data))==2 & dim(data)[1]==dim(data)[2]){
		if(max(abs(data-t(data)),na.rm=T)<1e-12 & (min(data,na.rm=T)>=-1 || max(data, na.rm=TRUE)<=1)){
    		dataIsExpr <- F
    	}
	}

	################################################
	if(dataIsExpr){
		if(verbose){
			message(sprintf("Determine soft-thresholding power beta from input data (expression-like matrix) (%s) ...", as.character(Sys.time())), appendLF=T)
		}
		invisible(utils::capture.output(beta_data <- WGCNA::pickSoftThreshold(t(data), networkType=networkType, RsquaredCut=RsquaredCut, powerVector=powerVector, moreNetworkConcepts=F, corFnc="cor", corOptions=list(use='p'), verbose=0)))
		
	}else{
		if(verbose){
			message(sprintf("Determine soft-thresholding power beta from input data (similarity matrix) (%s) ...", as.character(Sys.time())), appendLF=T)
		}
		beta_data <- WGCNA::pickSoftThreshold.fromSimilarity(data, RsquaredCut=RsquaredCut, powerVector=powerVector, verbose=0)
	}
	
	powerEstimate <- beta_data$powerEstimate
	# fitIndices: power, adjusted R^2 for the linear fit, the linear coefficient, adjusted R^2 for a more complicated fit models, mean connectivity, median connectivity and maximum connectivity; network density, centralization, and heterogeneity
	fitIndices <- beta_data$fitIndices

	# beta and r2
	if(!is.null(setBeta)){
		if(tolower(setBeta)=="wgcna"){
			beta <- beta_data$powerEstimate
			r2 <- fitIndices[fitIndices$Power==beta, "SFT.R.sq"]
		}else if(is.numeric(setBeta)){
			beta <- as.integer(setBeta)
			if(beta>=1 & beta<=nrow(fitIndices)){
				r2 <- fitIndices[fitIndices$Power == beta, 2]
			}else{
				setBeta <- NULL
			}
		}
	}
	# based on r2 changes (< 0.05)
	if(is.null(setBeta)){
		k <- fitIndices[,"mean.k."]
		powers <- fitIndices$Power
		fit <- -sign(fitIndices[,"slope"])*fitIndices[,"SFT.R.sq"]
		for(i in (1:(length(fit)-2))) {
			if(fit[i] >= 0.8) {
				d <- c(abs(fit[i]-fit[i+1]), abs(fit[i]-fit[i+2]), abs(fit[i+1]-fit[i+2]))
				if(max(d) < 0.05) {
					j <- which.max(k[i:(i+2)]) + i - 1
					beta <- powers[j]
					r2 <- fit[j]
					break
				}
			}
		}
	}
	if(verbose){
		message(sprintf("beta: %d", beta), appendLF=T)
		message(sprintf("r2: %f", r2), appendLF=T)
	}
	
	################################################
	if(dataIsExpr){
		if(verbose){
			message(sprintf("Calculate adjacency matrix from input data (expression-like matrix) typed as '%s' powered by %d (%s) ...", networkType, beta, as.character(Sys.time())), appendLF=T)
		}
		adj <- WGCNA::adjacency(t(data), power=beta, type=networkType, corFnc="cor", corOptions=list(use='p'))
	}else{
		if(verbose){
			message(sprintf("Calculate adjacency matrix from input data (similarity matrix) typed as '%s' powered by %d (%s) ...", networkType, beta, as.character(Sys.time())), appendLF=T)
		}
		adj <- WGCNA::adjacency(data, power=beta, type=networkType)
	}

	if(verbose){
		message(sprintf("Calculate Topological Overlap Matrix (TOM) typed as '%s' (%s) ...", TOMType, as.character(Sys.time())), appendLF=T)
	}
	adjMat <- adj
	if(TOMType=='signed'){
		if(dataIsExpr){
			adjMat <- adjMat * sign(WGCNA::cor(t(data), use="p", method="pearson", verbose=0))
		}else{
			adjMat <- adjMat * sign(data)
		}
	}
	tom <- WGCNA::TOMsimilarity(adjMat, TOMType=TOMType, verbose=F)

	if(verbose){
		message(sprintf("Cut tree to determine modules (each having >= %d nodes and '%d' deepSplit) (%s) ...", minClusterSize, deepSplit, as.character(Sys.time())), appendLF=T)
	}
	# TOM based distance measure
	diss <- 1-tom
	# tree
	tree <- stats::hclust(stats::as.dist(diss), method='average')
	# dynamic tree cut to determine modules
	mods <- dynamicTreeCut::cutreeDynamic(dendro=tree, minClusterSize=minClusterSize, method="hybrid", distM=diss, deepSplit=deepSplit, pamRespectsDendro=T, verbose=0)
	
	# merges similar modules
	if(merge & dataIsExpr){
		if(verbose){
			message(sprintf("Also merge close modules based on maximum dissimilarity (%1.2e) (%s) ...", cutHeight, as.character(Sys.time())), appendLF=T)
		}
		if(0){
			# Calculates eigengenes
			me_list <- WGCNA::moduleEigengenes(exprData=t(data), colors=mods, grey=0)
			me_eigen <- me_list$eigengenes
			# Calculates dissimilarity of module eigengenes
			me_diss <- 1 - stats::cor(me_eigen)
			# Clustering module eigengenes
			me_tree <- stats::hclust(as.dist(me_diss), method='average')                      
		}

		# Merging modules                    
		merged_mods <- WGCNA::mergeCloseModules(exprData=t(data), colors=mods, cutHeight=cutHeight, verbose=0)
		mods <- merged_mods$colors
	}
	
	# rename modules as: 1, 2, 3, ..., 0
	# rename modules (by size; in a decreasing order) as: 1, 2, 3, ..., 0
	onames <- setdiff(names(sort(table(mods),decreasing=T)), 0)
	nnames <- 1:length(onames)
	names(nnames) <- onames
	nnames["0"] <- 0
	mods <- nnames[as.character(mods)]
	
	# construct 'df_nodes_modules'
	df_nodes_modules <- data.frame(nodes=rownames(data), modules=mods, stringsAsFactors=F)
	if(verbose){
		message(sprintf("modules: %d", length(table(df_nodes_modules$modules))), appendLF=T)
	}
	
	# return dendrogram of the tree (if similar modules not merged)
	ls_tree <- NULL
	if(!merge){
		# convert to a dendrogram object
		dend <- stats::as.dendrogram(tree)
		# clusters/mods sorted by dend
		ind <- stats::order.dendrogram(dend)
		df_tips <- df_nodes_modules[ind,]
		if(dataIsExpr){
			mat <- data[ind,]
		}else{
			mat <- data[ind,ind]
		}
		ls_tree <- list(dend=dend, tips=df_tips, mat=mat)
		
		if(0){
			#dend <- cModule$io$tree$dend
			#df_tips <- cModule$io$tree$tips
			#mat <- cModule$io$tree$mat

			## visualisation using dendextend: color branches by clusters
			#cols <- xColormap("spectral")(length(unique(df_tips$modules)) - 1)
			#dend2 <- dend %>% dendextend::branches_attr_by_clusters(df_tips$modules, values=cols[unique(df_tips$modules)]) %>% dendextend::set("labels_cex", c(.5,0.5))
			#plot(dend2, type="rectangle", leaflab='none')
			# 1) using ggplot
			#ggd <- dendextend::as.ggdend(dend2)
			#gp <- ggplot(ggd, horiz=T)
			#gp <- ggplot(ggd, labels=F) + scale_y_reverse(expand=c(0.2,0)) + coord_polar(theta="x")
			# 2) using circlize
			#dendextend::circlize_dendrogram(dend2)
			# 3) using gplots
			#allcols <- c('black',cols)
			#RowSideColors <- allcols[df_tips$modules[match(rownames(data),df_tips$nodes)] + 1]
			#gplots::heatmap.2(data, dendrogram="row", Rowv=dend2, Colv=T,col=xColormap('spectral'),trace="none", labRow=NA, labCol=NA, srtCol=20, RowSideColors=RowSideColors)
		}
	}
	
	################################################
	if(verbose){
		message(sprintf("Extract module membership data info (%s) ...", as.character(Sys.time())), appendLF=T)
	}
	# df_nodes_modules_rank: containing nodes (grouped by modules and ranked by connectivity to other nodes within each module)
	modules <- NULL
	df_nodes_modules <- df_nodes_modules %>% dplyr::filter(modules!=0)
	ls_modules <- split(x=df_nodes_modules$nodes, f=df_nodes_modules$modules)
	ls_res <- lapply(1:length(ls_modules), function(i){
		x <- ls_modules[[i]]
		if(length(x)>1){
			mod_adj <- adj[x, x]			
			nodes <- names(sort(apply(mod_adj, 1, sum), decreasing=TRUE))
			data.frame(nodes=nodes, modules=as.numeric(names(ls_modules)[i]), rank=1:length(nodes), stringsAsFactors=F)
		}else{
			NULL
		}
	})
	df_nodes_modules_rank <- do.call(rbind, ls_res)
    
	################################################
	# nodes ordered according to df_nodes_modules_rank
	ind <- match(df_nodes_modules_rank$nodes, rownames(data))
	if(dataIsExpr){
		expr <- data[ind,]
		inputType <- 'expression-like'
	}else{
		expr <- data[ind,ind]
		inputType <- 'similarity'
	}
	adj <- adj[ind,ind]
	io <- list(inputType=inputType, networkType=networkType, TOMType=TOMType, minClusterSize=minClusterSize, beta=beta, r2=r2, fit=fitIndices[,c(1,2,3,5)], num_modules=length(table(df_nodes_modules_rank$modules)))
	
	# summarise modules by looking at intramodule pattern across samples
	io$summary_modules <- NULL
	if(dataIsExpr){
		if(verbose){
			message(sprintf("Summarise modules by pattern (%s) ...", as.character(Sys.time())), appendLF=T)
		}
		samples <- values <- modules <- rank <- NULL
		
		df <-  data.frame(modules=df_nodes_modules_rank$modules, expr, stringsAsFactors=F)
		df_tmp <- df %>% tidyr::gather(key=samples, value=values, -modules)
		## by mean
		modules_mean <- as.data.frame(df_tmp %>% dplyr::group_by(modules,samples) %>% dplyr::summarise(values=mean(values)) %>% tidyr::spread(key="samples", value="values"))
		## by median
		modules_median <- as.data.frame(df_tmp %>% dplyr::group_by(modules,samples) %>% dplyr::summarise(values=median(values)) %>% tidyr::spread(key="samples", value="values"))
		## by eigengene
		me_list <- WGCNA::moduleEigengenes(t(expr), colors=df_nodes_modules_rank$modules)
		mat_tmp <- t(me_list$eigengenes)
		colnames(mat_tmp) <- colnames(expr)
		rownames(mat_tmp) <- 1:nrow(mat_tmp)
		modules_eigengene <- data.frame(modules=rownames(mat_tmp), mat_tmp, stringsAsFactors=F)
		
		## by hub
		df <- data.frame(rank=df_nodes_modules_rank$rank, modules=df_nodes_modules_rank$modules, expr, stringsAsFactors=F)
		rownames(df) <- df_nodes_modules_rank$nodes
		modules_hub <- subset(df, rank==1)[,-1]
		
		io$summary_modules <- list(hub=modules_hub, eigen=modules_eigengene, mean=modules_mean, median=modules_median)
	}
	
	# gp_fit: plot of beta vs r2 and k
	if(1){
		if(verbose){
			message(sprintf("Plot of r2 versus k (%s) ...", as.character(Sys.time())), appendLF=T)
		}
		Power <- new_fit <- mean.k. <- NULL
		
		beta_power <- io$beta
		df <- io$fit
		df$new_fit <- -sign(df[,3])*df[,2]
		
		a <- ceiling(max(df$mean.k.)/5)
		b <- 10^floor(log10(a))
		c <- ceiling(a/b) * b * 5
		
		gp <- ggplot(df, aes(x=Power))
  		gp <- gp + geom_point(aes(y=new_fit, colour="r2")) + geom_line(aes(y=new_fit), color="blue")
  		gp <- gp + geom_point(aes(y=mean.k. /c, colour="k")) + geom_line(aes(y=mean.k. /c), color="lightblue")
  		gp <- gp + scale_y_continuous(sec.axis=sec_axis(~.*c, name="Mean network connectivity (k)"), breaks=seq(0,1,by=0.2))
		gp <- gp + scale_colour_manual(values=c("lightblue", "blue"))
  		gp <- gp + labs(y="Scale-free topology fit (r2)", x="Soft-threshold power (beta)", colour="Criteria")
  		gp <- gp + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())
  		if(!is.null(beta_power)){
			gp <- gp + geom_vline(xintercept=beta_power, color='red', linetype='dashed')
		}
		io$gp_fit <- gp
    }
	
	# tree
	io$tree <- ls_tree
	
####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
	
    cModule <- list(mem		= df_nodes_modules_rank,
                  	data 	= expr,
                  	adj 	= adj,
                  	io	    = io,
                  	call    = match.call()
                 	)
    class(cModule) <- "cModule"
 
    invisible(cModule)   
}