#' Function to identify a pathway crosstalk
#'
#' \code{xCrosstalkGenes} is supposed to identify maximum-scoring pathway crosstalk from an input graph with the node information on the significance (measured as p-values or fdr). It returns an object of class "cPath". 
#'
#' @param data a named input vector containing the significance level for genes (gene symbols) or genomic regions (GR). For this named vector, the element names are gene symbols or GR (in the format of 'chrN:start-end', where N is either 1-22 or X, start/end is genomic positional number; for example, 'chr1:13-20'), the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for gene symbols or GR, 2nd column for the significance level
#' @param entity the entity. It can be either "Gene" or "GR"
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level into scores. If given, those below this are considered significant and thus scored positively. Instead, those above this are considered insignificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to NULL, meaning that no capping is applied
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. It can be one of 'genehancer' (enhancer genes; PMID:28605766) or 'nearby' (nearby genes; if so, please also specify the relevant parameters 'nearby.distance.max', 'nearby.decay.kernel' and 'nearby.decay.exponent' below)
#' @param crosslink.customised the crosslink info with a score quantifying the link of a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if nor provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param cdf.function a character specifying how to transform the input crosslink score. It can be one of 'original' (no such transformation), and 'empirical'  for looking at empirical Cumulative Distribution Function (cdf; as such it is converted into pvalue-like values [0,1])
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR (also over Contexts if many). It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param nearby.distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param nearby.decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param nearby.decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param networks the built-in network. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways
#' @param seed.genes logical to indicate whether the identified network is restricted to seed genes (ie input genes with the signficant level). By default, it sets to true
#' @param subnet.significance the given significance threshold. By default, it is set to NULL, meaning there is no constraint on nodes/genes. If given, those nodes/genes with p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively
#' @param subnet.size the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option
#' @param ontologies the ontologies supported currently. It can be 'AA' for AA-curated pathways, KEGG pathways (including 'KEGG' for all, 'KEGGmetabolism' for 'Metabolism' pathways, 'KEGGgenetic' for 'Genetic Information Processing' pathways, 'KEGGenvironmental' for 'Environmental Information Processing' pathways, 'KEGGcellular' for 'Cellular Processes' pathways, 'KEGGorganismal' for 'Organismal Systems' pathways, and 'KEGGdisease' for 'Human Diseases' pathways), 'REACTOME' for REACTOME pathways or 'REACTOME_x' for its sub-ontologies (where x can be 'CellCellCommunication', 'CellCycle', 'CellularResponsesToExternalStimuli', 'ChromatinOrganization', 'CircadianClock', 'DevelopmentalBiology', 'DigestionAndAbsorption', 'Disease', 'DNARepair', 'DNAReplication', 'ExtracellularMatrixOrganization', 'GeneExpression(Transcription)', 'Hemostasis', 'ImmuneSystem', 'Metabolism', 'MetabolismOfProteins', 'MetabolismOfRNA', 'Mitophagy', 'MuscleContraction', 'NeuronalSystem', 'OrganelleBiogenesisAndMaintenance', 'ProgrammedCellDeath', 'Reproduction', 'SignalTransduction', 'TransportOfSmallMolecules', 'VesicleMediatedTransport')
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 2000
#' @param min.overlap the minimum number of overlaps. Only those terms with members that overlap with input data at least min.overlap (3 by default) will be processed
#' @param fdr.cutoff fdr cutoff used to declare the significant terms. By default, it is set to 0.05
#' @param glayout either a function or a numeric matrix configuring how the vertices will be placed on the plot. If layout is a function, this function will be called with the graph as the single parameter to determine the actual coordinates. This function can be one of "layout_nicely" (previously "layout.auto"), "layout_randomly" (previously "layout.random"), "layout_in_circle" (previously "layout.circle"), "layout_on_sphere" (previously "layout.sphere"), "layout_with_fr" (previously "layout.fruchterman.reingold"), "layout_with_kk" (previously "layout.kamada.kawai"), "layout_as_tree" (previously "layout.reingold.tilford"), "layout_with_lgl" (previously "layout.lgl"), "layout_with_graphopt" (previously "layout.graphopt"), "layout_with_sugiyama" (previously "layout.kamada.kawai"), "layout_with_dh" (previously "layout.davidson.harel"), "layout_with_drl" (previously "layout.drl"), "layout_with_gem" (previously "layout.gem"), "layout_with_mds", and "layout_as_bipartite". A full explanation of these layouts can be found in \url{http://igraph.org/r/doc/layout_nicely.html}
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "cPath", a list with following components:
#' \itemize{
#'  \item{\code{g}: an object of class "igraph". It has graph attribute (enrichment), ndoe attributes (crosstalk)}
#'  \item{\code{gp_paths}: a 'ggplot' object for pathway crosstalk visualisation}
#'  \item{\code{gp_heatmap}: a 'ggplot' object for pathway member gene visualisation}
#' }
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xDefineNet}}
#' @include xCrosstalk.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#'
#' }

xCrosstalk <- function(data, entity=c("Gene","GR"), significance.threshold=NULL, score.cap=NULL, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("genehancer","nearby"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring.scheme=c("max","sum","sequential"), nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, networks=c("KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease"), seed.genes=T, subnet.significance=0.01, subnet.size=NULL, ontologies=c("KEGGenvironmental","KEGG","KEGGmetabolism","KEGGgenetic","KEGGcellular","KEGGorganismal","KEGGdisease"), size.range=c(10,2000), min.overlap=10, fdr.cutoff=0.05, glayout=layout_with_kk, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    entity <- match.arg(entity)
    build.conversion <- match.arg(build.conversion)
    crosslink <- match.arg(crosslink)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
    #networks <- match.arg(networks)
    #ontologies <- match.arg(ontologies)
	
	default.networks <- c("KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease")
	ind <- match(default.networks, networks)
	networks <- default.networks[!is.na(ind)]
	
	network.customised <- NULL
	if(length(networks) == 0){
		network <- "KEGG"
	}else{
		if(length(networks)>1){
			if(any(networks %in% "KEGG")){
				network <- "KEGG"
			}else{
				ls_ig <- lapply(networks, function(network){
					g <- xDefineNet(network=network, weighted=FALSE, verbose=FALSE, RData.location=RData.location)
				})
				network.customised <- xCombineNet(ls_ig, combineBy='union', attrBy="intersect", verbose=TRUE)
			}
		}else{
			network <- networks
		}
	}
	
	if(entity=='Gene'){
		subg <- xSubneterGenes(data=data, network=network, network.customised=network.customised, seed.genes=seed.genes, subnet.significance=subnet.significance, subnet.size=subnet.size, verbose=verbose, RData.location=RData.location)
		
	}else if(entity=='GR'){
		subg <- xGR2xNet(data=data, significance.threshold=significance.threshold, score.cap=score.cap, build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function=cdf.function, scoring.scheme=scoring.scheme, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, network=network, network.customised=network.customised, seed.genes=seed.genes, subnet.significance=subnet.significance, subnet.size=subnet.size, verbose=verbose, RData.location=RData.location)
	}
	
	####################################################
	if(is.null(network.customised)){
		subg_bg <- xDefineNet(network=network, verbose=FALSE, RData.location=RData.location)
	}else{
		subg_bg <- network.customised
	}
	
	if(1){
		ls_eTerm <- xEnricherGenesAdv(list_vec=V(subg)$name, background=V(subg_bg)$name, ontologies=ontologies, size.range=size.range, min.overlap=min.overlap, test="fisher", verbose=F, RData.location=RData.location)
		adjp <- CIl <- NULL
		df_enrichment <- subset(ls_eTerm$df, adjp<fdr.cutoff & CIl>1)
		
	}else{
		eTerm <- xEnricherGenes(data=V(subg)$name, background=V(subg_bg)$name, ontology=ontologies, size.range=size.range, min.overlap=min.overlap, test="fisher", verbose=F, RData.location=RData.location)
		
		eTerm_concise <- xEnrichConciser(eTerm, cutoff=c(0.9,0.5), verbose=T)
		df <- xEnrichViewer(eTerm_concise, top_num="all", sortBy="or", details=TRUE)
			
		if(is.null(df)){
			return(NULL)
		}else{
			df_enrichment <- cbind(group=rep('G1',nrow(df)), ontology=rep(ontologies,nrow(df)), id=rownames(df), df, stringsAsFactors=F)
		}
	}
	
	#df_enrichment <- df_enrichment[,-1]
	####################################################
	
	## list of individual paths
	ls_path <- lapply(1:nrow(df_enrichment), function(j){
		scores <- rep(-1, vcount(subg))
		names(scores) <- V(subg)$name
		x <- df_enrichment$members[j]
		query <- unlist(strsplit(x, ", "))
		scores[query] <- 1
		path <- dnet::dNetFind(subg, scores)
	})
	names(ls_path) <- df_enrichment$name
	
	## remove redundant path
	vec_Redundant <- rep(0, length(ls_path))
	flag_found <- V(ls_path[[1]])$name
	for(j in 2:length(ls_path)){
		path <- ls_path[[j]]
		frac <- length(intersect(V(path)$name,flag_found)) / vcount(path)
		if(frac >=0.5){
			vec_Redundant[j] <- 1
		}else{
			flag_found <- union(V(path)$name,flag_found)
		}
	}
	### make sure each path has no less than 'min.overlap' members
	vec_Redundant[sapply(ls_path,vcount)<min.overlap] <- 1
	### update ls_path and df_enrichment
	ls_path <- ls_path[vec_Redundant==0]
	df_enrichment <- df_enrichment[vec_Redundant==0,]
	
	## combination of individual paths
	ls_path_tmp <- lapply(ls_path, function(path){
		path <- igraph::delete_vertex_attr(path, "score")
		path <- igraph::delete_vertex_attr(path, "type")
	})
	paths_tmp <- xCombineNet(ls_path_tmp, combineBy="union", attrBy="intersect", verbose=TRUE)
	query <- V(paths_tmp)$name
	scores <- rep(-1, vcount(subg))
	names(scores) <- V(subg)$name
	scores[query] <- 1
	paths <- dnet::dNetFind(subg, scores)
	
	## matrix of genes X paths
	ls_vec <- lapply(1:length(ls_path), function(j){
		path <- ls_path[[j]]
		ind <- match(V(paths)$name, V(path)$name)
		vec <- rep(NA, vcount(paths))
		vec[!is.na(ind)] <- names(ls_path)[j]
		vec
	})
	df_res <- do.call(cbind, ls_vec)
	colnames(df_res) <- names(ls_path)
	rownames(df_res) <- V(paths)$name
	vec_sum <- apply(!is.na(df_res), 1, sum)
	## order by vec_sum
	if(1){
		vec_sum <- rev(sort(vec_sum))
		ind <- match(names(vec_sum), rownames(df_res))
		df_res <- df_res[ind,]
	}
	
	## add node attribute 'crosstalk'
	### crosstalk
	crosstalk <- rep('Not assigned', length(vec_sum))
	names(crosstalk) <- names(vec_sum)
	crosstalk[vec_sum>1] <- 'Crosstalk points'
	for(j in 1:ncol(df_res)){
		ind <- !is.na(df_res[,j]) & vec_sum==1
		crosstalk[ind] <- df_res[ind,j]
	}
	### vec_crosstalk
	vec_crosstalk <- rep('Not assigned', vcount(paths))
	names(vec_crosstalk) <- V(paths)$name
	ind <- match(names(vec_crosstalk), names(crosstalk))
	vec_crosstalk[!is.na(ind)] <- crosstalk[ind[!is.na(ind)]]
	V(paths)$crosstalk <- vec_crosstalk
	
	## add graph attribute 'enrichment'
	ls_vec <- lapply(ls_path, function(x) V(x)$name)
	df_enrichment$path_members <- sapply(ls_vec, function(x) paste(x,collapse=", "))
	df_enrichment$nPath <- sapply(ls_vec, length)
	#df_enrichment$label <- paste0(df_enrichment$name, "\nOR=", df_enrichment$or, ",CI=[", df_enrichment$CIl, ",", df_enrichment$CIu, "], FDR=", df_enrichment$adjp, ", n=", df_enrichment$nPath)
	df_enrichment$label <- paste0(df_enrichment$name, "\n[OR=", df_enrichment$or, ", FDR=", df_enrichment$adjp, ", n=", df_enrichment$nPath, "]")
	paths$enrichment <- df_enrichment
	
	###############
	## visualisation
	###############
	if(1){
		#glayout <- layout_with_kk
		set.seed(825)
		glayout <- lapply(list(paths), glayout)[[1]]
		#glayout <- igraph::layout_with_kk(paths)
		#glayout <- igraph::layout_as_tree(paths,root=dnet::dDAGroot(paths),circular=TRUE,flip.y=TRUE)
		V(paths)$xcoord <- glayout[,1]
		V(paths)$ycoord <- glayout[,2]
	}
	V(paths)$color <- -log10(as.numeric(V(paths)$significance))
	
	gp_paths <- xGGnetwork(g=paths, node.label="name", node.label.size=2, node.label.color="black", node.label.alpha=0.8, node.label.padding=0.1, node.label.arrow=0, node.label.force=0.001, node.shape="crosstalk", node.xcoord="xcoord", node.ycoord="ycoord", node.color="color", node.color.title=expression(-log[10]("input")), colormap="jet.top", ncolors=64, node.size.range=5, edge.color="orange",edge.color.alpha=0.3,edge.curve=0,edge.arrow.gap=0.02, title=paste0("Pathway crosstalk involving ",vcount(paths)," genes"), zlim=NULL)
	#gp_paths
	
	df_heatmap <- 0 + !is.na(df_res)
	df_heatmap[df_heatmap==0] <- NA
	for(i in 1:nrow(df_heatmap)){
		x <- df_heatmap[i,]
		df_heatmap[i,!is.na(x)] <- vec_sum[i]
	}
	mat <- t(df_heatmap)
	ind <- match(rownames(mat), df_enrichment$name)
	rownames(mat) <- df_enrichment$label[ind]
	gp_heatmap <- xHeatmap(mat, reorder="none", colormap="grey-black", ncolors=64, barwidth=0.4, x.rotate=90, shape=19, size=2, x.text.size=6,y.text.size=6, na.color='transparent')
	gp_heatmap <- gp_heatmap + theme(legend.title=element_text(size=8), legend.position="none")
	#gp_heatmap
	
    cPath <- list(g = paths,
    			   	gp_paths = gp_paths,
    			   	gp_heatmap = gp_heatmap
                 )
    class(cPath) <- "cPath"
	
    return(cPath)
}
