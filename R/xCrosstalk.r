#' Function to identify a pathway crosstalk
#'
#' \code{xCrosstalkGenes} is supposed to identify maximum-scoring pathway crosstalk from an input graph with the node information on the significance (measured as p-values or fdr). It returns an object of class "cPath". 
#'
#' @param data a named input vector containing the significance level for genes (gene symbols) or genomic regions (GR). For this named vector, the element names are gene symbols or GR (in the format of 'chrN:start-end', where N is either 1-22 or X, start/end is genomic positional number; for example, 'chr1:13-20'), the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for gene symbols or GR, 2nd column for the significance level. Also supported is the input with GR only (without the significance level)
#' @param entity the entity. It can be either "Gene" or "GR"
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level into scores. If given, those below this are considered significant and thus scored positively. Instead, those above this are considered insignificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to NULL, meaning that no capping is applied
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. See \code{\link{xGR2xGenes}} for details
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
#' @param crosstalk.top the number of the top paths will be returned. By default, it is NULL meaning no such restrictions
#' @param glayout either a function or a numeric matrix configuring how the vertices will be placed on the plot. If layout is a function, this function will be called with the graph as the single parameter to determine the actual coordinates. This function can be one of "layout_nicely" (previously "layout.auto"), "layout_randomly" (previously "layout.random"), "layout_in_circle" (previously "layout.circle"), "layout_on_sphere" (previously "layout.sphere"), "layout_with_fr" (previously "layout.fruchterman.reingold"), "layout_with_kk" (previously "layout.kamada.kawai"), "layout_as_tree" (previously "layout.reingold.tilford"), "layout_with_lgl" (previously "layout.lgl"), "layout_with_graphopt" (previously "layout.graphopt"), "layout_with_sugiyama" (previously "layout.kamada.kawai"), "layout_with_dh" (previously "layout.davidson.harel"), "layout_with_drl" (previously "layout.drl"), "layout_with_gem" (previously "layout.gem"), "layout_with_mds", and "layout_as_bipartite". A full explanation of these layouts can be found in \url{http://igraph.org/r/doc/layout_nicely.html}
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "cPath", a list with following components:
#' \itemize{
#'  \item{\code{ig_paths}: an object of class "igraph". It has graph attribute (enrichment, and/or evidence, gp_evidence and membership if entity is 'GR'), ndoe attributes (crosstalk)}
#'  \item{\code{gp_paths}: a 'ggplot' object for pathway crosstalk visualisation}
#'  \item{\code{gp_heatmap}: a 'ggplot' object for pathway member gene visualisation}
#'  \item{\code{ig_subg}: an object of class "igraph".}
#' }
#' @export
#' @seealso \code{\link{xDefineNet}}, \code{\link{xCombineNet}}, \code{\link{xSubneterGenes}}, \code{\link{xGR2xNet}}, \code{\link{xEnricherGenesAdv}}, \code{\link{xGGnetwork}}, \code{\link{xHeatmap}}
#' @include xCrosstalk.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#'
#' # 1) at the gene level
#' data(Haploid_regulators)
#' ## only PD-L1 regulators and thir significance info (FDR)
#' data <- subset(Haploid_regulators, Phenotype=='PDL1')[,c('Gene','FDR')]
#' ## pathway crosstalk
#' cPath <- xCrosstalk(data, entity="Gene", network="KEGG", subnet.significance=0.05, subnet.size=NULL, ontologies="KEGGenvironmental", RData.location=RData.location)
#' cPath
#' ## visualisation
#' pdf("xCrosstalk_Gene.pdf", width=7, height=8)
#' gp_both <- gridExtra::grid.arrange(grobs=list(cPath$gp_paths,cPath$gp_heatmap), layout_matrix=cbind(c(1,1,1,1,2)))
#' dev.off()
#' 
#' # 2) at the genomic region (SNP) level
#' data(ImmunoBase)
#' ## all ImmunoBase GWAS SNPs and their significance info (p-values)
#' ls_df <- lapply(ImmunoBase, function(x) as.data.frame(x$variant))
#' df <- do.call(rbind, ls_df)
#' data <- unique(cbind(GR=paste0(df$seqnames,':',df$start,'-',df$end), Sig=df$Pvalue))
#' ## pathway crosstalk
#' df_xGenes <- xGR2xGenes(data[as.numeric(data[,2])<5e-8,1], format="chr:start-end", crosslink="PCHiC_combined", scoring=T, RData.location=RData.location)
#' mSeed <- xGR2xGeneScores(data, significance.threshold=5e-8, crosslink="PCHiC_combined", RData.location=RData.location)
#' subg <- xGR2xNet(data, significance.threshold=5e-8, crosslink="PCHiC_combined", network="KEGG", subnet.significance=0.1, RData.location=RData.location)
#' cPath <- xCrosstalk(data, entity="GR", significance.threshold=5e-8, crosslink="PCHiC_combined", networks="KEGG", subnet.significance=0.1, ontologies="KEGGenvironmental", RData.location=RData.location)
#' cPath
#' ## visualisation
#' pdf("xCrosstalk_SNP.pdf", width=7, height=8)
#' gp_both <- gridExtra::grid.arrange(grobs=list(cPath$gp_paths,cPath$gp_heatmap), layout_matrix=cbind(c(1,1,1,1,2)))
#' dev.off()
#' 
#' # 3) at the genomic region (without the significance info) level
#' Age_CpG <- xRDataLoader(RData.customised='Age_CpG', RData.location=RData.location)[-1,1]
#' CgProbes <- xRDataLoader(RData.customised='CgProbes', RData.location=RData.location)
#' ind <- match(Age_CpG, names(CgProbes))
#' gr_CpG <- CgProbes[ind[!is.na(ind)]]
#' df <- as.data.frame(gr_CpG)
#' data <- paste0(df$seqnames,':',df$start,'-',df$end)
#' ## pathway crosstalk
#' df_xGenes <- xGR2xGenes(data, format="chr:start-end", crosslink="PCHiC_combined", scoring=T, RData.location=RData.location)
#' subg <- xGR2xNet(data, crosslink="PCHiC_combined", network="KEGG", subnet.significance=0.1, RData.location=RData.location)
#' cPath <- xCrosstalk(data, entity="GR", crosslink="PCHiC_combined", networks="KEGG", subnet.significance=0.1, ontologies="KEGGenvironmental", RData.location=RData.location)
#' cPath
#' }

xCrosstalk <- function(data, entity=c("Gene","GR"), significance.threshold=NULL, score.cap=NULL, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("genehancer","PCHiC_combined","GTEx_V6p_combined","nearby"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring.scheme=c("max","sum","sequential"), nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, networks=c("KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease"), seed.genes=T, subnet.significance=0.01, subnet.size=NULL, ontologies=c("KEGGenvironmental","KEGG","KEGGmetabolism","KEGGgenetic","KEGGcellular","KEGGorganismal","KEGGdisease"), size.range=c(10,2000), min.overlap=10, fdr.cutoff=0.05, crosstalk.top=NULL, glayout=layout_with_kk, verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    entity <- match.arg(entity)
    build.conversion <- match.arg(build.conversion)
    #crosslink <- match.arg(crosslink)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
    #networks <- match.arg(networks)
    #ontologies <- match.arg(ontologies)
	
	default.networks <- c("KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease")
	ind <- match(default.networks, networks)
	networks <- default.networks[!is.na(ind)]
	
	network.customised <- NULL
	network <- "KEGG"
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
	
	ig_subg <- subg
	paths <- NULL
	gp_paths <- NULL
	gp_heatmap <- NULL
	
	if(!is.null(subg) && vcount(subg) >= min.overlap){

		####################################################
		if(is.null(network.customised)){
			subg_bg <- xDefineNet(network=network, verbose=FALSE, RData.location=RData.location)
		}else{
			subg_bg <- network.customised
		}
	
		if(1){
			ls_eTerm <- xEnricherGenesAdv(list_vec=V(subg)$name, background=V(subg_bg)$name, ontologies=ontologies, size.range=size.range, min.overlap=min.overlap, test="fisher", verbose=F, RData.location=RData.location)
			df_enrichment <- ls_eTerm$df
		
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
		
		if(!is.null(df_enrichment)){
	
			adjp <- CIl <- NULL
			df_enrichment <- subset(ls_eTerm$df, adjp<fdr.cutoff & CIl>1)
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
			### only keep the path whose members are unique in number (>50%)
			vec_Redundant <- rep(0, length(ls_path))
			flag_found <- V(ls_path[[1]])$name
			if(length(ls_path) > 2){
				for(j in 2:length(ls_path)){
					path <- ls_path[[j]]
					frac <- length(intersect(V(path)$name,flag_found)) / vcount(path)
					if(frac >=0.5){
						vec_Redundant[j] <- 1
					}else{
						flag_found <- union(V(path)$name,flag_found)
					}
				}
			}
			### make sure each path has no less than 'min.overlap' members
			vec_Redundant[sapply(ls_path,vcount) < 0.8*min.overlap] <- 1
			### update ls_path and df_enrichment
			ls_path <- ls_path[vec_Redundant==0]
			df_enrichment <- df_enrichment[vec_Redundant==0,]
			
			## crosstalk.top
			if(is.null(crosstalk.top)){
				crosstalk.top <- nrow(df_enrichment)
			}
			if(crosstalk.top > nrow(df_enrichment)){
				crosstalk.top <- nrow(df_enrichment)
			}
			crosstalk.top <- as.integer(crosstalk.top)
			crosstalk.cutoff <- df_enrichment[crosstalk.top,'or']
			ind <- which(df_enrichment$or >= crosstalk.cutoff)
			### update ls_path and df_enrichment
			ls_path <- ls_path[ind]
			df_enrichment <- df_enrichment[ind,]
			
			## merge into paths
			if(length(ls_path)>=1){
			
				if(verbose){
					message(sprintf("combination of individual paths (%s) ...", as.character(Sys.time())), appendLF=TRUE)
				}
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
				
				if(verbose){
					message(sprintf("matrix of genes X paths (%s) ...", as.character(Sys.time())), appendLF=TRUE)
				}
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
				### construct data.frame 'df_tmp'
				df_tmp <- data.frame(num=vec_sum, df_res, gene=rownames(df_res), stringsAsFactors=F)
				### arrange_at, 'num' in ascending order, and other variables except 'num' (select using -one_of("num")) in descending order
				#df_tmp <- df_tmp %>% dplyr::arrange_at(dplyr::vars(dplyr::desc("num"), -dplyr::one_of("num")))
				df_tmp <- df_tmp %>% dplyr::arrange_all()
				colnames(df_tmp)[2:(ncol(df_tmp)-1)] <- colnames(df_res)
				### update 'vec_sum' and 'df_res'
				vec_sum <- df_tmp$num
				names(vec_sum) <- df_tmp$gene
				df_res <- df_tmp[,2:(ncol(df_tmp)-1)]
				rownames(df_res) <- df_tmp$gene
				
				if(verbose){
					message(sprintf("add node attribute 'crosstalk' (%s) ...", as.character(Sys.time())), appendLF=TRUE)
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
	
				if(verbose){
					message(sprintf("add graph attribute 'enrichment' (%s) ...", as.character(Sys.time())), appendLF=TRUE)
				}
				## add graph attribute 'enrichment'
				ls_vec <- lapply(ls_path, function(x) V(x)$name)
				df_enrichment$path_members <- sapply(ls_vec, function(x) paste(x,collapse=", "))
				df_enrichment$nPath <- sapply(ls_vec, length)
				#df_enrichment$label <- paste0(df_enrichment$name, "\nOR=", df_enrichment$or, ",CI=[", df_enrichment$CIl, ",", df_enrichment$CIu, "], FDR=", df_enrichment$adjp, ", n=", df_enrichment$nPath)
				df_enrichment$label <- paste0(df_enrichment$name, "\n[OR=", df_enrichment$or, ", FDR=", df_enrichment$adjp, ", n=", df_enrichment$nPath, "]")
				paths$enrichment <- df_enrichment

				if(verbose){
					message(sprintf("visualisation (%s) ...", as.character(Sys.time())), appendLF=TRUE)
				}
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
				
				############
				## vec_crosstalk for gp_paths
				path_names <- c('Crosstalk points', names(ls_path))
				names(path_names) <- paste0(0:(length(path_names)-1), '. ', path_names)
				ind <- match(V(paths)$crosstalk, path_names)
				vec_crosstalk <- names(path_names)[ind]
				############

				if(verbose){
					message(sprintf("gp_paths (%s) ...", as.character(Sys.time())), appendLF=TRUE)
				}			
				gp_paths <- xGGnetwork(g=paths, node.label="name", node.label.size=2, node.label.color="black", node.label.alpha=0.8, node.label.padding=0.1, node.label.arrow=0, node.label.force=0.001, node.shape=vec_crosstalk, node.xcoord="xcoord", node.ycoord="ycoord", node.color="color", node.color.title=expression(-log[10]("input")), colormap="jet.top", ncolors=64, node.size.range=5, edge.color="orange",edge.color.alpha=0.3,edge.curve=0,edge.arrow.gap=0.02, title=paste0("Pathway crosstalk involving ",vcount(paths)," genes"), zlim=NULL)
				#gp_paths
			
				if(length(ls_path)>1){
				
					if(verbose){
						message(sprintf("gp_heatmap (%s) ...", as.character(Sys.time())), appendLF=TRUE)
					}
				
					df_heatmap <- 0 + !is.na(df_res)
					df_heatmap[df_heatmap==0] <- NA
					for(i in 1:nrow(df_heatmap)){
						x <- df_heatmap[i,]
						df_heatmap[i,!is.na(x)] <- vec_sum[i]
					}
					mat_heatmap <- t(df_heatmap)
					ind <- match(rownames(mat_heatmap), df_enrichment$name)
					rownames(mat_heatmap) <- df_enrichment$label[ind]
					gp_heatmap <- xHeatmap(mat_heatmap, reorder="none", colormap="white-skyblue-darkblue", zlim=c(0,max(mat_heatmap,na.rm=T)), ncolors=64, barwidth=0.4, x.rotate=90, shape=19, size=2, x.text.size=6,y.text.size=6, na.color='transparent')
					gp_heatmap <- gp_heatmap + theme(legend.title=element_text(size=8), legend.position="none")
					colsep <- cumsum(table(vec_sum))
					colsep <- colsep[-length(colsep)]
					gp_heatmap <- gp_heatmap + geom_vline(xintercept=colsep+0.5,color="grey90",size=0.5)
					#gp_heatmap
					
				
				}else{
					gp_heatmap <- ggplot() + theme_void()
				}
				
				## update graph attributes 'evidence' and 'gp_evidence'
				if(!is.null(paths$evidence)){
				
					if(verbose){
						message(sprintf("update graph attribute 'evidence' and 'gp_evidence' (%s) ...", as.character(Sys.time())), appendLF=TRUE)
					}
				
					### update graph attribute 'evidence'
					ind <- match(V(paths)$name, paths$evidence$Gene)
					evidence <- paths$evidence[ind[!is.na(ind)], c('GR','Gene','Score')]
					paths$evidence <- evidence
					### update graph attribute 'gp_evidence'
					Gene <- Score <- NULL
					mat_evidence <- tidyr::spread(evidence, key=Gene, value=Score)
					mat <- mat_evidence[,-1]
					rownames(mat) <- mat_evidence[,1]
					#### sort by chromosome, start and end
					ind <- xGRsort(rownames(mat))
					mat <- mat[ind,]
					####
					
					################
					## obtain rowsep
					rowsep <- xGRsep(rownames(mat))
					rowsep <- nrow(mat) - rowsep
					################
										
					if(verbose){
						message(sprintf("keep the same order columns in mat_heatmap (%s) ...", as.character(Sys.time())), appendLF=TRUE)
					}
					#### keep the same order columns in mat_heatmap
					if(!is.null(mat_heatmap)){
						ind <- match(colnames(mat_heatmap), colnames(mat))
						mat <- mat[,ind]
						df_membership <- rbind(mat_heatmap, mat)
						
						########################
						## replace '\n' with ' '
						rownames(df_membership) <- gsub('\n',' ', rownames(df_membership))
						########################
						
						## add graph attribute 'membership'
						paths$membership <- df_membership
						
						## add graph attribute 'gp_membership'
						zlim_max <- max(mat_heatmap,na.rm=T)
						zlim_min <- -zlim_max
						x <- mat
						x[!is.na(x)] <- zlim_min
						y <- rbind(mat_heatmap, x)
						rownames(y) <- gsub('\n.*',' ', rownames(y))
						gp_membership <- xHeatmap(y, reorder="none", colormap="grey-grey-white-skyblue-darkblue", zlim=c(zlim_min, zlim_max), ncolors=64, barwidth=0.4, x.rotate=90, shape=19, size=2, x.text.size=6,y.text.size=6, na.color='transparent')
						gp_membership <- gp_membership + theme(legend.title=element_text(size=8), legend.position="none") + scale_y_discrete(position="right")
						gp_membership <- gp_membership + geom_hline(yintercept=rowsep+0.5,color="grey90",size=0.5) + geom_vline(xintercept=colsep+0.5,color="grey90",size=0.5) + geom_hline(yintercept=nrow(mat)+0.5,color="grey50",size=0.5)
						
						paths$gp_membership <- gp_membership
					}
					
					if(verbose){
						message(sprintf("gp_evidence (%s) ...", as.character(Sys.time())), appendLF=TRUE)
					}
					
					## add graph attribute 'gp_evidence'
					gp_evidence <- xHeatmap(mat, reorder="none", colormap="spectral", ncolors=64, barwidth=0.4, x.rotate=90, shape=19, size=2, x.text.size=6,y.text.size=6, na.color='transparent')
					gp_evidence <- gp_evidence + theme(legend.title=element_text(size=8), legend.position="left") + scale_y_discrete(position="right")
					gp_evidence <- gp_evidence + geom_hline(yintercept=rowsep+0.5,color="grey90",size=0.5) + geom_vline(xintercept=colsep+0.5,color="grey90",size=0.5)
						
					paths$gp_evidence <- gp_evidence
				}
				
			}
		}
	}
	
    cPath <- list(ig_paths = paths,
    			  gp_paths = gp_paths,
    			  gp_heatmap = gp_heatmap,
    			  ig_subg =ig_subg
                 )
    class(cPath) <- "cPath"
	
    return(cPath)
}
