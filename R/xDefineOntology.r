#' Function to define ontology and its annotations
#'
#' \code{xDefineOntology} is supposed to define ontology and its annotations. It returns an object of class "aOnto".
#'
#' @param ontology the ontology supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "PSG" for phylostratigraphy (phylostratific age), "PS" for sTOL-based phylostratific age information, "PS2" for the collapsed PS version (inferred ancestors being collapsed into one with the known taxonomy information), "SF" for SCOP domain superfamilies, "Pfam" for Pfam domain families, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPCM" for Human Phenotype Clinical Modifier, "HPMA" for Human Phenotype Mortality Aging, "MP" for Mammalian Phenotype, "EF" for Experimental Factor Ontology (used to annotate GWAS Catalog genes), Drug-Gene Interaction database ("DGIdb") for druggable categories, tissue-specific eQTL-containing genes from GTEx ("GTExV4", "GTExV6p" and "GTExV7"), crowd extracted expression of differential signatures from CREEDS ("CreedsDisease", "CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP", "CreedsDrugDN", "CreedsGene", "CreedsGeneUP" and "CreedsGeneDN"), KEGG pathways (including 'KEGG' for all, 'KEGGmetabolism' for 'Metabolism' pathways, 'KEGGgenetic' for 'Genetic Information Processing' pathways, 'KEGGenvironmental' for 'Environmental Information Processing' pathways, 'KEGGcellular' for 'Cellular Processes' pathways, 'KEGGorganismal' for 'Organismal Systems' pathways, and 'KEGGdisease' for 'Human Diseases' pathways), 'REACTOME' for REACTOME pathways or 'REACTOME_x' for its sub-ontologies (where x can be 'CellCellCommunication', 'CellCycle', 'CellularResponsesToExternalStimuli', 'ChromatinOrganization', 'CircadianClock', 'DevelopmentalBiology', 'DigestionAndAbsorption', 'Disease', 'DNARepair', 'DNAReplication', 'ExtracellularMatrixOrganization', 'GeneExpression(Transcription)', 'Hemostasis', 'ImmuneSystem', 'Metabolism', 'MetabolismOfProteins', 'MetabolismOfRNA', 'Mitophagy', 'MuscleContraction', 'NeuronalSystem', 'OrganelleBiogenesisAndMaintenance', 'ProgrammedCellDeath', 'Reproduction', 'SignalTransduction', 'TransportOfSmallMolecules', 'VesicleMediatedTransport'), and the molecular signatures database (Msigdb, including "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7"), and the SIFTS database ("SIFTS2GOBP" for Gene Ontology Biological Process, "SIFTS2GOMF" for Gene Ontology Molecular Function, "SIFTS2GOCC" for Gene Ontology Cellular Component), and 'EnrichrX' for Enrichr libraries  (where X can be "AchillesFitnessD","AchillesFitnessI","ARCHS4Cells","ARCHS4Tissues","BioCarta","ChEA","ConsensusTFs","dbGaP","DSigDB","EncodeTF","HumanGeneAtlas","huMAP","JensenDiseases","JensenTissues","KEGG","NCIpathway","OMIM","OMIMexpanded","Panther","Reactome","TFlof","TFpert","TissueHumanProteomeMap","TissueProteomicsDB","WikiPathways")
#' @param anno.identity identity for gene annotations. It can be "GeneID" for Gene ID and "Symbol" for gene symbol
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' an object of class "aOnto", a list with two components (an igraph object 'g' and a list 'anno')
#' @note none
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xDefineOntology.r
#' @examples
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#'
#' \dontrun{
#' aOnto <- xDefineOntology("HPPA", RData.location=RData.location)
#' aOnto <- xDefineOntology("REACTOME_ImmuneSystem", RData.location=RData.location)
#' aOnto <- xDefineOntology("CGL", RData.location=RData.location)
#' }

xDefineOntology <- function(ontology=c(NA,"GOBP","GOMF","GOCC","PSG","PS","PS2","SF","Pfam","DO","HPPA","HPMI","HPCM","HPMA","MP", "EF", "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7", "DGIdb", "GTExV4", "GTExV6p", "GTExV7", "CreedsDisease", "CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP", "CreedsDrugDN", "CreedsGene", "CreedsGeneUP", "CreedsGeneDN", "KEGG","KEGGmetabolism","KEGGgenetic","KEGGenvironmental","KEGGcellular","KEGGorganismal","KEGGdisease", "REACTOME", "REACTOME_ImmuneSystem", "REACTOME_SignalTransduction", "CGL", "SIFTS2GOBP","SIFTS2GOMF","SIFTS2GOCC", "EnrichrARCHS4Cells","EnrichrARCHS4Tissues","EnrichrHumanGeneAtlas","EnrichrTissueHumanProteomeMap","EnrichrTissueProteomicsDB", "EnrichrAchillesFitnessD","EnrichrAchillesFitnessI","EnrichrDSigDB","EnrichrOMIM","EnrichrOMIMexpanded","EnrichrdbGaP", "EnrichrJensenDiseases","EnrichrJensenTissues", "EnrichrBioCarta","EnrichrKEGG","EnrichrNCIpathway","EnrichrPanther","EnrichrReactome","EnrichrWikiPathways","EnrichrhuMAP", "EnrichrChEA","EnrichrConsensusTFs","EnrichrEncodeTF","EnrichrTFlof","EnrichrTFpert"), anno.identity=c("GeneID","Symbol"), verbose=T, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    anno.identity <- match.arg(anno.identity)
    
    ontology <- ontology[1]
    
    g <- NULL
    anno <- NULL
    
    if(!is.na(ontology)){
	
		if(verbose){
			message(sprintf("Load the ontology %s and its gene annotations (%s) ...", ontology, as.character(Sys.time())), appendLF=T)
		}

		#########
		## load GS information
		## flag the simplified version of PS
		flag_PS2 <- FALSE
		if(ontology=="PS2"){
			flag_PS2 <- TRUE
			ontology <- "PS"
		}
		
		## flag the simplified version of REACTOME
		flag_REACTOME <- FALSE
		if(grepl('REACTOME_', ontology)){
			flag_REACTOME <- TRUE
			ontology_REACTOME <- ontology
			ontology <- "REACTOME"
		}
		
		GS <- xRDataLoader(RData.customised=paste('org.Hs.eg', ontology, sep=''), RData.location=RData.location, verbose=verbose)
		
		################
		if(flag_PS2){
			tmp <- as.character(unique(GS$set_info$name))
			inds <- sapply(tmp,function(x) which(GS$set_info$name==x))
		
			## new set_info
			set_info <- data.frame()
			for(i in 1:length(inds)){
				set_info<- rbind(set_info,as.matrix(GS$set_info[max(inds[[i]]),]))
			}
			## new gs
			gs <- list()
			for(i in 1:length(inds)){
				gs[[i]] <- unlist(GS$gs[inds[[i]]], use.names=F)
			}
			names(gs) <- rownames(set_info)
		
			## new GS
			GS$set_info <- set_info
			GS$gs <- gs
		}
		
		if(flag_REACTOME){
			flag <- unlist(strsplit(ontology_REACTOME, '_'))[2]
			if(flag %in% GS$set_info$namespace){
				## new GS
				GS$set_info <- GS$set_info[GS$set_info$namespace==flag, ]
				ind <- match(names(GS$gs), GS$set_info$setID)
				GS$gs <- GS$gs[!is.na(ind)]
			}
		}
		################
		
		#########
		## get annotation information
		anno <- GS$gs

		# replace EntrezGenes with gene symbols	
		if(anno.identity=="Symbol"){
			## load Enterz Gene information
			EG <- xRDataLoader(RData.customised=paste('org.Hs.eg', sep=''), RData.location=RData.location, verbose=verbose)
			allGeneID <- EG$gene_info$GeneID
			allSymbol <- as.vector(EG$gene_info$Symbol)
	
			## anno_symbols
			anno_symbols <- lapply(anno,function(x){
				ind <- match(x, allGeneID)
				allSymbol[ind]
			})
			anno <- anno_symbols	
		}


		#########
		## get ontology information
		## check the eligibility for the ontology
		all.ontologies <- c("GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPCM","HPMA","MP","EF","SIFTS2GOBP","SIFTS2GOMF","SIFTS2GOCC","REACTOME")
		flag_ontology <- ontology %in% all.ontologies
    	
    	if(flag_ontology){
    		#######################################
    		ontology <- gsub('^SIFTS2','',ontology)
    		#######################################
    		
			g <- xRDataLoader(RData.customised=paste('ig.', ontology, sep=''), RData.location=RData.location, verbose=verbose)
			if(is.null(V(g)$term_namespace)){
				V(g)$term_namespace <- ontology
			}
			
			#######################################
			## restricted to those with annotations for the simplified version of REACTOME
			if(flag_REACTOME){
				if(0){
					g <- dnet::dDAGinduce(g, nodes_query=GS$set_info$setID, path.mode="all_paths")
				}else{
					ind <- match(V(g)$name, GS$set_info$setID)
					vids <- V(g)$name[!is.na(ind)]
					neighs.in <- igraph::neighborhood(g, order=vcount(g), nodes=vids, mode="in")
					neighbors <- unique(names(unlist(neighs.in)))
					g <- dnet::dNetInduce(g, nodes_query=neighbors, knn=0, remove.loops=TRUE, largest.comp=T)
					#g <- delete_vertices(g, V(g)$term_name!=unlist(strsplit(ontology_REACTOME, '_'))[2])
				}
			}
			#######################################
						
		}else{
			# force ontology.algorithm to be 'none'
			ontology.algorithm <- 'none'
		
			nodes <- data.frame(name=as.character(GS$set_info$setID), term_id=as.character(GS$set_info$setID), term_name=as.character(GS$set_info$name), term_distance=as.character(GS$set_info$distance), term_namespace=as.character(GS$set_info$namespace), stringsAsFactors=F)
			nodes <- rbind(nodes, c('root','root','root','root','root'))
			relations <- data.frame(from='root', to=nodes$name)
			g <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
		}
	
	}
	
    aOnto <- list(g = g,
    			anno = anno
                 )
    class(aOnto) <- "aOnto"
	
	invisible(aOnto)
}
