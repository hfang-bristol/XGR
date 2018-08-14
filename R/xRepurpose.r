#' Function to obtain repurposing matrix
#'
#' \code{xRepurpose} is supposed to obtain repurposing matrix given a query list of genes. It returns an object of the class 'DR'.
#'
#' @param data an input vector containing gene symbols
#' @param phase.min the minumum phase of drugs allowed. By default it is 3 defining target genes of drugs reaching development phase 3 and above
#' @param target.max the maximum number of targets per drug allowed. By default it is 5. It is used to define non-promoscuous drug target genes
#' @param plot logical to indicate whether heatmap plot is drawn
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param DTT the drug therapeutic targets. It can be "ChEMBL_v24" for the version 24 (by default), and the version 23. Note: you can also load your customised object directly with columns ('target_number','efo_term','phase','pref_name_drug','Symbol')
#' @param restricted the disease areas restricted to. By default it is NULL
#' @param excluded the disease areas that are excluded. By default it is NULL
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param ... additional graphic parameters for xHeatmap
#' @return 
#' an object of class "DR", a list with following components:
#' \itemize{
#'  \item{\code{df}: a data frame of n x 5, where the 5 columns are "Target", "Disease", "Phase", "Drug", "Drug_index"}
#'  \item{\code{index}: a data frame of n x 2, where the 2 columns are "Drug_index", "Drug"}
#'  \item{\code{gp}: NULL if the plot is not drawn; otherwise, a 'ggplot' object}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xHeatmap}}
#' @include xRepurpose.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' # a) provide the input Genes of interest (eg 1000 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' set.seed(825)
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 1000))
#' 
#' # b) obtain repurposing matrix
#' DR <- xRepurpose(data, RData.location=RData.location, reorder="none", colormap="ggplot2.top", zlim=c(1,4), na.color='transparent', label.size=1.5, label.color="white")
#' DR$gp
#' write.table(DR$df, file="xRepurpose.txt", sep="\t", row.names=F, quote=F)
#' write.table(DR$index, file="xRepurpose_index.txt", sep="\t", row.names=F, quote=F)
#' }

xRepurpose <- function(data, phase.min=3, target.max=5, plot=TRUE, verbose=T, DTT=c("ChEMBL_v24","ChEMBL_v23"), restricted=NULL, excluded=NULL, RData.location="http://galahad.well.ox.ac.uk/bigdata", ...)
{
    
	if(class(DTT) == "data.frame"){
			ChEMBL <- DTT
	}else{
		ChEMBL <- xRDataLoader(RData.customised=DTT[1], verbose=verbose, RData.location=RData.location)
		if(is.null(ChEMBL)){
			ChEMBL <- "ChEMBL_v24"
			if(verbose){
				message(sprintf("Instead, %s will be used", ChEMBL), appendLF=T)
			}
			ChEMBL <- xRDataLoader(RData.customised=ChEMBL, verbose=verbose, RData.location=RData.location)
		}
    }
    
	if(!(all(c('target_number','efo_term','phase','pref_name_drug','Symbol') %in% colnames(ChEMBL)))){
		warnings("The input data.frame does not contain required columns: c('target_number','efo_term','phase','pref_name_drug','Symbol').\n")
		return(NULL)
	}
    
	phase.min <- as.integer(phase.min)
	if(phase.min>4 | phase.min<0){
		phase.min <- 3
	}
	target.max <- as.integer(target.max)
	if(target.max<0){
		target.max <- 5
	}
	
	target_number <- efo_term <- phase <- pref_name_drug <- Symbol <- NULL    
	df_chembl <- unique(ChEMBL %>% dplyr::filter(target_number<=target.max) %>% dplyr::select(efo_term,phase,pref_name_drug,Symbol))
	
	################################################################
	### UFTword
	UFTword <- function(x){
	  	#y <- strsplit(x, " ")
	  	y <- x
	  	sapply(y, function(s){
	  		paste(toupper(substring(s,1,1)), substring(s,2), sep="", collapse=" ")
	  	})
	}
	################################################################
	if(class(DTT) != "data.frame"){
    	df_chembl$pref_name_drug <- UFTword(tolower(df_chembl$pref_name_drug))
		df_chembl$efo_term <- UFTword(tolower(df_chembl$efo_term))
	}
	
	####################
	## restrict to data
	ind <- match(df_chembl$Symbol, data)
	df_chembl <- df_chembl[!is.na(ind),]

    if(nrow(df_chembl)==0){
    	return(NULL)
    }
	####################
    
    ## maximum phased drug per target (per disease)
    df <- df_chembl %>% dplyr::filter(phase>=phase.min) %>% dplyr::arrange(-phase)
    #df_target_efo_phase_drug <- df[!duplicated(df[,c("efo_term","Symbol")]),]
    df_target_efo_phase_drug <- df[!duplicated(df[,c("efo_term","Symbol","pref_name_drug")]),]
    
	####################
	if(!is.null(restricted)){
		#restricted <- c("Crohn's disease","Multiple sclerosis","Osteoarthritis","Psoriasis","Rheumatoid arthritis","Systemic lupus erythematosus","Type i diabetes mellitus","Ulcerative colitis")
		ind <- match(df_target_efo_phase_drug$efo_term, restricted)
		df_target_efo_phase_drug <- df_target_efo_phase_drug[!is.na(ind), ]
	}
	if(!is.null(excluded)){
		#excluded <- c("Cancer","Immune system disease")
		ind <- match(df_target_efo_phase_drug$efo_term, excluded)
		df_target_efo_phase_drug <- df_target_efo_phase_drug[is.na(ind), ]
	}
	####################

	###########################
	## df_uid_drug_index
	###########################
	df_target_efo_phase_drug$uid <- paste0(df_target_efo_phase_drug$Symbol,':',df_target_efo_phase_drug$efo_term)

	ls_uid <- split(x=df_target_efo_phase_drug$pref_name_drug, f=df_target_efo_phase_drug$uid)
	ls_uid <- lapply(ls_uid,function(x){
		paste(sort(x),collapse=',')
	})
	vec_uid <- unlist(ls_uid)
	vec_tmp <- sort(unique(vec_uid))
	ind <- match(vec_uid, vec_tmp)
	df_uid_drug_index <- data.frame(uid=names(vec_uid), drug=vec_uid, index=ind, stringsAsFactors=F)
	###########################
    
    ## 1st: data_matrix
	df <- unique(df_target_efo_phase_drug[,c('Symbol','efo_term','phase')])
	data_matrix <- df %>% tidyr::spread(efo_term, phase)
	rownames(data_matrix) <- data_matrix$Symbol
	data_matrix <- data_matrix[,-1]
	gp <- xHeatmap(data_matrix, reorder="none")
	df_res <- gp$data
	ind <- match(df_res$uid, df_uid_drug_index$uid)
	df_res$drug <- df_uid_drug_index$drug[ind]
	df_res$index <- df_uid_drug_index$index[ind]

    #########################################################
    sample <- val <- index <- NULL
    
    ## 2nd: data_matrix
	df_data <- df_res[,c('gene','sample','val')]
	data_matrix <- df_data %>% tidyr::spread(sample, val)
	rownames(data_matrix) <- data_matrix$gene
	data_matrix <- data_matrix[,-1]
	## data_label
	df_label <- df_res[,c('gene','sample','index')]
	data_label <- df_label %>% tidyr::spread(sample, index)
	rownames(data_label) <- data_label$gene
	data_label <- data_label[,-1]
    
    ############
    # by default, ordered by the input data
    ind <- match(data, rownames(data_matrix))
    data_matrix <- data_matrix[ind[!is.na(ind)], ]
    data_label <- data_label[ind[!is.na(ind)], ]
    ############
    
    
    ## heatmap view
    if(plot){
		gp_new <- xHeatmap(data_matrix, data.label=data_label, ...)
    }else{
    	gp_new <- NULL
    }
    
    #########################################################
    ## df_heatmap & df_index
    Target <- Disease <- Drug_index <- NULL
    df_heatmap <- df_res[!is.na(df_res$val), c('gene','sample','val','drug','index')]
    colnames(df_heatmap) <- c('Target','Disease','Phase','Drug','Drug_index')
    df_heatmap <- df_heatmap %>% dplyr::arrange(Target, Disease)
    rownames(df_heatmap) <- NULL
    df_index <- unique(df_heatmap[, c('Drug_index','Drug')]) %>% dplyr::arrange(Drug_index)
    rownames(df_index) <- NULL
    #########################################################
    
	if(verbose){
		if(phase.min==4){
			message(sprintf("%d rows X %d columns for approved drugs", length(unique(df_heatmap$Target)), length(unique(df_heatmap$Disease))), appendLF=TRUE)
		}else{
			message(sprintf("%d rows X %d columns for phase %d and above drugs", length(unique(df_heatmap$Target)), length(unique(df_heatmap$Disease)), phase.min), appendLF=TRUE)
		}
	}
    
    DR <- list(df = df_heatmap,
    		   index = df_index,
    		   gp = gp_new
               )
    class(DR) <- "DR"
    
    invisible(DR)
}
