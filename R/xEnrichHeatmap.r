#' Function to visualise enrichment results using heatmap
#'
#' \code{xEnrichHeatmap} is supposed to visualise enrichment results using heatmap. It returns an object of class "ggplot".
#'
#' @param list_eTerm an object of class "ls_eTerm". Alterntively, it can be a data frame having all these columns (named as 'group','ontology','name','adjp') and one of these columns ("zscore","fdr","pvalue","fc","or"). Note, the column 'fdr' can be inferred from the column 'adjp'
#' @param fdr.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05
#' @param displayBy which statistics will be used for comparison. It can be "fc" for enrichment fold change (by default), "adjp" for adjusted p value (or FDR), "pvalue" for p value, "zscore" for enrichment z-score, "or" for enrichment odd ratio
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param zlim the minimum and maximum z values for which colors should be plotted, defaulting to the range of the -log10(FDR)
#' @param reorder how to reorder rows and columns. It can be "none" for no reordering, "row" for reordering rows according to number of sharings (by default), "col" for reordering columns, and "both" for reordering rows and columns
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xHeatmap}}
#' @include xEnrichHeatmap.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' # provide the input Genes of interest (eg 100 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' set.seed(825)
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 100))
#' data
#' 
#' # optionally, provide the test background (if not provided, all human genes)
#' #background <- as.character(org.Hs.eg$gene_info$Symbol)
#' 
#' # 2) Gene-based enrichment analysis using ontologies (REACTOME and GOMF)
#' # perform enrichment analysis
#' ls_eTerm <- xEnricherGenesAdv(data, ontologies=c("REACTOME","GOMF"), RData.location=RData.location)
#' ## heatmap plot of enrichment results
#' gp <- xEnrichHeatmap(ls_eTerm, fdr.cutoff=0.1, displayBy="zscore")
#' }

xEnrichHeatmap <- function(list_eTerm, fdr.cutoff=0.05, displayBy=c("zscore","fdr","pvalue","fc","or"), colormap=NULL, zlim=NULL, reorder=c("none","row","col","both"))
{
    
    displayBy <- match.arg(displayBy)
    reorder <- match.arg(reorder)
    
    if(is.null(list_eTerm)){
        warnings("There is no enrichment in the 'eTerm' object.\n")
        return(NULL)
    }
    
    if(class(list_eTerm)=='ls_eTerm' | class(list_eTerm)=='data.frame'){
	
		if(class(list_eTerm)=='ls_eTerm'){
			## when 'auto', will keep the significant terms
			df <- list_eTerm$df
			
		}else if(class(list_eTerm)=='data.frame'){
			## when displayBy='fdr', always append the column 'fdr' if not provided
			if(displayBy=='fdr'){
				if(is.null(list_eTerm$fdr)){
					list_eTerm$fdr <- list_eTerm$adjp
				}
			}
			
			##################
			## always append the column 'id' if not provided
			if(is.null(list_eTerm$id) & !is.null(list_eTerm$name)){
				list_eTerm$id <- list_eTerm$name
			}
			##################
						
			if(all(c('group','ontology','id','name','adjp',displayBy) %in% colnames(list_eTerm))){
				df <- list_eTerm[,c('group','ontology','id','name','adjp',displayBy)]
			
			}else if(all(c('group','id','name','adjp',displayBy) %in% colnames(list_eTerm))){
				df <- list_eTerm[,c('group','id','name','adjp',displayBy)]
				df$ontology <- 'ontology'
			
			}else if(all(c('ontology','id','name','adjp',displayBy) %in% colnames(list_eTerm))){
				df <- list_eTerm[,c('ontology','id','name','adjp',displayBy)]
				df$group <- 'group'
			
			}else if(all(c('id','name','adjp',displayBy) %in% colnames(list_eTerm))){
				df <- list_eTerm[,c('id','name','adjp',displayBy)]
				df$group <- 'group'
				df$ontology <- 'ontology'
			
			}else{
				warnings("The input data.frame does not contain required columns: c('group','ontology','name','adjp').\n")
				return(NULL)
			}
			
		}
		
		df_all <- df
		
		## heatmap view
		if(!is.null(df_all)){

			adjp <- NULL
		
			gp <- NULL
			mat <- NULL
		
			ls_df <- split(x=df_all, f=df_all$ontology)
			#######
			## keep the same order for ontologies as input
			ls_df <- ls_df[unique(df_all$ontology)]
			#######
			ls_mat <- lapply(1:length(ls_df), function(i){
		
				df <- ls_df[[i]]
				ind <- which(df$adjp < fdr.cutoff)
				if(length(ind)>=1){
					df <- as.data.frame(df %>% dplyr::filter(adjp < fdr.cutoff))
				
					if(displayBy=='fdr'){
						mat <- as.matrix(xSparseMatrix(df[,c('name','group','adjp')], rows=unique(df$name), columns=unique(df_all$group)))
						mat[mat==0] <- NA
						mat <- -log10(mat)
					}else if(displayBy=='pvalue'){
						mat <- as.matrix(xSparseMatrix(df[,c('name','group','pvalue')], rows=unique(df$name), columns=unique(df_all$group)))
						mat[mat==0] <- NA
						mat <- -log10(mat)
					}else if(displayBy=='zscore'){
						mat <- as.matrix(xSparseMatrix(df[,c('name','group','zscore')], rows=unique(df$name), columns=unique(df_all$group)))
						mat[mat==0] <- NA
					}else if(displayBy=='fc'){
						mat <- as.matrix(xSparseMatrix(df[,c('name','group','fc')], rows=unique(df$name), columns=unique(df_all$group)))
						mat[mat==0] <- NA
						mat <- log2(mat)
					}else if(displayBy=='or'){
						mat <- as.matrix(xSparseMatrix(df[,c('name','group','or')], rows=unique(df$name), columns=unique(df_all$group)))
						mat[mat==0] <- NA
						mat <- log2(mat)
					}
				
					if(nrow(mat)==1){
						df_mat <- mat
					}else{
					
						## order by the length of names
						rname_ordered <- rownames(mat)[order(-nchar(rownames(mat)))]
						## order by the evolutionary ages
						if(names(ls_df)[i]=='PS2' || names(ls_df)[i]=='PSG'){
							df_tmp <- unique(df[,c('id','name')])
							df_tmp <- df_tmp[with(df_tmp, order(as.numeric(df_tmp$id))),]
							rname_ordered <- df_tmp$name
						}
					
						ind <- match(rname_ordered, rownames(mat))
						df_mat <- as.matrix(mat[ind,], ncol=ncol(mat))
						colnames(df_mat) <- colnames(mat)
					
						colnames(df_mat) <- colnames(mat)
					}
					return(df_mat)
					
				}else{
					return(NULL)
				}
			
			})
			mat <- do.call(rbind, ls_mat)
		
			if(!is.null(mat)){
				if(displayBy=='fdr' | displayBy=='pvalue'){
					if(is.null(colormap)){
						colormap <- 'grey100-darkorange'
					}
					if(is.null(zlim)){
						zlim <- c(0, ceiling(max(mat[!is.na(mat)])))
					}
					
					if(displayBy=='fdr'){
						legend.title <- expression(-log[10]("FDR"))
					}else if(displayBy=='pvalue'){
						legend.title <- expression(-log[10]("p-value"))
					}
				
				}else if(displayBy=='fc' | displayBy=='zscore' | displayBy=='or'){
					tmp_max <- ceiling(max(mat[!is.na(mat)]))
					tmp_min <- floor(min(mat[!is.na(mat)]))
					if(tmp_max>0 & tmp_min<0){
						if(is.null(colormap)){
							colormap <- 'deepskyblue-grey100-darkorange'
						}
						if(is.null(zlim)){
							tmp <- max(tmp_max, abs(tmp_min))
							zlim <- c(-tmp, tmp)
						}
					}else if(tmp_max<=0){
						if(is.null(colormap)){
							colormap <- 'deepskyblue-grey100'
						}
						if(is.null(zlim)){
							zlim <- c(tmp_min, 0)
						}
					}else if(tmp_min>=0){
						if(is.null(colormap)){
							colormap <- 'grey100-darkorange'
						}
						if(is.null(zlim)){
							zlim <- c(0, tmp_max)
						}
					}
				
					if(displayBy=='fc'){
						legend.title <- expression(log[2]("FC"))
					}else if(displayBy=='zscore'){	
						legend.title <- ("Z-score")
					}else if(displayBy=='or'){	
						legend.title <- expression(log[2]("OR"))
					}
				}
		
				gp <- xHeatmap(mat, reorder=reorder, colormap=colormap, ncolors=64, zlim=zlim, legend.title=legend.title, barwidth=0.4, x.rotate=60, shape=19, size=2, x.text.size=6,y.text.size=6, na.color='transparent',barheight=max(3,min(5,nrow(mat))))
				gp <- gp + theme(legend.title=element_text(size=8))
				gp$mat <- mat
			}
		
		}else{
			mat <- NULL
			gp <- NULL
		}
	
	}
	
	invisible(gp)
}

