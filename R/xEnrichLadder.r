#' Function to visualise enrichment results using ladder-like plot
#'
#' \code{xEnrichLadder} is supposed to visualise enrichment results using ladder-like plot in which rows for terms and columns for its members. The members are sorted first by sharings and then by individual terms. It returns an object of class "ggplot". 
#'
#' @param eTerm an object of class "eTerm"
#' @param sortBy which statistics will be used for sorting and viewing gene sets (terms). It can be "adjp" or "fdr" for adjusted p value (FDR), "pvalue" for p value, "zscore" for enrichment z-score, "fc" for enrichment fold change, "nAnno" for the number of sets (terms), "nOverlap" for the number in overlaps, "or" for the odds ratio, and "none" for ordering according to ID of terms
#' @param top_num the number of the top terms (sorted according to FDR or adjusted p-values). If it is 'auto', only the significant terms (see below FDR.cutoff) will be displayed
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05
#' @param CI.one logical to indicate whether to allow the inclusion of one in CI. By default, it is TURE (allowed)
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param x.rotate the angle to rotate the x tick labelings. By default, it is 60
#' @param x.text.size the text size of the x tick labelings. By default, it is 6
#' @param y.text.size the text size of the y tick labelings. By default, it is 6
#' @param shape the number specifying the shape. By default, it is 19
#' @param size the number specifying the shape size. By default, it is 2
#' @param verbose logical to indicate whether the messages will be displayed in the screen
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xEnrichViewer}}, \code{\link{xHeatmap}}
#' @include xEnrichLadder.r
#' @examples
#' \dontrun{
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#'
#' data(Haploid_regulators)
#' ## only IRF1 positive regulators
#' data <- subset(Haploid_regulators, Phenotype=='IRF1' & MI<0)[,c('Gene')]
#'
#' # 1) KEGGenvironmental
#' eTerm <- xEnricherGenes(data, ontology="KEGGenvironmental", size.range=c(10,2000), min.overlap=5, RData.location=RData.location)
#' gp_ladder <- xEnrichLadder(eTerm)
#' 
#' # 2) PSG
#' eTerm <- xEnricherGenes(data, ontology=c("PSG","Approved","GWAS","CGL")[1], size.range=c(1,20000), min.overlap=0, RData.location=RData.location)
#' gp_ladder <- xEnrichLadder(eTerm, sortBy="none", top_num="auto", FDR.cutoff=1)
#' gp_ladder+ coord_flip()
#' 
#' # 3) save into the file "xEnrichLadder.pdf"
#' mat <- xSparseMatrix(gp_ladder$data)
#' pdf("xEnrichLadder.pdf", width=2+ncol(mat)*0.075, height=2+nrow(mat)*0.1, compress=T)
#' print(gp_ladder)
#' dev.off()
#' }

xEnrichLadder <- function(eTerm, sortBy=c("or","adjp","fdr","pvalue","zscore","fc","nAnno","nOverlap","none"), top_num=10, FDR.cutoff=0.05, CI.one=T, colormap="skyblue-darkblue", x.rotate=60, x.text.size=6, y.text.size=6, shape=19, size=2, verbose=T)
{

    sortBy <- match.arg(sortBy)
    
    gp_heatmap <- NULL
    
    if(class(eTerm)=='eTerm'){
		## when 'auto', will keep the significant terms
		df <- xEnrichViewer(eTerm, top_num="all")
		
		############
		if(!CI.one){
			ind <- which(df$CIl>1 | df$CIu<1)
			df <- df[ind,]
		}
		############
		
		if(FDR.cutoff==1){
			FDR.cutoff <- FDR.cutoff + 0.01
		}
		
		if(top_num=='auto'){
			top_num <- sum(df$adjp<FDR.cutoff)
			if(top_num <= 1){
				top_num <- 10
			}
		}
		df_enrichment <- xEnrichViewer(eTerm, top_num=top_num, sortBy=sortBy, details=T)

		if(!is.null(df_enrichment)){

			df_enrichment$label <- paste0(df_enrichment$name, "\n[OR=", df_enrichment$or, ", P=", df_enrichment$pvalue, ", FDR=", df_enrichment$adjp, ", n=", df_enrichment$nOverlap, "/", df_enrichment$nAnno, "]")
	
			## list of individual paths
			ls_path <- lapply(1:nrow(df_enrichment), function(j){
				x <- df_enrichment$members[j]
				query <- unlist(strsplit(x, ", "))
			})
			names(ls_path) <- df_enrichment$name
			
			## only works if there are no less than 2 terms
			if(length(ls_path)>=1){
				
				all_genes <- unique(unlist(ls_path))
				
				## matrix of genes X paths
				ls_vec <- lapply(1:length(ls_path), function(j){
					ind <- match(all_genes, ls_path[[j]])
					vec <- rep(NA, length(all_genes))
					vec[!is.na(ind)] <- names(ls_path)[j]
					vec
				})
				df_res <- do.call(cbind, ls_vec)
				colnames(df_res) <- names(ls_path)
				rownames(df_res) <- all_genes
				vec_sum <- apply(!is.na(df_res), 1, sum)
				### construct data.frame 'df_tmp'
				df_tmp <- data.frame(num=vec_sum, df_res, gene=rownames(df_res), stringsAsFactors=F)
				df_tmp <- df_tmp %>% dplyr::arrange_all()
				#######
				## reverse
				df_tmp <- df_tmp[nrow(df_tmp):1, ]
				#######	
				colnames(df_tmp)[2:(ncol(df_tmp)-1)] <- colnames(df_res)
				### update 'vec_sum' and 'df_res'
				vec_sum <- df_tmp$num
				names(vec_sum) <- df_tmp$gene
				df_res <- df_tmp[,2:(ncol(df_tmp)-1)]
				rownames(df_res) <- df_tmp$gene

				if(verbose){
					message(sprintf("heatmap of %d rows X %d columns (%s) ...", nrow(df_res), ncol(df_res), as.character(Sys.time())), appendLF=TRUE)
				}
				
				###############
				## visualisation
				###############
				if(1){
					df_heatmap <- 0 + !is.na(df_res)
					df_heatmap[df_heatmap==0] <- NA
					for(i in 1:nrow(df_heatmap)){
						x <- df_heatmap[i,]
						df_heatmap[i,!is.na(x)] <- vec_sum[i]
					}
					mat_heatmap <- t(df_heatmap)
					ind <- match(rownames(mat_heatmap), df_enrichment$name)
					rownames(mat_heatmap) <- df_enrichment$label[ind]
					gp_heatmap <- xHeatmap(mat_heatmap, reorder="none", colormap=colormap, zlim=c(0,max(mat_heatmap,na.rm=T)), ncolors=64, barwidth=0.4, x.rotate=x.rotate, x.text.size=x.text.size, y.text.size=y.text.size, shape=shape, size=size, na.color='transparent')
					gp_heatmap <- gp_heatmap + theme(legend.title=element_text(size=8), legend.position="none") + scale_y_discrete(position="right")
					colsep <- cumsum(table(vec_sum))
					colsep <- length(vec_sum) - colsep[-length(colsep)]
					gp_heatmap <- gp_heatmap + geom_vline(xintercept=colsep+0.5,color="grey90",size=0.5)
					#gp_heatmap
				}
			}
		}
	}

    return(gp_heatmap)
}
