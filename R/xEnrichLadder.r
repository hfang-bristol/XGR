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
#' @param label how to label gene sets (terms). It can be "concise" or "full"
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
#' 
#' # 4) SIFTS2GOMF
#' ## df_fpocket
#' SIFTS_fpocket <- xRDataLoader(RData='SIFTS_fpocket',RData.location=RData.location)
#' df_fpocket <- as.data.frame(SIFTS_fpocket %>% dplyr::filter(druggable=='Y') %>% dplyr::group_by(Symbol,PDB_code) %>% dplyr::summarise(num_pockets=n()) %>% dplyr::arrange(Symbol,desc(num_pockets),PDB_code))
#' df_fpocket <- df_fpocket[!duplicated(df_fpocket$Symbol), ]
#' ## mat_fpocket
#' mat_fpocket <- df_fpocket %>% tidyr::spread(Symbol, num_pockets)
#' rownames(mat_fpocket) <- mat_fpocket[,1]
#' mat_fpocket <- mat_fpocket[,-1]
#' ## gp_ladder
#' set.seed(825)
#' data <- as.character(sample(unique(df_fpocket$Symbol), 100))
#' eTerm <- xEnricherGenes(data=data, ontology="SIFTS2GOMF", RData.location=RData.location)
#' gp_ladder <- xEnrichLadder(eTerm, sortBy="none", top_num=5, FDR.cutoff=0.01, x.rotate=90)
#' #gp_ladder + coord_flip()
#' ## data_matrix
#' ind <- match(colnames(gp_ladder$matrix), colnames(mat_fpocket))
#' data_matrix <- mat_fpocket[,ind[!is.na(ind)]]
#' ind <- which(apply(!is.na(data_matrix), 1, sum)!=0)
#' data_matrix <- data_matrix[ind,]
#' ind <- match(data, colnames(data_matrix))
#' data_matrix <- data_matrix[,ind[!is.na(ind)]]
#' ## gp_pdb
#' gp_pdb <- xHeatmap(t(data_matrix), reorder="row", colormap="jet.top", x.rotate=90, shape=19, size=1, x.text.size=6,y.text.size=5, na.color='transparent', legend.title='# pockets')
#' #gp_pdb + coord_flip()
#' ## plot_combined
#' #plot_combined <- cowplot::plot_grid(gp_ladder, gp_pdb, align="h", ncol=1, rel_heights=c(2,3))
#' 
#' ## enrichment analysis
#' SIFTS_fpocket <- xRDataLoader(RData='SIFTS_fpocket',RData.location=RData.location)
#' annotation.file <- SIFTS_fpocket[!duplicated(SIFTS_fpocket$Symbol), c('Symbol','druggable')]
#' ### 100 randomly chosen human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg', RData.location=RData.location)
#' set.seed(825)
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 100))
#' ### optionally, provide the test background (if not provided, all human genes)
#' background <- as.character(org.Hs.eg$gene_info$Symbol)
#' ### perform enrichment analysis
#' eTerm <- xEnricherYours(data.file=data, annotation.file=annotation.file, background.file=background, size.range=c(10,20000))
#' }

xEnrichLadder <- function(eTerm, sortBy=c("or","adjp","fdr","pvalue","zscore","fc","nAnno","nOverlap","none"), top_num=10, FDR.cutoff=0.05, CI.one=T, colormap="skyblue-darkblue", x.rotate=60, x.text.size=6, y.text.size=6, shape=19, size=2, label=c('concise','full'), verbose=T)
{

    sortBy <- match.arg(sortBy)
    label <- match.arg(label)
    
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
	
	}else if(class(eTerm)=='data.frame'){
		df_enrichment <- eTerm
	}
	
		if(!is.null(df_enrichment)){
			
			if(label=='concise'){
				#df_enrichment$label <- paste0(df_enrichment$name, " [OR=", df_enrichment$or, ", FDR=", df_enrichment$adjp, ", n=", df_enrichment$nOverlap, "]")
				df_enrichment$label <- paste0(df_enrichment$name, " [OR=", df_enrichment$or, ", FDR=", df_enrichment$adjp, ", n=", df_enrichment$nOverlap, "/", df_enrichment$nAnno, "]")
			}else{
				df_enrichment$label <- paste0(df_enrichment$name, "\n[OR=", df_enrichment$or, ", P=", df_enrichment$pvalue, ", FDR=", df_enrichment$adjp, ", n=", df_enrichment$nOverlap, "/", df_enrichment$nAnno, "]")
			}
			
			###############
			## remove those rows with equal name
			df_enrichment <- df_enrichment[!duplicated(df_enrichment$name),]
			###############
						
			## list of individual paths
			ls_path <- lapply(1:nrow(df_enrichment), function(j){
				x <- df_enrichment$members[j]
				query <- unlist(strsplit(x, ", "))
			})
			names(ls_path) <- df_enrichment$name
			
			## works if 1 term or more
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
				df_res <- df_tmp %>% dplyr::select(2:(ncol(df_tmp)-1))
				#df_res <- df_tmp[,2:(ncol(df_tmp)-1)]
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
					
					## append 'matrix'
					gene <- sample <- val <- NULL
					data_matrix <- gp_heatmap$data %>% dplyr::select(gene,sample,val) %>% tidyr::spread(sample, val)
					rownames(data_matrix) <- data_matrix$gene
					data_matrix <- data_matrix[,-1]
					gp_heatmap$matrix <- data_matrix
					#gp_heatmap
				}
			}
		}

    return(gp_heatmap)
}
