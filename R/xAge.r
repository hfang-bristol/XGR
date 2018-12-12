#' Function to compare evolutionary history age of genes based on their Phylostratigraphy
#'
#' \code{xAge} is supposed to compare evolutionary history age of genes based on their Phylostratigraphy. 
#
#' @param list_vec a named vector containing human genes (gene symbol). Alternatively it can be a list of named vectors
#' @param p.type the type of the p-value calcualted. It can be 'nominal' for nominal p-value (calculated using two-sample Wilcoxon rank sum tests) or 'empirical' for empirical p-value
#' @param seed an integer specifying the seed
#' @param nperm the number of random permutations
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param plot logical to indicate whether violin plot is drawn
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a list with three componets:
#' \itemize{
#'  \item{\code{df_summary}: a data frame of n x 7, where n is the number of named vectors, and the 5 columns are "name", "num" (i.e. number of genes used for calculation), "observed" (i.e. average gene age for input/observed genes), "other" (i.e. average gene age for the other genes), "ratio" (i.e. the observed divided by the other),  "pval" (i.e. p-value), "fdr" (if multiple vectors provided)}
#'  \item{\code{ls_gp_violin}: NULL if the plot is not drawn; otherwise, a list of 'ggplot' objects for violet plot}
#'  \item{\code{ls_gp_pdf}: NULL if the plot is not drawn; otherwise, a list of 'ggplot' objects for pdf plot for null distribution of age together with a vertical line for observed age}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xAge}}
#' @include xAge.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' # ExAC LoF intolerant
#' ExAC <- xRDataLoader('ExAC', RData.location=RData.location)
#' data <- subset(ExAC, pLI>0.9)$gene
#' res <- xAge(data, plot=TRUE, RData.location=RData.location)
#' res$ls_gp_violin[[1]]
#'
#' # Genetic regulators 
#' Haploid <- xRDataLoader('Haploid_regulators', RData.location=RData.location)
#' # positive regulators
#' df <- subset(Haploid, MI<0 & Phenotype!='HMGCR')
#' ls_group <- split(x=df$Gene, f=df$Phenotype)
#' ls_res <- xAge(ls_group, plot=TRUE, RData.location=RData.location)
#' ls_res$df_summary
#' gridExtra::grid.arrange(grobs=ls_res$ls_gp_violin[1:4], ncol=2)
#' 
#' ####################################################################
#' ## advanced further analysis comparing multiple groups
#' ls_df <- lapply(ls_res$ls_gp_violin, function(x) x$data)
#' res_df <- do.call(rbind, ls_df)
#' 
#' #################
#' res_df$data <- factor(res_df$data, levels=c("observed","other"))
#' # 1) facet by group
#' p <- ggpubr::ggviolin(res_df, x="data", y="Phylostrata", color="transparent", fill="group", facet.by="group", short.panel.labs=T, add="boxplot", add.params=list(width=0.05, fill='transparent', color="black")) + theme(legend.position='none') + ggpubr::stat_compare_means(label=c("p.format","p.signif")[1],label.y=18,label.x=1.3)
#' # 2) all group
#' yintercept <- as.data.frame(res_df %>% dplyr::group_by(data) %>% dplyr::summarise(median(Phylostrata)))[,2]
#' p <- ggpubr::ggviolin(res_df, x="data", y="Phylostrata", color="transparent", fill="data", add="boxplot", add.params=list(width=0.02, fill='transparent', color="black")) + theme(legend.position='none') + ggpubr::stat_compare_means(aes(group=data),label=c("p.format","p.signif")[1],label.y=18,label.x=1.4) + geom_hline(yintercept=yintercept, linetype=2)
#' 
#' #################
#' # 3) Pairwise comparison against all
#' df <- res_df %>% dplyr::filter(data!="other")
#' ##ggpubr::compare_means(Phylostrata ~ group,  data=df, ref.group=".all.", method="wilcox.test", p.adjust.method="BH")
#' p <- ggpubr::ggviolin(df, x="group", y="Phylostrata", color="transparent", fill="group", add="boxplot", add.params=list(width=0.05, fill='transparent', color="black")) + theme(legend.position='none') + ggpubr::stat_compare_means(method="anova", label.y=18) + ggpubr::stat_compare_means(label="p.signif", label.y=18, method="wilcox.test", ref.group=".all.", hide.ns=T) + geom_hline(yintercept=median(df$Phylostrata), linetype=2)
#' 
#' #################
#' # 4) Pairwise comparison against each other
#' df <- res_df %>% dplyr::filter(data!="other")
#' df_pair <- ggpubr::compare_means(Phylostrata ~ group,  data=df, method="wilcox.test", p.adjust.method="BH")
#' df_pair <- as.data.frame(df_pair %>% dplyr::filter(p.adj < 0.05) %>% dplyr::select(group1,group2))
#' my_comparisons <- lapply(1:nrow(df_pair), function(i) paste0(df_pair[i,]))
#' ggpubr::ggviolin(df, x="group", y="Phylostrata", color="transparent", fill="group", add="boxplot", add.params=list(width=0.05, fill='transparent', color="black")) + theme(legend.position='none') + ggpubr::stat_compare_means(label="p.format", method="wilcox.test", comparisons=my_comparisons) + geom_hline(yintercept=median(df$Phylostrata), linetype=2)
#' 
#' #################
#' # 5) Multiple pairwise tests against a reference group
#' df1 <- res_df %>% dplyr::filter(data!="other")
#' df2 <- res_df %>% dplyr::filter(data=="other") %>% dplyr::mutate(group=data)
#' ## genes not found in any group
#' df2 <- as.data.frame(df2 %>% dplyr::group_by(Symbol) %>% dplyr::group_by(num=n(),add=T) %>% dplyr::filter(num==12))
#' df2 <- df2 %>% dplyr::select(-num)
#' ## remove duplicated genes
#' df2 <- df2[!duplicated(df2$Symbol),]
#' df <- rbind(df1, df2)
#' ##ggpubr::compare_means(Phylostrata ~ group,  data=df, ref.group="other", method="wilcox.test", p.adjust.method="BH")
#' p <- ggpubr::ggviolin(df, x="group", y="Phylostrata", fill="group", color="transparent", add="boxplot", add.params=list(width=0.05, color="black")) + theme(legend.position='none') + ggpubr::stat_compare_means(method="anova", label.y=18) + ggpubr::stat_compare_means(label="p.signif", label.y=17, method="wilcox.test", ref.group="other", hide.ns=T) + geom_hline(yintercept=median(df$Phylostrata[df$group=='other']), linetype=2) + geom_hline(yintercept=median(df$Phylostrata[df$group!='other']), linetype=2)
#' 
#' ### 6) labels
#' df_PSG <- xRDataLoader('Phylostratigraphy', RData.location=RData.location)
#' df_tmp <- unique(df_PSG[,c("Phylostrata","Tax_name")])
#' lables <- paste0(df_tmp$Phylostrata, ' (', df_tmp$Tax_name, ')')
#' p + scale_y_continuous(breaks=0:18, limits=c(0,18), labels=c('',lables,rep('',2))) + theme(axis.text.y=element_text(size=7,color="black"), axis.title.x=element_blank()) 
#' 
#' }

xAge <- function(list_vec, p.type=c("nominal","empirical"), seed=825, nperm=2000, p.adjust.method=c("BH","BY","bonferroni","holm","hochberg","hommel"), plot=FALSE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    p.type <- match.arg(p.type)
    p.adjust.method <- match.arg(p.adjust.method)
    
    ############
    if(length(list_vec)==0){
    	return(NULL)
    }
    ############
    if(is.vector(list_vec) & class(list_vec)!="list"){
    	# assume a vector
		list_vec <- list(list_vec)
	}else if(class(list_vec)=="list"){
		## Remove null elements in a list
		list_vec <- base::Filter(base::Negate(is.null), list_vec)
		if(length(list_vec)==0){
			return(NULL)
		}
    }else{
        stop("The input data must be a named vector or a list of named vectors.\n")
    }
    
	list_names <- names(list_vec)
	if(is.null(list_names)){
		list_names <- paste0('V', 1:length(list_vec))
		names(list_vec) <- list_names
	}
    
    df_PSG <- xRDataLoader('Phylostratigraphy', RData.location=RData.location)
    
    ### labels
    df_tmp <- unique(df_PSG[,c("Phylostrata","Tax_name")])
    lables <- paste0(df_tmp$Phylostrata, ' (', df_tmp$Tax_name, ')')
    ### 
    
    ls_df_gp <- lapply(1:length(list_vec), function(i){
		
		data <- list_vec[[i]]
		
		ind <- match(df_PSG$Symbol, data)
		x <- df_PSG$Phylostrata[!is.na(ind)]
		y <- df_PSG$Phylostrata[is.na(ind)]
		
		res <- stats::wilcox.test(x, y, alternative=c("two.sided","less")[1])
		obs <- signif(mean(x), 2)
		other <- signif(mean(y), 2)
		
		if(p.type == 'nominal'){
			pval_obs <- res$p.value
			
			if(pval_obs < 0.05){
				pval_obs <- as.numeric(format(signif(res$p.value,2),scientific=TRUE))
			}else{
				pval_obs <- signif(res$p.value,3)
			}
			
			gp_pdf <- NULL
			
		}else if(p.type == 'empirical'){
			xy <- c(x,y)
			
			B <- nperm
			set.seed(seed)
			vec_exp <- sapply(1:B, function(i){
				exp <- mean(sample(xy, length(x)))
			})
			other <- signif(mean(vec_exp), 2)
			
			if(mean(vec_exp) > obs){
				## the observed is less than the expected
				pval_obs <- sum(vec_exp <= obs)/B
			}else{
				## the observed is greater than the expected
				pval_obs <- sum(vec_exp >= obs)/B			
			}
			
			if(pval_obs < 0.05){
				pval_obs <- as.numeric(format(signif(res$p.value,2),scientific=TRUE))
			}else{
				pval_obs <- signif(res$p.value,3)
			}
			
			##################################
			## gp_pdf
			age <- NULL
			gp <- ggplot(data.frame(age=vec_exp), aes(age))
			gp <- gp + geom_density(adjust=0.8,fill='grey',colour='grey',alpha=0.5)
			gp <- gp + geom_vline(xintercept=obs,color='red')

			gp <- gp + scale_x_continuous(breaks=1:16, limits=c(1,16), labels=lables)
			gp <- gp + theme_bw() + theme(axis.title.y=element_text(size=10,color="black"), axis.text.y=element_text(size=8,color="black"), axis.title.x=element_text(size=10,color="black"), axis.text.x=element_text(size=8,color="black"), panel.background=element_rect(fill="transparent"))
			gp <- gp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_flip()
			gp_pdf <- gp + labs(x="Phylostratigraphy", y=paste0("Probability density (null distribution)\nestimated based on ",B," permutations"), title=names(list_vec)[i], subtitle=paste0("Observed age: ",obs,' (empirical p: ',pval_obs,')')) + theme(plot.title=element_text(hjust=0.5, size=10), plot.subtitle=element_text(hjust=0.5, size=8))
			##################################
		}

		df <- data.frame(name=names(list_vec)[i], num=length(x), observed=obs, other=other, ratio=obs/other, pval=pval_obs, stringsAsFactors=FALSE)
		
		list(df=df, gp=gp_pdf)
    })
    
    ## ls_gp_pdf
    ls_gp_pdf <- lapply(ls_df_gp, function(x){
    	if(!is.null(x)){
    		x$gp
    	}
    })
    names(ls_gp_pdf) <- names(list_vec)
    ## dff
    ls_df <- lapply(ls_df_gp, function(x){
    	if(!is.null(x)){
    		x$df
    	}
    })
    dff <- do.call(rbind, ls_df)
    
    ###########
    if(is.null(dff)){
    	return(NULL)
    }
    ###########
    if(nrow(dff)>1){
		fdr <- stats::p.adjust(dff$pval, method=p.adjust.method)
		dff$fdr <- ifelse(fdr<0.05, as.numeric(format(signif(fdr,2),scientific=TRUE)), signif(fdr,3))
    }
    rownames(dff) <- 1:nrow(dff)
    
    if(plot){
		ls_gp_violin <- lapply(1:length(list_vec), function(i){
		
			data <- list_vec[[i]]
			
			df_PSG_data <- df_PSG
			df_PSG_data$data <- 'other'
			ind <- match(df_PSG_data$Symbol, data)
			df_PSG_data$data[!is.na(ind)] <- 'observed'
			
			## add the column 'group'
			df_PSG_data$group <- names(list_vec)[i]
			
			ind <- match(names(list_vec)[i], dff$name)
			if(is.na(ind)){
				return(NULL)
			}else{
				name_obs <- dff$name[ind]
				pval_obs <- dff$pval[ind]
				if(!is.null(dff$fdr)){
					fdr_obs <- dff$fdr[ind]
				}else{
					fdr_obs <- NULL
				}
				
			}
			
			Phylostrata <- NULL
			m <- ggplot(df_PSG_data, aes(x=data, y=Phylostrata))
			m <- m + scale_y_continuous(breaks=1:16, limits=c(1,16), labels=lables)
			m <- m + geom_violin(trim=F, fill='transparent', color="#56B4E9") + stat_summary(fun.y=mean, geom="point", shape=23, size=2) + geom_boxplot(width=0.05, outlier.shape=NA, fill='transparent', color="#E69F00")

			m <- m + theme_bw() + theme(legend.position="top", axis.title.y=element_text(size=10,color="black"), axis.text.y=element_text(size=8,color="black"), axis.title.x=element_blank(), axis.text.x=element_text(size=8,color="black"))
			m <- m + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
			
			if(is.null(fdr_obs)){
				subtitle <- paste0(p.type,' p: ',pval_obs)
			}else{
				subtitle <- paste0(p.type,' p: ',pval_obs,', fdr: ',fdr_obs)
			}
			
			gp_violin <- m + labs(x="data frame", title=paste0(name_obs, " (n=",sum(df_PSG_data$data=='observed'),")"), subtitle=subtitle) + theme(plot.title=element_text(hjust=0.5, size=10), plot.subtitle=element_text(hjust=0.5, size=8))
			
			return(gp_violin)
    	})
    	names(ls_gp_violin) <- names(list_vec)
    }else{
    	ls_gp_violin <- NULL
    }
    
    ls_res <- list(df_summary = dff,
    			  ls_gp_violin = ls_gp_violin,
    			  ls_gp_pdf	= ls_gp_pdf
                 )
    
    invisible(ls_res)
}
