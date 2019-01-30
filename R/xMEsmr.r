#' Function to conduct SMR analysis
#'
#' \code{xMEsmr} is supposed to conduct conduct SMR analysis.
#'
#' @param data.file an input data file, containing a list of entities (e.g. genes or SNPs) to test. The entities can be anything, for example, in this file \url{http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt}, the entities are InterPro domains (InterPro). As seen in this example, entries in the first column must be domains. If the file also contains other columns, these additional columns will be ignored. Alternatively, the data.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param annotation.file an input annotation file containing annotations between entities and ontology terms. For example, a file containing annotations between InterPro domains and GO Molecular Function (GOMF) terms can be found in \url{http://dcgor.r-forge.r-project.org/data/InterPro/Domain2GOMF.txt}. As seen in this example, the input file must contain two columns: 1st column for domains, 2nd column for ontology terms. If there are additional columns, these columns will be ignored. Alternatively, the annotation.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param background.file an input background file containing a list of entities as the test background. The file format is the same as 'data.file'. By default, it is NULL meaning all annotatable entities (i.g. those entities in 'annotation.file') are used as background
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 2000
#' @param min.overlap the minimum number of overlaps. Only those terms with members that overlap with input data at least min.overlap (3 by default) will be processed
#' @param test the test statistic used. It can be "fisher" for using fisher's exact test, "hypergeo" for using hypergeometric test, or "binomial" for using binomial test. Fisher's exact test is to test the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement). Hypergeometric test is to sample at random (without replacement) from the background containing annotated and non-annotated genes, and thus compare sampling to background. Unlike hypergeometric test, binomial test is to sample at random (with replacement) from the background with the constant probability. In terms of the ease of finding the significance, they are in order: hypergeometric test > fisher's exact test > binomial test. In other words, in terms of the calculated p-value, hypergeometric test < fisher's exact test < binomial test
#' @param background.annotatable.only logical to indicate whether the background is further restricted to the annotatable (covered by 'annotation.file'). By default, it is NULL: if the background not provided, it will be TRUE; otherwise FALSE. Surely, it can be explicitly stated. Notably, if only one annotation is provided in 'annotation.file', it should be false (also the background.file should be provided)
#' @param p.tail the tail used to calculate p-values. It can be either "two-tails" for the significance based on two-tails (ie both over- and under-overrepresentation)  or "one-tail" (by default) for the significance based on one tail (ie only over-representation)
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @return 
#' an object of class "eTerm", a list with following components:
#' \itemize{
#'  \item{\code{term_info}: a matrix of nTerm X 4 containing snp/gene set information, where nTerm is the number of terms, and the 4 columns are "id" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
#'  \item{\code{annotation}: a list of terms containing annotations, each term storing its annotations. Always, terms are identified by "id"}
#'  \item{\code{g}: an igraph object to represent DAG}
#'  \item{\code{data}: a vector containing input data in consideration. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{background}: a vector containing the background data. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{overlap}: a list of overlapped snp/gene sets, each storing snps overlapped between a snp/gene set and the given input data (i.e. the snps of interest). Always, gene sets are identified by "id"}
#'  \item{\code{fc}: a vector containing fold changes}
#'  \item{\code{zscore}: a vector containing z-scores}
#'  \item{\code{pvalue}: a vector containing p-values}
#'  \item{\code{adjp}: a vector containing adjusted p-values. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{or}: a vector containing odds ratio}
#'  \item{\code{CIl}: a vector containing lower bound confidence interval for the odds ratio}
#'  \item{\code{CIu}: a vector containing upper bound confidence interval for the odds ratio}
#'  \item{\code{cross}: a matrix of nTerm X nTerm, with an on-diagnal cell for the overlapped-members observed in an individaul term, and off-diagnal cell for the overlapped-members shared betwene two terms}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note None
#' @export
#' @seealso \code{\link{xMEsmr}}
#' @include xMEsmr.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' 
#' bfile <- "~/Sites/SVN/github/bigdata_dev/Merged_EUR"
#' beqtl.summary <- "~/Sites/SVN/github/bigdata_dev/Pi_eQTL_hg19/Neutrophil"
#' gwas.summary <- "summary_gwas.txt"
#' df_res <- xMEsmr(gwas.summary, beqtl.summary, bfile)
#' }

xMEsmr <- function(gwas.summary, beqtl.summary, bfile, test=c("both","cis","trans"), window.cis=1000, window.trans=1000, peqtl=5e-8, ld=c(0.05,0.9), heidi=T, heidi.peqtl=1.57e-3, heidi.num=c(3,20), freq.check=F, thread.num=1, verbose=T, silent=FALSE)
{
    startT <- Sys.time()
    if(!silent){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }else{
    	verbose <- FALSE
    }
    
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    test <- match.arg(test)
    
    ############
    ## bfile
    vec <- paste0(bfile, c('.bed','.fam','.bim'))
    if(any(!file.exists(vec))){
		if(verbose){
			message(sprintf("The bfile '%s' not found (%s)!", bfile, as.character(Sys.time())), appendLF=T)
		}
		return(NULL)
    }
    
    ############
    ## beqtl.summary
    vec <- paste0(beqtl.summary, c('.besd','.esi','.epi'))
    if(any(!file.exists(vec))){
		if(verbose){
			message(sprintf("The beqtl.summary '%s' not found (%s)!", beqtl.summary, as.character(Sys.time())), appendLF=T)
		}
		return(NULL)
    }
    
    ############
    # random numbers derived from time stamps
	tempnum <- gsub('\\.*', '', as.character(as.numeric(Sys.time())))
	my_cis <- paste0(tempnum,'_cis')
	my_trans <- paste0(tempnum,'_trans')
    ############
        
    cmd <- paste0("smr --bfile ", bfile, " --gwas-summary ", gwas.summary, " --beqtl-summary ", beqtl.summary, " --peqtl-smr ", peqtl, " --ld-lower-limit ", ld[1], " --ld-upper-limit ", ld[2], " --peqtl-heidi ", heidi.peqtl, " --heidi-min-m ", heidi.num[1], " --heidi-max-m ", heidi.num[2], " --thread-num ", thread.num)
    
    if(test %in% c('cis','both')){
    	cmd_cis <- paste0(cmd, " --cis-wind ", window.cis)
    	if(!heidi){
    		cmd_cis <- paste0(cmd_cis, " --heidi-off ")
    	}
    	if(!freq.check){
    		cmd_cis <- paste0(cmd_cis, " --disable-freq-ck ")
    	}
		cmd_cis <- paste0(cmd_cis, " --out ", my_cis, " > my_cis.log")
		
		if(verbose){
			message(sprintf("Analysing (cis) ... (%s)", as.character(Sys.time())), appendLF=T)
		}
		flag_cis <- try(system(cmd_cis), silent=TRUE)
		if(flag_cis==0){
			if(verbose){
				message(sprintf("'%s' executed successfully (%s)", cmd_cis, as.character(Sys.time())), appendLF=TRUE)
			}
		}else{
			if(verbose){
				message(sprintf("'%s' failed (%s)", cmd_cis, as.character(Sys.time())), appendLF=TRUE)
			}						
		}

    }
    
    if(test %in% c('trans','both')){
    	cmd_trans <- paste0(cmd, " --trans --trans-wind ", window.trans)
    	if(!heidi){
    		cmd_trans <- paste0(cmd_trans, " --heidi-off ")
    	}
    	if(!freq.check){
    		cmd_trans <- paste0(cmd_trans, " --disable-freq-ck ")
    	}
		cmd_trans <- paste0(cmd_trans, " --out ", my_trans, " > my_trans.log")
		
		if(verbose){
			message(sprintf("Analysing (trans) ... (%s)", as.character(Sys.time())), appendLF=T)
		}
		flag_trans <- try(system(cmd_trans), silent=TRUE)
		if(flag_trans==0){
			if(verbose){
				message(sprintf("'%s' executed successfully (%s)", cmd_trans, as.character(Sys.time())), appendLF=TRUE)
			}
		}else{
			if(verbose){
				message(sprintf("'%s' failed (%s)", cmd_trans, as.character(Sys.time())), appendLF=TRUE)
			}						
		}
    	
    }
    
    ###################
    df_cis <- NULL
    file_cis <- paste0(my_cis,".smr")
    if(file.exists(file_cis)){
    	df_cis <- utils::read.delim(file=file_cis, header=T, row.names=NULL, stringsAsFactors=F)
    	df_cis$mode <- 'cis'
    }
    df_trans <- NULL
    file_trans <- paste0(my_trans,".smr")
    if(file.exists(file_trans)){
    	df_trans <- utils::read.delim(file=file_trans, header=T, row.names=NULL, stringsAsFactors=F)
    	df_trans <- df_trans[,c(-5,-6,-7)]
    	df_trans$mode <- 'trans'
    }
    df_res <- do.call(rbind, list(df_cis, df_trans))
    ###################
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xMEsmr): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(df_res)
}
