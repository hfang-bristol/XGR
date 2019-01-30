#' Function to conduct colocalisation analysis through SMR
#'
#' \code{xMEsmr} is supposed to conduct conduct Summary-data-based Mendelian Randomisation (SMR) integrating GWAS and eQTL summary data.
#'
#' @param gwas.summary an input GWAS summary data file, containing columns 'snp', 'effect allele' (the allele assessed), 'other allele', 'freq' (frequency of the effect allele; not essential unless 'freq.check' is true), 'b' (effect size for the allele assessed; log(odds ratio) for a case-control study), 'se' (standard error), 'p' (p-value) adn 'n' (sample size; not required)
#' @param beqtl.summary a character specifying where to find eQTL summary data files in the BESD format containing three files (.esi for SNP information, .epi for probe information, and .besd for eQTL summary statistics)
#' @param bfile a character specifying where to find the LD reference data containing three files (.bed, .bim, and .fam)
#' @param mode a character specifying the SMR and HEIDI test mode. It can be 'cis' for the test in cis regions, 'trans' for the test in trans regions, and 'both' for both regions
#' @param peqtl eQTL p-value threshold for selecting a probe (with the top associated eQTL passing a p-value threshold) for the SMR test. In other words, a probe with the top associated eQTL not passing this threshold will be removed for the test. By default, it is 5e-8
#' @param window.cis an integer specifying a window centred around the probe to select cis-eQTLs (passing a p-value threshold) for the SMR test. The default value is 1000Kb
#' @param window.trans an integer specifying a window centred around the top associated trans-eQTL to select trans-eQTLs (passing a p-value threshold) for the SMR and HEIDI test. The default value is 1000Kb
#' @param heidi logical to indicate whether the HEIDI test is enabled. By default it is true
#' @param heidi.peqtl eQTL p-value threshold for selecting eQTLs per probe for the HEIDI test. The default value is 1.57e-3 (a chi-squared value 10 with df=1)
#' @param heidi.ld LD r2 threshold used to further prune SNPs (eQTLs) in the HEIDI test. By default only those between 0.05 and 0.9 will be used for the test
#' @param heidi.num the number of SNPs (eQTLs) left per probe in the HEIDI test. By default, the test skipped if the number of remaining SNPs is less than 3; and top 20 SNPs (ranked by eQTL p-value) will be used in the test
#' @param freq.check logical to indicate whether to remove SNPs withe discrepant allele frequencies between data sets. By default it is disabled
#' @param thread.num an integer specifying the number of OpenMP threads for parallel computing. By default it is 1
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @return 
#' a data frame with following columns ("mode", "probeID", "Gene", "ProbeChr", "Probe_bp", "topSNP", "topSNP_chr", "topSNP_bp", "A1", "A2", "b_GWAS", "b_eQTL", "b_SMR", "p_GWAS", "p_eQTL", "p_SMR", "fdr_SMR", "p_HEIDI", "fdr_HEIDI", "nsnp_HEIDI")
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
#' df_output <- xMEsmr(gwas.summary, beqtl.summary, bfile)
#' utils::write.table(df_output, file="df_output.txt", row.names=F, col.names=T, quote=F, sep="\t")
#' }

xMEsmr <- function(gwas.summary, beqtl.summary, bfile, mode=c("both","cis","trans"), peqtl=5e-8, window.cis=1000, window.trans=1000, heidi=T, heidi.peqtl=1.57e-3, heidi.ld=c(0.05,0.9), heidi.num=c(3,20), freq.check=F, thread.num=1, p.adjust.method=c("BH","BY","bonferroni","holm","hochberg","hommel"), verbose=T, silent=FALSE)
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
    mode <- match.arg(mode)
    p.adjust.method <- match.arg(p.adjust.method)
    
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
        
    cmd <- paste0("smr --bfile ", bfile, " --gwas-summary ", gwas.summary, " --beqtl-summary ", beqtl.summary, " --peqtl-smr ", peqtl, " --ld-lower-limit ", heidi.ld[1], " --ld-upper-limit ", heidi.ld[2], " --peqtl-heidi ", heidi.peqtl, " --heidi-min-m ", heidi.num[1], " --heidi-max-m ", heidi.num[2], " --thread-num ", thread.num)
    
    if(mode %in% c('cis','both')){
    	cmd_cis <- paste0(cmd, " --cis-wind ", window.cis)
    	if(!heidi){
    		cmd_cis <- paste0(cmd_cis, " --heidi-off ")
    	}
    	if(!freq.check){
    		cmd_cis <- paste0(cmd_cis, " --disable-freq-ck ")
    	}
		cmd_cis <- paste0(cmd_cis, " --out ", my_cis, " > my_cis.log")
		
		if(verbose){
			message(sprintf("cis mode ... (%s)", as.character(Sys.time())), appendLF=TRUE)
			message(sprintf("%s", cmd_cis), appendLF=TRUE)
		}
		flag_cis <- try(system(cmd_cis), silent=TRUE)
		if(flag_cis==1){
			if(verbose){
				message(sprintf("\tFailed (%s)", as.character(Sys.time())), appendLF=TRUE)
			}						
		}

    }
    
    if(mode %in% c('trans','both')){
    	cmd_trans <- paste0(cmd, " --trans --trans-wind ", window.trans)
    	if(!heidi){
    		cmd_trans <- paste0(cmd_trans, " --heidi-off ")
    	}
    	if(!freq.check){
    		cmd_trans <- paste0(cmd_trans, " --disable-freq-ck ")
    	}
		cmd_trans <- paste0(cmd_trans, " --out ", my_trans, " > my_trans.log")
		
		if(verbose){
			message(sprintf("trans mode ... (%s)", as.character(Sys.time())), appendLF=TRUE)
			message(sprintf("%s", cmd_trans), appendLF=TRUE)
		}
		flag_trans <- try(system(cmd_trans), silent=TRUE)
		if(flag_trans==1){
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
    
    ls_mode <- split(x=df_res, f=df_res$mode)
    ls_df <- lapply(ls_mode, function(x){
    	x$fdr_SMR <- stats::p.adjust(x$p_SMR, method=p.adjust.method)
    	x$fdr_HEIDI <- stats::p.adjust(x$p_HEIDI, method=p.adjust.method)
    	return(x)
    })
    df_output <- do.call(rbind, ls_df)
    
    df_output <- df_output[,c('mode','probeID','Gene','ProbeChr','Probe_bp','topSNP','topSNP_chr','topSNP_bp','A1','A2','b_GWAS','b_eQTL','b_SMR','p_GWAS','p_eQTL','p_SMR','fdr_SMR','p_HEIDI','fdr_HEIDI','nsnp_HEIDI')]
    rownames(df_output) <- NULL
    
    p_SMR <- p_GWAS <- NULL
    df_output <- df_output %>% dplyr::arrange(mode, p_SMR, p_GWAS)
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xMEsmr): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(df_output)
}
