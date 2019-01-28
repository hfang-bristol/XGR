#' Function to store eQTL summary data in the BESD format
#'
#' \code{xMEbesd} is supposed to store eQTL summary data in the BESD format ready for the use in SMR.
#
#' @param summary a data frame storing eQTL summary statistics. It must contain columns 'context', and columns ('snp_cse','snps','effect_allele','other_allele','effect_maf') for the esi file, and columns ('gene_cse','gene','Symbol') for the epi file, and columns ('snps','gene','beta','statistic','pvalue','FDR') for the besd file
#' @param contexts a vector specifying contexts included
#' @param outdir the output directory. By default it is '.'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' a logical vector
#' @note This function requires the software 'smr' at \url{http://cnsgenomics.com/software/smr}. Shell command lines in Terminal (Mac and Linux) are:
#' \itemize{
#' \item{1a) Mac: \code{wget http://cnsgenomics.com/software/smr/download/smr_Mac.zip && unzip smr_Mac.zip && mv smr_Mac ~/smr}}
#' \item{1b) Linux: \code{wget https://cnsgenomics.com/software/smr/download/smr_Linux.zip && unzip smr_Linux.zip && mv smr_Linux ~/smr}}
#' \item{2a) # Assuming a ROOT (sudo) privilege: \cr\code{sudo cp ~/smr /usr/local/bin}}
#' \item{2b) # Assuming without ROOT (sudo) privilege and adding the system PATH variable to your ~/.bash_profile file: \cr\code{export PATH=$HOME:$PATH}}
#' }
#' @export
#' @seealso \code{\link{xMEbesd}}
#' @include xMEbesd.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' summary <- xRDataLoader('JK_cohort_xMEdb', RData.location=RData.location)
#' xMEbesd(summary, contexts="CD14", outdir="Pi_eQTL_hg19")
#' }

xMEbesd <- function(summary, contexts=c("CD14","LPS2","LPS24","IFN","Bcell","NK","Neutrophil","Monocyte","Blood","CD4","CD8","shared_CD14","shared_IFN","shared_LPS2","shared_LPS24"), outdir=".", verbose=T)
{
	
	## check the folder and create it if not exist
	if(!dir.exists(outdir)){
		dir.create(outdir)
		if(verbose){
			message(sprintf("The folder '%s' created (%s)", outdir, as.character(Sys.time())), appendLF=TRUE)
		}
	}
	
	## df_summary_all
	if(class(summary)=='data.frame'){
		# remove snp with 'chrXY'
		ind <- grepl('^chrXY', summary$snp_cse)
		df_summary_all <- summary[!ind,]
		
		# remove gene with 'chrUn'
		ind <- grepl('^chrUn', df_summary_all$gene_cse)
		df_summary_all <- df_summary_all[!ind,]
		
		if(nrow(df_summary_all)==0){
			return(NULL)
		}
	}else{
		return(NULL)
	}
	
	## keep only contexts specified
	if(is.null(contexts)){
		vec_context <- unique(df_summary_all$context)
	}else{
		vec_context <- contexts
	}
	ind <- match(df_summary_all$context, vec_context)
	df_summary_all <- df_summary_all[!is.na(ind),]
	
	res_vec <- sapply(1:length(vec_context), function(j){
		if(verbose){
			message(sprintf("%d %s (%s) ...", j, vec_context[j], as.character(Sys.time())), appendLF=TRUE)
		}
		
		context <- vec_context[j]
		ind <- match(df_summary_all$context, context)
		df_summary <- df_summary_all[!is.na(ind), ]
		
		####################
		# df_esi -> output_esi
		####################
		# BESD format: an efficient format to store eQTL summary data
		# .esi (SNP information): chromosome, SNP, genetic distance (not required), basepair position, the effect (coded) allele, the other allele, frequency of the effect allele (not essential)
		ls_vec <- strsplit(df_summary$snp_cse, ':|-')
		df <- as.data.frame(do.call(rbind, ls_vec), stringsAsFactors=F)
		colnames(df) <- c('chr','start','end')
		df$chr <- gsub('chr','',df$chr)
		df_esi <- data.frame(chr=df$chr, snp=df_summary$snps, distance=0, pos=df$start, effect=df_summary$effect_allele, other=df_summary$other_allele, maf=df_summary$effect_maf, stringsAsFactors=F)
		output_esi <- file.path(outdir, paste0(context,".esi.txt"))
		utils::write.table(df_esi, file=output_esi, row.names=F, col.names=F, quote=F, sep="\t")
		if(verbose){
			message(sprintf("%s (%s) ...", output_esi, as.character(Sys.time())), appendLF=TRUE)
		}
		
		####################
		# df_epi -> output_epi
		####################
		# .epi (probe information): chromosome, probe, genetic distance (not required), physical position, gene, gene orientation (not essential)
		ls_vec <- strsplit(df_summary$gene_cse, ':')
		df <- as.data.frame(do.call(rbind, ls_vec), stringsAsFactors=F)
		colnames(df) <- c('chr','pos','orientation')
		df$chr <- gsub('chr','',df$chr)
		df$pos <- as.integer(apply(apply(do.call(rbind,strsplit(df$pos,'-')), 1:2, as.integer), 1, mean))
		df_epi <- data.frame(chr=df$chr, probe=df_summary$gene, distance=0, pos=df$pos, gene=df_summary$Symbol, orientation=df$orientation, stringsAsFactors=F)
		output_epi <- file.path(outdir, paste0(context,".epi.txt"))
		utils::write.table(df_epi, file=output_epi, row.names=F, col.names=F, quote=F, sep="\t")
		if(verbose){
			message(sprintf("%s (%s) ...", output_epi, as.character(Sys.time())), appendLF=TRUE)
		}
		
		####################
		# df_me -> output_me
		####################
		# .besd (a binary file to store the summary statistics)
		# store SNPs within 2Mb of a probe in either direction, SNPs within 1Mb of any trans-eQTL in either direction, and SNPs with p < 1e-5 in the rest of the genome
		# store effect size (b) and SE in the BESD file, and re-calculate p-value for analysis when necessary, assuming b / SE follows a standard normal distribution, i.e. N(0, 1). Strictly speaking, b / SE follows a t-distribution which is approximately N(0, 1) if sample size is large. For data sets with small sample sizes (e.g. GTEx), this might lead to a bias in p-value. In this scenario, we recommend users to compute z* based on the original p-value from a standard normal distribution, and adjust the standard error as SE = b / z*. This adjustment guarantees that the re-computed p-value from b and SE being exact the same as the original p-value.
		# Matrix eQTL output: SNP, gene, beta (SNP effect on gene expression), statistic, pvalue, fdr
		df_me <- data.frame(snp=df_summary$snps, probe=df_summary$gene, beta=df_summary$beta, statistic=df_summary$statistic, pvalue=df_summary$pvalue, fdr=df_summary$FDR, stringsAsFactors=F)
		output_me <- file.path(outdir, paste0(context,".me.txt"))
		utils::write.table(df_me, file=output_me, row.names=F, col.names=T, quote=F, sep="\t")
		if(verbose){
			message(sprintf("%s (%s) ...", output_me, as.character(Sys.time())), appendLF=TRUE)
		}
		
		###########################################################
		
		########
		## --out
		out <- file.path(outdir, context)
		########
				
		####################
		# generate .besd file
		####################
		#smr --eqtl-summary CD14.me.txt --matrix-eqtl-format --make-besd --out CD14
		cmd1 <- paste0("smr --eqtl-summary ", output_me, " --matrix-eqtl-format --make-besd --out ", out, " > besd.log")
		cmd2 <- paste0("$HOME/SMR/smr --eqtl-summary ", output_me, " --matrix-eqtl-format --make-besd --out ", out, " > besd.log")
		cmds <- c(cmd1, cmd2)
		besd_flag <- 1
		for(i in 1:length(cmds)){
			invisible(utils::capture.output(cmd <- try(system(cmds[i]), silent=TRUE), type='message'))
			if(cmd==0){
				besd_flag <- 0
				if(verbose){
					message(sprintf("'%s' executed successfully (%s)", cmds[i], as.character(Sys.time())), appendLF=TRUE)
				}
				break
			}
		}
		
		####################
		# update .esi and .epi files
		####################
		if(besd_flag==0){
			vec_files <- paste0(out, c('.besd','.esi','.epi'))
			if(all(file.exists(vec_files))){
				## update .esi file
				cmd1 <- paste0("smr --beqtl-summary ", out, " --update-esi ", output_esi, " > esi.log")
				cmd2 <- paste0("$HOME/SMR/smr --beqtl-summary ", out, " --update-esi ", output_esi, " > esi.log")
				cmds <- c(cmd1, cmd2)
				esi_flag <- 1
				for(i in 1:length(cmds)){
					invisible(utils::capture.output(cmd <- try(system(cmds[i]), silent=TRUE), type='message'))
					if(cmd==0){
						esi_flag <- 0
						if(verbose){
							message(sprintf("'%s' executed successfully (%s)", cmds[i], as.character(Sys.time())), appendLF=TRUE)
						}
						break
					}else{
						if(verbose){
							message(sprintf("'%s' failed (%s)", cmds[i], as.character(Sys.time())), appendLF=TRUE)
						}						
					}
				}

				## update .epi file
				cmd1 <- paste0("smr --beqtl-summary ", out, " --update-epi ", output_epi, " > epi.log")
				cmd2 <- paste0("$HOME/SMR/smr --beqtl-summary ", out, " --update-epi ", output_epi, " > epi.log")
				cmds <- c(cmd1, cmd2)
				epi_flag <- 1
				for(i in 1:length(cmds)){
					invisible(utils::capture.output(cmd <- try(system(cmds[i]), silent=TRUE), type='message'))
					if(cmd==0){
						epi_flag <- 0
						if(verbose){
							message(sprintf("'%s' executed successfully (%s)", cmds[i], as.character(Sys.time())), appendLF=TRUE)
						}					
						break
					}else{
						if(verbose){
							message(sprintf("'%s' failed (%s)", cmds[i], as.character(Sys.time())), appendLF=TRUE)
						}						
					}
				}
				
				## check files generated
				if(all(esi_flag==0, epi_flag==0)){
					if(verbose){
						message(sprintf("Congratulations! files ('%s', '%s' and '%s') created in the directory '%s'", vec_files[1], vec_files[2], vec_files[3], outdir), appendLF=T)
					}
					
					##################################
					## remove *.bak.esi
					unlink(file.path(outdir, paste0(context,'.bak.esi')))
					## remove *.bak.epi
					unlink(file.path(outdir, paste0(context,'.bak.epi')))
					##################################
					if(1){
						## remove *.esi.txt
						unlink(output_esi)
						## remove *.epi.txt
						unlink(output_epi)
						## remove *.me.txt
						unlink(output_me)
					}
					
					if(1){
						## remove esi.log
						unlink('esi.log')
						## remove epi.log
						unlink('epi.log')
						## remove besd.log
						unlink('besd.log')
					}
					
				}else{
					warning("Unfortunately, fail to produce the files. Please install smr first and then put it into the system PATH variable.\n")
					return(NULL)
				}
			
			}
		
		}else{
			warning("Unfortunately, fail to produce the files. Please install smr first and then put it into the system PATH variable.\n")
			return(NULL)
		}
		
		return(all(besd_flag==0, esi_flag==0, epi_flag==0))
		
	})
	names(res_vec) <- vec_context
	
    invisible(res_vec)
}
