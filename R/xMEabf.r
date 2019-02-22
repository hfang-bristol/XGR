#' Function to conduct colocalisation analysis through Wakefield's Approximate Bayes Factor approach integrating GWAS and eQTL summary data
#'
#' \code{xMEabf} is supposed to conduct colocalisation analysis integrating GWAS and eQTL summary data through Wakefield's Approximate Bayes Factor (ABF).
#'
#' @param eqtl.summary an input eQTL summary data for a region (eg the eQTLs for a gene), a list with mandatory components 'beta' (a vector for eQTL effect size), 'varbeta' (a vector for beta variance), 'N' (an integer specifying number of samples), 'MAF' (minor allele frequency, eg effect allele frequency), 'snp' (a vector for dbSNP identity)
#' @param gwas.summary an input GWAS summary data, a list with mandatory components 'beta' (a vector for GWAS SNP effect size), 'varbeta' (a vector for beta variance), 'snp' (a vector for dbSNP identity)
#' @param prior.eqtl the prior probability an eQTL associated with the eQTL trait. The default value is 1e-4
#' @param prior.gwas the prior probability an SNP associated with the GWAS trait. The default value is 1e-4
#' @param prior.both the prior probability an eQTL/SNP associated with both eQTL/GWAS traits. The default value is 1e-5
#' @return 
#' a list with two compenents (1) the component 'summary', a vector of 'nsnps' (number of SNPs analysed), 'PP.H0.abf' (posterior probabilities of H0 - no causal variant), 'PP.H1.abf' (posterior probabilities of H1 - causal variant for eQTL trait only), 'PP.H2.abf' (posterior probabilities of H2 - causal variant for GWAS trait only), 'PP.H3.abf' (posterior probabilities of H3 - two distinct causal variants), and 'PP.H4.abf' (posterior probabilities of H4 - one shared causal variant), and (2) the component 'results', a data frame with a column 'snp' (SNPs analysed), columns for eQTL statistics calcualted ('eqtl.V', 'eqtl.z', 'eqtl.r' and 'eqtl.lABF'), columns for GWAS statistics calculated ('gwas.V', 'gwas.z', 'gwas.r' and 'gwas.lABF'), a column 'both.sum.lABF' (the sum of 'eqtl.lABF' and 'gwas.lABF') and a column 'SNP.PP.H4' (the posterior probability of the SNP being causal for both traits).
#' @export
#' @seealso \code{\link{xMEabf}}
#' @include xMEabf.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' 
#' res <- xMEabf(eqtl.summary, gwas.summary)
#' utils::write.table(res$results, file="df_abf.txt", row.names=F, col.names=T, quote=F, sep="\t")
#' }

xMEabf <- function(eqtl.summary, gwas.summary, prior.eqtl=1e-4, prior.gwas=1e-4, prior.both=1e-5)
{

    ###################################################
    # eqtl.summary -> df_eQTL
    ###################################################
    if(all(c('beta','varbeta','N','MAF','snp') %in% names(eqtl.summary))){
    
    	## 'sdY' estimated from maf and varbeta
    	oneover <- 1/eqtl.summary$MAF
    	nvx <- 2 * eqtl.summary$N * eqtl.summary$MAF * (1 - eqtl.summary$MAF)
    	m <- stats::lm(nvx ~ oneover - 1)
    	cf <- stats::coef(m)[['oneover']]
    	
    	if(cf > 0){
    		eqtl.summary$sdY <- sqrt(cf)
    		
    		## Approximate Bayes Factors (ABF) calcuated using variance of the regression coefficients
    		## Wakefield method (PMID:18642345) yields a Bayes facto measuring relative support for a model in which the SNP is associated with the trait compared to the null model of no association
    		### Z statistic
    		z <- eqtl.summary$beta / sqrt(eqtl.summary$varbeta)
    		### variance of the effect estimate
    		V <- eqtl.summary$varbeta
    		### standard deviation of the normal prior
    		sd.prior <- 0.15 * eqtl.summary$sdY
    		### the ratio of the variance of the prior and the total variance
    		r <- sd.prior^2 / (sd.prior^2 + V)
    		### the log-transformed ABF (lABF)
    		lABF <- 0.5 * (log(1 - r) + (r * z^2))
    		
    		## output: df_eQTL
    		df_eQTL <- data.frame(eqtl.V=V, eqtl.z=z, eqtl.r=r, eqtl.lABF=lABF, snp=eqtl.summary$snp, stringsAsFactors=F)
    		
    	}else{
    		return(NULL)
    	}
    	
    }else{
    	return(NULL)
    }
    
    
    ###################################################
    # gwas.summary -> df_gwas
    ###################################################
    if(all(c('beta','varbeta','snp') %in% names(gwas.summary))){
    	
    	## Approximate Bayes Factors (ABF)
    	## Wakefield method (PMID:18642345) yields a Bayes facto measuring relative support for a model in which the SNP is associated with the trait compared to the null model of no association
    	### Z statistic
    	z <- gwas.summary$beta / sqrt(gwas.summary$varbeta)
    	### variance of the effect estimate
    	V <- gwas.summary$varbeta
    	### standard deviation of the normal prior
    	sd.prior <- 0.2
    	### the ratio of the variance of the prior and the total variance
    	r <- sd.prior^2 / (sd.prior^2 + V)
    	### the log-transformed ABF (lABF)
    	lABF <- 0.5 * (log(1 - r) + (r * z^2))
    	
    	## output: df_gwas
    	df_gwas <- data.frame(gwas.V=V, gwas.z=z, gwas.r=r, gwas.lABF=lABF, snp=gwas.summary$snp, stringsAsFactors=F)

    }else{
    	return(NULL)
    }
    
    
    ###################################################
    # Bayesian colocalisation analysis -> merged.df
    ###################################################
    merged.df <- merge(df_eQTL, df_gwas)
    if(nrow(merged.df)==0){
		warnings("fail to merge df_eQTL and df_gwas. Both should contain the common snps")
		return(NULL)
		
	}else{
		merged.df$both.sum.lABF <- merged.df$eqtl.lABF + merged.df$gwas.lABF
		
		####################################
		## logsum function: the log of the sum of the exponentiated logs taking out the max to ensure the sum is not Inf
		## max(x) + log(sum(exp(x - max(x))))
		logsum <- function(x) {
  			max(x) + log(sum(exp(x - max(x))))
		}
		
		## logdiff function: the log of the difference of the exponentiated logs taking out the max to ensure the difference is not negative
		## max(x,y) + log(exp(x - max(x,y)) - exp(y-max(x,y)))
		logdiff <- function(x,y) {
			max(x,y) + log(exp(x - max(x,y)) - exp(y - max(x,y)))
		}
		####################################
		
		## add SNP.PP.H4: posterior probability that each SNP is the causal variant for a shared signal
		x <- merged.df$both.sum.lABF
		my.denom.log.abf <- logsum(x)
		merged.df$SNP.PP.H4 <- exp(x - my.denom.log.abf)
		## sort by SNP.PP.H4
		SNP.PP.H4 <- NULL
		merged.df <- merged.df %>% dplyr::arrange(-SNP.PP.H4)
		
		## calculate posterior probabilities for 5 different configurations, given lABF for each SNP and prior probabilities
		p1 <- prior.eqtl
		p2 <- prior.gwas
		p12 <- prior.both
		l1 <- merged.df$eqtl.lABF
		l2 <- merged.df$gwas.lABF
		lsum <- l1 + l2
		lH0.abf <- 0
		lH1.abf <- log(p1) + logsum(l1)
		lH2.abf <- log(p2) + logsum(l2)
		lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
		lH4.abf <- log(p12) + logsum(lsum)
		all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
		pp.abf <- exp(all.abf - logsum(all.abf))
		names(pp.abf) <- paste0("PP.H", (1:length(pp.abf)) - 1, ".abf")
		
		## output
  		output <- list(summary=c(nsnps=nrow(merged.df),pp.abf), results=merged.df)
  		
		invisible(output)
	}

}
