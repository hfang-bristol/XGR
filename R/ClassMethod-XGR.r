######################################################################
# eTerm
######################################################################
#' @title Definition for S3 class \code{eTerm}
#' @description \code{eTerm} mush have following components: term_info, annotation, g, data, background, overlap, fc, zscore, pvalue, adjp, cross.
#' @param term_info a data frame
#' @param annotation a list
#' @param g an 'igraph' object
#' @param data a vector
#' @param background a vector
#' @param overlap a vector
#' @param fc a vector
#' @param zscore a vector
#' @param pvalue a vector
#' @param adjp a vector
#' @param cross a matrix
#' @return an object of S3 class \code{eTerm}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' eTerm(term_info, annotation, g, data, background, overlap, fc, zscore, pvalue, adjp, cross)
#' }
eTerm <- function(term_info, annotation, g, data, background, overlap, fc, zscore, pvalue, adjp, cross){
	## integrity checks
	if(class(term_info)!='data.frame' | class(g)!='igraph'){
		stop("The S3 class 'eTerm' object failed to pass integrity checks!\n")
	}
	value <- list(term_info=term_info, annotation=annotation, g=g, data=data, background=background, overlap=overlap, fc=fc, zscore=zscore, pvalue=pvalue, adjp=adjp, cross=cross)
	class(value) <- "eTerm"
	return(value)
}
#' @param x an object of class \code{eTerm}
#' @param ... other parameters
#' @rdname eTerm
#' @export
print.eTerm <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components including:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $term_info: a data frame of %d rows X %d columns", dim(x$term_info)[1],dim(x$term_info)[2]), "\n", sep="")
	cat(sprintf("  $data: a vector (%d in total)", length(x$data)), "\n", sep="")
	cat(sprintf("  $background: a vector (%d in total)", length(x$background)), "\n", sep="")
	cat(sprintf("  $adjp: a vector (%d in total)", length(x$adjp)), "\n", sep="")
	cat(sprintf("  $cross: a matrix of %d X %d", dim(x$cross)[1], dim(x$cross)[2]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("xEnrichViewer(eTerm):\n")
	print(xEnrichViewer(x), row.names=TRUE)
	cat("......\n")
}

######################################################################
# mSeed
######################################################################
#' @title Definition for S3 class \code{mSeed}
#' @description \code{cTarget} has 2 components: GR and Gene.
#' @param GR a data frame
#' @param Gene a data frame
#' @return an object of S3 class \code{mSeed}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' mSeed(priority, predictor)
#' }
mSeed <- function(GR, Gene){
	## integrity checks
	if(class(GR)!='data.frame' | class(Gene)!='data.frame'){
		stop("The S3 class 'mSeed' object failed to pass integrity checks!\n")
	}
	value <- list(GR=GR, Gene=Gene)
	class(value) <- "mSeed"
	return(value)
}
#' @param x an object of class \code{mSeed}
#' @param ... other parameters
#' @rdname mSeed
#' @export
print.mSeed <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $GR: a data frame of %d rows X %d columns", dim(x$GR)[1],dim(x$GR)[2]), "\n", sep="")
	cat(sprintf("  $Gene: a data frame of %d rows X %d columns", dim(x$Gene)[1],dim(x$Gene)[2]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$GR:\n")
	print(x$GR[1:2,], row.names=FALSE)
	cat("......\n")
	cat("$Gene:\n")
	print(x$Gene[1:2,], row.names=FALSE)
	cat("......\n")
}
