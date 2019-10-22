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
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
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
#' @method print eTerm
print.eTerm <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components including:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $term_info: a data frame of %d rows X %d columns", dim(x$term_info)[1],dim(x$term_info)[2]), "\n", sep="")
	cat(sprintf("  $data: a vector (%d in total)", length(x$data)), "\n", sep="")
	cat(sprintf("  $background: a vector (%d in total)", length(x$background)), "\n", sep="")
	cat(sprintf("  $adjp: a vector (%d in total)", length(x$adjp)), "\n", sep="")
	cat(sprintf("  $cross: a matrix of %d X %d", dim(x$cross)[1], dim(x$cross)[2]), "\n", sep="")
	cat(sprintf("  $g: an 'igraph' object"), "\n", sep="")
	cat(sprintf("  $g$ontology: '%s'", x$g$ontology), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("xEnrichViewer(eTerm):\n")
	print(xEnrichViewer(x), row.names=TRUE)
	cat("......\n")
}

######################################################################
# ls_eTerm
######################################################################
#' @title Definition for S3 class \code{ls_eTerm}
#' @description \code{ls_eTerm} has 3 components: df, mat and gp.
#' @param df a data frame
#' @param mat a matrix
#' @param gp a ggplot object
#' @return an object of S3 class \code{ls_eTerm}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' ls_eTerm(df, mat, gp)
#' }
ls_eTerm <- function(df, mat, gp){
	## integrity checks
	if(class(df)!='data.frame' | class(mat)!='matrix' | all(class(gp) %in% c('ggplot','gg'))){
		stop("The S3 class 'ls_eTerm' object failed to pass integrity checks!\n")
	}
	value <- list(df=df, mat=mat, gp=gp)
	class(value) <- "ls_eTerm"
	return(value)
}
#' @param x an object of class \code{ls_eTerm}
#' @param ... other parameters
#' @rdname ls_eTerm
#' @export
#' @method print ls_eTerm
print.ls_eTerm <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $df: a data frame of %d rows X %d columns", dim(x$df)[1],dim(x$df)[2]), "\n", sep="")
	cat(sprintf("  $mat: a data matrix of %d rows X %d columns", dim(x$mat)[1],dim(x$mat)[2]), "\n", sep="")
	cat(sprintf("  $gp: a ggplot object"), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$df:\n")
	print(x$df[1:min(2,nrow(x$df)),1:13], row.names=FALSE)
	cat("......\n")
}

######################################################################
# aOnto
######################################################################
#' @title Definition for S3 class \code{aOnto}
#' @description \code{aOnto} has 2 components: g, anno.
#' @param g an igraph object
#' @param anno a list
#' @return an object of S3 class \code{aOnto}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' aOnto(g, anno)
#' }
aOnto <- function(g, anno){
	## integrity checks
	if(class(g)!='igraph' | class(anno)!='list'){
		stop("The S3 class 'aOnto' object failed to pass integrity checks!\n")
	}
	value <- list(g=g, anno=anno)
	class(value) <- "aOnto"
	return(value)
}
#' @param x an object of class \code{aOnto}
#' @param ... other parameters
#' @rdname aOnto
#' @export
#' @method print aOnto
print.aOnto <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $g: an igraph object or NULL"), "\n", sep="")
	cat(sprintf("  $anno: a list with %d length or NULL", length(x$anno)), "\n", sep="")
	
	cat("\n--------------------------------------------------\n")
	cat("$g:\n")
	print(x$g)
	cat("......\n")
}

######################################################################
# iSubg
######################################################################
#' @title Definition for S3 class \code{iSubg}
#' @description \code{iSubg} has 2 components: g, ls_subg.
#' @param g an igraph object
#' @param ls_subg a list of igraph objects
#' @return an object of S3 class \code{iSubg}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' iSubg(g, ls_subg)
#' }
iSubg <- function(g, ls_subg){
	## integrity checks
	if(class(g)!='igraph' | class(ls_subg)!='list'){
		stop("The S3 class 'iSubg' object failed to pass integrity checks!\n")
	}
	value <- list(g=g, ls_subg=ls_subg)
	class(value) <- "iSubg"
	return(value)
}
#' @param x an object of class \code{iSubg}
#' @param ... other parameters
#' @rdname iSubg
#' @export
#' @method print iSubg
print.iSubg <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $g: an igraph object or NULL"), "\n", sep="")
	cat(sprintf("  $ls_subg: a list with %d length or NULL", length(x$ls_subg)), "\n", sep="")
	
	cat("\n--------------------------------------------------\n")
	cat("$g:\n")
	print(x$g)
	cat("......\n")
}
