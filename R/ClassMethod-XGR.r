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
#' @description \code{cTarget} has 3 components: GR, Gene, Link.
#' @param GR a data frame
#' @param Gene a data frame
#' @param Link a data frame
#' @return an object of S3 class \code{mSeed}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' mSeed(GR, Gene, Link)
#' }
mSeed <- function(GR, Gene, Link){
	## integrity checks
	if(class(GR)!='data.frame' | class(Gene)!='data.frame' | class(Link)!='data.frame'){
		stop("The S3 class 'mSeed' object failed to pass integrity checks!\n")
	}
	value <- list(GR=GR, Gene=Gene, Link=Link)
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
	cat(sprintf("  $Link: a data frame of %d rows X %d columns", dim(x$Link)[1],dim(x$Link)[2]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$GR:\n")
	print(x$GR[1:min(2,nrow(x$GR)),], row.names=FALSE)
	cat("......\n")
	cat("$Gene:\n")
	print(x$Gene[1:min(2,nrow(x$Gene)),], row.names=FALSE)
	cat("......\n")
	cat("$Link:\n")
	print(x$Link[1:min(2,nrow(x$Link)),], row.names=FALSE)
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
# cPath
######################################################################
#' @title Definition for S3 class \code{cPath}
#' @description \code{cPath} has 4 components: ig_paths, gp_paths, gp_heatmap, ig_subg.
#' @param ig_paths an igraph object
#' @param gp_paths a ggplot object
#' @param gp_heatmap a ggplot object
#' @param ig_subg an igraph object
#' @return an object of S3 class \code{cPath}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' cPath(ig_paths, gp_paths, gp_heatmap, ig_subg)
#' }
cPath <- function(ig_paths, gp_paths, gp_heatmap, ig_subg){
	## integrity checks
	if(class(ig_paths)!='igraph' | all(class(gp_paths) %in% c('ggplot','gg')) | all(class(gp_heatmap) %in% c('ggplot','gg')) | class(ig_subg)!='igraph'){
		stop("The S3 class 'cPath' object failed to pass integrity checks!\n")
	}
	value <- list(ig_paths=ig_paths, gp_paths=gp_paths, gp_heatmap=gp_heatmap, ig_subg=ig_subg)
	class(value) <- "cPath"
	return(value)
}
#' @param x an object of class \code{cPath}
#' @param ... other parameters
#' @rdname cPath
#' @export
print.cPath <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $ig_paths: an igraph object or NULL"), "\n", sep="")
	cat(sprintf("  $gp_paths: a ggplot object or NULL"), "\n", sep="")
	cat(sprintf("  $gp_heatmap: a ggplot object or NULL"), "\n", sep="")
	cat(sprintf("  $ig_subg: ab igraph object or NULL"), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$ig_paths$enrichment:\n")
	print(x$ig_paths$enrichment[1:min(2,nrow(x$ig_paths$enrichment)),2:13], row.names=FALSE)
	cat("......\n")
}

######################################################################
# bLD
######################################################################
#' @title Definition for S3 class \code{bLD}
#' @description \code{bLD} has 2 components: best, block.
#' @param best a GR object
#' @param block a GRL object
#' @return an object of S3 class \code{bLD}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' bLD(best, block)
#' }
bLD <- function(best, block){
	## integrity checks
	if(class(best)!='GRanges' | class(block)!='GRangesList'){
		stop("The S3 class 'bLD' object failed to pass integrity checks!\n")
	}
	value <- list(best=best, block=block)
	class(value) <- "bLD"
	return(value)
}
#' @param x an object of class \code{bLD}
#' @param ... other parameters
#' @rdname bLD
#' @export
print.bLD <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $best: a GR object or NULL"), "\n", sep="")
	cat(sprintf("  $block: a GRL object or NULL"), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$best:\n")
	print(head(x$best,1))
	cat("......\n")
	cat("$block:\n")
	print(x$block)
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
# DR
######################################################################
#' @title Definition for S3 class \code{DR}
#' @description \code{DR} has 3 components: df, index, gp.
#' @param df a data frame
#' @param index a data frame
#' @param gp a ggplot object
#' @return an object of S3 class \code{DR}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' DR(df, index, gp)
#' }
DR <- function(df, index, gp){
	## integrity checks
	if(class(df)!='data.frame' | class(index)!='data.frame' | any(class(gp)!='ggplot')){
		stop("The S3 class 'DR' object failed to pass integrity checks!\n")
	}
	value <- list(df=df, index=index, gp=gp)
	class(value) <- "DR"
	return(value)
}
#' @param x an object of class \code{DR}
#' @param ... other parameters
#' @rdname DR
#' @export
print.DR <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $df: a data frame of %d rows X %d columns", dim(x$df)[1],dim(x$df)[2]), "\n", sep="")
	cat(sprintf("  $index: a data frame of %d rows X %d columns", dim(x$index)[1],dim(x$index)[2]), "\n", sep="")
	cat(sprintf("  $gp: a ggplot object"), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$df:\n")
	print(head(x$df,2))
	cat("......\n")
	cat("$index:\n")
	print(head(x$index,2))
	cat("......\n")
}