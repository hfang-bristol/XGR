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
#' @method print mSeed
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
#' @method print cPath
print.cPath <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $ig_paths: an igraph object or NULL"), "\n", sep="")
	cat(sprintf("  $gp_paths: a ggplot object or NULL"), "\n", sep="")
	cat(sprintf("  $gp_heatmap: a ggplot object or NULL"), "\n", sep="")
	cat(sprintf("  $ig_subg: an igraph object or NULL"), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$ig_paths$enrichment:\n")
	print(x$ig_paths$enrichment[1:min(2,nrow(x$ig_paths$enrichment)),c(2:13,18:19)], row.names=FALSE)
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
#' @method print bLD
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
#' @method print DR
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

######################################################################
# cModule
######################################################################
#' @title Definition for S3 class \code{cModule}
#' @description \code{cModule} has 4 components: mem, data, adj, io.
#' @param mem a data frame
#' @param data a matrix
#' @param adj an adjacency matrix
#' @param io a named list
#' @return an object of S3 class \code{cModule}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' cModule(mem, data, adj, io)
#' }
cModule <- function(mem, data, adj, io){
	## integrity checks
	if(class(mem)!='data.frame' | class(data)!='matrix' | class(adj)!='matrix' | class(io)!='list'){
		stop("The S3 class 'cModule' object failed to pass integrity checks!\n")
	}
	value <- list(mem=mem, data=data, adj=adj, io=io)
	class(value) <- "cModule"
	return(value)
}
#' @param x an object of class \code{cModule}
#' @param ... other parameters
#' @rdname cModule
#' @export
#' @method print cModule
print.cModule <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $mem: a data frame of %d rows X %d columns", dim(x$mem)[1],dim(x$mem)[2]), "\n", sep="")
	cat(sprintf("  $data: a data matrix of %d rows X %d columns", dim(x$data)[1],dim(x$data)[2]), "\n", sep="")
	cat(sprintf("  $adj: an adjacency matrix of %d rows X %d columns", dim(x$adj)[1],dim(x$adj)[2]), "\n", sep="")
	cat(sprintf("  $io: a named list with %d components", length(x$io)), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$mem:\n")
	print(head(x$mem,3))
	cat("......\n")
	cat("$io$inputType:\n")
	print(x$io$inputType)
	cat("$io$fit:\n")
	print(head(x$io$fit,3))
	cat("......\n")
	cat("$io$beta:\n")
	print(x$io$beta)
	cat("$io$num_modules:\n")
	print(x$io$num_modules)

}

######################################################################
# pPerf
######################################################################
#' @title Definition for S3 class \code{pPerf}
#' @description \code{pPerf} mush have following components: data, auroc, fmax, amax, direction, gp, Pred_obj.
#' @param data a data frame
#' @param auroc a scalar
#' @param fmax a scalar
#' @param amax a scalar
#' @param direction a character
#' @param gp a 'ggplot' object or NULL
#' @param Pred_obj a ROCR prediction-class object
#' @return an object of S3 class \code{pPerf}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' pPerf(data, auroc, fmax, amax, direction, gp, Pred_obj)
#' }
pPerf <- function(data, auroc, fmax, amax, direction, gp, Pred_obj){
	## integrity checks
	if(class(data)!='data.frame' | class(Pred_obj)!='prediction'){
		stop("The S3 class 'pPerf' object failed to pass integrity checks!\n")
	}
	value <- list(data=data, auroc=auroc, fmax=fmax, amax=amax, direction=direction, gp=gp, Pred_obj=Pred_obj)
	class(value) <- "pPerf"
	return(value)
}
#' @param x an object of class \code{pPerf}
#' @param ... other parameters
#' @rdname pPerf
#' @export
#' @method print pPerf
print.pPerf <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components including:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $data: a data frame of %d rows X %d columns", dim(x$data)[1],dim(x$data)[2]), "\n", sep="")
	cat(sprintf("  $auroc: %.3f", x$auroc), "\n", sep="")
	cat(sprintf("  $fmax: %.3f", x$fmax), "\n", sep="")
	cat(sprintf("  $amax: %.3f", x$amax), "\n", sep="")
	cat(sprintf("  $direction: %s", x$direction), "\n", sep="")
	cat(sprintf("  $Pred_obj: an object of S3 class '%s'", class(x$Pred_obj)), "\n", sep="")
	cat("\n--------------------------------------------------\n")
}

######################################################################
# sClass
######################################################################
#' @title Definition for S3 class \code{sClass}
#' @description \code{sClass} mush have following components: prediction, predictor, performance.
#' @param prediction a data frame
#' @param predictor a data frame
#' @param performance a data frame
#' @return an object of S3 class \code{sClass}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' sClass(prediction, predictor, performance)
#' }
sClass <- function(prediction, predictor, performance){
	## integrity checks
	if(class(prediction)!='data.frame'){
		stop("The S3 class 'sClass' object failed to pass integrity checks!\n")
	}
	value <- list(prediction=prediction, predictor=predictor, performance=performance)
	class(value) <- "sClass"
	return(value)
}
#' @param x an object of class \code{sClass}
#' @param ... other parameters
#' @rdname sClass
#' @export
#' @method print sClass
print.sClass <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', resulted from %d cross-validation, containing %d components:", class(x), length(x$cv_model), length(names(x))), "\n", sep="")
	cat(sprintf("  $prediction: a data frame of %d rows X %d columns", dim(x$prediction)[1], dim(x$prediction)[2], dim(x$prediction)[2]-2), "\n", sep="")
	cat(sprintf("  $predictor: a data frame of %d rows X %d columns (%d predictors)", dim(x$predictor)[1], dim(x$predictor)[2], dim(x$predictor)[2]-2), "\n", sep="")
	cat(sprintf("  $performance: a data frame of %d rows X %d columns", dim(x$performance)[1],dim(x$performance)[2]), "\n", sep="")
	cat(sprintf("  $cv_auroc[1,1]: cross-validated AUC = %.3f", x$cv_auroc[1,1]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$prediction:\n")
	print(x$prediction[1:5,], row.names=FALSE)
	cat("......\n")
	cat("$predictor:\n")
	print(x$predictor[1:2,], row.names=FALSE)
	cat("......\n")
	cat("$performance:\n")
	print(x$performance[1:5,], row.names=TRUE)
	cat("......\n")
	cat("$call:\n")
	print(x$call)	
}

######################################################################
# sMap
######################################################################
#' @title Definition for S3 class \code{sMap}
#' @description \code{sClass} mush have following components: nHex, xdim, ydim, r, lattice, shape, coord, polygon, init, codebook, hits
#' @param nHex an integer
#' @param xdim an integer
#' @param ydim an integer
#' @param r an integer
#' @param lattice a character
#' @param shape a character
#' @param coord a data frame
#' @param polygon a data frame
#' @param init a character
#' @param codebook a data frame
#' @param hits a vector
#' @return an object of S3 class \code{sMap}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' sMap(nHex=nHex, xdim=xdim, ydim=ydim, r=r, lattice=lattice, shape=shape, coord=coord, polygon=polygon, init=init, codebook=codebook, hits=hits)
#' }
sMap <- function(nHex, xdim, ydim, r, lattice, shape, coord, polygon, init, codebook, hits){
	## integrity checks
	if(class(codebook)!='data.frame' | class(coord)!='data.frame' | class(polygon)!='data.frame'){
		stop("The S3 class 'sMap' object failed to pass integrity checks!\n")
	}
	value <- list(nHex=nHex, xdim=xdim, ydim=ydim, r=r, lattice=lattice, shape=shape, coord=coord, polygon=polygon, init=init, codebook=codebook, hits=hits)
	class(value) <- "sMap"
	return(value)
}
#' @param x an object of class \code{sMap}
#' @param ... other parameters
#' @rdname sMap
#' @export
#' @method print sMap
print.sMap <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', containing %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $shape: %s (grid shape)", x$shape), "\n", sep="")
	cat(sprintf("  $lattice: %s (grid lattice)", x$lattice), "\n", sep="")
	cat(sprintf("  $r: %d (grid radius)", x$r), "\n", sep="")
	cat(sprintf("  $xdim: %d hexagons (x-dimension)", x$xdim), "\n", sep="")
	cat(sprintf("  $xdim: %d hexagons (y-dimension)", x$ydim), "\n", sep="")
	cat(sprintf("  $nHex: %d hexagons (in total)", x$nHex), "\n", sep="")
	cat(sprintf("  $hits: a vector of %d in length", length(x$hits)), "\n", sep="")
	cat(sprintf("  $codebook: a data frame of %d rows X %d columns", dim(x$codebook)[1], dim(x$codebook)[2]), "\n", sep="")
	cat(sprintf("  $coord: a data frame of %d rows X %d columns", dim(x$coord)[1], dim(x$coord)[2]), "\n", sep="")
	cat(sprintf("  $polygon: a data frame of %d rows X %d columns", dim(x$polygon)[1], dim(x$polygon)[2]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$codebook:\n")
	print(x$codebook[1:2,1:min(5,ncol(x$codebook))], row.names=FALSE)
	cat("......\n")
	cat("$coord:\n")
	print(x$coord[1:5,], row.names=FALSE)
	cat("......\n")
	cat("$polygon:\n")
	print(x$polygon[1:6,], row.names=TRUE)
	cat("......\n")
	cat("$call:\n")
	print(x$call)	
}

