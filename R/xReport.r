#' Function to generate a graphml file from a pathway upon query
#'
#' \code{xReport} is supposed to generate a graphml file from a pathway upon query. If data is provided, pathway gene members are color-coded.
#'
#' @param obj an R object. Usually a S3-class object storing results such as an 'eTerm' object 
#' @param rmd the R markdown file. If NULL, the pre-prepared one in the directory 'inst/DynamicReport' of the XGR package will be used
#' @param output_format the output format rendered from the R markdown file. If NULL, the output format is the first one defined within the R markdown file. The advanced use is to pass an output format object via rmarkdown::html_document()
#' @param output_file the name of the output file. If NULL, the output filename will be based on the filename of R markdown file (extension replaced)
#' @param output_dir the directory of the output file
#' @param quiet the logic specifying whether to suppress printing of the pandoc command line
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param ... additional parameters used in rmarkdown::render
#' @return
#' the message on the rendered output file and directory.
#' @note none
#' @export
#' @seealso \code{\link{xReport}}
#' @include xReport.r
#' @examples
#' # Load the XGR package and specify the location of built-in data
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' \donttest{
#' res <- xReport(eTerm)
#' 
#' # advanced use
#' output_format <- rmarkdown::html_document(number_sections=T,theme="journal", hightlight="espresso",code_folding="hide")
#' res <- xReport(eTerm, output_format=output_format)
#' }

xReport <- function(obj, rmd=NULL, output_format=NULL, output_file=NULL, output_dir=NULL, quiet=T, verbose=T, ...)
{
	
	res <- NULL
	
	if(is.null(output_dir)){
		output_dir <- getwd()
	}
	
    if(class(obj)=='eTerm'){
    	if(is.null(rmd)){
    		rmd <- system.file("DynamicReport", "eTerm.Rmd", package="XGR")
    	}
    	
    	if(file.exists(rmd)){
			rmarkdown::render(rmd, output_format=output_format, output_file=output_file, output_dir=output_dir, quiet=quiet, ...)

			if(is.null(output_file)){
				output_file <- gsub("\\.rmd$","\\.html", basename(rmd), ignore.case=T)
			}
		
			res <- sprintf("Congratulations! A html file '%s' created in the directory '%s'", output_file, output_dir)
		
			if(verbose){
				message(res, appendLF=TRUE)
			}
		}
		
	}
    
    invisible(res)
}






