#' Function to encode a file as a base64 string for embedding
#'
#' \code{xAuxEmbed} is supposed to encode a file as a base64 string for embedding such as into the R markdown rendering html file. It returns a hyperlink.
#'
#' @param file a file to encode
#' @param download.attribute the download attribute specifying the target to be downloaded instead of being explored. By default, it is the filename of the input file. The filename of the downloaded file can be different from the input file if provided differently
#' @param link.text the link text (the visible part of the hyperlink)
#' @return 
#' a hyperlink in the form of: <a href="encoded base64 string" download="download.attribute">"link.text"</a>
#' @note This auxiliary function helps embed a file into the R markdown rendering html file for the download.
#' @export
#' @seealso \code{\link{xAuxEmbed}}
#' @include xAuxEmbed.r
#' @examples
#' # file <- system.file("DESCRIPTION",package="XGR")
#' # res <- xAuxEmbed(file)

xAuxEmbed <- function(file, download.attribute=basename(file), link.text=paste0("Download ",download.attribute))
{
	
	download.attribute <- basename(file)
	
	ext <- tools::file_ext(file)
	
	if(ext %in% c("txt")){
		type <- "text/plain"
		
	}else if(ext %in% c("csv")){
		type <- "text/csv"
		
	}else if(ext %in% c("html")){
		type <- "text/html"
		
	}else if(ext %in% c("pdf")){
		type <- "application/pdf"
		
	}else if(ext %in% c("Rmd","rmd")){
		type <- "text/x-markdown"
		
	}else{
		# by default
		type <- "text/plain"
		
	}
	
	# file encoded as a base64 string
	encoded <- sprintf('data:%s;base64,%s', type, base64enc::base64encode(file))
	
	res <- sprintf("<a href='%s' download='%s'>%s</a>", encoded, download.attribute, link.text)
	
	invisible(res)
}
