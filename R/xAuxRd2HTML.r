#' Function to convert Rd files to HTML files
#'
#' \code{xAuxRd2HTML} is supposed to convert Rd files to HTML files.
#'
#' @param path.from a directory containing Rd files converted from
#' @param path.to a directory containing HTML files converted to
#' @return 
#' none
#' @note This auxiliary function helps create a new package.
#' @export
#' @seealso \code{\link{xAuxRd2HTML}}
#' @include xAuxRd2HTML.r
#' @examples
#' # xAuxRd2HTML(path.from="./XGR/man", path.to="./XGR/vignettes")

xAuxRd2HTML <- function(path.from="./XGR/man", path.to="./XGR/vignettes")
{
	for(nm in list.files(path.from, pattern="\\.Rd$")) {
		nm.out <- gsub(pattern="\\.Rd$", "\\.html", nm, ignore.case=T, perl=T)
		
		input_file <- file.path(path.from, nm)
		output_file <- file.path(path.to, nm.out)
		message(sprintf("%s => %s", input_file, output_file), appendLF=T)
		tools::Rd2HTML(Rd=input_file, out=output_file, package="", stages=c("install", "render"))
	}
}
