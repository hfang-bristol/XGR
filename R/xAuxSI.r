#' Function to collect information on the R session
#'
#' \code{xAuxSI} is supposed to collect information on the R session. Only attached 
#'
#' @param si an object of the class "sessionInfo"
#' @param file a file name. By defaul it is NULL
#' @return 
#' If the file is not NULL, it will write into the file. Otherwise it will return a tibble with columns 'kind' ('other' other attached packages, 'loaded' for loaded but not attached packages), 'pkg' (the package name), 'version' (the package version), 'title' (the package title) and 'description' (the package description)
#' @note none.
#' @export
#' @seealso \code{\link{xAuxSI}}
#' @include xAuxSI.r
#' @examples
#' # res <- xAuxSI()

xAuxSI <- function(si=utils::sessionInfo(), file=NULL)
{
	
	i <- pkg <- NULL
	
	## other attached packages
	ls_tl <- lapply(si$otherPkgs, function(x){
		unlist(x) %>% tibble::enframe(name='key',value='value')
	})
	tibble::tibble(i=1:length(ls_tl), kind="other", pkg=names(ls_tl)) %>% 
	dplyr::mutate(
		version=purrr::map_chr(i, ~ls_tl[[.x]] %>% dplyr::filter(key=='Version') %>% dplyr::pull(value)),
		title=purrr::map_chr(i, ~ls_tl[[.x]] %>% dplyr::filter(key=='Title') %>% dplyr::pull(value)),
		description=purrr::map_chr(i, ~ls_tl[[.x]] %>% dplyr::filter(key=='Description') %>% dplyr::pull(value))
	) %>% dplyr::arrange(pkg) %>% dplyr::select(-i) -> res_other
	
	## loaded packages (not attached)
	ls_tl <- lapply(si$loadedOnly, function(x){
		unlist(x) %>% tibble::enframe(name='key',value='value')
	})
	tibble::tibble(i=1:length(ls_tl), kind="loaded", pkg=names(ls_tl)) %>% 
	dplyr::mutate(
		version=purrr::map_chr(i, ~ls_tl[[.x]] %>% dplyr::filter(key=='Version') %>% dplyr::pull(value)),
		title=purrr::map_chr(i, ~ls_tl[[.x]] %>% dplyr::filter(key=='Title') %>% dplyr::pull(value)),
		description=purrr::map_chr(i, ~ls_tl[[.x]] %>% dplyr::filter(key=='Description') %>% dplyr::pull(value))
	) %>% dplyr::arrange(pkg) %>% dplyr::select(-i) -> res_loaded
	
	## bind by row
	res <- dplyr::bind_rows(list(res_other, res_loaded))
	
	## output
	if(!is.null(file)){
		res %>% readr::write_delim(file, delim="\t")
	}else{
		invisible(res)
	}
	
}
