#' Function to assign (and evaluate) arguments with default values for a given function
#'
#' \code{xFunArgs} is supposed to assign (and evaluate) arguments with default values for a given function.
#'
#' @param fun character specifying the name of the function
#' @param action logical to indicate whether the function will act as it should be (with assigned values in the current environment). By default, it sets to NULL, return a string specifying the assignment to be evalated
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' If action is logical, a list containing arguments and their default values. If action is NULL, a string specifying the assignment to be evalated.
#' @note
#' This function is potentially useful when debugging as it frees developers from specifying default values for all arguments except those arguments of interest
#' @export
#' @seealso \code{\link{xFunArgs}}
#' @include xFunArgs.r
#' @examples
#' xFunArgs(fun="xRDataLoader")

xFunArgs <- function(fun, action=NULL, verbose=TRUE)
{
    
    args_list <- base::formals(fun)
    args_names <- names(args_list)
    
    if(is.null(action)){
    	
    	ls_vec <- lapply(1:length(args_list), function(i){
    		lft <- args_names[[i]]
    		rgt <- paste(base::deparse(args_list[[i]]),collapse='')
    		if(rgt!=''){
    			paste0(lft, '=', rgt, '[1];')
    		}else{
    			NULL	
    		}
    	})
    	res <- noquote(paste(unlist(ls_vec), collapse=''))
    	return(res)
    	
    }else{
		for(i in 1:length(args_list)){
			lft <- args_names[[i]]
			rgt <- paste(base::deparse(args_list[[i]]),collapse='')
			if(rgt!=''){
				tmp <- paste(lft, '<-', rgt, sep=' ')
				if(action==T){
					base::eval(base::parse(text=tmp), envir=parent.frame())
				}else if(action==F){
					base::eval(base::parse(text=tmp))
				}
			}
		}
    }
    
    if(verbose){
        if(action==T){
            message(sprintf("For the function '%s', %d arguments have been assigned with default values in the current environment.", fun, length(args_names)), appendLF=T)
        }else{
            message(sprintf("For the function '%s', there are %d arguments.", fun, length(args_names)), appendLF=T)
        }
    }
    
    invisible(args_list)
}
