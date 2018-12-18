#' Function to estimate memory allocated for an R variable or a file
#'
#' \code{xObjSize} is supposed to estimate memory allocated for an R variable or a file.
#'
#' @param obj character specifying an R variable or a local file
#' @param type the object type specifying an R object or a local file. It can be 'variable', 'file' or 'auto' (by default; automatically determined)
#' @param units the units to be used in formatting the size. It can be 'auto', 'Kb', 'Mb', 'Gb', 'Tb', 'Pb'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' If action is logical, a list containing arguments and their default values. If action is NULL, a string specifying the assignment to be evalated.
#' @note
#' This function is potentially useful when debugging as it frees developers from specifying default values for all arguments except those arguments of interest
#' @export
#' @seealso \code{\link{xObjSize}}
#' @include xObjSize.r
#' @examples
#' xObjSize(ls()[1])
#' #res <- lapply(ls(), xObjSize)
#' #res <- lapply(list.files(), xObjSize)

xObjSize <- function(obj, type=c("auto","variable","file"), units="auto", verbose=TRUE)
{
    
    type <- match.arg(type)
    
    res <- NULL
    if(type=='auto'){
    	if(base::file.exists(obj)){
    		type <- "file"
    	}else{
    		type <- "variable"
    	}
    }
    
    if(type=='variable'){
		if(class(suppressWarnings(try(base::get(obj), T)))=="try-error"){
    		warnings(sprintf("The R variable '%s' NOT found!", obj))
        	return(NULL)
		}else{
			res <- base::format(utils::object.size(base::get(obj)), units)
			if(verbose){
				message(sprintf("The R variable '%s' in size: %s", obj, res), appendLF=T)
			}
		}
    }
    
    if(type=='file'){
    	if(base::file.exists(obj)){
    		#res <- utils::format.object_size(file.info(obj)$size, units)
    		
			if(verbose){
				#message(sprintf("The file '%s' in size: %s", obj, res), appendLF=T)
				message(sprintf("The file '%s' in size: %d", obj, file.info(obj)$size), appendLF=T)
			}
    	}else{
    		warnings(sprintf("The file '%s' NOT found!", obj))
        	return(NULL)
    	}
    }
    
    invisible(res)
}
