#' Function to create a data frame (with three columns) from a (sparse) matrix
#'
#' \code{xSM2DF} is supposed to create a data frame (with three columns) from a (sparse) matrix. Only nonzero entries from the matrix will be kept in the resulting data frame.
#'
#' @param data a matrix or an object of the dgCMatrix class (a sparse matrix)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return a data frame containing three columns: 1st column for row names, 2nd for column names, and 3rd for numeric values
#' @note none
#' None
#' @export
#' @seealso \code{\link{xSM2DF}}
#' @include xSM2DF.r
#' @examples
#' # create a sparse matrix of 4 X 2
#' input.file <- rbind(c('R1','C1',1), c('R2','C1',1), c('R2','C2',1), c('R3','C2',2), c('R4','C1',1))
#' data <- xSparseMatrix(input.file)
#' # convert into a data frame a full matrix
#' res_df <- xSM2DF(data)
#' res_df

xSM2DF <- function(data, verbose=TRUE)
{
    
    if(class(data) == 'dgCMatrix' | is.data.frame(data)){
    	#data <- as.matrix(data)
	}
    
    ## row names
    if(is.null(rownames(data))){
    	names_row <- 1:nrow(data)
    }else{
    	names_row <- rownames(data)
    }
    
    ## column names
    if(is.null(colnames(data))){
    	names_col <- 1:ncol(data)
    }else{
    	names_col <- colnames(data)
    }
	
	if(is.data.frame(data) | class(data) == 'matrix'){
		data <- as.matrix(data)
		
		## values
		xy <- which(data!=0, arr.ind=TRUE)
		ind <- which(data!=0, arr.ind=FALSE)

		if(nrow(xy) > 0){
			res_df <- data.frame(rownames=names_row[xy[,1]], colnames=names_col[xy[,2]], values=data[ind], stringsAsFactors=FALSE)
		
			res_df <- res_df[order(res_df$rownames,res_df$values,decreasing=FALSE),]
		
			if(verbose){
				message(sprintf("A matrix of %d X %d has been converted into a data frame of %d X 3.", dim(data)[1], dim(data)[2], nrow(res_df)), appendLF=T)
			}
		
		}else{
			res_df <- NULL
		}

	}else if(class(data) == 'dgCMatrix'){
		ijx <- summary(data)
		if(nrow(ijx)>0){
			res_df <- data.frame(rownames=names_row[ijx[,1]], colnames=names_col[ijx[,2]], values=ijx[,3], stringsAsFactors=FALSE)
		
			res_df <- res_df[order(res_df$rownames,res_df$values,decreasing=FALSE),]
		
			if(verbose){
				message(sprintf("A matrix of %d X %d has been converted into a data frame of %d X 3.", dim(data)[1], dim(data)[2], nrow(res_df)), appendLF=T)
			}
		}else{
			res_df <- NULL
		}
		
	}
    
    invisible(res_df)
}
