#' Function to visualise distance kernel functions
#'
#' \code{xVisKernels} is supposed to visualise distance kernels, each of which is a decaying function of: i) the relative distance \eqn{d_{gs}} between the gene \eqn{g} and the SNP \eqn{s}, and ii) the decay exponent \eqn{\lambda}.
#'
#' @param exponent an integer specifying decay exponent. By default, it sets to 2
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @return 
#' invisible
#' @note There are five kernels that are currently supported:
#' \itemize{
#' \item{For "slow decay" kernel, \eqn{h_{ds}(t)={1-{d_{gs}/D}\lambda}*(d_{gs} \le D)}}
#' \item{For "linear decay" kernel, \eqn{h_{ds}(t)={1-d_{gs}/D}*(d_{gs} \le D)}}
#' \item{For "rapid decay" kernel, \eqn{h_{ds}(t)={{1-d_{gs}/D}^\lambda}*(d_{gs} \le D)}}
#' }
#' @export
#' @seealso \code{\link{xVisKernels}}
#' @include xVisKernels.r
#' @examples
#' # visualise distance kernels
#' xVisKernels(exponent=2)
#' xVisKernels(exponent=3)

xVisKernels <-function(exponent=2, newpage=TRUE) 
{

    if (newpage){
        dev.new(width=7, height=7)
    }
    
    cl <- rainbow(3)
    ph <- c(21,22,17)

    x <- seq(0,1,by=0.01)
    
    ## slow decay
    y1 <- 1-(x)^exponent
    ## linear decay
    y2<- 1-x
    ## rapid decay
    y3 <- (1-x)^exponent
    
    plot(0,0, xlim = c(0,1),ylim = c(0,1), type = "n",
        xlab=expression(paste("Relative distance of ", d[gs], " (within window D) between gene g and SNP s", collapse=" ")), 
        ylab=expression(paste("Distance kernel function ", h[gs], collapse=" ")),
        main=bquote("Decay exponent is" ~ lambda == .(exponent))
    )
    lines(x, y1, type = "b", pch=ph[1], col = cl[1])
    lines(x, y2, type = "b", pch=ph[2], col = cl[2])
    lines(x, y3, type = "b", pch=ph[3], col = cl[3])
    
    leg.txt <- c(
    	expression(paste("slow decay:", h[gs], "=", (1-(frac(d[gs],D)) ^ lambda)*(d[gs]<=D), collapse=" ")),
    	expression(paste("linear decay:", h[gs], "=", (1-frac(d[gs],D))*(d[gs]<=D), collapse=" ")), 
        expression(paste("rapid decay: ", , h[gs], "=", ((1-frac(d[gs],D)) ^ lambda)*(d[gs]<=D), collapse=" "))
    )
            
    legend("topright", legend=leg.txt, pch=ph, col=cl, border="transparent", box.col="transparent", cex=0.6)
    
    invisible()
}
