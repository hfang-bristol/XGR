#' Function to define a colormap
#'
#' \code{xColormap} is supposed to define a colormap. It returns a function, which will take an integer argument specifying how many colors interpolate the given colormap.
#'
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @return 
#' \itemize{
#'  \item{\code{palette.name}: a function that takes an integer argument for generating that number of colors interpolating the given sequence}
#' }
#' @note The input colormap includes: 
#' \itemize{
#' \item{"jet": jet colormap}
#' \item{"bwr": blue-white-red}
#' \item{"gbr": green-black-red}
#' \item{"wyr": white-yellow-red}
#' \item{"br": black-red}
#' \item{"yr": yellow-red}
#' \item{"wb": white-black}
#' \item{"rainbow": rainbow colormap, that is, red-yellow-green-cyan-blue-magenta}
#' \item{"ggplot2": emulating ggplot2 default color palette}
#' \item{Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkblue-lightblue-lightyellow-darkorange", "darkgreen-white-darkviolet", "darkgreen-lightgreen-lightpink-darkred". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}}
#' }
#' @export
#' @include xColormap.r
#' @examples
#' # 1) define "blue-white-red" colormap
#' palette.name <- xColormap(colormap="bwr")
#' # use the return function "palette.name" to generate 10 colors spanning "bwr"
#' palette.name(10)
#' 
#' # 2) define default colormap from ggplot2
#' palette.name <- xColormap(colormap="ggplot2")
#' # use the return function "palette.name" to generate 3 default colors used by ggplot2
#' palette.name(3)

xColormap <- function(colormap=c("bwr","jet","gbr","wyr","br","yr","rainbow","wb","heat","terrain","topo","cm","ggplot2"))
{
    if(length(colormap)>1){
        colormap <- colormap[1]
    }
    
    if(colormap=='ggplot2'){
    
		my_hue_pal <- function (h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, direction = 1){
			function(n) {
				if ((diff(h)%%360) < 1) {
					h[2] <- h[2] - 360/n
				}
				rotate <- function(x)(x + h.start)%%360 * direction
				hues <- rotate(seq(h[1], h[2], length.out = n))
				grDevices::hcl(hues, c, l)
			}
		}
    	palette.name <- my_hue_pal(h=c(0,360)+15, c=100, l=65, h.start=0, direction=1)
    	
	}else if(colormap == "heat"){
		palette.name <- grDevices::heat.colors
	}else if(colormap == "terrain"){
		palette.name <- grDevices::terrain.colors
	}else if(colormap == "topo"){
		palette.name <- grDevices::topo.colors
	}else if(colormap == "cm"){
		palette.name <- grDevices::cm.colors
    	
    }else{
    	palette.name <- supraHex::visColormap(colormap=colormap)
    }

    invisible(palette.name)
}