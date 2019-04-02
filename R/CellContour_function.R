#' Determine Cell Contour
#'
#' This function allows you to determine the cell contour based upon a black and white image 
#' @param array of x and y values of the image
#' @keywords contour
#' @export
#' @examples
#' CellContour()

#' @useDynLib cellshape CellContour_
CellContour <- function(x) {
  .C(CellContour_,as.integer(x),dim(x),as.integer(array(0,c(100000,4))))[[3]]
}
