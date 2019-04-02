#' Determine Cell Contour
#'
#' This function allows you to determine the cell contour based upon a black and white image array.
#' The function returns an array containing the number of contour points; x-position; y-position; total area.
#' @param x a 2D array indicating all pixels assigned to the cell shape.
#' 
#' @keywords contour
#' @export
#' @examples
#' img <- readPNG("cell.png")
#' tmp<-array(CellContour(as.matrix(img)),dim=c(100000,4))
#' contour<-array(tmp,dim=c(tmp[1][1],4))

#' @useDynLib cellshape CellContour_
CellContour <- function(x) {
  .C(CellContour_,as.integer(x),dim(x),as.integer(array(0,c(100000,4))))[[3]]
}
