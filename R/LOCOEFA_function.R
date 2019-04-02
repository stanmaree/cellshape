#' Calculate LOCO-EFA Coefficients
#'
#' This function allows you to calculate the LOCO-EFA coefficients based upon a cell outline 
#' @param x an array containing the number of contour points; x-position; y-position; total area.
#' The input array is most straightforwardly generated using the CellContour function.  
#' @keywords LOCO-EFA
#' @export
#' @examples
#' img <- readPNG("cell.png")
#' tmp<-array(CellContour(as.matrix(img)),dim=c(100000,4))
#' contour<-array(tmp,dim=c(tmp[1][1],4))
#' shapedescriptors<-array(LOCOEFA(as.matrix(contour)),dim=c(52,7))

#' @useDynLib cellshape LOCOEFA_
LOCOEFA <- function(x) {
  .C(LOCOEFA_,as.numeric(x),dim(x),as.numeric(array(0,c(52,7))))[[3]]
}
