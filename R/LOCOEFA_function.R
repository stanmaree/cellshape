#' Calculate LOCO-EFA Coefficients
#'
#' This function allows you to calculate the LOCO-EFA coefficients based upon a cell outline 
#' @param array of x and y coordinates of the outline
#' @keywords LOCO-EFA
#' @export
#' @examples
#' LOCOEFA()

#' @useDynLib cellshape LOCOEFA_
LOCOEFA <- function(x) {
  .C(LOCOEFA_,as.numeric(x),dim(x),as.numeric(array(0,c(52,6))))[[3]]
}
