#' Calculate effect size
#' @author Matthijs Knigge
#'
#'
#' @param se Vector of SE on genetic effect
#' @param Z Vector of Z scores from genetic effects
#'
#' @keywords se
#' @export
#' @examples
#' mr.calculate.B()
#'
#' @return List with the following elements:
#'         B: genetic effect size
mr.calculate.B<- function(se, Z){
  # se
  B <- Z*se
  # return B
  return(list(B = B))
}
