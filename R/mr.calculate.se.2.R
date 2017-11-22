#' Calculate standard error of effect size
#' @author Matthijs Knigge
#'
#'
#' @param p effect allele frequence
#' @param n Sample size
#' @param Z Vector of Z scores from genetic effects
#'
#' @keywords se
#' @export
#' @examples
#' mr.calculate.se.2()
#'
#' @return List with the following elements:
#'         se: standard error of genetic effect
mr.calculate.se.2<- function(p, n, Z){
  # se
  se <- 1/sqrt(2*p*(1-p)*(n + Z^2))
  # return se
  return(list(se = se))
}
