#' Calculate standard error of effect size
#' @author Matthijs Knigge
#'
#'
#' @param B Vector of genetic effects
#' @param pval Vector of p values of genetic effects
#'
#' @keywords se
#' @export
#' @examples
#' mr.calculate.se()
#'
#' @return List with the following elements:
#'         se: standard error of genetic effect
mr.calculate.se <- function(B, pval){
  # se
  se <-sqrt(B^2/(qchisq(pval, 1, lower.tail = FALSE)))

  return(list(se = se))
}
