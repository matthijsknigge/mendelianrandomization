#' Perform wald-ratio
#' @author Matthijs Knigge
#'
#'
#' @param By Vector of genetic effects of outcome
#' @param Bx Vector of genetic effects of exposure
#' @param By.se Standard errors of genetic effects on outcome
#' @param Bx.se Standard errors of genetic effects on exposure
#'
#' @keywords wald-ratio
#' @export
#' @examples
#' mr.wald.ratio()
#'
#'
#' @return List with the following elements:
#'         iv: causal effect estimate
#'         iv.se: standard error
#'         iv.p: p-value
mr.wald.ratio <- function(By, Bx, By.se, Bx.se){
  # iv
  iv <- By / Bx
  # z-scores
  By.z <- By / By.se
  Bx.z <- Bx / Bx.se
  # iv.se
  iv.se <-sqrt (iv^2 / ((By.z^2 * Bx.z^2) /
                        (By.z^2 + Bx.z^2)))
  # iv.p
  iv.p <-  pnorm((iv / iv.se), lower.tail = FALSE)

  return(list(iv = iv, iv.se = iv.se, iv.p = iv.p))
}
