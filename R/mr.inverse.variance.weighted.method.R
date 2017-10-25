#' Perform inverse variance weighted method
#' @author Matthijs Knigge
#'
#' @param By Vector of genetic effects on exposure
#' @param Bx Vector of genetic effects on outcome
#' @param By.se Standard errors of genetic effects on exposure
#' @param Bx.se Standard errors of genetic effects on outcome
#'
#' @keywords inverse variance weighted method
#' @export
#' @examples
#' mr.inverse.variance.weighted.method()
#'
#' @return List with the following elements:
#'         ivw: MR estimate
#'         ivw.se: standard error
#'         ivw.p: p-value
mr.inverse.variance.weighted.method <- function(By, Bx, By.se, Bx.se){
  # beta.ivw
  ivw     <- sum(By*Bx*By.se^-2)/sum(Bx^2*By.se^-2)
  # beta.ivw.se
  ivw.se  <- 1/sqrt(sum(Bx^2*By.se^-2))
  # beta.ivw.p
  ivw.p   <- pchisq((ivw / ivw.se)^2, 1, lower.tail = FALSE)


  return(list(ivw = ivw, ivw.se = ivw.se, ivw.p = ivw.p))
}
