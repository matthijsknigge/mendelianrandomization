#' Perform inverse variance weighted method
#' @author Matthijs Knigge
#'
#' @param By Vector of genetic effects of outcome
#' @param Bx Vector of genetic effects of exposure
#' @param By.se Standard errors of genetic effects on outcome
#' @param Bx.se Standard errors of genetic effects on exposure
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
  # test if vector is empty
  if(length(By) == 0){
    return(list(ivw = NA, ivw.se = NA, ivw.p = NA))
  }
  # ivw
  ivw     <- sum(By*Bx*By.se^-2)/sum(Bx^2*By.se^-2)
  # ivw.se
  ivw.se  <- 1/sqrt(sum(Bx^2*By.se^-2))
  # ivw.p
  ivw.p   <- pchisq((ivw / ivw.se)^2, 1, lower.tail = FALSE)
  # return ivw
  return(list(ivw = ivw, ivw.se = ivw.se, ivw.p = ivw.p))
}
