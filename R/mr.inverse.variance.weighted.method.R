#' Perform inverse variance weighted method
#' @author Matthijs Knigge
#'
#' @param By Vector of genetic effects of outcome
#' @param Bx Vector of genetic effects of exposure
#' @param By.se Standard errors of genetic effects on outcome
#' @param Bx.se Standard errors of genetic effects on exposure
#' @param subset vector containing 1-0 indicating to use it in the calculation or not. Default is NULL.
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
mr.inverse.variance.weighted.method <- function(By, Bx, By.se, Bx.se, subset = NULL){
  # if sum subset == 0, do not apply
  if(!is.null(subset) & sum(subset) == 0){
    return(list(ivw = NA, ivw.se = NA, ivw.p = NA))
  }
  # if subset is used, filter data
  if(!is.null(subset)){
    By <- By[which(subset == 1)]
    Bx <- Bx[which(subset == 1)]
    By.se <- By.se[which(subset == 1)]
    Bx.se <- Bx.se[which(subset == 1)]
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
