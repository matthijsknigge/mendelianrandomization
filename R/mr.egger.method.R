#' Perform mendelian randomization-egger method
#' @author Matthijs Knigge
#'
#'
#' @param By Vector of genetic effects of outcome
#' @param Bx Vector of genetic effects of exposure
#' @param By.se Standard errors of genetic effects on outcome
#' @param Bx.se Standard errors of genetic effects on exposure
#'
#' @keywords mendelian randomization egger method
#' @export
#' @examples
#' mr.egger.method()
#'
#' @return List of with the following elements:
#'         egger: MR estimate slope
#'         egger.se: Standard error of MR estimate
#'         egger.p: p-value of MR estimate
#'         egger.i: Estimate of horizontal pleiotropy (intercept)
#'         egger.i.se: Standard error of intercept
#'         egger.i.p: p-value of intercept
mr.egger.method <- function(By, Bx, By.se, Bx.se){
  # test if vector is empty
  if(length(By) == 0){
    return(list(egger = NA, egger.se = NA, egger.p = NA, egger.i = NA, egger.i.se = NA, egger.i.p = NA))
  }
  # if not enough data
  if(length(By) < 3){
    return(list(egger = NA, egger.se = NA, egger.p = NA, egger.i = NA, egger.i.se = NA, egger.i.p = NA))
  }
  # if slope can not be calculated
  if (nrow(summary(lm(By ~ Bx, weights = By.se^-2))$coef) <= 1){
    return(list(egger = NA, egger.se = NA, egger.p = NA, egger.i = NA, egger.i.se = NA, egger.i.p = NA))
  }
  # egger
  egger <- summary(lm(By~Bx, weights=By.se^-2))$coef[2,1]
  # egger.se
  egger.se <- summary(lm(By ~ Bx, weights = By.se^-2))$coef[2,2] / summary(lm(By ~ Bx, weights = By.se^-2))$sigma
  # egger.p
  egger.p <- 2 * pt(abs(egger/egger.se), df = length(Bx) - 2, lower.tail = F)

  # egger.i
  egger.i <- summary(lm(By ~ Bx, weights = By.se^-2))$coef[1,1]
  # egger.i.se
  egger.i.se <- summary(lm(By ~ Bx, weights = By.se^-2))$coef[1,2] / min(summary(lm(By ~ Bx, weights = By.se^-2))$sigma, 1)
  # egger.i.p
  egger.i.p <- 2 * pt(abs(egger.i / egger.i.se), df = length(Bx) - 2, lower.tail = F)
  # return egger
  return(list(egger = egger, egger.se = egger.se, egger.p = egger.p, egger.i = egger.i, egger.i.se = egger.i.se, egger.i.p = egger.i.p))

}
