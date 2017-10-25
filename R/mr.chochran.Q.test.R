#' Perform Chochran's Q test
#' @author Matthijs Knigge
#'
#'
#' @return calculation phenotype, SNP, beta.exposure, se.exposure, beta.outcome, se.outcome, beta.iv, beta.iv.se, beta.iv.p
#' @return method beta.ivw, beta.ivw.se, beta.ivw.p, beta.egger, beta.egger.se, beta.egger.p, beta.egger.i, beta.egger.i.se, beta.egger.i.p
#' @param p.threshold p value threshold
#' @keywords chochrans Q
#' @export
#' @examples
#' mr.chochran.Q.test()
#'
#' @return calculation with removed pleiotrpoic effect, and re-calculated egger, and ivw method
mr.chochran.Q.test <- function(calculation, method, p.threshold){
  # add data slots for chochran Q
  method$beta.egger.Q <- 0; method$beta.egger.p.Q <- 0; method$beta.egger.i.Q <- 0; method$beta.ivw.Q <- 0; method$beta.ivw.p.Q <- 0; method$beta.ivw.se.Q <- 0; method$beta.egger.se.Q <- 0
  # save copy of data
  calculation.copy <<- calculation
  method.copy <<- method
  # delete snp if beta.iv.se is not available
  if(length(which(is.na(calculation.copy$beta.iv.se)) > 0)){
    calculation.copy <<- calculation.copy[-which(is.na(calculation.copy$beta.iv.se)), ]
  }
  # store chochrans.q.p
  chochrans.q.p <<- 0.0

  while(p.threshold > chochrans.q.p & length(calculation.copy$SNP) >= 3){
    # perform chochrans.q
    calculation.copy$chochrans.q <<- 1/calculation.copy$beta.iv.se * (calculation.copy$beta.iv - method.copy$beta.ivw)^2
    # determine maximum term
    chochrans.q.max <- calculation.copy[which.max(calculation.copy$chochrans.q), ]$SNP
    # calcule the p of the sum of chochrans.q terms
    chochrans.q.p <<-  pchisq(sum(calculation.copy$chochrans.q), df = length(calculation.copy$chochrans.q-1), lower.tail = FALSE)
    # remove maximum term
    calculation.copy <<- calculation.copy[-which(calculation.copy$SNP == chochrans.q.max), ]
    # re-calculate ivw
    method.copy$beta.ivw <- mr.beta.ivw(By = calculation.copy$beta.outcome, Bx = calculation.copy$beta.exposure, By.se = calculation.copy$se.outcome, Bx.se = calculation.copy$se.exposure)$beta.ivw

  }
  # if not able to perform chochrans.q
  if(length(calculation.copy$SNP) < 3){
    return(list(calculation = calculation, method = method))
  }
  # re-calculate ivw, and egger
  ivw <- mr.beta.ivw(By = calculation.copy$beta.outcome, Bx = calculation.copy$beta.exposure, By.se = calculation.copy$se.outcome, Bx.se = calculation.copy$se.exposure)
  egger <- mr.egger(By = calculation.copy$beta.outcome, Bx = calculation.copy$beta.exposure, By.se = calculation.copy$se.outcome, Bx.se = calculation.copy$se.exposure)
  # re-assign ivw, and egger
  method$beta.ivw.Q <- ivw$beta.ivw; method$beta.ivw.p.Q <- ivw$beta.ivw.p
  method$beta.egger.Q <- egger$beta.egger; method$beta.egger.p.Q <- egger$beta.egger.p; method$beta.egger.i.Q <- egger$beta.egger.i
  method$beta.egger.se.Q <- egger$beta.egger.se; method$beta.ivw.se.Q <- ivw$beta.ivw.se
  # intercept between and determine which snps did not participate in the calculation
  int <- Reduce(intersect, list(calculation.copy$SNP, calculation$SNP))
  calculation[calculation$SNP %in% int, ]$chochrans.q <- 1

  return(list(calculation = calculation, method = method))
}
