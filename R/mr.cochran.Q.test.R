#' Perform Cochran's Q test
#' @author Matthijs Knigge
#' @description Cochran's Q is a iterative method that tries to find pleiotropic SNPs, and filters those SNPs out.
#'
#' @param data data.frame containen SNP, By, Bx, By.se, Bx.se, iv, iv.se
#' @param pval p value threshold
#' @keywords chochrans Q
#' @export
#' @examples
#' mr.cochran.Q.test()
#'
#' @return vector containing  chochran's Q which says if the SNPs is used or not in the method.
mr.cochran.Q.test <- function(data, pval){
  # no data
  if(is.data.frame(h) && nrow(h) == 0){
    return(NULL)
  }
  # if not enough data
  if(length(data$SNP) < 3){
    data$cochran.Q <- 0
    return(list(cochran.Q = data$cochran.Q))
  }
  # add slot
  data$cochran.Q <- 0
  # deep copy of ivw
  IVW <<- mr.inverse.variance.weighted.method(By = data$By, Bx = data$Bx, By.se = data$By.se, Bx.se = data$Bx.se)$ivw
  # save copy of data
  data.copy <<- data
  # delete snp if beta.iv.se is not available
  if(length(which(is.na(data.copy$iv.se)) > 0)){
    data.copy <<- data.copy[-which(is.na(data.copy$iv.se)), ]
  }
  # store chochrans.q.p
  cochran.Q.p <<- 0.0

  # before entering loop, find out if it is necessary
  Q <- 1/data.copy$iv.se * (data.copy$iv - IVW)^2
  Q.p <-  pchisq(sum(Q), df = length(Q-1), lower.tail = FALSE)
  # if q-term does not fall below threshold
  if(pval < Q.p){
    return(list(cochran.Q = data$cochran.Q))
  }

  while(pval > cochran.Q.p & length(data.copy$SNP) >= 3){
    # perform chochrans.q
    data.copy$cochran.Q <<- 1/data.copy$iv.se * (data.copy$iv - IVW)^2
    # determine maximum term
    cochran.Q.max <- data.copy[which.max(data.copy$cochran.Q), ]$SNP
    # calcule the p of the sum of chochrans.q terms
    cochran.Q.p <<-  pchisq(sum(data.copy$cochran.Q), df = length(data.copy$cochran.Q-1), lower.tail = FALSE)
    # remove maximum term
    data.copy <<- data.copy[-which(data.copy$SNP == cochran.Q.max), ]
    # re-calculate ivw
    IVW <<- mr.inverse.variance.weighted.method(By = data.copy$By, Bx = data.copy$Bx, By.se = data.copy$By.se, Bx.se = data.copy$Bx.se)$ivw

  }
  # if not able to perform chochrans.q
  if(length(data.copy$SNP) < 3){
    return(data$cochran.Q)
  }
  # intercept between and determine which snps did not participate in the calculation
  int <- Reduce(intersect, list(data.copy$SNP, data$SNP))
  # 0 means not used in Chochrans Q test thus pleiotropic
  data[data$SNP %in% int, ]$cochran.Q <- 1
  # clear workspace
  rm(data.copy, envir = .GlobalEnv); rm(IVW, envir = .GlobalEnv); rm(cochran.Q.p, envir = .GlobalEnv);
  return(list(cochran.Q = data$cochran.Q))
}
