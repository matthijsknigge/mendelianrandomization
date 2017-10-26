#' Perform Chochran's Q test
#' @author Matthijs Knigge
#' @description Chochran's Q is a iterative method that tries to find pleiotropic SNPs, and filters those SNPs out.
#'
#' @param data data.frame containen SNP, By, Bx, By.se, Bx.se, iv, iv.se
#' @param pval p value threshold
#' @keywords chochrans Q
#' @export
#' @examples
#' mr.chochran.Q.test()
#'
#' @return data with column chochran's Q column added which says if the SNPs is used or not in the method.
mr.chochran.Q.test <- function(data, pval){
  # add slot
  data$chochrans.q <- 0
  # deep copy of ivw
  ivw <<- 0
  # save copy of data
  data.copy <<- data
  # delete snp if beta.iv.se is not available
  if(length(which(is.na(data.copy$iv.se)) > 0)){
    data.copy <<- data.copy[-which(is.na(data.copy$iv.se)), ]
  }
  # store chochrans.q.p
  chochrans.q.p <<- 0.0

  while(pval > chochrans.q.p & length(data.copy$SNP) >= 3){
    # perform chochrans.q
    data.copy$chochrans.q <<- 1/data.copy$iv.se * (data.copy$iv - ivw)^2
    # determine maximum term
    chochrans.q.max <- data.copy[which.max(data.copy$chochrans.q), ]$SNP
    # calcule the p of the sum of chochrans.q terms
    chochrans.q.p <<-  pchisq(sum(data.copy$chochrans.q), df = length(data.copy$chochrans.q-1), lower.tail = FALSE)
    # remove maximum term
    data.copy <<- data.copy[-which(data.copy$SNP == chochrans.q.max), ]
    # re-calculate ivw
    ivw <<- mr.inverse.variance.weighted.method(By = data.copy$By, Bx = data.copy$Bx, By.se = data.copy$By.se, Bx.se = data.copy$Bx.se)$ivw

  }
  # if not able to perform chochrans.q
  if(length(data.copy$SNP) < 3){
    return(data)
  }
  # intercept between and determine which snps did not participate in the calculation
  int <- Reduce(intersect, list(data.copy$SNP, data$SNP))
  # 0 means not used in Chochrans Q test thus pleiotropic
  data[data$SNP %in% int, ]$chochrans.q <- 1
  # clear workspace
  rm(data.copy, envir = .GlobalEnv); rm(ivw, envir = .GlobalEnv); rm(chochrans.q.p, envir = .GlobalEnv);
  return(data)
}
