#' Count number of SNPs for methods
#' @author Matthijs Knigge
#'
#'
#' @param SNPs.frame data frame of SNPs that are used
#' @param method.frame summary table of methods used on SNPs
#' @param column.identifier.phenotype column name of phenotype identifier
#' @param column.identifier.SNP column name of SNP identifier
#'
#' @keywords p adjust
#' @export
#' @examples
#' mr.adjust.p()
#'
#' @return Adjusted data frame with number of SNPs added to summary table
#'
mr.count.n.SNP <- function(SNPs.frame, method.frame, column.identifier.phenotype, column.identifier.SNP){
  # initiate column for nSNP
  method.frame$nSNP <- 0

  for(p in method.frame[, get(column.identifier.phenotype)]){
    # get method for given phenotype
    # m <- method.frame[which(method.frame[, get(column.identifier.phenotype)] == p), ]
    # get SNPs that belong to given phenotype
    snp <- SNPs.frame[which(SNPs.frame[, get(column.identifier.phenotype)] == p), ]
    # amount of SNP
    amount.SNP <- length(snp[, get(column.identifier.SNP)])
    # assign number of SNPs
    method.frame[which(method.frame[, get(column.identifier.phenotype)] == p), ]$nSNP <- amount.SNP
  }

  # return method.frame
  return(method.frame)
}
