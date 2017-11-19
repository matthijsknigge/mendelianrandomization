#' Perform pre-process
#' @author Matthijs Knigge
#' @description pre-process filters on genome-wide significance, SNPs without effectsize, no allelic information, removal of duplicates, removes alleles from which it is not possible to measure the direction, flips alleles and beta when needed,
#'
#' @param B Vector of genetic effects
#' @param B.se Standard errors of genetic effects
#' @param pval p-value
#' @param effect_allele vector containing the effect_alleles
#' @param other_allele vector containing the effect_alleles
#' @param SNP vector containing rs_ids
#' @keywords pre process
#' @export
#' @examples
#' mr.pre.process()
#'
#' @return cleaned data
mr.pre.process <- function(B, B.se, pval, effect_allele, other_allele, SNP){
  # bind to frame
  data <- data.frame(SNP = SNP, beta = B, se = B.se, effect_allele = effect_allele, other_allele = other_allele, pval = pval, stringsAsFactors=FALSE)
  # filter on p-value
  data <- data[which(data$pval < 5*10^-8),]
  # delete exposures without effectsize or OR
  if(length(which(is.na(data$beta)) > 0)){
    data <- data[-which(is.na(data$beta)), ]
  }
  # delete exposures without ellelic information
  if(length(which(data$effect_allele == "")) > 0){
    data <- data[-which(data$effect_allele == ""), ]
  }
  if(length(which(data$other_allele == "")) > 0){
    data <- data[-which(data$other_allele == ""), ]
  }
  # delete exposures without ellelic information
  if(length(which(is.na(data$effect_allele))) > 0){
    data <- data[-which(is.na(data$effect_allele)), ]
  }
  if(length(which(is.na(data$other_allele))) > 0){
    data <- data[-which(is.na(data$other_allele)), ]
  }
  # remove duplicates
  data <- data[!duplicated(data$SNP),]
  # determine alleles to flip
  effect_allele <- data[which(data$beta < 0),]$effect_allele;  other_allele  <- data[which(data$beta < 0),]$other_allele
  # flip alleles
  data[which(data$beta < 0),]$effect_allele <- other_allele; data[which(data$beta < 0),]$other_allele <- effect_allele
  # flip beta
  data[which(data$beta < 0),]$beta <- data[which(data$beta < 0),]$beta *-1
  # remove alleles from which it is not possible to measure direction
  remove.alleles <- data[- which(data$effect_allele == "C" & data$other_allele == "G" | data$effect_allele == "G" & data$other_allele == "C" | data$effect_allele == "A" & data$other_allele == "T" | data$effect_allele == "T" & data$other_allele == "A" ),]
  # test
  if (length(remove.alleles$SNP) != 0){
    data <- remove.alleles
  }
  return(data)
}
