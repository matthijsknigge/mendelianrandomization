#' Perform mr harmonization
#' @author Matthijs Knigge
#' @description allelic harmonization between outcome and exposure. The intersect is taken from both outcome and exposure. The alleles and effect sizes are aligned. Mismatch alles are removed. And duplicates are removed.
#'
#' @param By Vector of genetic effects  from outcome
#' @param By.se Standard errors of genetic effects from outcome
#' @param outcome.pval p-value from outcome
#' @param outcome.effect_allele vector containing the effect_alleles from outcome
#' @param outcome.other_allele vector containing the effect_alleles from outcome
#' @param outcome.SNP vector containing rs_ids from outcome
#' @param Bx Vector of genetic effects from exposure
#' @param Bx.se Standard errors of genetic effects from exposure
#' @param exposure.pval p-value from exposure
#' @param exposure.effect_allele vector containing the effect_alleles from exposure
#' @param exposure.other_allele vector containing the effect_alleles from exposure
#' @param exposure.SNP vector containing SNPs from exposure
#' @keywords harmonize
#' @export
#' @examples
#' mr.harmonize()
#'
#' @return harmonized data
mr.harmonize <- function(By, Bx, By.se, Bx.se, outcome.pval, exposure.pval, outcome.effect_allele, exposure.effect_allele, exposure.other_allele, outcome.SNP, exposure.SNP){
  # assign data slots
  o <- data.frame(beta = By, se = By.se, pval = outcome.pval, effect_allele = outcome.effect_allele, SNP = outcome.SNP, stringsAsFactors=FALSE)
  e <- data.frame(beta = Bx, se = Bx.se, pval = exposure.pval, effect_allele = exposure.effect_allele, other_allele = exposure.other_allele, SNP = exposure.SNP, stringsAsFactors=FALSE)
  # delete exposures without effectsize or OR
  if(length(which(is.na(e$beta)) > 0)){
    e <- e[-which(is.na(e$beta)), ]
  }
  # delete outcome without effectsize or OR
  if(length(which(is.na(o$beta)) > 0)){
    o <- o[-which(is.na(o$beta)), ]
  }
  # find intersect between exposures and Celiac
  int <- Reduce(intersect, list(o$SNP, e$SNP))
  # intersect between data sets
  e <- e[e$SNP %in% int, ]; o <- o[o$SNP %in% int, ]
  # remove duplicates
  e <- e[!duplicated(e$SNP),]; o <- o[!duplicated(o$SNP),]
  # order data
  o <- o[order(o$SNP), ]; e <- e[order(e$SNP), ];
  # always assume that the alignment is incorrect
  o[which(as.character(o$effect_allele) != as.character(e$effect_allele)), ]$beta <-  o[which(as.character(o$effect_allele) != as.character(e$effect_allele)), ]$beta * -1
  # order on snps
  o <- o[order(o$SNP), ]; e <- e[order(e$SNP), ];
  # find mismatch between outcome effect_allele != exposure effect_allele
  remove.effect.allele.snp <- o[which(as.character(o$effect_allele) != as.character(e$effect_allele)), ]$SNP
  # find mismatch between outcome effect_allele != exposure other_allele
  remove.other.allele.snp <- o[which(as.character(o$effect_allele) != as.character(e$other_allele)), ]$SNP
  # intercept between this finds true mismatches
  remove.snp <- Reduce(intersect, list(remove.effect.allele.snp, remove.other.allele.snp))
  # adjust data.frame
  e <- e[!e$SNP %in% remove.snp, ]; o <- o[!o$SNP %in% remove.snp, ]
  # assign return data slots
  data <- data.frame(SNP = o$SNP, By = o$beta, Bx = e$beta, By.se = o$se, Bx.se = e$se, pval = e$pval, effect_allele = e$effect_allele, other_allele = e$other_allele)
  return(data)
}
