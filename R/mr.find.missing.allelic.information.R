#' Find missing allele
#' @author Matthijs Knigge
#'
#'
#' @param refdata absolute location to 1000Genomes file
#' @param data data.frame that must contain effect_allele, and SNP columns
#' @keywords missing allele
#' @export
#' @examples
#' mr.find.missing.allele(data = outcome, thousand.G = '/path/to/thousand/genomes/file.bim')
#'
#' @return data with both alleles effect_allele and other_allele
mr.find.missing.allele <- function(data, thousand.G){
  require(R.utils)
  # read thousand G
  thousand.G <- fread(thousand.G)
  # intersect between refdata and exposure
  int <- Reduce(intersect, list(data$SNP, thousand.G$SNP))
  # apply intersect
  data <- data[data$SNP %in% int, ]; subset.thousand.G <- thousand.G[thousand.G$SNP %in% int, ]
  # order data for merging
  data <- data[order(data$SNP), ]; subset.thousand.G <- subset.thousand.G[order(subset.thousand.G$SNP), ]
  # remove duplicates
  data <- data[!duplicated(data$SNP),]
  # merge
  data$A1 <- subset.thousand.G$A1; data$A2 <- subset.thousand.G$A2
  # determine which is effect and other allele and apply
  data[which(data$effect_allele == data$A1), ]$other_allele <- data[which(data$effect_allele == data$A1), ]$A2
  data[which(data$effect_allele == data$A2), ]$other_allele <- data[which(data$effect_allele == data$A2), ]$A1
  return(data)
}
