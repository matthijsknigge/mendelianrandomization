#' Find missing allele
#' @author Matthijs Knigge
#' @description If the outcome or exposure is missing allelic information, this function will sort out what is the missing allele, and what is the effect and other allele.
#'
#' @param thousand.G absolute location to 1000Genomes file
#' @param data data.frame that must contain effect_allele, and SNP columns
#' @keywords missing allele
#' @export
#' @examples
#' mr.find.missing.allelic.information(data = outcome, thousand.G = '/path/to/thousand/genomes/file.bim')
#'
#' @return data with both alleles effect_allele and other_allele
mr.find.missing.allelic.information <- function(data, thousand.G){
  require(data.table)
  # read thousand G
  thousand.G <- fread(thousand.G); colnames(thousand.G) <- c("chr_name", "SNP", "pos", "start", "A1", "A2")
  # intersect between refdata and exposure
  int <- Reduce(intersect, list(data$SNP, thousand.G$SNP))
  # apply intersect
  data <- data[data$SNP %in% int, ]; subset.thousand.G <- thousand.G[thousand.G$SNP %in% int, ]
  # order data for merging
  data <- data[order(data$SNP), ]; subset.thousand.G <- subset.thousand.G[order(subset.thousand.G$SNP), ]
  # remove duplicates
  data <- data[!duplicated(data$SNP),]
  # merge
  data$A1 <- subset.thousand.G$A1; data$A2 <- subset.thousand.G$A2; data$other_allele <- ""
  # determine which is effect and other allele and apply
  data[which(data$effect_allele == data$A1), ]$other_allele <- data[which(data$effect_allele == data$A1), ]$A2
  data[which(data$effect_allele == data$A2), ]$other_allele <- data[which(data$effect_allele == data$A2), ]$A1
  # clear workspace
  data$A1 <- NULL; data$A2 <- NULL
  return(data)
}
