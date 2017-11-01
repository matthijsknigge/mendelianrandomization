#' Remove region from file
#' @author Matthijs Knigge
#'
#'
#' @param data data.frame must contain column pos, SNP, and chr
#' @param chr chrosome number from wich the region must be romved
#' @param lower.limit lower limit of the region to be removed in kb
#' @param upper.limit upper limit of the region to be removed in kb
#'
#' @keywords remove region
#' @export
#' @examples
#' mr.remove.region()
#'
#' @return data without region selected
mr.remove.region <- function(chr, lower.limit, upper.limit, data){
  require(data.table)
  # get chromosome from data
  data.chr <- data[which(data$chr_name == chr), ]
  # remove chromosome from data
  data <- data[which(data$chr_name != chr),]
  # get SNP from region
  int <- data.chr[which(data.chr$pos > lower.limit & data.chr$pos < upper.limit), ]$SNP
  # intersect and thereby remove the region
  data.chr <- data.chr[!data.chr$SNP %in% int, ]
  # bind together data with chromosome removed with chromosome with removed region
  data <- rbind(data, data.chr)
  return(data)
}
