#' Get chromosome and position based on reference file
#' @author Matthijs Knigge
#'
#'
#' @param data data.frame must contain column SNP
#' @param refdata absolute path to reference data
#'
#' @keywords get chromosome position
#' @export
#' @examples
#' mr.get.chr.pos()
#'
#' @return data with chr and pos added
mr.get.chr.pos <- function(data, refdata){
  require(data.table)
  # read reference data
  ref <- fread(refdata); colnames(ref) <- c("chr_name", "SNP", "morgan", "pos", "A1", "A2")
  # remove columns for merging
  ref$morgan <- NULL; ref$A1 <- NULL; ref$A2 <- NULL
  # merge on SNP column
  data <- merge(data, ref, by = "SNP")

  return(data)
}
