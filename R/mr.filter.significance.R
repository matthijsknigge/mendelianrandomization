#' Filter on pvalue
#' @author Matthijs Knigge
#'
#'
#' @param summary.table data frame summarizing all methods used
#' @param column.ivw.fdr column name of pvalues to be filtered for ivw
#' @param column.egger.fdr column name of pvalues to be filtered for egger
#' @param threshold threshold for filtering pvalues
#'
#' @keywords p adjust
#' @export
#' @examples
#' mr.filter.significance()
#'
#' @return Adjusted data frame with filtered on threshold
#'
mr.filter.significance <- function(summary.table, column.ivw.fdr, column.egger.fdr, threshold){

  # filter on given columns and return which are below threshold
  summary.table <- summary.table[which(summary.table[, get(column.ivw.fdr)] < threshold | summary.table[, get(column.egger.fdr)] < threshold), ]
  # return the table
  return(summary.table)

}
