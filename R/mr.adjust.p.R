#' Do p value correction on data set
#' @author Matthijs Knigge
#'
#'
#' @param data Data frame containing column to do p value correction on.
#' @param column.name The name of the column to do p value correction on, This must be a list, for it can adjust multiple columns at once
#' @param method The method to use for p value correction.
#'
#' @keywords p adjust
#' @export
#' @examples
#' mr.adjust.p()
#'
#' @return Adjusted data frame with adjusted pvalues
#'
mr.adjust.p <- function(data, column.names, method){
  # find column to adjust, do pvalue correction with given method
  for(i in 1:length(column.names)){
    # assign column name slot
    column.name <- column.names[i]

    # adjust column
    data[, paste0(column.name, ".", method)] <- p.adjust(p = data[, get(as.character(column.name))], method = "fdr", n = length(data[, get(as.character(column.name))]))

  }
  # return frame
  return(data)
}
