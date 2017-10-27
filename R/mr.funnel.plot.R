#' MR funnel plot
#' @author Matthijs Knigge
#'
#' @param iv vector of causal effect estimate
#' @param iv.se vector of standard error of causal effect estimate
#'
#' @keywords funnel
#' @export
#' @examples
#' mr.funnel.plot()
#'
#' @return funnel plot
mr.funnel.plot <- function(iv, iv.se){
  require(ggplot2); require(latex2exp); require(gridExtra); require("RColorBrewer");

  p <- ggplot(data = NULL, aes(x = iv, y = iv.se))
  p <- p + geom_point()


  return(p)
}
