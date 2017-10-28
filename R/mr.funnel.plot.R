#' MR funnel plot
#' @author Matthijs Knigge
#'
#' @param iv vector of causal effect estimate
#' @param iv.se vector of standard error of causal effect estimate
#' @param method.estimate.before.Q the beta estimate of the chosen method to plot before Chochren's Q test
#' @param method.estimate.before.Q the beta estimate of the chosen method to plot after Chochren's Q test. Default NULL.
#' @param method.name name of the method, used for legend.
#' @param linetype what linetype. Default dashed.
#'
#' @keywords funnel
#' @export
#' @examples
#' mr.funnel.plot()
#'
#' @return funnel plot
mr.funnel.plot <- function(iv, iv.se, method.estimate.before.Q, method.estimate.after.Q = NULL, method.name, linetype = "dashed"){
  require(ggplot2); require(latex2exp); require("RColorBrewer"); require(ggExtra)
  # calculate z-score for coloring points
  iv.z <- iv / iv.se

  # before Chochran's Q, calculate confidence interval of the chosen method
  # vector of values that spans the range from 0 to the max value of impression
  calc.95.99.CL <- function(method.estimate, iv.se){
    se.seq<- seq(0, max(iv.se), 0.001)
    # 95% CI region
    ll95 = method.estimate-(1.96*se.seq)
    ul95 = method.estimate+(1.96*se.seq)
    # 99% CI region
    ll99 = method.estimate-(3.29*se.seq)
    ul99 = method.estimate+(3.29*se.seq)
    # return
    return(data.frame(ll95, ul95, ll99, ul99, se.seq, method.estimate))
  }
  # frame
  dfCI.before.Q = calc.95.99.CL(method.estimate = method.estimate.before.Q, iv.se = iv.se)

  # initiate plot
  p <- ggplot(data = NULL, aes(x = iv, y = iv.se))
  # points
  p <- p + geom_point(aes(colour = iv.z), alpha = 0.8, size=4)

  # method estimate before Chochran's Q test
  p <- p +  geom_vline(aes(xintercept=method.estimate.before.Q, linetype = method.name), color = "darkblue")
  # method estimate after Chochran's Q test
  if(!is.null(method.estimate.after.Q)){
    p <- p +  geom_vline(aes(xintercept=method.estimate.after.Q, linetype = paste0(method.name, ".Q")), color = "darkblue", alpha = .3)
  }
  # add confidence interval
  p <- p + geom_line(aes(y = se.seq, x = ll95), linetype = 'dotted', data = dfCI.before.Q, color = "grey")
  p <- p + geom_line(aes(y = se.seq, x = ul95), linetype = 'dotted', data = dfCI.before.Q, color = "grey")
  p <- p + geom_line(aes(y = se.seq, x = ll99), linetype = 'dashed', data = dfCI.before.Q, color = "grey")
  p <- p + geom_line(aes(y = se.seq, x = ul99), linetype = 'dashed', data = dfCI.before.Q, color = "grey")

  # flip y axis
  p <- p + scale_y_reverse()
  # labs for the axes
  p <- p +  labs(x = TeX("$\\beta_{IV}$"), y = TeX("$\\sigma_{\\beta_{IV}}$"))
  # font size and type
  p <- p + theme(axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=14,face="bold"))
  # legend for scale gradient
  p <- p + scale_colour_gradientn(colours = brewer.pal(9, "Blues")[2:9], name = TeX('$\\frac{\\hat{\\beta}_{IV}}{\\sigma_{\\hat{\\beta}_{IV}}}$'))
  # manually change linetye or override the existing one
  p <- p + scale_linetype_manual(name = "Method", values = c(assign(method.name, linetype), assign(paste0(method.name, ".Q"), linetype)))
  # manually change legend, ugly code
  if(!is.null(method.estimate.after.Q)){lty <- c(assign(method.name, linetype), assign(paste0(method.name, ".Q"), linetype)); co <- c("darkblue", "darkblue"); alp <- c(0.1,1)}
  if(is.null(method.estimate.after.Q)){lty <- c(assign(method.name, linetype)); co <- c("darkblue"); alp <- c(1)}
  p <- p + guides( linetype = guide_legend(override.aes = list(linetype=lty, color=co, alpha=alp)))


  #  add rug
  p <- p + geom_rug(aes(color = iv.z))

  # stolen theme for plot
  p <- p + theme_minimal()

  return(p)
}
