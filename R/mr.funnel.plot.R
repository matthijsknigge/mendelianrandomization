#' MR funnel plot
#' @author Matthijs Knigge
#'
#' @param iv vector of causal effect estimate
#' @param iv.se vector of standard error of causal effect estimate
#' @param method.estimate.before.Q the beta estimate of the chosen method to plot before Chochren's Q test
#' @param method.estimate.before.Q the beta estimate of the chosen method to plot after Chochren's Q test. Default NULL.
#' @param method.name name of the method, used for legend.
#' @param linetype what linetype. Default dashed.
#' @param legend boolean. Use legend, or remove. Default is TRUE
#' @param position position of the legend. Default bottom
#' @param cochran.Q vector indicating if the SNPs is removed by Chochran's Q test due to pleiotropic effects. Default is NULL.
#' @param outcome.name string of outcome name
#' @param exposure.name string of exposure name
#'
#' @keywords funnel
#' @export
#' @examples
#' mr.funnel.plot()
#'
#' @return funnel plot
mr.funnel.plot <- function(iv, iv.se, method.estimate.before.Q, method.estimate.after.Q = NULL, method.name, linetype = "dashed", legend = TRUE, position = "bottom", cochran.Q = NULL,
                           outcome.name, exposure.name){
  require(ggplot2); require(latex2exp); require("RColorBrewer"); require(ggExtra); require(cowplot); require(gridExtra)
  # calculate z-score for coloring points
  iv.z <- iv / iv.se
  # check if Chochran's Q is used on data set
  if(sum(cochran.Q) == 0){
    cochran.Q <- NULL
    method.estimate.after.Q <- NULL
  }

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

  # Chochran's Q points
  if(!is.null(cochran.Q)){
    p <- p + geom_point(aes(colour = iv.z, shape = factor(cochran.Q)), alpha = 0.8, size=4)
    p <- p + geom_point(data=NULL, aes(x=iv[which(cochran.Q == 0)], y=iv.se[which(cochran.Q == 0)]), colour="red", size=4.1)
    p <- p + guides(shape = guide_legend(override.aes = list(size = 5, shape=c(21,21), colour="white", fill=c("red", brewer.pal(9, "Blues")[7]))))
    p <- p + scale_shape_manual(values=c(19, 19), name = TeX("$\\chi^2_{homogeneity}$"), labels = c("p > 0.05", "p < 0.05"))
  }

  # method estimate before Chochran's Q test
  p <- p +  geom_vline(aes(xintercept=method.estimate.before.Q, linetype = method.name), color = "darkblue", alpha = .3)
  # method estimate after Chochran's Q test
  if(!is.null(cochran.Q)){
    # method after Q
    p <- p +  geom_vline(aes(xintercept=method.estimate.after.Q, linetype = paste0(method.name, ".Q")), color = "darkblue")
    # calc confidence intervals after Q test
    dfCI.after.Q = calc.95.99.CL(method.estimate = method.estimate.after.Q, iv.se = iv.se)
    # add confidence interval
    p <- p + geom_line(aes(y = se.seq, x = ll95), linetype = 'dotted', data = dfCI.after.Q, color = "darkblue", alpha = .6)
    p <- p + geom_line(aes(y = se.seq, x = ul95), linetype = 'dotted', data = dfCI.after.Q, color = "darkblue", alpha = .6)
    p <- p + geom_line(aes(y = se.seq, x = ll99), linetype = 'twodash', data = dfCI.after.Q, color = "darkblue", alpha = .6)
    p <- p + geom_line(aes(y = se.seq, x = ul99), linetype = 'twodash', data = dfCI.after.Q, color = "darkblue", alpha = .6)
    # shade area
    p <- p + geom_area(aes(x = dfCI.after.Q$ll99, y = dfCI.after.Q$se.seq), color = "grey", alpha = .1, size = 0)
    p <- p + geom_area(aes(x = dfCI.after.Q$ul99, y = dfCI.after.Q$se.seq), color = "grey", alpha = .1, size = 0)
  }

  # add confidence interval
  p <- p + geom_line(aes(y = se.seq, x = ll95), linetype = 'dotted', data = dfCI.before.Q, color = "grey")
  p <- p + geom_line(aes(y = se.seq, x = ul95), linetype = 'dotted', data = dfCI.before.Q, color = "grey")
  p <- p + geom_line(aes(y = se.seq, x = ll99), linetype = 'twodash', data = dfCI.before.Q, color = "grey")
  p <- p + geom_line(aes(y = se.seq, x = ul99), linetype = 'twodash', data = dfCI.before.Q, color = "grey")

  # shade region
  if(is.null(cochran.Q)){
    p <- p + geom_area(aes(x = dfCI.before.Q$ll99, y = dfCI.before.Q$se.seq), color = "grey", alpha = .1, size = 0)
    p <- p + geom_area(aes(x = dfCI.before.Q$ul99, y = dfCI.before.Q$se.seq), color = "grey", alpha = .1, size = 0)
  }

  # flip y axis
  p <- p + scale_y_reverse()
  # labs for the axes
  p <- p +  labs(x = TeX("$\\beta_{IV}$"), y = TeX("$\\sigma_{\\beta_{IV}}$"),
                 title = paste0(outcome.name, " ~ ", exposure.name))
  # font size and type
  p <- p + theme(axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=14,face="bold"))
  # legend for scale gradient
  p <- p + scale_colour_gradientn(colours = brewer.pal(9, "Blues")[2:9], name = TeX('$\\frac{\\hat{\\beta}_{IV}}{\\sigma_{\\hat{\\beta}_{IV}}}$'))
  # manually change linetye or override the existing one
  p <- p + scale_linetype_manual(name = "Method", values = c(assign(method.name, linetype), assign(paste0(method.name, ".Q"), linetype)))
  # manually change legend, ugly code
  if(!is.null(cochran.Q)){lty <- c(assign(method.name, linetype), assign(paste0(method.name, ".Q"), linetype)); co <- c("darkblue", "darkblue"); alp <- c(0.1,1)}
  if(is.null(cochran.Q)){lty <- c(assign(method.name, linetype)); co <- c("darkblue"); alp <- c(.1)}
  p <- p + guides( linetype = guide_legend(override.aes = list(linetype=lty, color=co, alpha=alp)))

  #  add rug
  p <- p + geom_rug(aes(color = iv.z))
  # add rug for Chochran's Q
  if(!is.null(cochran.Q)){
    p <- p + geom_rug(data=NULL, aes(x=iv[which(cochran.Q == 0)], y=iv.se[which(cochran.Q == 0)]), colour="red")
  }

  # stolen theme for plot
  p <- p + theme_minimal()

  # legend position
  p <- p + theme(legend.position=position)

  # legend
  if(legend == FALSE){
    p <- p + theme(legend.position="none")
  }

  # shade part in the end, otherwise can't query max/min limits of x/y axis
  if(!is.null(cochran.Q)){
    # shade blocks that fall out confidence interval for lower part
    lowest.value.x.axis <- layer_scales(p)$x$range$range[1]
    min.term.ll99 <- min(dfCI.after.Q$ll99)
    lowest.value.y.axis <- layer_scales(p)$y$range$range[1]
    highest.value.y.axis <- layer_scales(p)$y$range$range[2]
    p <- p + annotate("rect", xmin=lowest.value.x.axis, xmax=min.term.ll99, ymin=lowest.value.y.axis*-1, ymax=highest.value.y.axis, alpha=0.3, fill="grey")

    # shade blocks that fall out confidence interval for lower part
    highest.value.x.axis <- layer_scales(p)$x$range$range[2]
    max.term.ul99 <- max(dfCI.after.Q$ul99)
    p <- p + annotate("rect", xmin=max.term.ul99, xmax=highest.value.x.axis, ymin=lowest.value.y.axis*-1, ymax=highest.value.y.axis, alpha=0.3, fill="grey")
  }
  if(is.null(cochran.Q)){
    # shade blocks that fall out confidence interval for lower part
    lowest.value.x.axis <- layer_scales(p)$x$range$range[1]
    min.term.ll99 <- min(dfCI.before.Q$ll99)
    lowest.value.y.axis <- layer_scales(p)$y$range$range[1]
    highest.value.y.axis <- layer_scales(p)$y$range$range[2]
    p <- p + annotate("rect", xmin=lowest.value.x.axis, xmax=min.term.ll99, ymin=lowest.value.y.axis*-1, ymax=highest.value.y.axis, alpha=0.3, fill="grey")

    # shade blocks that fall out confidence interval for lower part
    highest.value.x.axis <- layer_scales(p)$x$range$range[2]
    max.term.ul99 <- max(dfCI.before.Q$ul99)
    p <- p + annotate("rect", xmin=max.term.ul99, xmax=highest.value.x.axis, ymin=lowest.value.y.axis*-1, ymax=highest.value.y.axis, alpha=0.3, fill="grey")
  }
  # for some #*$%! reason the variable must be in the global scope of R to be inserted into a grid??????????
  iv <<- iv; cochran.Q <<- cochran.Q;
  # x-axis
  xdens <- axis_canvas(p, axis = "x")
  if(is.null(cochran.Q)){
    xdens <- xdens + geom_density(data = NULL, aes(x = iv), fill = "steelblue", alpha = 1, size = 0.2)
    xdens <- xdens + geom_vline(xintercept = method.estimate.before.Q, color = "darkblue", linetype = linetype)
  }
  if(!is.null(cochran.Q)){
    xdens <- xdens + geom_density(data = NULL, aes(x = iv[which(cochran.Q == 1)]), fill = "steelblue", alpha = 1, size = 0.2)
    xdens <- xdens + geom_density(data = NULL, aes(x = iv[which(cochran.Q == 0)]), fill = "red", alpha = .3, size = 0.2)
    xdens <- xdens + geom_vline(xintercept = method.estimate.before.Q, color = "grey", linetype = linetype)
    xdens <- xdens + geom_vline(xintercept = method.estimate.after.Q, color = "darkblue", linetype = linetype)
  }
  # insert into grid
  p <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  # remove from global environment
  rm(iv, envir = .GlobalEnv); rm(cochran.Q, envir = .GlobalEnv)
  return(p)
}
