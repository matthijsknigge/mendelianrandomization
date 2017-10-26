#' Plot MR results
#' @author Matthijs Knigge
#'
#'
#' @param By Vector of genetic effects on exposure
#' @param Bx Vector of genetic effects on outcome
#' @param By.se Standard errors of genetic effects on exposure
#' @param Bx.se Standard errors of genetic effects on outcome
#' @param iv causal effect estimate
#' @param iv.se standard error of causal effect estimate
#'
#' @export
#' @return plot object
mr.plot <- function(By, Bx, By.se, Bx.se, iv, iv.se, ivw = NULL, egger = NULL, egger.i = NULL){
  require(ggplot2); require("RColorBrewer"); require(latex2exp)
  # calculate z-score of Instrumental Variable
  iv.z <- iv / iv.se
  # plot exposure effectsize against outcome effectsize
  p <- ggplot(data = NULL, aes(x = Bx, y = By))
  # color scale the points based on z-score, absolute value of se different from zero
  p <- p + geom_point(aes(colour = iv.z), alpha = 0.8, size=4)
  # ablines
  p <- p + geom_hline(yintercept=0, color="darkgrey")
  p <- p + geom_vline(xintercept=0, color="darkgrey")
  # errorbars
  p <- p + geom_errorbar(aes(ymin = By - qnorm(0.975)*By.se, ymax = By + qnorm(0.975)*By.se), colour = "darkblue", alpha = 0.3)
  p <- p + geom_errorbarh(aes(xmin = Bx - qnorm(0.975)*Bx.se, xmax = Bx + qnorm(0.975)*Bx.se), colour = "darkblue", alpha = 0.3)
  # draw inverse varaince weighted
  if(!is.null(ivw)){
  p <- p + geom_abline(aes(intercept = 0, slope = ivw, linetype = "IVW"), color = "darkblue", alpha = .2)
  }
  # draw egger method
  if(!is.null(egger) & !is.null(egger.i)){
  p <- p + geom_abline(aes(intercept = egger.i, slope = egger, linetype = "MRegger"), color = "darkblue", alpha = .2)
  }
  # manually change legend for lines
  p <- p + scale_linetype_manual(name = "Method", values = c(IVW = "dashed", MRegger = "solid", IVW.Q = "dashed", MRegger.Q = "solid"))


  # stolen theme for plot
  p <- p + theme_minimal()
  # labs and titles
  p <- p + labs(y= TeX('$\\hat{\\beta}_{outcome}$'),
                x= TeX('$\\hat{\\beta}_{exposure}$'),
                title = "Something ~ something")
  p <- p + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=14,face="bold"))




  return(p)
}
