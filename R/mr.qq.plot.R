#' Make quantile-quantile plot for inverse variance weighted method and MR-egger method
#' @author Matthijs Knigge
#'
#'
#' @param By Vector of genetic effects from outcome
#' @param Bx Vector of genetic effects from exposure
#' @param By.se Vector of standard errors from genetic effect of outcome
#' @param Bx.se Vector of standard errors from genetic effect of exposure
#' @param egger boolean indicating if the residuals from this method should be plotted. Default FALSE
#' @param ivw boolean indicating if the residuals from this method should be plotted. Default FALSE
#' @param chochran.Q vector containing 0-1 indicating which observations where not used in recalculating the effect. Default NULL.
#' @param outcome.name name of the outcome
#' @param exposure.name name of the exposure
#'
#' @keywords qqplot
#' @export
#' @examples
#' mr.qq.plot()
#'
#' @return qqplot
mr.qq.plot <- function(By, Bx, By.se, Bx.se, egger = FALSE, ivw = FALSE, chochran.Q = NULL, outcome.name, exposure.name) {
  require(ggplot2); require(latex2exp);
  # check Cochran's Q
  if(sum(chochran.Q) == 0){
    chochran.Q <- NULL
  }
  prepare.qqplot <- function(lm.model){
    # define model
    fitted.lm <- lm.model
    # extract residuals from lm object
    res = residuals(fitted.lm)

    # calculate slope and interncept for qqline
    slope = (quantile(res, .75) - quantile(res, .25)) / (qnorm(.75) - qnorm(.25))
    intercept = quantile(res,.25) - slope*qnorm(.25)
    qq_line = data.frame(intercept = intercept, slope = slope)
    return(list(fitted.lm = fitted.lm, res = res, slope = slope, intercept = intercept, qq_line = qq_line))
  }
  # generate ggplot for qqplot
  p <- ggplot(data = NULL)

  # add qq for ivw
  if(egger == TRUE){
    # model for egger
    egger.lm <- prepare.qqplot(lm(By ~ Bx, weights = By.se^2))
    # qq stats
    p <- p + stat_qq(aes(sample = egger.lm$res), size = 1, color = "black")
    # qq line
    p <- p + geom_abline(data = egger.lm$qq_line, aes(intercept = intercept, slope = slope, linetype = "MRegger"), color = "darkblue", size = 1, alpha = .3)
  }
  # add qq for egger
  if(ivw == TRUE){
    # model for ivw
    ivw.lm <- prepare.qqplot(lm(By ~ Bx-1, weights = By.se^2))
    # qq stats
    p <- p + stat_qq(aes(sample = ivw.lm$res), size = 1, color = "black")
    # qq line
    p <- p + geom_abline(data = ivw.lm$qq_line, aes(intercept = intercept, slope = slope, linetype = "IVW"), color = "darkblue", size = 1, alpha = .3)
  }
  # add qq for Cochran's Q
  if(!is.null(chochran.Q)){
    # Cochran's Q on IVW
    if(ivw == TRUE){
      # model for ivw for red points
      ivw.lm.Q <- prepare.qqplot(lm(By[which(chochran.Q == 0)] ~ Bx[which(chochran.Q == 0)]-1, weights = By.se[which(chochran.Q == 0)]^2))
      # model for ivw for stats
      ivw.lm.Q.s <- prepare.qqplot(lm(By[which(chochran.Q == 1)] ~ Bx[which(chochran.Q == 1)]-1, weights = By.se[which(chochran.Q == 1)]^2))
      # qq stats
      p <- p + stat_qq(aes(sample = ivw.lm.Q$res), size = 1, color = "red")
      # qq line
      p <- p + geom_abline(data = ivw.lm.Q.s$qq_line, aes(intercept = intercept, slope = slope, linetype = "IVW.Q"), color = "darkblue", size = 1)
    }
    # Cochran's Q on egger
    if(egger == TRUE){
      # model for egger for red points
      egger.lm.Q <- prepare.qqplot(lm(By[which(chochran.Q == 0)] ~ Bx[which(chochran.Q == 0)], weights = By.se[which(chochran.Q == 0)]^2))
      # model for egger for stats
      egger.lm.Q.s <- prepare.qqplot(lm(By[which(chochran.Q == 1)] ~ Bx[which(chochran.Q == 1)], weights = By.se[which(chochran.Q == 1)]^2))
      # qq stats
      p <- p + stat_qq(aes(sample = egger.lm.Q$res), size = 1, color = "red")
      # qq line
      p <- p + geom_abline(data = egger.lm.Q.s$qq_line, aes(intercept = intercept, slope = slope, linetype = "MRegger.Q"), color = "darkblue", size = 1)
    }
  }
  # add labs
  p <- p + labs(x = "Theoretical Quantile", y = "Standardized Residual", title = paste0(outcome.name, " ~ ", exposure.name))
  # stolen theme
  p <- p + theme_minimal()
  # font size and type
  p <- p + theme(axis.text=element_text(size=14, face="bold"),
                 axis.title=element_text(size=14,face="bold"))

  # override linetype
  p <- p + scale_linetype_manual(name = "Method", values = c(IVW = "dashed", MRegger = "solid", IVW.Q = "dashed", MRegger.Q = "solid"))
  # manually override legende. EXTREMELY SMELLY CODE
  if(egger == TRUE & ivw == FALSE){
    lty <- c(MRegger = "solid"); co <- c("darkblue"); alp <- c(.1); p <- p + guides( linetype = guide_legend(override.aes = list(linetype=lty, color=co, alpha=alp)))
  }
  if(egger == FALSE & ivw == TRUE){
    lty <- c(IVW = "dashed"); co <- c("darkblue"); alp <- c(.1); p <- p + guides( linetype = guide_legend(override.aes = list(linetype=lty, color=co, alpha=alp)))
  }
  if(egger == TRUE & ivw == TRUE){
    lty <- c(IVW = "dashed", MRegger = "solid"); co <- c("darkblue", "darkblue"); alp <- c(.1, .1); p <- p + guides( linetype = guide_legend(override.aes = list(linetype=lty, color=co, alpha=alp)))
  }
  if(!is.null(chochran.Q)){
    # add fake points so that the legend can be manipulated
    p <- p + geom_point(aes(x = c(0.01, 0.01), y = c(0.01, 0.01), shape = factor(c(1,0))), color = "white")
    if(ivw == TRUE & egger == FALSE){
      lty <- c(IVW = "dashed", IVW.Q = "dashed"); co <- c("darkblue", "darkblue"); alp <- c(.1, 1); p <- p + guides( linetype = guide_legend(override.aes = list(linetype=lty, color=co, alpha=alp)))
    }
    if(ivw == FALSE & egger == TRUE){
      lty <- c(MRegger = "solid", MRegger.Q = "solid"); co <- c("darkblue", "darkblue"); alp <- c(.1, 1); p <- p + guides( linetype = guide_legend(override.aes = list(linetype=lty, color=co, alpha=alp)))
    }
    if(ivw == TRUE & egger == TRUE){
      lty <- c(MRegger = "solid", MRegger.Q = "solid", IVW = "dashed", IVW.Q = "dashed"); co <- c("darkblue", "darkblue", "darkblue", "darkblue"); alp <- c(.1, 1, .1, 1); p <- p + guides( linetype = guide_legend(override.aes = list(linetype=lty, color=co, alpha=alp)))
    }

    p <- p + guides(shape = guide_legend(override.aes = list(size = 5, shape=c(21,21), colour=c("red", "black"), fill=c("red", "black"))))
    p <- p + scale_shape_manual(values=c(19, 19), name = TeX("$\\chi^2_{homogeneity}$"), labels = c("p > 0.05", "p < 0.05"))
  }


  return(p)
}
