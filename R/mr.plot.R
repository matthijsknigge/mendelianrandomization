#' Plot MR results
#' @author Matthijs Knigge
#'
#'
#' @param By Vector of genetic effects of outcome
#' @param Bx Vector of genetic effects of exposure
#' @param By.se Standard errors of genetic effects on outcome
#' @param Bx.se Standard errors of genetic effects on exposure
#' @param iv vector of causal effect estimate
#' @param iv.se vector of standard error of causal effect estimate
#' @param ivw numeric causal estimate inverse variance weighted method. Default is NULL.
#' @param egger numeric causal estimate of egger method, this is the slope. Default is NULL.
#' @param egger.i numeric pleiotropic estimate of egger, this is the intersept. Default is NULL.
#' @param chochran.Q vector indicating if the SNPs is removed by Chochran's Q test due to pleiotropic effects. Default is NULL.
#' @param outcome.name string of the outcome name. Default is NULL.
#' @param exposure.name string of the exposure name. Default is NULL.
#' @param ivw.Q numeric after correcting for pleotropic effect this is the re-calculated inverse variance weighted estimate. Default is NULL.
#' @param egger.Q numeric after correcting for pleotropic effect this is the re-calculated egger estimate. Default is NULL.
#' @param egger.i.Q numeric after correcting for pleotropic effect this is the re-calculated egger pleiotropic estimate. Default is NULL.
#' @param legend boolean use show legend, or do not. Default is TRUE
#' @param position position of legend. Default bottom.
#' @param show.stats boolean show statistics summary. Default is TRUE.
#' @param egger.p.fdr value indicating slope significant different from zero after FDR correction.
#' @param ivw.p.fdr value indicating slope significant different from zero after FDR correction.
#' @param egger.i.p value indicating intercept significant different from zero after FDR correction.
#'
#' @keywords mr.plot
#' @export
#' @examples
#' mr.plot <- function(By, Bx, By.se, Bx.se, iv, iv.se, ivw = NULL, egger = NULL, egger.i = NULL, chochran.Q = NULL,
#'                     ivw.Q = NULL, egger.Q = NULL, egger.i.Q = NULL, outcome.name, exposure.name, legend = TRUE)
#' @return plot object
mr.plot <- function(By, Bx, By.se, Bx.se, iv, iv.se, ivw = NULL, egger = NULL, egger.i = NULL, chochran.Q = NULL,
                    ivw.Q = NULL, egger.Q = NULL, egger.i.Q = NULL, egger.p.fdr = NULL, ivw.p.fdr = NULL, egger.i.p = NULL,
                    outcome.name, exposure.name, legend = TRUE, position = "bottom", show.stats = TRUE){
  require(ggplot2); require("RColorBrewer"); require(latex2exp); require(cowplot); require(gridExtra)
  # check if Chochran's Q is used on data set
  if(sum(chochran.Q) == length(By)){
    chochran.Q <- NULL
  }
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
  # apply Chochran's Q test
  if(!is.null(chochran.Q)){
  p <- p + geom_point(aes(colour = iv.z, shape = factor(chochran.Q)), alpha = 0.8, size=4)
  p <- p + geom_point(data=NULL, aes(x=Bx[which(chochran.Q == 0)], y=By[which(chochran.Q == 0)]), colour="red", size=4.1)
  p <- p + guides(shape = guide_legend(override.aes = list(size = 5, shape=c(21,21), colour="white", fill=c("red", brewer.pal(9, "Blues")[7]))))
  p <- p + scale_shape_manual(values=c(19, 19), name = TeX("$\\chi^2_{homogeneity}$"), labels = c("p > 0.05", "p < 0.05"))
  }

  # draw inverse variance weighted effect after Chochran's Q test
  if(!is.null(ivw.Q)){
    p <- p + geom_abline(aes(intercept = 0, slope = ivw.Q, linetype = "IVW.Q"), color = "darkblue")
  }
  # draw egger method effect after Chochran's Q test
  if(!is.null(egger.Q) & !is.null(egger.i.Q)){
    p <- p + geom_abline(aes(intercept = egger.i.Q, slope = egger.Q, linetype = "MRegger.Q"), color = "darkblue")
  }
  # manually change legend or override the existing one
  p <- p + scale_linetype_manual(name = "Method", values = c(IVW = "dashed", MRegger = "solid", IVW.Q = "dashed", MRegger.Q = "solid"))

  # waring! extremely ugly code. Override legends is not flexible, this is a dirty workaround.
  if(!is.null(ivw.Q)){lty <- c(IVW = "dashed", IVW.Q = "dashed", MRegger = "solid"); alp <- c(0.1,1,0.1); co  <- c("darkblue", "darkblue", "darkblue")}
  if(!is.null(egger.Q) & !is.null(egger.i.Q)){lty <- c(IVW = "dashed", MRegger.Q = "solid", MRegger = "solid"); alp <- c(0.1,1,0.1); co  <- c("darkblue", "darkblue", "darkblue")}
  if(!is.null(ivw.Q) & !is.null(egger.Q) & !is.null(egger.i.Q)){lty <- c(IVW = "dashed", IVW.Q = "dashed", MRegger = "solid", MRegger.Q = "solid"); alp <- c(0.1,1,0.1, 1); co  <- c("darkblue", "darkblue", "darkblue", "darkblue")}
  if(!is.null(ivw.Q) | !is.null(egger.Q) & !is.null(egger.i.Q)){p <- p + guides( linetype = guide_legend(override.aes = list(linetype=lty, color=co, alpha=alp)))}
  # legend for z-score
  p <- p + scale_colour_gradientn(colours = brewer.pal(9, "Blues")[2:9], name = TeX('$\\frac{\\hat{\\beta}_{IV}}{\\sigma_{\\hat{\\beta}_{IV}}}$'))

  # stolen theme for plot
  p <- p + theme_minimal()
  # labs and titles
  p <- p + labs(y= TeX('$\\hat{\\beta}_{outcome}$'),
                x= TeX('$\\hat{\\beta}_{exposure}$'),
                title = paste0(outcome.name, " ~ ", exposure.name))
  # font size and type
  p <- p + theme(axis.text=element_text(size=14, face="bold"),
                 axis.title=element_text(size=14,face="bold"))

  # add rug
  p <- p + geom_rug(aes(color = iv.z))
  # add rug for Chochran's Q
  if(!is.null(chochran.Q)){
   p <- p + geom_rug(data=NULL, aes(x=Bx[which(chochran.Q == 0)], y=By[which(chochran.Q == 0)]), colour="red")
  }
  # position
  p <- p + theme(legend.position = position)
  # remove legends
  if(legend == FALSE){
    p <- p + theme(legend.position="none")
  }
  # show stats
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(fill = blues9[1:4], col=NA),
              fg_params=list(fontface=3)),
    colhead=list(fg_params=list(col="navyblue", fontface=4L)),
    rowhead=list(fg_params=list(col="orange", fontface=3L)))

  if(show.stats == TRUE){
    mytable <- cbind(c("nSNP"," Egger FDR","IVW FDR","Horizontal pleiotropy"), c(length(By) , signif(egger.p.fdr, digits = 3), signif(ivw.p.fdr, digits = 3), signif(egger.i.p, digits = 3)))
    p <- p + annotation_custom(tableGrob(mytable, theme = tt3),
                               xmin = layer_scales(p)$x$range$range[2]*.8,
                               xmax = layer_scales(p)$x$range$range[2]*.8,
                               ymax = layer_scales(p)$y$range$range[1]*.8,
                               ymin = layer_scales(p)$y$range$range[1]*.8)
  }

  # for some #*$%! reason the variable must be in the global scope of R to be inserted into a grid??????????
  Bx <<- Bx; By <<- By; chochran.Q <<- chochran.Q
  # x-axis
  xdens <- axis_canvas(p, axis = "x")
  if(is.null(chochran.Q)){
   xdens <- xdens + geom_density(data = NULL, aes(x = Bx), fill = "steelblue", alpha = 1, size = 0.2)
  }
  if(!is.null(chochran.Q)){
    xdens <- xdens + geom_density(data = NULL, aes(x = Bx[which(chochran.Q == 1)]), fill = "steelblue", alpha = 1, size = 0.2)
    xdens <- xdens + geom_density(data = NULL, aes(x = Bx[which(chochran.Q == 0)]), fill = "red", alpha = .3, size = 0.2)
  }
  # y-axis
  ydens <- axis_canvas(p, axis = "y", coord_flip = T)
  if(is.null(chochran.Q)){
   ydens <- ydens + geom_density(data = NULL, aes(x = By), fill = "steelblue",  alpha = 1, size = 0.2)
  }
   if(!is.null(chochran.Q)){
   ydens <- ydens + geom_density(data = NULL, aes(x = By[which(chochran.Q == 1)]), fill = "steelblue",  alpha = 1, size = 0.2)
   ydens <- ydens + geom_density(data = NULL, aes(x = By[which(chochran.Q == 0)]), fill = "red",  alpha = .3, size = 0.2)
  }
  ydens <- ydens + coord_flip()
  # insert into grid
  p <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  p <- insert_yaxis_grob(p, ydens, grid::unit(.2, "null"), position = "right")

  # remove from global environment
  rm(Bx, envir = .GlobalEnv); rm(By, envir = .GlobalEnv); rm(chochran.Q, envir = .GlobalEnv)
  return(p)
}
