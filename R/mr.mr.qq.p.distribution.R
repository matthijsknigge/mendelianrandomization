#' MR funnel plot
#' @author Matthijs Knigge
#'
#' @param egger vector of pvalues from egger estimate. Default NULL.
#' @param ivw vector of pvalues from ivw estimate. Default NULL.
#' @param egger.Q vector of pvalues from egger estimate after Chochren's Q test. Default NULL.
#' @param ivw.Q vector of pvalues from ivw estimate after Chochren's Q test. Default NULL.
#' @param CI confidence interval
#'
#' @keywords qqplot
#' @export
#' @examples
#' mr.qq.p.distribution()
#'
#' @return funnel plot
mr.qq.p.distribution <- function(egger = NULL, ivw = NULL, egger.Q = NULL, ivw.Q = NULL, CI = 0.95) {
  require(ggplot2)
  # remove NA's
  if(!is.null(egger)){egger <- egger[!is.na(egger)]};
  # ivw
  if(!is.null(ivw)){ivw <- ivw[!is.na(ivw)]};
  # egger after Q
  if(!is.null(egger.Q)){egger.Q <- egger.Q[!is.na(egger.Q)]};
  # ivw after Q
  if(!is.null(ivw.Q)){ivw.Q <- ivw.Q[!is.na(ivw.Q)]}
  # get max length
  len <- c("egger", "ivw", "egger.Q", "ivw.Q")[which.max(c(length(egger), length(ivw.Q), length(egger.Q), length(ivw.Q)))]
  # poolsize
  N  <- length(get(len))
  # calculate observed vs expected
  p.obs.exp <- function(pval, N){
    # log and sort p-values for observations
    observed <- -log10(sort(pval))
    # the expected values
    expected <- -log10(1:N / N)
    return(list(observed = observed, expected = expected))
  }
  # upper part of confidence interval
  clower   <- -log10(qbeta(CI, 1:N, N - 1:N + 1))
  # lower part of confidence interval
  cupper   <- -log10(qbeta(1 - CI, 1:N, N - 1:N + 1))
  # expression of expected
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  # expression of observed
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  # initiate plot
  p <- ggplot()
  # to shade, first get the maximum term. And draw it once
  max.method <- c("egger", "ivw", "egger.Q", "ivw.Q")[which.max(c(length(egger), length(ivw.Q), length(egger.Q), length(ivw.Q)))]
  # shade
  max.method <- p.obs.exp(pval = get(max.method), N = length(get(max.method)))
  p <- p + geom_ribbon(aes(ymin = clower, ymax = cupper, x = max.method$expected), fill = "grey", alpha = .5)
  # egger
  if(!is.null(egger)){
    # p distribution for egger
    egger.p.obs.exp <- p.obs.exp(pval = egger, N = length(egger))
    # add line
    p <- p + geom_line(aes(egger.p.obs.exp$expected, egger.p.obs.exp$observed), size = 1, color = "darkblue", linetype = "solid", alpha = .3)
  }
  # egger.Q
  if(!is.null(egger.Q)){
    # p distribution for ivw
    egger.Q.p.obs.exp <- p.obs.exp(pval = egger.Q, N = length(egger.Q))
    # add line
    p <- p + geom_line(aes(egger.Q.p.obs.exp$expected, egger.Q.p.obs.exp$observed), size = 1, color = "darkblue", linetype = "solid")
  }
  # IVW
  if(!is.null(ivw)){
    # p distribution for ivw
    ivw.p.obs.exp <- p.obs.exp(pval = ivw, N = length(ivw))
    # add line
    p <- p + geom_line(aes(ivw.p.obs.exp$expected, ivw.p.obs.exp$observed), size = 1, color = "darkblue", linetype = "dashed", alpha = .3)
  }
  # IVW.Q
  if(!is.null(ivw.Q)){
    # p distribution for ivw
    ivw.Q.p.obs.exp <- p.obs.exp(pval = ivw.Q, N = length(ivw.Q))
    # add line
    p <- p + geom_line(aes(ivw.Q.p.obs.exp$expected, ivw.Q.p.obs.exp$observed), size = 1, color = "darkblue", linetype = "dashed")
  }
  # add slope
  p <- p + geom_abline(intercept = 0, slope = 1, alpha = 0.5)
  # add upper CI
  # p <- p + geom_line(aes(expected, cupper), linetype = 2, color = "black")
  # add lower CI
  # p <- p + geom_line(aes(expected, clower), linetype = 2, color = 'black')
  # add labs
  p <- p + xlab(log10Pe)
  p <- p + ylab(log10Po)
  # stolen theme for plot
  p <- p + theme_minimal()
  # font size and type
  p <- p + theme(axis.text=element_text(size=14, face="bold"),
                 axis.title=element_text(size=14,face="bold"))

  return(p)

}
