#' Plot forest
#' @author Matthijs Knigge
#'
#' @param iv vector of causal effect estimate
#' @param iv.se vector of standard error of causal effect estimate
#' @param chochran.Q vector indicating if the SNPs is removed by Chochran's Q test due to pleiotropic effects. Default is NULL.
#' @param ivw the beta estimate of inverse variance weighted method
#' @param ivw.se  the standard error of the inverse variance weighted method
#' @param egger the beta estimate of egger method
#' @param egger.se Standard error of MR egger estimate
#' @param ivw.Q numeric after correcting for pleotropic effect this is the re-calculated inverse variance weighted estimate. Default is NULL.
#' @param egger.Q numeric after correcting for pleotropic effect this is the re-calculated egger estimate. Default is NULL.
#' @param ivw.se.Q numeric after correcting for pleotropic effect this is the re-calculated standard error of the inverse variance weighted estimate. Default is NULL.
#' @param egger.se.Q numeric after correcting for pleotropic effect this is the re-calculated standard error of the egger estimate. Default is NULL.
#'
#' @keywords forest
#' @export
#' @examples
#' mr.forest.plot()
#'
#' @return forest plot
mr.forest.plot <- function(SNP, iv, iv.se, chochran.Q = NULL, ivw = NULL, ivw.se = NULL, egger = NULL, egger.se = NULL, ivw.Q = NULL, ivw.se.Q = NULL, egger.Q = NULL, egger.se.Q = NULL,
                           outcome.name, exposure.name){
  require(ggplot2); require(latex2exp);
  # check if Chochran's Q is used on data set
  if(sum(chochran.Q) == 0){
    chochran.Q <- NULL
    ivw.se.Q <- NULL
    egger.Q <- NULL
    ivw.Q <- NULL
    egger.se.Q <- NULL
  }

  if(!is.null(chochran.Q)){
    # create frame
    df <- data.frame(SNP = SNP, iv = iv, iv.se = iv.se, chochran.Q = chochran.Q, color = "black", stringsAsFactors=FALSE);
    # different color Chochran's Q test
    df[which(df$chochran.Q == 0), ]$color <- "red"
    # add methods to data.frame
    new.row <- data.frame(SNP = "IVW",                iv = ivw,     iv.se = ivw.se,     chochran.Q = NA,  color = "grey"); df <- rbind(df, new.row)
    new.row <- data.frame(SNP = "Egger",              iv = egger,   iv.se = egger.se,   chochran.Q = NA,  color = "grey"); df <- rbind(df, new.row)
    new.row <- data.frame(SNP = "IVW Chocharn's Q",   iv = ivw.Q,   iv.se = ivw.se.Q,   chochran.Q = NA, color = "darkblue"); df <- rbind(df, new.row)
    new.row <- data.frame(SNP = "Egger Chocharn's Q", iv = egger.Q, iv.se = egger.se.Q, chochran.Q = NA, color = "purple"); df <- rbind(df, new.row)
  }
  if(is.null(chochran.Q)){
    # create data.frame
    df <- data.frame(SNP = SNP, iv = iv, iv.se = iv.se, color = "black", stringsAsFactors=FALSE);
    # add methods to data.frame
    new.row <- data.frame(SNP = "IVW", iv = ivw, iv.se = ivw.se, color = "darkblue"); df <- rbind(df, new.row)
    new.row <- data.frame(SNP = "Egger", iv = egger, iv.se = egger.se, color = "purple"); df <- rbind(df, new.row)
  }

  # initiate plot
  p <- ggplot(data = NULL, aes(x=df$iv, y=df$SNP))
  # points
  p <- p + geom_point(col = df$color)
  # center line
  p <- p + geom_vline(xintercept=0, lty=3)
  # error bar 95 confidence interval
  p <- p + geom_errorbarh(aes(xmin=df$iv - qnorm(0.975)*df$iv.se, xmax=df$iv + qnorm(0.975)*df$iv.se), alpha = 0.5, color = df$color)

  # horizontal line
  p <- p +   geom_hline(aes(yintercept = 0), colour="grey")
  # add labs
  p <- p +  labs(x = TeX("$\\beta_{IV}$"), y = "", title = paste0(outcome.name, " ~ ", exposure.name))
  # stolen theme
  p <- p +  theme_minimal()
  # font size and type
  p <- p + theme(axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=14,face="bold"))

  if(!is.null(chochran.Q)){
    p <- p +  annotate("rect",
                       xmin = egger - qnorm(0.975)*egger.se,
                       xmax = egger + qnorm(0.975)*egger.se,
                       ymin = -Inf,
                       ymax = df$SNP[length(df$SNP)-4],
                       fill = "grey", alpha = .6, color = NA)

    p <- p +  annotate("rect",
                       xmin = ivw - qnorm(0.975)*ivw.se,
                       xmax = ivw + qnorm(0.975)*ivw.se,
                       ymin = -Inf,
                       ymax = df$SNP[length(df$SNP)-4],
                       fill = "grey", alpha = .4, color = NA)

    p <- p +   annotate("rect",
                        xmin = egger.Q - qnorm(0.975)*egger.se.Q,
                        xmax = egger.Q + qnorm(0.975)*egger.se.Q,
                        ymin = -Inf,
                        ymax = df$SNP[length(df$SNP)-4],
                        fill = "steelblue", alpha = .6, color = NA)

    p <- p +  annotate("rect",
                       xmin = ivw.Q - qnorm(0.975)*ivw.se.Q,
                       xmax = ivw.Q + qnorm(0.975)*ivw.se.Q,
                       ymin = -Inf,
                       ymax = df$SNP[length(df$SNP)-4],
                       fill = "purple", alpha = .4, color = NA)
  }
  if(is.null(chochran.Q)){
    p <- p +  annotate("rect",
                       xmin = egger - qnorm(0.975)*egger.se,
                       xmax = egger + qnorm(0.975)*egger.se,
                       ymin = -Inf,
                       ymax = df$SNP[length(df$SNP)-2],
                       fill = "purple", alpha = .6, color = NA)

    p <- p +  annotate("rect",
                       xmin = ivw - qnorm(0.975)*ivw.se,
                       xmax = ivw + qnorm(0.975)*ivw.se,
                       ymin = -Inf,
                       ymax = df$SNP[length(df$SNP)-2],
                       fill = "darkblue", alpha = .4, color = NA)
  }

  return(p)

}
