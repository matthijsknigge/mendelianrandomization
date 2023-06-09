% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr.wald.ratio.R
\name{mr.wald.ratio}
\alias{mr.wald.ratio}
\title{Perform wald-ratio}
\usage{
mr.wald.ratio(By, Bx, By.se, Bx.se)
}
\arguments{
\item{By}{Vector of genetic effects of outcome}

\item{Bx}{Vector of genetic effects of exposure}

\item{By.se}{Standard errors of genetic effects on outcome}

\item{Bx.se}{Standard errors of genetic effects on exposure}
}
\value{
List with the following elements:
        iv: causal effect estimate
        iv.se: standard error
        iv.p: p-value
}
\description{
Perform wald-ratio
}
\examples{
mr.wald.ratio()


}
\author{
Matthijs Knigge
}
\keyword{wald-ratio}
